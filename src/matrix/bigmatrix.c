/*
 * bigmatrix.c
 *
 * Matrix arithmetic over bignums.
 *
 * This file implements matrices with bignum entries: initialization,
 * freeing, copying, element/row/column access, structural check
 * predicates (square, zero, identity, diagonal, triangular, symmetric),
 * printing, component-wise addition, schoolbook multiplication,
 * matrix-vector products, the Hadamard bound on the determinant, and
 * the determinant (computed via the RNS path; a direct Bareiss
 * implementation is kept below it, currently disabled).
 *
 * A bigmatrix is a row-major dynamic array of bignums with row and
 * column counts.
 *
 * Copyright (C) 2026 Diego Strebel
 * SPDX-License-Identifier: GPL-3.0-or-later
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "../../include/bigmatrix.h"

#include <omp.h>
#include <stdio.h>
#include <string.h>

#include "../../include/bigmath.h"
#include "../../include/bigpoly.h"
#include "../../include/bigrns.h"
#include "../../include/primes.h"

void bigmatrix_init(bigmatrix* M, u64 r, u64 c)
<%
  M->c_size = c;
  M->r_size = r;

  M->data = malloc(r * c * sizeof(bignum));

  for (u64 i = 0; i < r * c; i++) {
    bn_init(&M->data[i]);
  }
}

void bigmatrix_free(bigmatrix* A)
{
  if (!A || !A->data) return;

  for (u64 i = 0; i < A->r_size * A->c_size; i++) {
    bn_free(&A->data[i]);
  }

  free(A->data);

  A->data = NULL;
}

void bigmatrix_print_python(const bigmatrix* A)
{
  printf("[");  // Start outer list
  for (u64 i = 0; i < A->r_size; i++) {
    printf("[");  // Start row
    for (u64 j = 0; j < A->c_size; j++) {
      bn_print(GET(A, i, j));
      if (j < A->c_size - 1) printf(", ");
    }
    printf("]");  // End row
    if (i < A->r_size - 1) printf(",\n ");
  }
  printf("]\n");  // End outer list
}

void bigmatrix_copy(bigmatrix* R, const bigmatrix* A)
{
  // check if sizes match
  if (A->c_size != R->c_size || A->r_size != R->r_size) return;

  for (u64 i = 0; i < A->c_size * A->r_size; i++) {
    bn_copy(&R->data[i], &A->data[i]);
  }
}

void bigmatrix_get(bignum* R, const bigmatrix* A, u64 r, u64 c)
{
  if (r >= A->r_size || c >= A->c_size) return;

  bn_copy(R, GET(A, r, c));
}

void bigmatrix_set(bigmatrix* A, const bignum* a, u64 r, u64 c)
{
  if (r >= A->r_size || c >= A->c_size) return;

  bn_copy(GET(A, r, c), a);
}

void bigmatrix_get_col(bigvector* c, const bigmatrix* A, u64 col)
{
  if (c->size != A->r_size) printf("Size does not match \n");

  for (u64 i = 0; i < A->r_size; i++) {
    bigvector_set(c, GET(A, i, col), i);
  }
}

void bigmatrix_get_row(bigvector* r, const bigmatrix* A, u64 row)
{
  if (r->size != A->c_size) printf("Size does not match \n");

  for (u64 i = 0; i < A->c_size; i++) {
    bigvector_set(r, GET(A, row, i), i);
  }
}

bool bigmatrix_is_square(const bigmatrix* A) { return A->r_size == A->c_size; }

bool bigmatrix_is_zero(const bigmatrix* A)
{
  // single sequential pass over the flat row-major storage; bn_is_zero is
  // an O(1) size/limb inspection, so no temporary bignums are needed
  u64 total = A->r_size * A->c_size;
  for (u64 i = 0; i < total; i++) {
    if (!bn_is_zero(&A->data[i])) {
      return false;
    }
  }
  return true;
}

bool bigmatrix_is_identity(const bigmatrix* A)
{
  // dimension short-circuit: a non-square matrix cannot be the identity,
  // bail out before touching any entry
  if (!bigmatrix_is_square(A)) {
    return false;
  }

  u64 n = A->r_size;
  for (u64 i = 0; i < n; i++) {
    bignum* row = &A->data[i * A->c_size];
    for (u64 j = 0; j < n; j++) {
      if (i == j) {
        // diagonal entry must be exactly 1 (bn_is_one rejects -1)
        if (!bn_is_one(&row[j])) {
          return false;
        }
      } else {
        // every off-diagonal entry must be 0
        if (!bn_is_zero(&row[j])) {
          return false;
        }
      }
    }
  }
  return true;
}

bool bigmatrix_is_diagonal(const bigmatrix* A)
{
  // only the off-diagonal sub-regions can violate: in row i those are the
  // columns left of the diagonal, j < i (capped at c_size), and the
  // columns right of it, j > i (empty once i >= c_size). The diagonal
  // entry (i, i) itself is never inspected
  for (u64 i = 0; i < A->r_size; i++) {
    bignum* row = &A->data[i * A->c_size];

    // left of the diagonal: j in [0, min(i, c_size))
    u64 left = MIN(i, A->c_size);
    for (u64 j = 0; j < left; j++) {
      if (!bn_is_zero(&row[j])) {
        return false;
      }
    }

    // right of the diagonal: j in (i, c_size)
    for (u64 j = i + 1; j < A->c_size; j++) {
      if (!bn_is_zero(&row[j])) {
        return false;
      }
    }
  }
  return true;
}

bool bigmatrix_is_upper_triangular(const bigmatrix* A)
{
  // only the strict lower triangle (row > col) can violate; row 0 has no
  // entries below the diagonal, so the scan starts at i = 1
  for (u64 i = 1; i < A->r_size; i++) {
    bignum* row = &A->data[i * A->c_size];

    // below the diagonal in row i: j in [0, min(i-1, c_size-1)], i.e.
    // j < min(i, c_size) — the cap handles tall (c < r) matrices
    u64 limit = MIN(i, A->c_size);
    for (u64 j = 0; j < limit; j++) {
      if (!bn_is_zero(&row[j])) {
        return false;
      }
    }
  }
  return true;
}

bool bigmatrix_is_lower_triangular(const bigmatrix* A)
{
  // only the strict upper triangle (col > row) can violate; in row i the
  // first candidate is column i+1, and the inner loop is empty once
  // i >= c_size, so rows below the diagonal region cost nothing
  for (u64 i = 0; i < A->r_size; i++) {
    bignum* row = &A->data[i * A->c_size];
    for (u64 j = i + 1; j < A->c_size; j++) {
      if (!bn_is_zero(&row[j])) {
        return false;
      }
    }
  }
  return true;
}

bool bigmatrix_is_symmetric(const bigmatrix* A)
{
  // dimension short-circuit: a non-square matrix cannot be symmetric,
  // bail out before touching any entry
  if (!bigmatrix_is_square(A)) {
    return false;
  }

  u64 n = A->r_size;
  // walk only the strict upper triangle: each mirror pair (i, j)/(j, i)
  // with i < j is compared exactly once; row i is swept sequentially,
  // and GET(A, j, i) is a constant offset into row j's contiguous block
  for (u64 i = 0; i < n; i++) {
    bignum* row_i = &A->data[i * A->c_size];
    for (u64 j = i + 1; j < n; j++) {
      if (bn_cmp(&row_i[j], GET(A, j, i)) != 0) {
        return false;
      }
    }
  }
  return true;
}

void bigmatrix_scalar(bigmatrix* R, const bigmatrix* A, const bignum* a)
{
  u64 total = A->r_size * A->c_size;

  for (u64 x = 0; x < total; x++) {
    bn_mul(&R->data[x], &A->data[x], a);
  }
}

void bigmatrix_div_exact_scalar(bigmatrix* R, const bigmatrix* A,
                                const bignum* a)
{
  u64 total = A->r_size * A->c_size;

  for (u64 x = 0; x < total; x++) {
    bn_div_exact(&R->data[x], &A->data[x], a);
  }
}

void bigmatrix_mv(bigvector* r, const bigmatrix* A, bigvector* v)
{
  if (r->size != A->r_size || v->size != A->c_size) {
    printf("Size does not match \n");
    return;
  }

  bigvector tmp;
  bigvector_init(&tmp, A->c_size);

  bignum dot;
  bn_init(&dot);

  for (u64 i = 0; i < A->r_size; i++) {
    bigmatrix_get_row(&tmp, A, i);
    bigvector_dot(&dot, v, &tmp);
    bigvector_set(r, &dot, i);
  }

  bn_free(&dot);
  bigvector_free(&tmp);
}

void bigmatrix_vm(bigvector* r, const bigmatrix* A, bigvector* v)
{
  if (r->size != A->c_size || v->size != A->r_size) {
    printf("Size does not match \n");
    return;
  }

  bigvector tmp;
  bigvector_init(&tmp, A->r_size);

  bignum dot;
  bn_init(&dot);

  for (u64 i = 0; i < A->c_size; i++) {
    bigmatrix_get_col(&tmp, A, i);
    bigvector_dot(&dot, v, &tmp);
    bigvector_set(r, &dot, i);
  }

  bn_free(&dot);
  bigvector_free(&tmp);
}

void bigmatrix_hadamard(bignum* r, const bigmatrix* A)
{
  bn_set_u64(r, 1);

  bigvector v;
  bigvector_init(&v, A->r_size);
  bignum tmp;
  bn_init(&tmp);

  for (u64 i = 0; i < A->c_size; i++) {
    bigmatrix_get_col(&v, A, i);

    bigvector_norm(&tmp, &v);
    // add 1 since norm uses isqrt
    bn_add_u64(&tmp, &tmp, 1);
    bn_mul(r, r, &tmp);
  }

  bn_free(&tmp);
  bigvector_free(&v);
}

void bigmatrix_add(bigmatrix* R, const bigmatrix* A, const bigmatrix* B)
{
  // check if sizes match
  if (A->c_size != B->c_size || A->r_size != B->r_size) return;

  for (u64 i = 0; i < A->c_size * A->r_size; i++) {
    bn_add(&R->data[i], &A->data[i], &B->data[i]);
  }

  return;
}

void bigmatrix_sub(bigmatrix* R, const bigmatrix* A, const bigmatrix* B)
{
  // check if sizes match
  if (A->c_size != B->c_size || A->r_size != B->r_size) return;

  for (u64 i = 0; i < A->c_size * A->r_size; i++) {
    bn_sub(&R->data[i], &A->data[i], &B->data[i]);
  }

  return;
}
// maximum printed width of a single entry in the truncated modes
#define BIGMATRIX_PRINT_MAX_WIDTH 16

// upper bound on the printed width (digits plus sign) of a single entry,
// computed without a decimal conversion: ceil(bits * log10(2)) using the
// rational over-approximation 30103/100000 > log10(2), so the result is
// always at least the true digit count
static size_t bigmatrix_entry_width(const bignum* x)
{
  if (bn_is_zero(x)) return 1;

  size_t bits = (size_t)bn_bit_length(x);
  size_t w = (bits * 30103 + 99999) / 100000;
  if (x->is_neg) w++;
  return w;
}

static void bigmatrix_print_impl(const bigmatrix* A, bool full, bool tail)
{
  if (A->r_size == 0 || A->c_size == 0) return;

  // per-column width: widest entry of the column, capped in the truncated
  // modes. A single outlier (e.g. the last invariant factor of a Smith form
  // sitting in a column of zeros) must not inflate the column, so in the
  // truncated modes outliers are excluded from the column width; they are
  // rendered at the capped width and simply overflow the column to the
  // right
  size_t cap = BIGMATRIX_PRINT_MAX_WIDTH;
  size_t* colw = calloc(A->c_size, sizeof(size_t));
  for (u64 i = 0; i < A->r_size; i++) {
    for (u64 j = 0; j < A->c_size; j++) {
      size_t xw = bigmatrix_entry_width(GET(A, i, j));
      if (!full && xw > cap) continue;
      if (xw > colw[j]) colw[j] = xw;
    }
  }

  for (u64 i = 0; i < A->r_size; i++) {
    for (u64 j = 0; j < A->c_size; j++) {
      char* s = bn_to_string(GET(A, i, j));
      size_t len = strlen(s);
      size_t w = colw[j];

      if (full || len <= w) {
        printf("%*s", (int)w, s);
      } else if (len <= cap) {
        // fits within the cap: print it whole, overflowing by its own width
        printf("%s", s);
      } else if (tail) {
        size_t sign = (s[0] == '-') ? 1 : 0;
        if (sign) fputc('-', stdout);
        size_t keep = cap - 1 - sign;
        printf("\u2026%.*s", (int)keep, s + len - keep);
      } else {
        printf("%.*s\u2026", (int)(cap - 1), s);
      }

      free(s);
      printf(" ");
    }
    printf("\n");
  }

  free(colw);
}

void bigmatrix_print(const bigmatrix* A)
{
  bigmatrix_print_impl(A, false, false);
}

void bigmatrix_print_full(const bigmatrix* A)
{
  bigmatrix_print_impl(A, true, false);
}

void bigmatrix_print_tail(const bigmatrix* A)
{
  bigmatrix_print_impl(A, false, true);
}

void bigmatrix_mul(bigmatrix* R, const bigmatrix* A, const bigmatrix* B)
{
  if (A->c_size != B->r_size) return;

  // Transpose B
  bigmatrix B_T;
  bigmatrix_init(&B_T, B->c_size, B->r_size);

  for (u64 i = 0; i < B->r_size; i++) {
    for (u64 j = 0; j < B->c_size; j++) {
      bn_copy(&B_T.data[j * B_T.c_size + i], &B->data[i * B->c_size + j]);
    }
  }
#pragma omp parallel
  {
    bignum sum, tmp;
    bn_init_multi(&sum, &tmp, NULL);

#pragma omp for collapse(2) schedule(dynamic)
    for (u64 i = 0; i < A->r_size; i++) {
      for (u64 j = 0; j < B_T.r_size; j++) {
        bn_set_i64(&sum, 0);

        for (u64 k = 0; k < A->c_size; k++) {
          // tmp = A[i][k] * B_T[j][k]
          bn_mul(&tmp, &A->data[i * A->c_size + k],
                 &B_T.data[j * B_T.c_size + k]);
          bn_add(&sum, &sum, &tmp);
        }

        bigmatrix_set(R, &sum, i, j);
      }
    }

    bn_free(&sum);
    bn_free(&tmp);
  }

  bigmatrix_free(&B_T);
}

void bigmatrix_det(bignum* d, const bigmatrix* A)
{
  if (A->c_size != A->r_size) {
    printf("Matrix must be square\n");
    return;
  }
  u64 n = A->c_size;

  bigmatrix T;
  bigmatrix_init(&T, A->c_size, A->r_size);
  bigmatrix_copy(&T, A);

  u64 k = rns_estimate_determinant(&T);

  ctx_rns ctx;
  rns_context_init(&ctx, RNS_PRIMES, k);

  bigmatrix_det_rns(d, &T, &ctx);

  rns_context_free(&ctx);
  bigmatrix_free(&T);

  return;
}

void bigmatrix_neg(const bigmatrix* A)
{
  for (u64 i = 0; i < A->r_size; i++) {
    for (u64 j = 0; j < A->c_size; j++) {
      bn_neg(GET(A, i, j), GET(A, i, j));
    }
  }
}

void bigmatrix_trace(bignum* t, const bigmatrix* A)
{
  assert(A->c_size == A->r_size);
  u64 n = A->c_size;
  bn_set_u64(t, 0);
  for (u64 i = 0; i < n; i++) {
    bn_add(t, t, GET(A, i, i));
  }
}

void bigmatrix_swap(bigmatrix* a, bigmatrix* b)
{
  bigmatrix t = *a;
  *a = *b;
  *b = t;
}

void bigmatrix_charpoly_adj(bigpoly* p, bigmatrix* J, const bigmatrix* A)
{
  assert(A->c_size == A->r_size);

  // 1. [Initialize]
  u64 n = A->c_size;

  bigmatrix C, T;
  bigmatrix_init(&C, n, n);
  bigmatrix_init(&T, n, n);
  bigmatrix_id(&C, n);

  printf("A: \n");
  bigmatrix_print(A);
  printf("\n");
  printf("C: \n");
  bigmatrix_print(&C);
  printf("\n");

  bignum temp, trace;
  bn_init_multi(&temp, &trace, NULL);
  bignum* a = malloc((n + 1) * sizeof(bignum));
  for (u64 i = 0; i <= n; i++) {
    bn_init(&a[i]);
  }
  bn_set_u64(&a[0], 1);

  // 2. [Finished?]
  for (u64 i = 1; i <= n; i++) {
    // 3. [Compute next a_i and C]
    bigmatrix_mul(&T, A, &C);

    bigmatrix_trace(&trace, &T);

    bn_set_u64(&temp, i);
    // bn_div(&a[i], &trace, &temp);
    bn_divmod_u64(&a[i], &trace, i);

    bn_neg(&a[i], &a[i]);

    if (J && i == n) {
      bigmatrix_copy(J, &C);
      if ((n - 1) & 1) {
        bigmatrix_neg(J);
      }
    }
    bigmatrix_copy(&C, &T);
    for (u64 j = 0; j < n; j++) {
      bn_add(&temp, GET(&C, j, j), &a[i]);
      bigmatrix_set(&C, &temp, j, j);
    }
  }

  bigpoly_alloc(p, n + 1);

  for (u64 i = 0; i <= n; i++) {
    bn_println(&a[i]);
    // bn_init(&p->coeff[n - i]);
    bn_copy(&p->coeff[n - i], &a[i]);
    bn_free(&a[i]);
  }

  p->deg = n;

  free(a);
  bn_free_multi(&temp, &trace, NULL);
  bigmatrix_free(&C);
  bigmatrix_free(&T);
}

void bigmatrix_gcd_all(bignum* g, const bigmatrix* A)
{
  u64 size = A->c_size * A->r_size;

  if (size == 0) {
    bn_set_u64(g, 1);
    return;
  }

  // start with first element
  bn_copy(g, &A->data[0]);

  for (u64 i = 1; i < size; i++) {
    if (bn_is_eq_i64(g, 1)) return;

    if (bn_is_zero(&A->data[i])) continue;

    bn_gcd(g, g, &A->data[i]);
  }
}

void bigmatrix_id(bigmatrix* I, const u64 n)
{
  bignum a, b;
  bn_init_multi(&a, &b, NULL);
  bn_set_u64(&a, 1);
  bn_set_u64(&b, 0);

  for (u64 i = 0; i < n; i++) {
    for (u64 j = 0; j < n; j++) {
      if (i == j) {
        bigmatrix_set(I, &a, i, j);
      } else {
        bigmatrix_set(I, &b, i, j);
      }
    }
  }

  bn_free_multi(&a, &b, NULL);
}

void bigmatrix_neg_col(bigmatrix* A, u64 j)
{
  for (u64 i = 0; i < A->r_size; i++) {
    bn_neg(GET(A, i, j), GET(A, i, j));
  }
}

void bigmatrix_neg_row(bigmatrix* A, u64 i)
{
  for (u64 j = 0; j < A->c_size; j++) {
    bn_neg(GET(A, i, j), GET(A, i, j));
  }
}

static bool bigmatrix_row_prefix_is_zero(const bigmatrix* A, u64 i, u64 k)
{
  for (u64 j = 0; j < k; j++) {
    if (!bn_is_zero(GET(A, i, j))) {
      return false;
    }
  }
  return true;
}

static u64 bigmatrix_min_abs_index(const bigmatrix* A, int i, int k)
{
  int j, j0 = -1;
  for (j = 0; j <= k; j++) {
    if (bn_is_zero(GET(A, i, j))) continue;

    if (j0 == -1 || bn_cmp_abs(GET(A, i, j), GET(A, i, j0)) < 0) j0 = j;
  }

  return j0;
}

static void bigmatrix_swap_col(bigmatrix* A, u64 k1, u64 k2)
{
  for (u64 r = 0; r < A->r_size; r++) {
    bn_swap(GET(A, r, k1), GET(A, r, k2));
  }
}

void bigmatrix_hermite(bigmatrix* W, const bigmatrix* A)
{
  // 1. [Initialize]
  u64 m = A->r_size;
  u64 n = A->c_size;
  if (m == 0 || n == 0) return;

  i64 i = m - 1;
  i64 k = n - 1;
  u64 l = 0;
  if (m > n) {
    l = m - n;
  }

  bigmatrix A_work;
  bigmatrix_init(&A_work, m, n);
  bigmatrix_copy(&A_work, A);

  bignum b, q, rem, temp, one, t2;
  bn_init_multi(&b, &q, &rem, &temp, &one, &t2, NULL);
  bn_set_u64(&one, 1);

// 2. [Row finished?]
step2:
  if (bigmatrix_row_prefix_is_zero(&A_work, i, k)) {
    if (GET(&A_work, i, k)->is_neg) {
      bigmatrix_neg_col(&A_work, k);
    }
    goto step5;
  }

  // 3. [Choose non-zero entry]
  u64 j_0 = bigmatrix_min_abs_index(&A_work, i, k);
  if (j_0 < k) {
    bigmatrix_swap_col(&A_work, k, j_0);
  }
  if (GET(&A_work, i, k)->is_neg) {
    bigmatrix_neg_col(&A_work, k);
  }

  bn_copy(&b, GET(&A_work, i, k));

  // 4. [Reduce]
  for (u64 j = 0; j < k; j++) {
    bn_div_euclid(&q, GET(&A_work, i, j), &b);
    for (u64 x = 0; x < m; x++) {
      bn_mul(&temp, &q, GET(&A_work, x, k));
      bn_sub(&t2, GET(&A_work, x, j), &temp);
      bigmatrix_set(&A_work, &t2, x, j);
    }
  }
  goto step2;

// 5. [Final reductions]
step5:
  bn_copy(&b, GET(&A_work, i, k));
  if (bn_is_zero(&b)) {
    k++;
    goto step6;
  }

  for (u64 j = k + 1; j < n; j++) {
    bn_div_euclid(&q, GET(&A_work, i, j), &b);
    for (u64 x = 0; x < m; x++) {
      bn_mul(&temp, &q, GET(&A_work, x, k));
      bn_sub(&t2, GET(&A_work, x, j), &temp);
      bigmatrix_set(&A_work, &t2, x, j);
    }
  }

// 6. [Finished?]
step6:
  if (i != l) {
    i--;
    k--;
    goto step2;
  }

  bigmatrix_free(W);
  bigmatrix_init(W, m, n - k);

  // copy result
  for (u64 j = 0; j < n - k; j++) {
    for (u64 x = 0; x < m; x++) {
      bn_copy(GET(W, x, j), GET(&A_work, x, j + k));
    }
  }

  bn_free_multi(&b, &q, &rem, &temp, &one, NULL);
  bigmatrix_free(&A_work);
}

static void bigmatrix_hnf_centered_mod(bignum* r, const bignum* x,
                                       const bignum* R, const bignum* half)
{
  if (bn_cmp_abs(x, half) <= 0) {
    bn_copy(r, x);
    return;
  }

  bn_mod_pos(r, x, R);
  if (bn_cmp(r, half) > 0) {
    bn_sub(r, r, R);
  }
}

// Cohen's "Important Remark" after Algorithm 2.4.5: the Euclidean steps of
// the HNF algorithms need the Bezout pair (u, v) of u*a + v*b = d with
// small coefficients. All solutions are (u0 + t*(b/d), v0 - t*(a/d)), so
// this rewrites the pair in place to the unique one with v in the centered
// half-period (-|a/d|/2, |a/d|/2] (which minimizes |v|); when a | b
// (d = |a|) the book's essential condition v = 0, u = sign(a) applies
// instead. The identity u*a + v*b = d is preserved in all cases.
static void hnf_bezout_minimal(bignum* u, bignum* v, const bignum* d,
                               const bignum* a, const bignum* b)
{
  bignum p, pmag, half, quot, r, t, two;
  bn_init_multi(&p, &pmag, &half, &quot, &r, &t, &two, NULL);

  if (bn_is_zero(a)) {
    // v*sign(b) = d/|b| = 1 is forced; u is free, take u = 0
    bn_set_u64(u, 0);
    if (bn_is_zero(b)) {
      bn_set_u64(v, 0);
    } else if (b->is_neg) {
      bn_set_i64(v, -1);
    } else {
      bn_set_u64(v, 1);
    }
    goto cleanup;
  }

  if (bn_is_zero(b)) {
    // u*sign(a) = 1 is forced; v = 0
    if (a->is_neg) {
      bn_set_i64(u, -1);
    } else {
      bn_set_u64(u, 1);
    }
    bn_set_u64(v, 0);
    goto cleanup;
  }

  bn_copy(&p, a);
  p.is_neg = false;

  // a | b  =>  d = |a|  =>  v = 0, u = sign(a)  (the book's essential
  // condition; it also keeps the Euclidean step from stalling)
  if (bn_cmp(&p, d) == 0) {
    if (a->is_neg) {
      bn_set_i64(u, -1);
    } else {
      bn_set_u64(u, 1);
    }
    bn_set_u64(v, 0);
    goto cleanup;
  }

  // centered half-period: v = v0 - t*(a/d) with |a/d| >= 2 and
  // v in (-|a/d|/2, |a/d|/2]; recompute u = (d - v*b)/a (exact)
  bn_div_exact(&p, a, d);
  bn_copy(&pmag, &p);
  pmag.is_neg = false;
  bn_set_u64(&two, 2);
  bn_div(&half, &pmag, &two);

  // r = v0 mod |p| in [0, |p|)
  bn_divmod(&quot, &r, v, &pmag);  // |r| < |p|, sign(r) = sign(v0)
  if (r.is_neg) {
    bn_add(&r, &r, &pmag);
  }
  if (bn_cmp(&r, &half) > 0) {
    bn_sub(v, &r, &pmag);  // v in (-|p|/2, 0)
  } else {
    bn_copy(v, &r);  // v in [0, |p|/2]
  }

  bn_mul(&t, v, b);
  bn_sub(&t, d, &t);
  bn_div_exact(u, &t, a);

cleanup:
  bn_free_multi(&p, &pmag, &half, &quot, &r, &t, &two, NULL);
}

void bigmatrix_hermite_mod_d(bigmatrix* W, const bigmatrix* A, const bignum* D)
{
  // 1. [Initialize]
  u64 m = A->r_size;
  u64 n = A->c_size;
  if (m == 0 || n == 0) return;

  if (m > n) {
    bigmatrix_hermite_gcd(W, A);
    return;
  }

  bigmatrix A_work;
  bigmatrix_init(&A_work, m, n);
  bigmatrix_copy(&A_work, A);

  bignum R, u, v, d, q, t1, t2, t3, half;
  bn_init_multi(&R, &u, &v, &d, &q, &t1, &t2, &t3, &half, NULL);
  bn_copy(&R, D);
  R.is_neg = false;
  bn_rshift(&half, &R, 1);

  // reduce the input entries into (-R/2, R/2]; all column operations are
  // taken modulo R, so only the classes modulo R matter
  for (u64 x = 0; x < m * n; x++) {
    bigmatrix_hnf_centered_mod(&A_work.data[x], &A_work.data[x], &R, &half);
  }

  bignum* B;
  B = malloc(m * sizeof(bignum));
  for (u64 x = 0; x < m; x++) {
    bn_init(&B[x]);
  }

  bigmatrix W_local;
  bigmatrix_init(&W_local, m, m);

  for (i64 i = (i64)m - 1; i >= 0; i--) {
    i64 k = (i64)n - (i64)m + i;

    // 2. [Check zero] / 3. [Euclidean step]
    for (i64 j = k - 1; j >= 0; j--) {
      if (bn_is_zero(GET(&A_work, i, j))) {
        continue;
      }

      bn_gcd_extended_lehmer(&u, &v, &d, GET(&A_work, i, k),
                             GET(&A_work, i, j));
      hnf_bezout_minimal(&u, &v, &d, GET(&A_work, i, k), GET(&A_work, i, j));

      // B = u*A_k + v*A_j, reduced into (-R/2, R/2]
      for (u64 x = 0; x < m; x++) {
        bn_mul(&B[x], &u, GET(&A_work, x, k));
        bn_mul(&t2, &v, GET(&A_work, x, j));
        bn_add(&B[x], &B[x], &t2);
        bigmatrix_hnf_centered_mod(&B[x], &B[x], &R, &half);
      }

      // A_j = (a_{i,k}/d)*A_j - (a_{i,j}/d)*A_k, reduced into (-R/2, R/2]
      bn_div_exact(&t2, GET(&A_work, i, k), &d);
      bn_div_exact(&q, GET(&A_work, i, j), &d);
      for (u64 x = 0; x < m; x++) {
        bn_mul(&t3, &t2, GET(&A_work, x, j));
        bn_mul(&t1, &q, GET(&A_work, x, k));
        bn_sub(&t1, &t3, &t1);
        bigmatrix_hnf_centered_mod(GET(&A_work, x, j), &t1, &R, &half);
        bn_copy(GET(&A_work, x, k), &B[x]);
      }

      // a_{i,k} = u*a_{i,k} + v*a_{i,j} = d exactly
      bn_copy(GET(&A_work, i, k), &d);
    }

    // 4. [Next row]
    // u*a_{i,k} + v*R = d = gcd(a_{i,k}, R)
    bn_gcd_extended_lehmer(&u, &v, &d, GET(&A_work, i, k), &R);
    hnf_bezout_minimal(&u, &v, &d, GET(&A_work, i, k), &R);

    // W_i = u*A_k mod R, taken in [0, R)
    for (u64 x = 0; x < m; x++) {
      bn_mul(&t1, &u, GET(&A_work, x, k));
      bn_mod_pos(GET(&W_local, x, i), &t1, &R);
    }

    // if d = R (i.e. R | a_{i,k}) the diagonal entry would be 0; set it to d
    if (bn_cmp(&d, &R) == 0) {
      bn_copy(GET(&W_local, i, i), &d);
    }

    // final reductions: W_j -= floor(W_{i,j}/W_{i,i}) * W_i for j > i
    for (u64 j = (u64)i + 1; j < m; j++) {
      if (bn_is_zero(GET(&W_local, i, j))) {
        continue;
      }
      bn_div_euclid(&q, GET(&W_local, i, j), GET(&W_local, i, i));
      for (u64 x = 0; x <= (u64)i; x++) {
        bn_mul(&t1, &q, GET(&W_local, x, i));
        bn_sub(GET(&W_local, x, j), GET(&W_local, x, j), &t1);
      }
    }

    R.is_neg = false;
    bn_div_exact(&R, &R, &d);
    bn_rshift(&half, &R, 1);

    if (i > 0) {
      // working modulo R, a_{i-1,k-1} may have reduced to zero; replace it
      // by any nonzero multiple of R
      if (bn_is_zero(GET(&A_work, i - 1, k - 1))) {
        bn_copy(GET(&A_work, i - 1, k - 1), &R);
      }
    }
  }

  bigmatrix_swap(W, &W_local);
  bigmatrix_free(&W_local);

  bn_free_multi(&R, &u, &v, &d, &q, &t1, &t2, &t3, &half, NULL);
  for (u64 x = 0; x < m; x++) {
    bn_free(&B[x]);
  }
  free(B);
  bigmatrix_free(&A_work);
}

void bigmatrix_hermite_gcd(bigmatrix* W, const bigmatrix* A)
{
  // 1. [Initialize]
  u64 m = A->r_size;
  u64 n = A->c_size;
  if (m == 0 || n == 0) return;

  u64 i = m - 1;
  u64 j = n - 1;
  u64 k = n - 1;
  u64 l = 0;
  if (m > n) {
    l = m - n;
  }

  bigmatrix A_work;
  bigmatrix_init(&A_work, m, n);
  bigmatrix_copy(&A_work, A);

  bignum* B;
  B = malloc(m * sizeof(bignum));
  for (u64 x = 0; x < m; x++) {
    bn_init(&B[x]);
  }

  bignum u, v, d, temp, temp2, temp3, one, b, q_k, q_j;
  bn_init_multi(&u, &v, &d, &temp, &temp2, &temp3, &one, &b, &q_k, &q_j, NULL);
  bn_set_u64(&one, 1);

// 2. [Check zero]
step2:
  if (j == 0) {
    goto step4;
  }
  j--;
  if (bn_is_zero(GET(&A_work, i, j))) {
    goto step2;
  }

  // 3. [Euclidean step]
  bn_gcd_extended_lehmer(&u, &v, &d, GET(&A_work, i, k), GET(&A_work, i, j));
  hnf_bezout_minimal(&u, &v, &d, GET(&A_work, i, k), GET(&A_work, i, j));

  for (u64 x = 0; x < m; x++) {
    bn_mul(&B[x], &u, GET(&A_work, x, k));
    bn_mul(&temp, &v, GET(&A_work, x, j));
    bn_add(&B[x], &B[x], &temp);
  }

  bn_div_euclid(&q_k, GET(&A_work, i, k), &d);
  bn_div_euclid(&q_j, GET(&A_work, i, j), &d);

  for (u64 x = 0; x < m; x++) {
    bn_mul(&temp, &q_k, GET(&A_work, x, j));
    bn_mul(&temp2, &q_j, GET(&A_work, x, k));
    bn_sub(GET(&A_work, x, j), &temp, &temp2);
    bn_copy(GET(&A_work, x, k), &B[x]);
  }
  goto step2;

  // 4. [Final reduction]
step4:
  bn_copy(&b, GET(&A_work, i, k));
  if (b.is_neg && !bn_is_zero(&b)) {
    for (u64 x = 0; x < m; x++) {
      bn_neg(GET(&A_work, x, k), GET(&A_work, x, k));
    }
    bn_neg(&b, &b);
  }

  if (bn_is_zero(&b)) {
    k++;
    goto step5;
  }

  for (u64 j_0 = k + 1; j_0 < n; j_0++) {
    bn_div_euclid(&temp2, GET(&A_work, i, j_0), &b);
    for (u64 x = 0; x < m; x++) {
      bn_mul(&temp, &temp2, GET(&A_work, x, k));
      bn_sub(GET(&A_work, x, j_0), GET(&A_work, x, j_0), &temp);
    }
  }

step5:
  // 5. [Finished]
  if (i == l) {
    bigmatrix_free(W);
    bigmatrix_init(W, m, n - k);
    for (u64 x = 0; x < n - k; x++) {
      for (u64 y = 0; y < m; y++) {
        bn_copy(GET(W, y, x), GET(&A_work, y, x + k));
      }
    }
    goto cleanup;
  }

  i--;
  k--;
  j = k;
  goto step2;

cleanup:

  bigmatrix_free(&A_work);
  bn_free_multi(&u, &v, &d, &temp, &temp2, &temp3, &one, &b, &q_k, &q_j, NULL);
  for (u64 x = 0; x < m; x++) {
    bn_free(&B[x]);
  }
}

void bigmatrix_smith(bigmatrix* S, const bigmatrix* A)
{
  if (A->c_size != A->r_size) {
    printf("Matrix must be square\n");
    return;
  }

  // 1. [Initialize i]
  i64 n = A->r_size;
  i64 i = n;
  bignum R;
  bn_init(&R);
  bigmatrix_det(&R, A);
  R.is_neg = false;

  if (n == 1) {
    bigmatrix_set(S, &R, 0, 0);
    bn_free(&R);
    return;
  }

  bignum u, v, d, q_i, q_j, t1, t2, b;
  bn_init_multi(&u, &v, &d, &q_i, &q_j, &t1, &t2, &b, NULL);

  bigmatrix A_work;
  bigmatrix_init(&A_work, n, n);
  bigmatrix_copy(&A_work, A);

  bignum* B = malloc(n * sizeof(bignum));
  for (u64 x = 0; x < n; x++) {
    bn_init(&B[x]);
  }

  bignum* a_i_i = NULL;
  bignum* a_i_j = NULL;

// 2. [Initialize j for row reduction]
step2:
  i64 j = i;
  i64 c = 0;

// 3. [Check zero]
step3:
  if (j == 1) {
    goto step5;
  }

  j--;
  if (bn_is_zero(GET(&A_work, i - 1, j - 1))) {
    goto step3;
  }

  // 4. [Euclidean step]
  a_i_i = GET(&A_work, i - 1, i - 1);
  a_i_j = GET(&A_work, i - 1, j - 1);

  bn_gcd_extended_lehmer(&u, &v, &d, a_i_i, a_i_j);

  // use Remark to find u,v minimal
  bignum a_abs;
  bn_init(&a_abs);
  bn_copy(&a_abs, GET(&A_work, i - 1, j - 1));
  a_abs.is_neg = false;

  bignum a_ii_abs;
  bn_init(&a_ii_abs);
  bn_copy(&a_ii_abs, GET(&A_work, i - 1, i - 1));
  a_ii_abs.is_neg = false;

  // if d == |a_{i,i}|, force u = sign(a_{i,i}), v = 0
  if (bn_cmp(&d, &a_ii_abs) == 0) {
    bn_set_u64(&u, 1);
    bn_set_u64(&v, 0);
    if (GET(&A_work, i - 1, i - 1)->is_neg) bn_set_i64(&u, -1);
  }
  // if d == |a_{i,j}|, force u = 0, v = sign(a_{i,j})
  else if (bn_cmp(&d, &a_abs) == 0) {
    bn_set_u64(&u, 0);
    bn_set_u64(&v, 1);
    if (GET(&A_work, i - 1, j - 1)->is_neg) bn_set_i64(&v, -1);
  }

  bn_free(&a_abs);
  bn_free(&a_ii_abs);

  // B = uA_i + vA_j
  for (u64 x = 0; x < n; x++) {
    bn_mul(&B[x], &u, GET(&A_work, x, i - 1));
    bn_mul(&t1, &v, GET(&A_work, x, j - 1));
    bn_add(&B[x], &B[x], &t1);
  }

  // A_j = ((a_ii / d)A_j - (a_ij / d)A_i) mod R
  bn_div(&q_i, GET(&A_work, i - 1, i - 1), &d);
  bn_div(&q_j, GET(&A_work, i - 1, j - 1), &d);
  for (u64 x = 0; x < n; x++) {
    bn_mul(&t1, &q_i, GET(&A_work, x, j - 1));
    bn_mul(&t2, &q_j, GET(&A_work, x, i - 1));
    bn_sub(&t1, &t1, &t2);
    bn_mod(&t1, &t1, &R);
    bn_copy(GET(&A_work, x, j - 1), &t1);
  }

  // A_i = B mod R
  for (u64 x = 0; x < n; x++) {
    bn_mod(GET(&A_work, x, i - 1), &B[x], &R);
  }

  goto step3;

// 5. [Initialize j for column reduction]
step5:
  j = i;

// 6. [Check zero]
step6:
  if (j == 1) {
    goto step8;
  }
  j--;
  if (bn_is_zero(GET(&A_work, j - 1, i - 1))) {
    goto step6;
  }

  // 7. [Euclidean step]
  a_i_i = GET(&A_work, i - 1, i - 1);
  a_i_j = GET(&A_work, j - 1, i - 1);

  bn_gcd_extended_lehmer(&u, &v, &d, a_i_i, a_i_j);

  // use Remark to find u,v minimal
  // use Remark to find u,v minimal safely
  bignum a_abs_col;
  bn_init(&a_abs_col);
  bn_copy(&a_abs_col, GET(&A_work, j - 1, i - 1));
  a_abs_col.is_neg = false;

  bignum a_ii_abs_col;
  bn_init(&a_ii_abs_col);
  bn_copy(&a_ii_abs_col, GET(&A_work, i - 1, i - 1));
  a_ii_abs_col.is_neg = false;

  // if d == |a_{i,i}|, force u = sign(a_{i,i}), v = 0
  if (bn_cmp(&d, &a_ii_abs_col) == 0) {
    bn_set_u64(&u, 1);
    bn_set_u64(&v, 0);
    if (GET(&A_work, i - 1, i - 1)->is_neg) bn_set_i64(&u, -1);
  }
  // if d == |a_{j,i}|, force u = 0, v = sign(a_{j,i})
  else if (bn_cmp(&d, &a_abs_col) == 0) {
    bn_set_u64(&u, 0);
    bn_set_u64(&v, 1);
    if (GET(&A_work, j - 1, i - 1)->is_neg) bn_set_i64(&v, -1);
  }

  bn_free(&a_abs_col);
  bn_free(&a_ii_abs_col);

  // B = uA'_i + vA'_j
  for (u64 x = 0; x < n; x++) {
    bn_mul(&B[x], &u, GET(&A_work, i - 1, x));
    bn_mul(&t1, &v, GET(&A_work, j - 1, x));
    bn_add(&B[x], &B[x], &t1);
  }

  // A'_j = ((a_ii / d)A'_j - (a_ij / d)A'_i) mod R
  bn_div(&q_i, GET(&A_work, i - 1, i - 1), &d);
  bn_div(&q_j, GET(&A_work, j - 1, i - 1), &d);
  for (u64 x = 0; x < n; x++) {
    bn_mul(&t1, &q_i, GET(&A_work, j - 1, x));
    bn_mul(&t2, &q_j, GET(&A_work, i - 1, x));
    bn_sub(&t1, &t1, &t2);
    bn_mod(&t1, &t1, &R);
    bn_copy(GET(&A_work, j - 1, x), &t1);
  }

  // A'_i = B mod R
  for (u64 x = 0; x < n; x++) {
    bn_mod(GET(&A_work, i - 1, x), &B[x], &R);
  }

  c++;
  goto step6;

// 8. [Repeat stage i?]
step8:
  if (c > 0) {
    goto step2;
  }

  // 9. [Check the rest of the matrix]
  bn_copy(&b, GET(&A_work, i - 1, i - 1));
  for (i64 k_0 = 0; k_0 < i - 1; k_0++) {
    for (i64 l_0 = 0; l_0 < i - 1; l_0++) {
      bool not_divisible = false;
      if (bn_is_zero(&b)) {
        if (!bn_is_zero(GET(&A_work, k_0, l_0))) {
          not_divisible = true;
        }
      } else {
        bn_mod(&t1, GET(&A_work, k_0, l_0), &b);
        if (!bn_is_zero(&t1)) {
          not_divisible = true;
        }
      }

      if (not_divisible) {
        for (i64 x = 0; x < n; x++) {
          bn_add(GET(&A_work, i - 1, x), GET(&A_work, i - 1, x),
                 GET(&A_work, k_0, x));
        }
        goto step2;
      }
    }
  }

  // 10. [Next stage]
  bn_gcd_lehmer(&t1, GET(&A_work, i - 1, i - 1), &R);
  bigmatrix_set(&A_work, &t1, i - 1, i - 1);
  bn_div(&R, &R, &t1);

  if (i == 2) {
    bn_gcd_lehmer(&t2, GET(&A_work, 0, 0), &R);
    bigmatrix_set(&A_work, &t2, 0, 0);

    // order diagonal entries in proper order
    bignum zero;
    bn_init(&zero);
    for (i64 r = 0; r < n; r++) {
      for (i64 c = 0; c < n; c++) {
        if (r == c) {
          bigmatrix_set(S, GET(&A_work, n - 1 - r, n - 1 - r), r, r);
        } else {
          bigmatrix_set(S, &zero, r, c);
        }
      }
    }
    bn_free(&zero);
    goto cleanup;
  }

  i--;
  goto step2;

cleanup:
  for (i64 x = 0; x < n; x++) {
    bn_free(&B[x]);
  }
  free(B);
  bn_free(&R);
  bn_free_multi(&u, &v, &d, &q_i, &q_j, &t1, &t2, &b, NULL);
  bigmatrix_free(&A_work);
}

bool bigmatrix_equal(const bigmatrix* A, const bigmatrix* B)
{
  if (A->r_size != B->r_size || A->c_size != B->c_size) return false;

  u64 total = A->r_size * A->c_size;
  for (u64 i = 0; i < total; i++) {
    const bignum* a = &A->data[i];
    const bignum* b = &B->data[i];

    // bn_cmp does not normalize the two representations of zero
    // (size 0 from bigmatrix_init vs size 1 with a zero limb), so
    // settle the zero cases explicitly
    bool za = bn_is_zero(a);
    bool zb = bn_is_zero(b);
    if (za || zb) {
      if (za != zb) return false;
      continue;
    }
    if (bn_cmp(a, b) != 0) return false;
  }
  return true;
}

bool bigmatrix_hnf_check_structure(const bigmatrix* H)
{
  u64 r = H->r_size;
  u64 c = H->c_size;

  // once a zero row is hit, every later row must be zero as well
  bool zero_row_seen = false;

  for (u64 i = 0; i < r; i++) {
    bignum* row = &H->data[i * c];

    // upper triangular: H[i][j] == 0 for all j < i (capped at c for
    // rows below the last column)
    u64 left = MIN(i, c);
    for (u64 j = 0; j < left; j++) {
      if (!bn_is_zero(&row[j])) {
        return false;
      }
    }

    // zero rows only in a trailing block: from here on everything must
    // be zero
    if (zero_row_seen) {
      for (u64 j = 0; j < c; j++) {
        if (!bn_is_zero(&row[j])) {
          return false;
        }
      }
      continue;
    }

    // a nonzero row needs a diagonal pivot, so rows at or beyond the
    // last column can only be zero rows
    if (i >= c) {
      for (u64 j = 0; j < c; j++) {
        if (!bn_is_zero(&row[j])) {
          return false;
        }
      }
      zero_row_seen = true;
      continue;
    }

    bignum* pivot = &row[i];

    if (bn_is_zero(pivot)) {
      // a zero pivot may only start the trailing zero block: everything
      // right of the pivot must also be zero
      for (u64 j = i + 1; j < c; j++) {
        if (!bn_is_zero(&row[j])) {
          return false;
        }
      }
      zero_row_seen = true;
      continue;
    }

    // positive pivot
    if (pivot->is_neg) {
      return false;
    }

    // reduced entries: the entries right of the pivot, H[i][j] for
    // j > i, must lie in [0, pivot) — this is the reduction invariant
    // established by the column-operation HNF algorithms in this
    // library (each row is reduced modulo its own pivot)
    for (u64 j = i + 1; j < c; j++) {
      if (row[j].is_neg || bn_cmp(&row[j], pivot) >= 0) {
        return false;
      }
    }
  }

  return true;
}

bool bigmatrix_hnf_check_transformation(const bigmatrix* A, const bigmatrix* H,
                                        const bigmatrix* U)
{
  // column-operation convention: H = A * U; bigmatrix_mul() silently
  // no-ops unless A->c_size == U->r_size, so reject shape mismatches
  // before allocating the product
  if (A->r_size != H->r_size || A->c_size != U->r_size ||
      U->c_size != H->c_size) {
    return false;
  }

  bigmatrix P;
  bigmatrix_init(&P, H->r_size, H->c_size);

  bigmatrix_mul(&P, A, U);

  bool ok = bigmatrix_equal(&P, H);

  bigmatrix_free(&P);
  return ok;
}

bool bigmatrix_hnf_check_unimodular(const bigmatrix* U)
{
  if (!bigmatrix_is_square(U)) {
    return false;
  }
  // the empty matrix acts as the identity transformation
  if (U->r_size == 0) {
    return true;
  }

  bignum det;
  bn_init(&det);

  bigmatrix_det(&det, U);

  // |det| == 1: clear the sign flag and use the O(1) one-check
  det.is_neg = false;
  bool ok = bn_is_one(&det);

  bn_free(&det);
  return ok;
}

bool bigmatrix_hnf_verify(const bigmatrix* A, const bigmatrix* H,
                          const bigmatrix* U)
{
  // order matters: the cheap structure check runs first, the
  // multiplication and the determinant only if it passes
  return bigmatrix_hnf_check_structure(H) &&
         bigmatrix_hnf_check_transformation(A, H, U) &&
         bigmatrix_hnf_check_unimodular(U);
}

// void bigmatrix_LLL(bigmatrix* B, u64 n, double delta, bigmatrix* H)
// {
//   u64 k = 2;
//   u64 k_max = 1;

//   // Allocate Gram-Schmidt orthogonalized vectors and squared norms
//   bigvector* b_star = malloc(n * sizeof(bigvector));
//   bignum* B_vals = malloc(n * sizeof(bignum));
//   for (u64 i = 0; i < n; i++) {
//     bigvector_init(&b_star[i], B->rows);
//     bn_init(&B_vals[i]);
//   }

//   // Allocate table for Gram-Schmidt coefficients mu[k][j] stored as
//   doubles double* mu = calloc(n * n, sizeof(double));

//   // Step 1 [Initialize]
//   bigvector b_k;
//   bigvector_init(&b_k, B->rows);
//   bigmatrix_get_col(&b_star[0], B, 0);
//   bigvector_copy(&b_k, &b_star[0]);
//   bigvector_dot(&B_vals[0], &b_k, &b_k);

//   bigmatrix_id(H, n);

//   bigvector temp_vec;
//   bigvector_init(&temp_vec, B->rows);

//   // Helper lambda or inline logic for RED(k, l)
//   // RED(k, l): makes |mu[k][l]| <= 0.5 by reducing b_k using b_l
//   #auto_red = ^(u64 ki, u64 li) {
//     double mukl = mu[ki * n + li];
//     if (fabs(mukl) > 0.5) {
//       long long x = (long long)round(mukl);
//       // b_k = b_k - x * b_l
//       bigvector_sub_mul(B, ki, li, x); // Assumes a function to subtract
//       scaled column
//       // H_k = H_k - x * H_l
//       bigmatrix_sub_mul(H, ki, li, x);
//       mu[ki * n + li] -= x;
//       for (u64 j = 0; j < li; j++) {
//         mu[ki * n + j] -= x * mu[li * n + j];
//       }
//     }
//   };

//   while (k <= n) {
//     // Step 2 [Incremental Gram-Schmidt]
//     if (k <= k_max) {
//       // Proceed to Step 3
//     } else {
//       k_max = k;
//       bigmatrix_get_col(&b_k, B, k - 1);
//       bigvector_copy(&b_star[k - 1], &b_k);

//       for (u64 j = 1; j < k; j++) {
//         bignum dot_val;
//         bn_init(&dot_val);
//         bigvector_dot(&dot_val, &b_k, &b_star[j - 1]);

//         double dot_d = bn_to_double(&dot_val);
//         double Bj_d = bn_to_double(&B_vals[j - 1]);
//         mu[(k - 1) * n + (j - 1)] = dot_d / Bj_d;

//         // b*_k = b*_k - mu_{k,j} * b*_j (approximate or exact vector
//         subtraction)
//         // Implementation depends on scalar-vector vector subtraction
//         helpers
//         // ...
//       }
//       bigvector_dot(&B_vals[k - 1], &b_star[k - 1], &b_star[k - 1]);

//       if (bn_is_zero(&B_vals[k - 1])) {
//         // Error: vectors did not form a basis (linearly dependent)
//         return;
//       }
//     }

//     // Step 3 [Test LLL condition]
//     // Execute RED(k, k - 1)
//     // (Inline or helper call for RED calculation)
//     double muk_k1 = mu[(k - 1) * n + (k - 2)];
//     if (fabs(muk_k1) > 0.5) {
//       long long x = (long long)round(muk_k1);
//       // Apply reduction to basis column k and transform matrix H
//       mu[(k - 1) * n + (k - 2)] -= x;
//       for (u64 j = 0; j < k - 2; j++) {
//         mu[(k - 1) * n + j] -= x * mu[(k - 2) * n + j];
//       }
//     }

//     bignum Bk, Bk_prev, threshold_term;
//     bn_init(&Bk); bn_init(&Bk_prev); bn_init(&threshold_term);

//     // Check LLL inequality: B_k < (delta - mu_{k,k-1}^2) * B_{k-1}
//     // Using floating point conversion for norm comparison or exact bignum
//     arithmetic double Bk_d = bn_to_double(&B_vals[k - 1]); double Bk_prev_d
//     = bn_to_double(&B_vals[k - 2]); double rhs = (delta - muk_k1 * muk_k1)
//     * Bk_prev_d;

//     if (Bk_d < rhs) {
//       // SWAP(k)
//       // Swap columns k and k-1 in B and H, update Gram-Schmidt data
//       accordingly
//       // ...
//       k = (k > 2) ? k - 1 : 2;
//     } else {
//       for (long long l = (long long)k - 3; l >= 0; l--) {
//         // Execute RED(k, l)
//         double mukl = mu[(k - 1) * n + l];
//         if (fabs(mukl) > 0.5) {
//           long long x = (long long)round(mukl);
//           mu[(k - 1) * n + l] -= x;
//           for (u64 j = 0; j < (u64)l; j++) {
//             mu[(k - 1) * n + j] -= x * mu[l * n + j];
//           }
//         }
//       }
//       k++;
//     }

//     // Step 4 [Finished?] handled by the while (k <= n) loop condition.
//   }

//   // Cleanup allocations
//   for (u64 i = 0; i < n; i++) {
//     bigvector_free(&b_star[i]);
//     bn_free(&B_vals[i]);
//   }
//   free(b_star);
//   free(B_vals);
//   free(mu);
//   bigvector_free(&b_k);
//   bigvector_free(&temp_vec);
// }
