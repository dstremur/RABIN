/*
 * bigmatrix.c
 *
 * Matrix arithmetic over bignums.
 *
 * This file implements matrices with bignum entries: initialization,
 * freeing, copying, element/row/column access, printing,
 * component-wise addition, schoolbook multiplication, matrix-vector
 * products, the Hadamard bound on the determinant, and the
 * determinant (computed via the RNS path; a direct Bareiss
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

#include "../../include/bigmath.h"
#include "../../include/bigpoly.h"
#include "../../include/bigrns.h"
#include "../../include/primes.h"

void bigmatrix_init(bigmatrix* M, u64 r, u64 c)
{
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

void bigmatrix_print(const bigmatrix* A)
{
  bignum temp;
  bn_init(&temp);
  for (u64 i = 0; i < A->r_size; i++) {
    for (u64 j = 0; j < A->c_size; j++) {
      bigmatrix_get(&temp, A, i, j);
      bn_print(&temp);
      printf(" ");
    }
    printf("\n");
  }

  bn_free(&temp);
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

void bigmatrix_hermite_mod_d(bigmatrix* W, const bigmatrix* A, const bignum* D)
{
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

  bignum a_abs;
  bn_init(&a_abs);
  bn_copy(&a_abs, GET(&A_work, i, k));
  a_abs.is_neg = false;

  // if d == |a_{i,k}|, force v == 0, u = sign(a_{i,k})
  if (bn_cmp(&d, &a_abs) == 0) {
    bn_set_u64(&v, 0);
    bn_set_u64(&u, 1);
    if (GET(&A_work, i, k)->is_neg) {
      bn_set_i64(&u, -1);
    }
  }
  bn_free(&a_abs);

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
step4:
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
step7:
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
step9:
  bn_copy(&b, GET(&A_work, i - 1, i - 1));
  for (u64 k_0 = 0; k_0 < (u64)(i - 1); k_0++) {
    for (u64 l_0 = 0; l_0 < (u64)(i - 1); l_0++) {
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
        for (u64 x = 0; x < n; x++) {
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
    for (u64 r = 0; r < n; r++) {
      for (u64 c = 0; c < n; c++) {
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
  for (u64 x = 0; x < n; x++) {
    bn_free(&B[x]);
  }
  free(B);
  bn_free(&R);
  bn_free_multi(&u, &v, &d, &q_i, &q_j, &t1, &t2, &b, NULL);
  bigmatrix_free(&A_work);
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
