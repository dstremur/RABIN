/*
 * rmat.c
 *
 * Matrix arithmetic over bignums.
 *
 * This file implements matrices with rz_t entries: initialization,
 * freeing, copying, element/row/column access, structural check
 * predicates (square, zero, identity, diagonal, triangular, symmetric),
 * printing, component-wise addition, schoolbook multiplication,
 * matrix-vector products, the Hadamard bound on the determinant, and
 * the determinant (computed via the RNS path; a direct Bareiss
 * implementation is kept below it, currently disabled).
 *
 * A rmat_t is a row-major dynamic array of bignums with row and
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

#include "../../include/rmat.h"

#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <string.h>

#include "../../include/rmatq.h"
#include "../../include/rns.h"
#include "../../include/rns_primes.h"
#include "../../include/rpol.h"
#include "../../include/rq.h"
#include "../../include/rzmath.h"

rabin_err_t rmat_init(rmat_t* M, u64 r, u64 c)
{
  if (M == NULL) return RABIN_ERR_NULL_PTR;

  // overflow-safe row * cols
  if (r != 0 && c > UINT64_MAX / r) return RABIN_ERR_OVERFLOW;

  // idempotent init: release any previous allocation
  if (M->data != NULL) {
    rmat_clear(M);
  }

  if (r * c == 0) {
    M->data = NULL;
    M->rows = r;
    M->cols = c;
    return RABIN_SUCCESS;
  }

  M->data = malloc(r * c * sizeof(rz_t));
  if (M->data == NULL) return RABIN_ERR_OUT_OF_MEMORY;

  M->rows = r;
  M->cols = c;

  for (u64 i = 0; i < r * c; i++) {
    rz_init(&M->data[i]);
  }
  return RABIN_SUCCESS;
}

rabin_err_t rmat_clear(rmat_t* A)
{
  if (A == NULL) return RABIN_ERR_NULL_PTR;
  if (!A->data) return RABIN_SUCCESS;

  for (u64 i = 0; i < A->rows * A->cols; i++) {
    rz_clear(&A->data[i]);
  }

  free(A->data);

  A->data = NULL;
  A->rows = 0;
  A->cols = 0;
  return RABIN_SUCCESS;
}

rabin_err_t rmat_print_python(const rmat_t* A)
{
  if (A == NULL) return RABIN_ERR_NULL_PTR;

  printf("[");  // Start outer list
  for (u64 i = 0; i < A->rows; i++) {
    printf("[");  // Start row
    for (u64 j = 0; j < A->cols; j++) {
      rz_print(RMAT_GET(A, i, j));
      if (j < A->cols - 1) printf(", ");
    }
    printf("]");  // End row
    if (i < A->rows - 1) printf(",\n ");
  }
  printf("]\n");  // End outer list
  return RABIN_SUCCESS;
}

rabin_err_t rmat_copy(rmat_t* R, const rmat_t* A)
{
  if (R == NULL || A == NULL) return RABIN_ERR_NULL_PTR;

  // check if sizes match
  if (A->cols != R->cols || A->rows != R->rows) return RABIN_ERR_MATRIX_DIM;

  for (u64 i = 0; i < A->cols * A->rows; i++) {
    rabin_err_t err = rz_copy(&R->data[i], &A->data[i]);
    if (err != RABIN_SUCCESS) return err;
  }
  return RABIN_SUCCESS;
}

rabin_err_t rmat_get(rz_t* R, const rmat_t* A, u64 r, u64 c)
{
  if (R == NULL || A == NULL) return RABIN_ERR_NULL_PTR;
  if (r >= A->rows || c >= A->cols) return RABIN_ERR_INVALID_ARG;

  return rz_copy(R, RMAT_GET(A, r, c));
}

rabin_err_t rmat_set(rmat_t* A, const rz_t* a, u64 r, u64 c)
{
  if (A == NULL || a == NULL) return RABIN_ERR_NULL_PTR;
  if (r >= A->rows || c >= A->cols) return RABIN_ERR_INVALID_ARG;

  return rz_copy(RMAT_GET(A, r, c), a);
}

rabin_err_t rmat_get_col(rvec_t* c, const rmat_t* A, u64 col)
{
  if (c == NULL || A == NULL) return RABIN_ERR_NULL_PTR;
  if (col >= A->cols) return RABIN_ERR_INVALID_ARG;
  if (c->size != A->rows) return RABIN_ERR_MATRIX_DIM;

  for (u64 i = 0; i < A->rows; i++) {
    rabin_err_t err = rvec_set(c, RMAT_GET(A, i, col), i);
    if (err != RABIN_SUCCESS) return err;
  }
  return RABIN_SUCCESS;
}

rabin_err_t rmat_get_row(rvec_t* r, const rmat_t* A, u64 row)
{
  if (r == NULL || A == NULL) return RABIN_ERR_NULL_PTR;
  if (row >= A->rows) return RABIN_ERR_INVALID_ARG;
  if (r->size != A->cols) return RABIN_ERR_MATRIX_DIM;

  for (u64 i = 0; i < A->cols; i++) {
    rabin_err_t err = rvec_set(r, RMAT_GET(A, row, i), i);
    if (err != RABIN_SUCCESS) return err;
  }
  return RABIN_SUCCESS;
}

bool rmat_is_square(const rmat_t* A)
{
  if (A == NULL) return false;
  return A->rows == A->cols;
}

bool rmat_is_zero(const rmat_t* A)
{
  if (A == NULL) return false;
  // single sequential pass over the flat row-major storage; rz_is_zero is
  // an O(1) size/limb inspection, so no temporary bignums are needed
  u64 total = A->rows * A->cols;
  for (u64 i = 0; i < total; i++) {
    if (!rz_is_zero(&A->data[i])) {
      return false;
    }
  }
  return true;
}

bool rmat_is_identity(const rmat_t* A)
{
  if (A == NULL) return false;
  // dimension short-circuit: a non-square matrix cannot be the identity,
  // bail out before touching any entry
  if (!rmat_is_square(A)) {
    return false;
  }

  u64 n = A->rows;
  for (u64 i = 0; i < n; i++) {
    rz_t* row = &A->data[i * A->cols];
    for (u64 j = 0; j < n; j++) {
      if (i == j) {
        // diagonal entry must be exactly 1 (rz_is_one rejects -1)
        if (!rz_is_one(&row[j])) {
          return false;
        }
      } else {
        // every off-diagonal entry must be 0
        if (!rz_is_zero(&row[j])) {
          return false;
        }
      }
    }
  }
  return true;
}

bool rmat_is_diagonal(const rmat_t* A)
{
  if (A == NULL) return false;
  // only the off-diagonal sub-regions can violate: in row i those are the
  // columns left of the diagonal, j < i (capped at cols), and the
  // columns right of it, j > i (empty once i >= cols). The diagonal
  // entry (i, i) itself is never inspected
  for (u64 i = 0; i < A->rows; i++) {
    rz_t* row = &A->data[i * A->cols];

    // left of the diagonal: j in [0, min(i, cols))
    u64 left = RZ_MIN(i, A->cols);
    for (u64 j = 0; j < left; j++) {
      if (!rz_is_zero(&row[j])) {
        return false;
      }
    }

    // right of the diagonal: j in (i, cols)
    for (u64 j = i + 1; j < A->cols; j++) {
      if (!rz_is_zero(&row[j])) {
        return false;
      }
    }
  }
  return true;
}

bool rmat_is_upper_triangular(const rmat_t* A)
{
  if (A == NULL) return false;
  // only the strict lower triangle (row > col) can violate; row 0 has no
  // entries below the diagonal, so the scan starts at i = 1
  for (u64 i = 1; i < A->rows; i++) {
    rz_t* row = &A->data[i * A->cols];

    // below the diagonal in row i: j in [0, min(i-1, cols-1)], i.e.
    // j < min(i, cols) — the cap handles tall (c < r) matrices
    u64 limit = RZ_MIN(i, A->cols);
    for (u64 j = 0; j < limit; j++) {
      if (!rz_is_zero(&row[j])) {
        return false;
      }
    }
  }
  return true;
}

bool rmat_is_lower_triangular(const rmat_t* A)
{
  if (A == NULL) return false;
  // only the strict upper triangle (col > row) can violate; in row i the
  // first candidate is column i+1, and the inner loop is empty once
  // i >= cols, so rows below the diagonal region cost nothing
  for (u64 i = 0; i < A->rows; i++) {
    rz_t* row = &A->data[i * A->cols];
    for (u64 j = i + 1; j < A->cols; j++) {
      if (!rz_is_zero(&row[j])) {
        return false;
      }
    }
  }
  return true;
}

bool rmat_is_symmetric(const rmat_t* A)
{
  if (A == NULL) return false;
  // dimension short-circuit: a non-square matrix cannot be symmetric,
  // bail out before touching any entry
  if (!rmat_is_square(A)) {
    return false;
  }

  u64 n = A->rows;
  // walk only the strict upper triangle: each mirror pair (i, j)/(j, i)
  // with i < j is compared exactly once; row i is swept sequentially,
  // and RMAT_GET(A, j, i) is a constant offset into row j's contiguous block
  for (u64 i = 0; i < n; i++) {
    rz_t* row_i = &A->data[i * A->cols];
    for (u64 j = i + 1; j < n; j++) {
      if (rz_cmp(&row_i[j], RMAT_GET(A, j, i)) != 0) {
        return false;
      }
    }
  }
  return true;
}

rabin_err_t rmat_scalar(rmat_t* R, const rmat_t* A, const rz_t* a)
{
  if (R == NULL || A == NULL || a == NULL) return RABIN_ERR_NULL_PTR;
  if (R->rows != A->rows || R->cols != A->cols) return RABIN_ERR_MATRIX_DIM;

  u64 total = A->rows * A->cols;

  for (u64 x = 0; x < total; x++) {
    rabin_err_t err = rz_mul(&R->data[x], &A->data[x], a);
    if (err != RABIN_SUCCESS) return err;
  }
  return RABIN_SUCCESS;
}

rabin_err_t rmat_div_exact_scalar(rmat_t* R, const rmat_t* A, const rz_t* a)
{
  if (R == NULL || A == NULL || a == NULL) return RABIN_ERR_NULL_PTR;
  if (R->rows != A->rows || R->cols != A->cols) return RABIN_ERR_MATRIX_DIM;

  u64 total = A->rows * A->cols;

  for (u64 x = 0; x < total; x++) {
    rabin_err_t err = rz_div_exact(&R->data[x], &A->data[x], a);
    if (err != RABIN_SUCCESS) return err;
  }
  return RABIN_SUCCESS;
}

rabin_err_t rmat_mv(rvec_t* r, const rmat_t* A, rvec_t* v)
{
  if (r == NULL || A == NULL || v == NULL) return RABIN_ERR_NULL_PTR;
  if (r->size != A->rows || v->size != A->cols) return RABIN_ERR_MATRIX_DIM;

  rabin_err_t err = RABIN_SUCCESS;
  rvec_t tmp = {0};
  rz_t dot;
  rz_init(&dot);

  TRY(rvec_init(&tmp, A->cols));

  for (u64 i = 0; i < A->rows; i++) {
    TRY(rmat_get_row(&tmp, A, i));
    TRY(rvec_dot(&dot, v, &tmp));
    TRY(rvec_set(r, &dot, i));
  }

cleanup:
  rz_clear(&dot);
  rvec_clear(&tmp);
  return err;
}

rabin_err_t rmat_vm(rvec_t* r, const rmat_t* A, rvec_t* v)
{
  if (r == NULL || A == NULL || v == NULL) return RABIN_ERR_NULL_PTR;
  if (r->size != A->cols || v->size != A->rows) return RABIN_ERR_MATRIX_DIM;

  rabin_err_t err = RABIN_SUCCESS;
  rvec_t tmp = {0};
  rz_t dot;
  rz_init(&dot);

  TRY(rvec_init(&tmp, A->rows));

  for (u64 i = 0; i < A->cols; i++) {
    TRY(rmat_get_col(&tmp, A, i));
    TRY(rvec_dot(&dot, v, &tmp));
    TRY(rvec_set(r, &dot, i));
  }

cleanup:
  rz_clear(&dot);
  rvec_clear(&tmp);
  return err;
}

rabin_err_t rmat_hadamard(rz_t* r, const rmat_t* A)
{
  if (r == NULL || A == NULL) return RABIN_ERR_NULL_PTR;

  rabin_err_t err = RABIN_SUCCESS;
  rvec_t v = {0};
  rz_t tmp;
  rz_init(&tmp);

  TRY(rz_set_u64(r, 1));
  TRY(rvec_init(&v, A->rows));

  for (u64 i = 0; i < A->cols; i++) {
    TRY(rmat_get_col(&v, A, i));
    TRY(rvec_norm(&tmp, &v));
    // add 1 since norm uses isqrt
    TRY(rz_add_u64(&tmp, &tmp, 1));
    TRY(rz_mul(r, r, &tmp));
  }

cleanup:
  rz_clear(&tmp);
  rvec_clear(&v);
  return err;
}

rabin_err_t rmat_add(rmat_t* R, const rmat_t* A, const rmat_t* B)
{
  if (R == NULL || A == NULL || B == NULL) return RABIN_ERR_NULL_PTR;

  // check if sizes match
  if (A->cols != B->cols || A->rows != B->rows) return RABIN_ERR_MATRIX_DIM;
  if (R->cols != A->cols || R->rows != A->rows) return RABIN_ERR_MATRIX_DIM;

  for (u64 i = 0; i < A->cols * A->rows; i++) {
    rabin_err_t err = rz_add(&R->data[i], &A->data[i], &B->data[i]);
    if (err != RABIN_SUCCESS) return err;
  }

  return RABIN_SUCCESS;
}

rabin_err_t rmat_sub(rmat_t* R, const rmat_t* A, const rmat_t* B)
{
  if (R == NULL || A == NULL || B == NULL) return RABIN_ERR_NULL_PTR;

  // check if sizes match
  if (A->cols != B->cols || A->rows != B->rows) return RABIN_ERR_MATRIX_DIM;
  if (R->cols != A->cols || R->rows != A->rows) return RABIN_ERR_MATRIX_DIM;

  for (u64 i = 0; i < A->cols * A->rows; i++) {
    rabin_err_t err = rz_sub(&R->data[i], &A->data[i], &B->data[i]);
    if (err != RABIN_SUCCESS) return err;
  }

  return RABIN_SUCCESS;
}
// maximum printed width of a single entry in the truncated modes
#define RMAT_PRINT_MAX_WIDTH 16

// upper bound on the printed width (digits plus sign) of a single entry,
// computed without a decimal conversion: ceil(bits * log10(2)) using the
// rational over-approximation 30103/100000 > log10(2), so the result is
// always at least the true digit count
static size_t rmat_entry_width(const rz_t* x)
{
  if (rz_is_zero(x)) return 1;

  size_t bits = (size_t)rz_bit_length(x);
  size_t w = (bits * 30103 + 99999) / 100000;
  if (x->is_neg) w++;
  return w;
}

static rabin_err_t rmat_print_impl(const rmat_t* A, bool full, bool tail)
{
  if (A == NULL) return RABIN_ERR_NULL_PTR;
  if (A->rows == 0 || A->cols == 0) return RABIN_SUCCESS;

  // per-column width: widest entry of the column, capped in the truncated
  // modes. A single outlier (e.g. the last invariant factor of a Smith form
  // sitting in a column of zeros) must not inflate the column, so in the
  // truncated modes outliers are excluded from the column width; they are
  // rendered at the capped width and simply overflow the column to the
  // right
  size_t cap = RMAT_PRINT_MAX_WIDTH;
  size_t* colw = calloc(A->cols, sizeof(size_t));
  if (colw == NULL) return RABIN_ERR_OUT_OF_MEMORY;
  for (u64 i = 0; i < A->rows; i++) {
    for (u64 j = 0; j < A->cols; j++) {
      size_t xw = rmat_entry_width(RMAT_GET(A, i, j));
      if (!full && xw > cap) continue;
      if (xw > colw[j]) colw[j] = xw;
    }
  }

  rabin_err_t err = RABIN_SUCCESS;
  for (u64 i = 0; i < A->rows && err == RABIN_SUCCESS; i++) {
    for (u64 j = 0; j < A->cols && err == RABIN_SUCCESS; j++) {
      char* s = rz_to_string(RMAT_GET(A, i, j));
      if (s == NULL) {
        err = RABIN_ERR_OUT_OF_MEMORY;
        break;
      }
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
  return err;
}

rabin_err_t rmat_print(const rmat_t* A)
{
  return rmat_print_impl(A, false, false);
}

rabin_err_t rmat_print_full(const rmat_t* A)
{
  return rmat_print_impl(A, true, false);
}

rabin_err_t rmat_print_tail(const rmat_t* A)
{
  return rmat_print_impl(A, false, true);
}

rabin_err_t rmat_mul(rmat_t* R, const rmat_t* A, const rmat_t* B)
{
  if (R == NULL || A == NULL || B == NULL) return RABIN_ERR_NULL_PTR;
  if (A->cols != B->rows) return RABIN_ERR_MATRIX_DIM;
  if (R->rows != A->rows || R->cols != B->cols) return RABIN_ERR_MATRIX_DIM;

  // Transpose B
  rabin_err_t err = RABIN_SUCCESS;
  rmat_t B_T = {0};
  TRY(rmat_init(&B_T, B->cols, B->rows));

  for (u64 i = 0; i < B->rows; i++) {
    for (u64 j = 0; j < B->cols; j++) {
      TRY(rz_copy(&B_T.data[j * B_T.cols + i], &B->data[i * B->cols + j]));
    }
  }
#pragma omp parallel
  {
    rz_t sum, tmp;
    rz_init_multi(&sum, &tmp, NULL);

#pragma omp for collapse(2) schedule(dynamic)
    for (u64 i = 0; i < A->rows; i++) {
      for (u64 j = 0; j < B_T.rows; j++) {
        rz_set_i64(&sum, 0);

        for (u64 k = 0; k < A->cols; k++) {
          // tmp = A[i][k] * B_T[j][k]
          rz_mul(&tmp, &A->data[i * A->cols + k], &B_T.data[j * B_T.cols + k]);
          rz_add(&sum, &sum, &tmp);
        }

        rmat_set(R, &sum, i, j);
      }
    }

    rz_clear(&sum);
    rz_clear(&tmp);
  }

cleanup:
  rmat_clear(&B_T);
  return err;
}

rabin_err_t rmat_det(rz_t* d, const rmat_t* A)
{
  if (d == NULL || A == NULL) return RABIN_ERR_NULL_PTR;
  if (A->cols != A->rows) return RABIN_ERR_MATRIX_DIM;

  rabin_err_t err = RABIN_SUCCESS;
  rmat_t T = {0};
  rns_ctx_t ctx = {0};
  TRY(rmat_init(&T, A->cols, A->rows));
  TRY(rmat_copy(&T, A));

  u64 k = rns_estimate_det(&T);

  // not enough primes in the table: fall back to the exact Bareiss route
  const u64 np = sizeof(RNS_PRIMES) / sizeof(RNS_PRIMES[0]);
  if (k == 0 || k > np) {
    rmat_clear(&T);
    return rmat_det_bareiss(d, A);
  }

  TRY(rns_ctx_init(&ctx, RNS_PRIMES, k));
  TRY(rmat_det_rns(d, &T, &ctx));

cleanup:
  rns_ctx_clear(&ctx);
  rmat_clear(&T);
  return err;
}

rabin_err_t rmat_det_bareiss(rz_t* det, const rmat_t* A)
{
  if (det == NULL || A == NULL) return RABIN_ERR_NULL_PTR;
  if (!rmat_is_square(A)) return RABIN_ERR_MATRIX_DIM;

  u64 n = A->rows;

  if (n == 0) {
    return rz_set_u64(det, 1);
  }

  if (n == 1) {
    return rz_copy(det, RMAT_GET(A, 0, 0));
  }

  // Work on a local copy of the matrix to preserve the input
  rabin_err_t err = RABIN_SUCCESS;
  rmat_t T = {0};
  rz_t p0;
  rz_init(&p0);

  TRY(rmat_init(&T, n, n));
  TRY(rmat_copy(&T, A));
  TRY(rz_set_u64(&p0, 1));

  int sign = 1;
  bool det_zero = false;

#pragma omp parallel
  {
    // Each thread gets its own private temporary bignums initialized once
    // to avoid heap allocation churn inside the nested loops.
    rz_t prod1, prod2, diff;
    rz_init_multi(&prod1, &prod2, &diff, NULL);

    for (u64 k = 0; k < n - 1; k++) {
      // 1. Sequential Pivot Selection & Row Swap (executed by one thread)
#pragma omp single
      {
        if (rz_is_zero(RMAT_GET(&T, k, k))) {
          u64 pivot_row = n;
          for (u64 i = k + 1; i < n; i++) {
            if (!rz_is_zero(RMAT_GET(&T, i, k))) {
              pivot_row = i;
              break;
            }
          }

          if (pivot_row == n) {
            det_zero = true;
          } else {
            for (u64 col = 0; col < n; col++) {
              rz_swap(RMAT_GET(&T, k, col), RMAT_GET(&T, pivot_row, col));
            }
            sign = -sign;
          }
        }
      }  // Implicit barrier here ensures row swap and det_zero flag are visible
         // to all threads

      if (det_zero) {
        // Fast-forward through remaining steps if pivot search failed
        continue;
      }

      // 2. Parallel Bareiss Update Step
      // Row k and Column k entries are read-only; each (i, j) cell is modified
      // by exactly one thread.
#pragma omp for collapse(2) schedule(dynamic)
      for (u64 i = k + 1; i < n; i++) {
        for (u64 j = k + 1; j < n; j++) {
          rz_mul(&prod1, RMAT_GET(&T, k, k), RMAT_GET(&T, i, j));
          rz_mul(&prod2, RMAT_GET(&T, i, k), RMAT_GET(&T, k, j));
          rz_sub(&diff, &prod1, &prod2);

          // Exact division guaranteed by Bareiss property
          rz_div_exact(RMAT_GET(&T, i, j), &diff, &p0);
        }
      }

      // 3. Update division factor p0 for step k+1
#pragma omp single
      {
        rz_copy(&p0, RMAT_GET(&T, k, k));
      }  // Implicit barrier ensures p0 is updated before step k+1 starts
    }

    rz_clear_multi(&prod1, &prod2, &diff, NULL);
  }

  // Set final result
  if (det_zero) {
    TRY(rz_set_u64(det, 0));
  } else {
    TRY(rz_copy(det, RMAT_GET(&T, n - 1, n - 1)));
    if (sign < 0) {
      TRY(rz_neg(det, det));
    }
  }

cleanup:
  rz_clear(&p0);
  rmat_clear(&T);
  return err;
}

static rabin_err_t rmat_neg(rmat_t* A)
{
  if (A == NULL) return RABIN_ERR_NULL_PTR;

  for (u64 i = 0; i < A->rows; i++) {
    for (u64 j = 0; j < A->cols; j++) {
      rabin_err_t err = rz_neg(RMAT_GET(A, i, j), RMAT_GET(A, i, j));
      if (err != RABIN_SUCCESS) return err;
    }
  }
  return RABIN_SUCCESS;
}

rabin_err_t rmat_trace(rz_t* t, const rmat_t* A)
{
  if (t == NULL || A == NULL) return RABIN_ERR_NULL_PTR;
  if (A->cols != A->rows) return RABIN_ERR_MATRIX_DIM;

  u64 n = A->cols;
  rabin_err_t err = rz_set_u64(t, 0);
  if (err != RABIN_SUCCESS) return err;

  for (u64 i = 0; i < n; i++) {
    if ((err = rz_add(t, t, RMAT_GET(A, i, i))) != RABIN_SUCCESS) return err;
  }
  return RABIN_SUCCESS;
}

rabin_err_t rmat_swap(rmat_t* a, rmat_t* b)
{
  if (a == NULL || b == NULL) return RABIN_ERR_NULL_PTR;

  rmat_t t = *a;
  *a = *b;
  *b = t;
  return RABIN_SUCCESS;
}

rabin_err_t rmat_charpoly_adj(rpol_t* p, rmat_t* J, const rmat_t* A)
{
  if (p == NULL || A == NULL) return RABIN_ERR_NULL_PTR;
  if (A->cols != A->rows) return RABIN_ERR_MATRIX_DIM;

  // 1. [Initialize]
  u64 n = A->cols;

  rabin_err_t err = RABIN_SUCCESS;
  rmat_t C = {0}, T = {0};
  rz_t temp, trace;
  rz_init_multi(&temp, &trace, NULL);
  rz_t* a = NULL;

  TRY(rmat_init(&C, n, n));
  TRY(rmat_init(&T, n, n));
  TRY(rmat_id(&C, n));

  a = malloc((n + 1) * sizeof(rz_t));
  if (a == NULL) {
    err = RABIN_ERR_OUT_OF_MEMORY;
    goto cleanup;
  }
  for (u64 i = 0; i <= n; i++) {
    rz_init(&a[i]);
  }
  TRY(rz_set_u64(&a[0], 1));

  // 2. [Finished?]
  for (u64 i = 1; i <= n; i++) {
    // 3. [Compute next a_i and C]
    TRY(rmat_mul(&T, A, &C));

    TRY(rmat_trace(&trace, &T));

    TRY(rz_set_u64(&temp, i));
    rz_divmod_u64(&a[i], &trace, i);

    TRY(rz_neg(&a[i], &a[i]));

    if (J && i == n) {
      TRY(rmat_copy(J, &C));
      if ((n - 1) & 1) {
        TRY(rmat_neg(J));
      }
    }
    TRY(rmat_copy(&C, &T));
    for (u64 j = 0; j < n; j++) {
      TRY(rz_add(&temp, RMAT_GET(&C, j, j), &a[i]));
      TRY(rmat_set(&C, &temp, j, j));
    }
  }

  TRY(rpol_alloc(p, n + 1));

  for (u64 i = 0; i <= n; i++) {
    TRY(rz_copy(&p->coeff[n - i], &a[i]));
    rz_clear(&a[i]);
  }

  p->deg = n;

cleanup:
  if (a != NULL) {
    for (u64 i = 0; i <= n; i++) {
      rz_clear(&a[i]);
    }
    free(a);
  }
  rmat_clear(&T);
  rmat_clear(&C);
  rz_clear_multi(&temp, &trace, NULL);
  return err;
}

rabin_err_t rmat_gcd_all(rz_t* g, const rmat_t* A)
{
  if (g == NULL || A == NULL) return RABIN_ERR_NULL_PTR;

  u64 size = A->cols * A->rows;

  if (size == 0) {
    return rz_set_u64(g, 1);
  }

  // start with first element
  rabin_err_t err = rz_copy(g, &A->data[0]);
  if (err != RABIN_SUCCESS) return err;

  for (u64 i = 1; i < size; i++) {
    if (rz_is_eq_i64(g, 1)) return RABIN_SUCCESS;

    if (rz_is_zero(&A->data[i])) continue;

    if ((err = rz_gcd(g, g, &A->data[i])) != RABIN_SUCCESS) return err;
  }
  return RABIN_SUCCESS;
}

rabin_err_t rmat_id(rmat_t* I, const u64 n)
{
  if (I == NULL) return RABIN_ERR_NULL_PTR;
  if (I->rows != n || I->cols != n) return RABIN_ERR_MATRIX_DIM;

  rabin_err_t err = RABIN_SUCCESS;
  rz_t a, b;
  rz_init_multi(&a, &b, NULL);
  TRY(rz_set_u64(&a, 1));
  TRY(rz_set_u64(&b, 0));

  for (u64 i = 0; i < n; i++) {
    for (u64 j = 0; j < n; j++) {
      if (i == j) {
        TRY(rmat_set(I, &a, i, j));
      } else {
        TRY(rmat_set(I, &b, i, j));
      }
    }
  }

cleanup:
  rz_clear_multi(&a, &b, NULL);
  return err;
}

rabin_err_t rmat_neg_col(rmat_t* A, u64 j)
{
  if (A == NULL) return RABIN_ERR_NULL_PTR;
  if (j >= A->cols) return RABIN_ERR_INVALID_ARG;

  for (u64 i = 0; i < A->rows; i++) {
    rabin_err_t err = rz_neg(RMAT_GET(A, i, j), RMAT_GET(A, i, j));
    if (err != RABIN_SUCCESS) return err;
  }
  return RABIN_SUCCESS;
}

static bool rmat_row_prefix_is_zero(const rmat_t* A, u64 i, u64 k)
{
  for (u64 j = 0; j < k; j++) {
    if (!rz_is_zero(RMAT_GET(A, i, j))) {
      return false;
    }
  }
  return true;
}

static u64 rmat_min_abs_index(const rmat_t* A, int i, int k)
{
  int j, j0 = -1;
  for (j = 0; j <= k; j++) {
    if (rz_is_zero(RMAT_GET(A, i, j))) continue;

    if (j0 == -1 || rz_cmp_abs(RMAT_GET(A, i, j), RMAT_GET(A, i, j0)) < 0)
      j0 = j;
  }

  return j0;
}

static rabin_err_t rmat_swap_col(rmat_t* A, u64 k1, u64 k2)
{
  if (A == NULL) return RABIN_ERR_NULL_PTR;
  if (k1 >= A->cols || k2 >= A->cols) return RABIN_ERR_INVALID_ARG;

  for (u64 r = 0; r < A->rows; r++) {
    rabin_err_t err = rz_swap(RMAT_GET(A, r, k1), RMAT_GET(A, r, k2));
    if (err != RABIN_SUCCESS) return err;
  }
  return RABIN_SUCCESS;
}

rabin_err_t rmat_hermite(rmat_t* W, const rmat_t* A)
{
  if (W == NULL || A == NULL) return RABIN_ERR_NULL_PTR;

  // 1. [Initialize]
  u64 m = A->rows;
  u64 n = A->cols;
  if (m == 0 || n == 0) return RABIN_SUCCESS;

  i64 i = m - 1;
  i64 k = n - 1;
  u64 l = 0;
  if (m > n) {
    l = m - n;
  }

  rabin_err_t err = RABIN_SUCCESS;
  rmat_t A_work = {0};
  TRY(rmat_init(&A_work, m, n));
  TRY(rmat_copy(&A_work, A));

  rz_t b, q, rem, temp, one, t2;
  rz_init_multi(&b, &q, &rem, &temp, &one, &t2, NULL);
  TRY(rz_set_u64(&one, 1));

// 2. [Row finished?]
step2:
  if (rmat_row_prefix_is_zero(&A_work, i, k)) {
    if (RMAT_GET(&A_work, i, k)->is_neg) {
      TRY(rmat_neg_col(&A_work, k));
    }
    goto step5;
  }

  // 3. [Choose non-zero entry]
  u64 j_0 = rmat_min_abs_index(&A_work, i, k);
  if ((i64)j_0 < k) {
    TRY(rmat_swap_col(&A_work, k, j_0));
  }
  if (RMAT_GET(&A_work, i, k)->is_neg) {
    TRY(rmat_neg_col(&A_work, k));
  }

  TRY(rz_copy(&b, RMAT_GET(&A_work, i, k)));

  // 4. [Reduce]
  for (u64 j = 0; j < (u64)k; j++) {
    TRY(rz_div_euclid(&q, RMAT_GET(&A_work, i, j), &b));
    for (u64 x = 0; x < m; x++) {
      TRY(rz_mul(&temp, &q, RMAT_GET(&A_work, x, k)));
      TRY(rz_sub(&t2, RMAT_GET(&A_work, x, j), &temp));
      TRY(rmat_set(&A_work, &t2, x, j));
    }
  }
  goto step2;

// 5. [Final reductions]
step5:
  TRY(rz_copy(&b, RMAT_GET(&A_work, i, k)));
  if (rz_is_zero(&b)) {
    k++;
    goto step6;
  }

  for (u64 j = k + 1; j < n; j++) {
    TRY(rz_div_euclid(&q, RMAT_GET(&A_work, i, j), &b));
    for (u64 x = 0; x < m; x++) {
      TRY(rz_mul(&temp, &q, RMAT_GET(&A_work, x, k)));
      TRY(rz_sub(&t2, RMAT_GET(&A_work, x, j), &temp));
      TRY(rmat_set(&A_work, &t2, x, j));
    }
  }

// 6. [Finished?]
step6:
  if (i != (i64)l) {
    i--;
    k--;
    goto step2;
  }

  TRY(rmat_init(W, m, n - k));

  // copy result
  for (u64 j = 0; j < n - k; j++) {
    for (u64 x = 0; x < m; x++) {
      TRY(rz_copy(RMAT_GET(W, x, j), RMAT_GET(&A_work, x, j + k)));
    }
  }

cleanup:
  rz_clear_multi(&b, &q, &rem, &temp, &one, &t2, NULL);
  rmat_clear(&A_work);
  return err;
}

static rabin_err_t rmat_hnf_centered_mod(rz_t* r, const rz_t* x, const rz_t* R,
                                         const rz_t* half)
{
  if (r == NULL || x == NULL || R == NULL || half == NULL)
    return RABIN_ERR_NULL_PTR;

  if (rz_cmp_abs(x, half) <= 0) {
    return rz_copy(r, x);
  }

  rabin_err_t err = rz_mod_pos(r, x, R);
  if (err != RABIN_SUCCESS) return err;
  if (rz_cmp(r, half) > 0) {
    if ((err = rz_sub(r, r, R)) != RABIN_SUCCESS) return err;
  }
  return RABIN_SUCCESS;
}

// Cohen's "Important Remark" after Algorithm 2.4.5: the Euclidean steps of
// the HNF algorithms need the Bezout pair (u, v) of u*a + v*b = d with
// small coefficients. All solutions are (u0 + t*(b/d), v0 - t*(a/d)), so
// this rewrites the pair in place to the unique one with v in the centered
// half-period (-|a/d|/2, |a/d|/2] (which minimizes |v|); when a | b
// (d = |a|) the book's essential condition v = 0, u = sign(a) applies
// instead. The identity u*a + v*b = d is preserved in all cases.
static rabin_err_t hnf_bezout_minimal(rz_t* u, rz_t* v, const rz_t* d,
                                      const rz_t* a, const rz_t* b)
{
  rz_t p, pmag, half, quot, r, t, two;
  rz_init_multi(&p, &pmag, &half, &quot, &r, &t, &two, NULL);
  rabin_err_t err = RABIN_SUCCESS;

  if (rz_is_zero(a)) {
    // v*sign(b) = d/|b| = 1 is forced; u is free, take u = 0
    TRY(rz_set_u64(u, 0));
    if (rz_is_zero(b)) {
      TRY(rz_set_u64(v, 0));
    } else if (b->is_neg) {
      TRY(rz_set_i64(v, -1));
    } else {
      TRY(rz_set_u64(v, 1));
    }
    goto cleanup;
  }

  if (rz_is_zero(b)) {
    // u*sign(a) = 1 is forced; v = 0
    if (a->is_neg) {
      TRY(rz_set_i64(u, -1));
    } else {
      TRY(rz_set_u64(u, 1));
    }
    TRY(rz_set_u64(v, 0));
    goto cleanup;
  }

  TRY(rz_copy(&p, a));
  p.is_neg = false;

  // a | b  =>  d = |a|  =>  v = 0, u = sign(a)  (the book's essential
  // condition; it also keeps the Euclidean step from stalling)
  if (rz_cmp(&p, d) == 0) {
    if (a->is_neg) {
      TRY(rz_set_i64(u, -1));
    } else {
      TRY(rz_set_u64(u, 1));
    }
    TRY(rz_set_u64(v, 0));
    goto cleanup;
  }

  // centered half-period: v = v0 - t*(a/d) with |a/d| >= 2 and
  // v in (-|a/d|/2, |a/d|/2]; recompute u = (d - v*b)/a (exact)
  TRY(rz_div_exact(&p, a, d));
  TRY(rz_copy(&pmag, &p));
  pmag.is_neg = false;
  TRY(rz_set_u64(&two, 2));
  TRY(rz_div(&half, &pmag, &two));

  // r = v0 mod |p| in [0, |p|)
  TRY(rz_divmod(&quot, &r, v, &pmag));
  if (r.is_neg) {
    TRY(rz_add(&r, &r, &pmag));
  }
  if (rz_cmp(&r, &half) > 0) {
    TRY(rz_sub(v, &r, &pmag));
  } else {
    TRY(rz_copy(v, &r));
  }

  TRY(rz_mul(&t, v, b));
  TRY(rz_sub(&t, d, &t));
  TRY(rz_div_exact(u, &t, a));

cleanup:
  rz_clear_multi(&p, &pmag, &half, &quot, &r, &t, &two, NULL);
  return err;
}

/*
 * One [Euclidean step] of Cohen 2.4.5 (shared by rmat_hermite_gcd and
 * hnf_with_transform): with pivot row i and columns j < k, compute the
 * minimal Bezout pair (u, v, d) of (a_{i,k}, a_{i,j}) and apply the
 * column operation
 *
 *   col_k <- u*col_k + v*col_j                       (places d at row i)
 *   col_j <- (a_{i,k}/d)*col_j - (a_{i,j}/d)*col_k  (places 0 at row i)
 *
 * buf must hold A_work->rows initialized scratch rz_t values. If
 * U != NULL (with scratch vector C of U->rows values, both non-NULL or
 * both NULL) the identical operations are applied to columns k and j of
 * U, maintaining the invariant orig * U = A_work.
 */
static rabin_err_t hnf_euclid_step(rmat_t* A_work, rmat_t* U, rz_t* buf,
                                   rz_t* C, u64 i, u64 j, u64 k, rz_t* u,
                                   rz_t* v, rz_t* d, rz_t* q_k, rz_t* q_j,
                                   rz_t* temp, rz_t* temp2)
{
  if (A_work == NULL || buf == NULL || u == NULL || v == NULL || d == NULL ||
      q_k == NULL || q_j == NULL || temp == NULL || temp2 == NULL ||
      (U == NULL) != (C == NULL))
    return RABIN_ERR_NULL_PTR;

  rabin_err_t err = rz_gcd_extended_lehmer(u, v, d, RMAT_GET(A_work, i, k),
                                           RMAT_GET(A_work, i, j));
  if (err != RABIN_SUCCESS) return err;
  if ((err = hnf_bezout_minimal(u, v, d, RMAT_GET(A_work, i, k),
                                RMAT_GET(A_work, i, j))) != RABIN_SUCCESS)
    return err;

  // buf[x] = u*A[x][k] + v*A[x][j]
  for (u64 x = 0; x < A_work->rows; x++) {
    if ((err = rz_mul(&buf[x], u, RMAT_GET(A_work, x, k))) != RABIN_SUCCESS)
      return err;
    if ((err = rz_mul(temp, v, RMAT_GET(A_work, x, j))) != RABIN_SUCCESS)
      return err;
    if ((err = rz_add(&buf[x], &buf[x], temp)) != RABIN_SUCCESS) return err;
  }

  if ((err = rz_div_euclid(q_k, RMAT_GET(A_work, i, k), d)) != RABIN_SUCCESS)
    return err;
  if ((err = rz_div_euclid(q_j, RMAT_GET(A_work, i, j), d)) != RABIN_SUCCESS)
    return err;

  // col_j <- q_k*col_j - q_j*col_k,  col_k <- buf
  for (u64 x = 0; x < A_work->rows; x++) {
    if ((err = rz_mul(temp, q_k, RMAT_GET(A_work, x, j))) != RABIN_SUCCESS)
      return err;
    if ((err = rz_mul(temp2, q_j, RMAT_GET(A_work, x, k))) != RABIN_SUCCESS)
      return err;
    if ((err = rz_sub(RMAT_GET(A_work, x, j), temp, temp2)) != RABIN_SUCCESS)
      return err;
    if ((err = rz_copy(RMAT_GET(A_work, x, k), &buf[x])) != RABIN_SUCCESS)
      return err;
  }

  if (U != NULL) {
    for (u64 x = 0; x < U->rows; x++) {
      if ((err = rz_mul(&C[x], u, RMAT_GET(U, x, k))) != RABIN_SUCCESS)
        return err;
      if ((err = rz_mul(temp, v, RMAT_GET(U, x, j))) != RABIN_SUCCESS)
        return err;
      if ((err = rz_add(&C[x], &C[x], temp)) != RABIN_SUCCESS) return err;
    }
    for (u64 x = 0; x < U->rows; x++) {
      if ((err = rz_mul(temp, q_k, RMAT_GET(U, x, j))) != RABIN_SUCCESS)
        return err;
      if ((err = rz_mul(temp2, q_j, RMAT_GET(U, x, k))) != RABIN_SUCCESS)
        return err;
      if ((err = rz_sub(RMAT_GET(U, x, j), temp, temp2)) != RABIN_SUCCESS)
        return err;
      if ((err = rz_copy(RMAT_GET(U, x, k), &C[x])) != RABIN_SUCCESS)
        return err;
    }
  }
  return RABIN_SUCCESS;
}

rabin_err_t rmat_hermite_mod_d(rmat_t* W, const rmat_t* A, const rz_t* D)
{
  if (W == NULL || A == NULL || D == NULL) return RABIN_ERR_NULL_PTR;
  if (rz_is_zero(D)) return RABIN_ERR_DIV_BY_ZERO;

  // 1. [Initialize]
  u64 m = A->rows;
  u64 n = A->cols;
  if (m == 0 || n == 0) return RABIN_SUCCESS;

  if (m > n) {
    return rmat_hermite_gcd(W, A);
  }

  rabin_err_t err = RABIN_SUCCESS;
  rmat_t A_work = {0}, W_local = {0};
  rz_t R, u, v, d, q, t1, t2, t3, half;
  rz_init_multi(&R, &u, &v, &d, &q, &t1, &t2, &t3, &half, NULL);
  rz_t* B = NULL;

  TRY(rmat_init(&A_work, m, n));
  TRY(rmat_copy(&A_work, A));
  TRY(rz_copy(&R, D));
  R.is_neg = false;
  TRY(rz_rshift(&half, &R, 1));

  // reduce the input entries into (-R/2, R/2]; all column operations are
  // taken modulo R, so only the classes modulo R matter
  for (u64 x = 0; x < m * n; x++) {
    TRY(rmat_hnf_centered_mod(&A_work.data[x], &A_work.data[x], &R, &half));
  }

  if (m > 0) {
    B = malloc(m * sizeof(rz_t));
    if (B == NULL) {
      err = RABIN_ERR_OUT_OF_MEMORY;
      goto cleanup;
    }
    for (u64 x = 0; x < m; x++) {
      rz_init(&B[x]);
    }
  }

  TRY(rmat_init(&W_local, m, m));

  for (i64 i = (i64)m - 1; i >= 0; i--) {
    i64 k = (i64)n - (i64)m + i;

    // 2. [Check zero] / 3. [Euclidean step]
    for (i64 j = k - 1; j >= 0; j--) {
      if (rz_is_zero(RMAT_GET(&A_work, i, j))) {
        continue;
      }

      TRY(rz_gcd_extended_lehmer(&u, &v, &d, RMAT_GET(&A_work, i, k),
                                 RMAT_GET(&A_work, i, j)));
      TRY(hnf_bezout_minimal(&u, &v, &d, RMAT_GET(&A_work, i, k),
                             RMAT_GET(&A_work, i, j)));

      // B = u*A_k + v*A_j, reduced into (-R/2, R/2]
      for (u64 x = 0; x < m; x++) {
        TRY(rz_mul(&B[x], &u, RMAT_GET(&A_work, x, k)));
        TRY(rz_mul(&t2, &v, RMAT_GET(&A_work, x, j)));
        TRY(rz_add(&B[x], &B[x], &t2));
        TRY(rmat_hnf_centered_mod(&B[x], &B[x], &R, &half));
      }

      // A_j = (a_{i,k}/d)*A_j - (a_{i,j}/d)*A_k, reduced into (-R/2, R/2]
      TRY(rz_div_exact(&t2, RMAT_GET(&A_work, i, k), &d));
      TRY(rz_div_exact(&q, RMAT_GET(&A_work, i, j), &d));
      for (u64 x = 0; x < m; x++) {
        TRY(rz_mul(&t3, &t2, RMAT_GET(&A_work, x, j)));
        TRY(rz_mul(&t1, &q, RMAT_GET(&A_work, x, k)));
        TRY(rz_sub(&t1, &t3, &t1));
        TRY(rmat_hnf_centered_mod(RMAT_GET(&A_work, x, j), &t1, &R, &half));
        TRY(rz_copy(RMAT_GET(&A_work, x, k), &B[x]));
      }

      // a_{i,k} = u*a_{i,k} + v*a_{i,j} = d exactly
      TRY(rz_copy(RMAT_GET(&A_work, i, k), &d));
    }

    // 4. [Next row]
    // u*a_{i,k} + v*R = d = gcd(a_{i,k}, R)
    TRY(rz_gcd_extended_lehmer(&u, &v, &d, RMAT_GET(&A_work, i, k), &R));
    TRY(hnf_bezout_minimal(&u, &v, &d, RMAT_GET(&A_work, i, k), &R));

    // W_i = u*A_k mod R, taken in [0, R)
    for (u64 x = 0; x < m; x++) {
      TRY(rz_mul(&t1, &u, RMAT_GET(&A_work, x, k)));
      TRY(rz_mod_pos(RMAT_GET(&W_local, x, i), &t1, &R));
    }

    // if d = R (i.e. R | a_{i,k}) the diagonal entry would be 0; set it to d
    if (rz_cmp(&d, &R) == 0) {
      TRY(rz_copy(RMAT_GET(&W_local, i, i), &d));
    }

    // final reductions: W_j -= floor(W_{i,j}/W_{i,i}) * W_i for j > i
    for (u64 j = (u64)i + 1; j < m; j++) {
      if (rz_is_zero(RMAT_GET(&W_local, i, j))) {
        continue;
      }
      TRY(rz_div_euclid(&q, RMAT_GET(&W_local, i, j),
                        RMAT_GET(&W_local, i, i)));
      for (u64 x = 0; x <= (u64)i; x++) {
        TRY(rz_mul(&t1, &q, RMAT_GET(&W_local, x, i)));
        TRY(rz_sub(RMAT_GET(&W_local, x, j), RMAT_GET(&W_local, x, j), &t1));
      }
    }

    R.is_neg = false;
    TRY(rz_div_exact(&R, &R, &d));
    TRY(rz_rshift(&half, &R, 1));

    if (i > 0) {
      // working modulo R, a_{i-1,k-1} may have reduced to zero; replace it
      // by any nonzero multiple of R
      if (rz_is_zero(RMAT_GET(&A_work, i - 1, k - 1))) {
        TRY(rz_copy(RMAT_GET(&A_work, i - 1, k - 1), &R));
      }
    }
  }

  TRY(rmat_swap(W, &W_local));
  rmat_clear(&W_local);

cleanup:
  if (B != NULL) {
    for (u64 x = 0; x < m; x++) {
      rz_clear(&B[x]);
    }
    free(B);
  }
  rz_clear_multi(&R, &u, &v, &d, &q, &t1, &t2, &t3, &half, NULL);
  rmat_clear(&A_work);
  rmat_clear(&W_local);
  return err;
}

rabin_err_t rmat_hermite_gcd(rmat_t* W, const rmat_t* A)
{
  if (W == NULL || A == NULL) return RABIN_ERR_NULL_PTR;

  // 1. [Initialize]
  u64 m = A->rows;
  u64 n = A->cols;
  if (m == 0 || n == 0) return RABIN_SUCCESS;

  u64 i = m - 1;
  u64 j = n - 1;
  u64 k = n - 1;
  u64 l = 0;
  if (m > n) {
    l = m - n;
  }

  rabin_err_t err = RABIN_SUCCESS;
  rmat_t A_work = {0};
  rz_t* B = NULL;

  TRY(rmat_init(&A_work, m, n));
  TRY(rmat_copy(&A_work, A));

  B = malloc(m * sizeof(rz_t));
  if (B == NULL) {
    err = RABIN_ERR_OUT_OF_MEMORY;
    goto cleanup;
  }
  for (u64 x = 0; x < m; x++) {
    rz_init(&B[x]);
  }

  rz_t u, v, d, temp, temp2, temp3, one, b, q_k, q_j;
  rz_init_multi(&u, &v, &d, &temp, &temp2, &temp3, &one, &b, &q_k, &q_j, NULL);
  TRY(rz_set_u64(&one, 1));

// 2. [Check zero]
step2:
  if (j == 0) {
    goto step4;
  }
  j--;
  if (rz_is_zero(RMAT_GET(&A_work, i, j))) {
    goto step2;
  }

  // 3. [Euclidean step]
  TRY(hnf_euclid_step(&A_work, NULL, B, NULL, i, j, k, &u, &v, &d, &q_k, &q_j,
                      &temp, &temp2));
  goto step2;

  // 4. [Final reduction]
step4:
  TRY(rz_copy(&b, RMAT_GET(&A_work, i, k)));
  if (b.is_neg) {  // zero is never negative, no extra rz_is_zero test needed
    for (u64 x = 0; x < m; x++) {
      TRY(rz_neg(RMAT_GET(&A_work, x, k), RMAT_GET(&A_work, x, k)));
    }
    TRY(rz_neg(&b, &b));
  }

  if (rz_is_zero(&b)) {
    k++;
    goto step5;
  }

  for (u64 j_0 = k + 1; j_0 < n; j_0++) {
    TRY(rz_div_euclid(&temp2, RMAT_GET(&A_work, i, j_0), &b));
    for (u64 x = 0; x < m; x++) {
      TRY(rz_mul(&temp, &temp2, RMAT_GET(&A_work, x, k)));
      TRY(rz_sub(RMAT_GET(&A_work, x, j_0), RMAT_GET(&A_work, x, j_0), &temp));
    }
  }

step5:
  // 5. [Finished]
  if (i == l) {
    TRY(rmat_init(W, m, n - k));
    for (u64 x = 0; x < n - k; x++) {
      for (u64 y = 0; y < m; y++) {
        TRY(rz_copy(RMAT_GET(W, y, x), RMAT_GET(&A_work, y, x + k)));
      }
    }
    goto cleanup;
  }

  i--;
  k--;
  j = k;
  goto step2;

cleanup:
  if (B != NULL) {
    for (u64 x = 0; x < m; x++) {
      rz_clear(&B[x]);
    }
    free(B);
  }
  rmat_clear(&A_work);
  rz_clear_multi(&u, &v, &d, &temp, &temp2, &temp3, &one, &b, &q_k, &q_j, NULL);
  return err;
}

rabin_err_t rmat_smith(rmat_t* S, const rmat_t* A)
{
  if (S == NULL || A == NULL) return RABIN_ERR_NULL_PTR;
  if (A->cols != A->rows) return RABIN_ERR_MATRIX_DIM;
  if (S->rows != A->rows || S->cols != A->cols) return RABIN_ERR_MATRIX_DIM;

  // 1. [Initialize i]
  i64 n = A->rows;
  i64 i = n;
  rabin_err_t err = RABIN_SUCCESS;
  rz_t R, u, v, d, q_i, q_j, t1, t2, b;
  rz_t a_abs, a_ii_abs, a_abs_col, a_ii_abs_col, zero;
  rz_init_multi(&R, &u, &v, &d, &q_i, &q_j, &t1, &t2, &b, NULL);
  rz_init_multi(&a_abs, &a_ii_abs, &a_abs_col, &a_ii_abs_col, &zero, NULL);
  TRY(rmat_det(&R, A));
  R.is_neg = false;

  if (n == 1) {
    err = rmat_set(S, &R, 0, 0);
    goto cleanup;
  }

  rmat_t A_work = {0};
  rz_t* B = NULL;
  TRY(rmat_init(&A_work, n, n));
  TRY(rmat_copy(&A_work, A));

  B = malloc(n * sizeof(rz_t));
  if (B == NULL) {
    err = RABIN_ERR_OUT_OF_MEMORY;
    goto cleanup;
  }
  for (u64 x = 0; x < (u64)n; x++) {
    rz_init(&B[x]);
  }

  rz_t* a_i_i = NULL;
  rz_t* a_i_j = NULL;

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
  if (rz_is_zero(RMAT_GET(&A_work, i - 1, j - 1))) {
    goto step3;
  }

  // 4. [Euclidean step]
  a_i_i = RMAT_GET(&A_work, i - 1, i - 1);
  a_i_j = RMAT_GET(&A_work, i - 1, j - 1);

  if ((err = rz_gcd_extended_lehmer(&u, &v, &d, a_i_i, a_i_j)) != RABIN_SUCCESS)
    goto cleanup;

  // use Remark to find u,v minimal
  TRY(rz_copy(&a_abs, RMAT_GET(&A_work, i - 1, j - 1)));
  a_abs.is_neg = false;
  TRY(rz_copy(&a_ii_abs, RMAT_GET(&A_work, i - 1, i - 1)));
  a_ii_abs.is_neg = false;

  // if d == |a_{i,i}|, force u = sign(a_{i,i}), v = 0
  if (rz_cmp(&d, &a_ii_abs) == 0) {
    TRY(rz_set_u64(&u, 1));
    TRY(rz_set_u64(&v, 0));
    if (RMAT_GET(&A_work, i - 1, i - 1)->is_neg) TRY(rz_set_i64(&u, -1));
  }
  // if d == |a_{i,j}|, force u = 0, v = sign(a_{i,j})
  else if (rz_cmp(&d, &a_abs) == 0) {
    TRY(rz_set_u64(&u, 0));
    TRY(rz_set_u64(&v, 1));
    if (RMAT_GET(&A_work, i - 1, j - 1)->is_neg) TRY(rz_set_i64(&v, -1));
  }

  rz_clear(&a_abs);
  rz_clear(&a_ii_abs);

  // B = uA_i + vA_j
  for (u64 x = 0; x < (u64)n; x++) {
    TRY(rz_mul(&B[x], &u, RMAT_GET(&A_work, x, i - 1)));
    TRY(rz_mul(&t1, &v, RMAT_GET(&A_work, x, j - 1)));
    TRY(rz_add(&B[x], &B[x], &t1));
  }

  // A_j = ((a_ii / d)A_j - (a_ij / d)A_i) mod R
  TRY(rz_div(&q_i, RMAT_GET(&A_work, i - 1, i - 1), &d));
  TRY(rz_div(&q_j, RMAT_GET(&A_work, i - 1, j - 1), &d));
  for (u64 x = 0; x < (u64)n; x++) {
    TRY(rz_mul(&t1, &q_i, RMAT_GET(&A_work, x, j - 1)));
    TRY(rz_mul(&t2, &q_j, RMAT_GET(&A_work, x, i - 1)));
    TRY(rz_sub(&t1, &t1, &t2));
    TRY(rz_mod(&t1, &t1, &R));
    TRY(rz_copy(RMAT_GET(&A_work, x, j - 1), &t1));
  }

  // A_i = B mod R
  for (u64 x = 0; x < (u64)n; x++) {
    TRY(rz_mod(RMAT_GET(&A_work, x, i - 1), &B[x], &R));
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
  if (rz_is_zero(RMAT_GET(&A_work, j - 1, i - 1))) {
    goto step6;
  }

  // 7. [Euclidean step]
  a_i_i = RMAT_GET(&A_work, i - 1, i - 1);
  a_i_j = RMAT_GET(&A_work, j - 1, i - 1);

  if ((err = rz_gcd_extended_lehmer(&u, &v, &d, a_i_i, a_i_j)) != RABIN_SUCCESS)
    goto cleanup;

  // use Remark to find u,v minimal safely
  TRY(rz_copy(&a_abs_col, RMAT_GET(&A_work, j - 1, i - 1)));
  a_abs_col.is_neg = false;
  TRY(rz_copy(&a_ii_abs_col, RMAT_GET(&A_work, i - 1, i - 1)));
  a_ii_abs_col.is_neg = false;

  // if d == |a_{i,i}|, force u = sign(a_{i,i}), v = 0
  if (rz_cmp(&d, &a_ii_abs_col) == 0) {
    TRY(rz_set_u64(&u, 1));
    TRY(rz_set_u64(&v, 0));
    if (RMAT_GET(&A_work, i - 1, i - 1)->is_neg) TRY(rz_set_i64(&u, -1));
  }
  // if d == |a_{j,i}|, force u = 0, v = sign(a_{j,i})
  else if (rz_cmp(&d, &a_abs_col) == 0) {
    TRY(rz_set_u64(&u, 0));
    TRY(rz_set_u64(&v, 1));
    if (RMAT_GET(&A_work, j - 1, i - 1)->is_neg) TRY(rz_set_i64(&v, -1));
  }

  rz_clear(&a_abs_col);
  rz_clear(&a_ii_abs_col);

  // B = uA'_i + vA'_j
  for (u64 x = 0; x < (u64)n; x++) {
    TRY(rz_mul(&B[x], &u, RMAT_GET(&A_work, i - 1, x)));
    TRY(rz_mul(&t1, &v, RMAT_GET(&A_work, j - 1, x)));
    TRY(rz_add(&B[x], &B[x], &t1));
  }

  // A'_j = ((a_ii / d)A'_j - (a_ij / d)A'_i) mod R
  TRY(rz_div(&q_i, RMAT_GET(&A_work, i - 1, i - 1), &d));
  TRY(rz_div(&q_j, RMAT_GET(&A_work, j - 1, i - 1), &d));
  for (u64 x = 0; x < (u64)n; x++) {
    TRY(rz_mul(&t1, &q_i, RMAT_GET(&A_work, j - 1, x)));
    TRY(rz_mul(&t2, &q_j, RMAT_GET(&A_work, i - 1, x)));
    TRY(rz_sub(&t1, &t1, &t2));
    TRY(rz_mod(&t1, &t1, &R));
    TRY(rz_copy(RMAT_GET(&A_work, j - 1, x), &t1));
  }

  // A'_i = B mod R
  for (u64 x = 0; x < (u64)n; x++) {
    TRY(rz_mod(RMAT_GET(&A_work, i - 1, x), &B[x], &R));
  }

  c++;
  goto step6;

// 8. [Repeat stage i?]
step8:
  if (c > 0) {
    goto step2;
  }

  // 9. [Check the rest of the matrix]
  TRY(rz_copy(&b, RMAT_GET(&A_work, i - 1, i - 1)));
  for (i64 k_0 = 0; k_0 < i - 1; k_0++) {
    for (i64 l_0 = 0; l_0 < i - 1; l_0++) {
      bool not_divisible = false;
      if (rz_is_zero(&b)) {
        if (!rz_is_zero(RMAT_GET(&A_work, k_0, l_0))) {
          not_divisible = true;
        }
      } else {
        TRY(rz_mod(&t1, RMAT_GET(&A_work, k_0, l_0), &b));
        if (!rz_is_zero(&t1)) {
          not_divisible = true;
        }
      }

      if (not_divisible) {
        for (i64 x = 0; x < n; x++) {
          TRY(rz_add(RMAT_GET(&A_work, i - 1, x), RMAT_GET(&A_work, i - 1, x),
                     RMAT_GET(&A_work, k_0, x)));
        }
        goto step2;
      }
    }
  }

  // 10. [Next stage]
  TRY(rz_gcd_lehmer(&t1, RMAT_GET(&A_work, i - 1, i - 1), &R));
  TRY(rmat_set(&A_work, &t1, i - 1, i - 1));
  TRY(rz_div(&R, &R, &t1));

  if (i == 2) {
    TRY(rz_gcd_lehmer(&t2, RMAT_GET(&A_work, 0, 0), &R));
    TRY(rmat_set(&A_work, &t2, 0, 0));

    // order diagonal entries in proper order
    for (i64 r = 0; r < n; r++) {
      for (i64 c = 0; c < n; c++) {
        if (r == c) {
          TRY(rmat_set(S, RMAT_GET(&A_work, n - 1 - r, n - 1 - r), r, r));
        } else {
          TRY(rmat_set(S, &zero, r, c));
        }
      }
    }
    rz_clear(&zero);
    goto cleanup;
  }

  i--;
  goto step2;

cleanup:
  if (B != NULL) {
    for (i64 x = 0; x < n; x++) {
      rz_clear(&B[x]);
    }
    free(B);
  }
  rz_clear(&a_abs_col);
  rz_clear(&a_ii_abs_col);
  rz_clear(&a_abs);
  rz_clear(&a_ii_abs);
  rz_clear(&zero);
  rmat_clear(&A_work);
  rz_clear(&R);
  rz_clear_multi(&u, &v, &d, &q_i, &q_j, &t1, &t2, &b, NULL);
  return err;
}

bool rmat_equal(const rmat_t* A, const rmat_t* B)
{
  if (A == NULL || B == NULL) return false;
  if (A->rows != B->rows || A->cols != B->cols) return false;

  u64 total = A->rows * A->cols;
  for (u64 i = 0; i < total; i++) {
    const rz_t* a = &A->data[i];
    const rz_t* b = &B->data[i];

    // rz_cmp does not normalize the two representations of zero
    // (size 0 from rmat_init vs size 1 with a zero limb), so
    // settle the zero cases explicitly
    bool za = rz_is_zero(a);
    bool zb = rz_is_zero(b);
    if (za || zb) {
      if (za != zb) return false;
      continue;
    }
    if (rz_cmp(a, b) != 0) return false;
  }
  return true;
}

bool rmat_hnf_check_structure(const rmat_t* H)
{
  if (H == NULL) return false;

  u64 r = H->rows;
  u64 c = H->cols;

  // once a zero row is hit, every later row must be zero as well
  bool zero_row_seen = false;

  for (u64 i = 0; i < r; i++) {
    rz_t* row = &H->data[i * c];

    // upper triangular: H[i][j] == 0 for all j < i (capped at c for
    // rows below the last column)
    u64 left = RZ_MIN(i, c);
    for (u64 j = 0; j < left; j++) {
      if (!rz_is_zero(&row[j])) {
        return false;
      }
    }

    // zero rows only in a trailing block: from here on everything must
    // be zero
    if (zero_row_seen) {
      for (u64 j = 0; j < c; j++) {
        if (!rz_is_zero(&row[j])) {
          return false;
        }
      }
      continue;
    }

    // a nonzero row needs a diagonal pivot, so rows at or beyond the
    // last column can only be zero rows
    if (i >= c) {
      for (u64 j = 0; j < c; j++) {
        if (!rz_is_zero(&row[j])) {
          return false;
        }
      }
      zero_row_seen = true;
      continue;
    }

    rz_t* pivot = &row[i];

    if (rz_is_zero(pivot)) {
      // a zero pivot may only start the trailing zero block: everything
      // right of the pivot must also be zero
      for (u64 j = i + 1; j < c; j++) {
        if (!rz_is_zero(&row[j])) {
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
      if (row[j].is_neg || rz_cmp(&row[j], pivot) >= 0) {
        return false;
      }
    }
  }

  return true;
}

bool rmat_hnf_check_transformation(const rmat_t* A, const rmat_t* H,
                                   const rmat_t* U)
{
  // column-operation convention: H = A * U; reject shape mismatches
  // before allocating the product
  if (A == NULL || H == NULL || U == NULL) return false;
  if (A->rows != H->rows || A->cols != U->rows || U->cols != H->cols) {
    return false;
  }

  rmat_t P = {0};
  if (rmat_init(&P, H->rows, H->cols) != RABIN_SUCCESS) return false;

  if (rmat_mul(&P, A, U) != RABIN_SUCCESS) {
    rmat_clear(&P);
    return false;
  }

  bool ok = rmat_equal(&P, H);

  rmat_clear(&P);
  return ok;
}

bool rmat_hnf_check_unimodular(const rmat_t* U)
{
  if (U == NULL) return false;
  if (!rmat_is_square(U)) {
    return false;
  }
  // the empty matrix acts as the identity transformation
  if (U->rows == 0) {
    return true;
  }

  rz_t det;
  rz_init(&det);

  if (rmat_det(&det, U) != RABIN_SUCCESS) {
    rz_clear(&det);
    return false;
  }

  // |det| == 1: clear the sign flag and use the O(1) one-check
  det.is_neg = false;
  bool ok = rz_is_one(&det);

  rz_clear(&det);
  return ok;
}

bool rmat_hnf_verify(const rmat_t* A, const rmat_t* H, const rmat_t* U)
{
  if (A == NULL || H == NULL || U == NULL) return false;

  // order matters: the cheap structure check runs first, the
  // multiplication and the determinant only if it passes
  return rmat_hnf_check_structure(H) &&
         rmat_hnf_check_transformation(A, H, U) && rmat_hnf_check_unimodular(U);
}

/* ------------------------------------------------------------------ */
/* LLL basis reduction                                                 */
/* ------------------------------------------------------------------ */
/*
 * rmat_lll()  implements Cohen's Algorithm 2.6.7 (all-integer LLL,
 *             full rank only).
 * rmat_mlll() handles zero / dependent columns as a two-phase rebuild:
 *             (1) a rank-revealing HNF column reduction (Cohen 2.4.5)
 *             tracking the unimodular transform U with
 *                 B_orig * U = [ 0_{m x (n-p)} | W ],
 *             so the first n-p columns of U are a Z-basis of the
 *             relation space and W (m x p) a Z-basis of the column
 *             lattice;
 *             (2) rmat_lll() on W, W <- W * Q, and
 *                 B <- [ 0 | W ],  H <- [ U_rel | U_lift * Q ].
 *
 * rmat_lll() maintains the current generator vectors V (initially a
 * copy of B) together with the unimodular transform H such that
 * V = B_orig * H holds at all times; on success B is overwritten with
 * V. It stores the scaled Gram-Schmidt data (1-based indices as in
 * Cohen; 0-based here)
 *
 *     lambda[i][j] = d[j+1] * mu[i+1][j+1]     (0 <= j < i < n)
 *     d[0] = 1,  d[i+1] = ||b_{i+1}||^2 * d[i]
 *
 * so that all quantities are integers. In rmat_lll(), RED(k, l)
 * size-reduces vector k against vector l and SWAP exchanges vectors
 * k-1 and k (0-based).
 */

/* scratch used by the RED/SWAP steps */
typedef struct {
  rz_t u, v, w;  // general temporaries
  rz_t one, two;
} lll_scr_t;

static rabin_err_t lll_scr_init(lll_scr_t* s)
{
  if (rz_init_multi(&s->u, &s->v, &s->w, &s->one, &s->two, NULL) !=
      RABIN_SUCCESS)
    return RABIN_ERR_OUT_OF_MEMORY;

  rabin_err_t err = rz_set_u64(&s->one, 1);
  if (err == RABIN_SUCCESS) err = rz_set_u64(&s->two, 2);
  return err;
}

static void lll_scr_clear(lll_scr_t* s)
{
  rz_clear_multi(&s->u, &s->v, &s->w, &s->one, &s->two, NULL);
}

/*
 * delta as an exact binary rational. Any double delta in (0, 1] is one:
 * frexp() gives delta = m * 2^e with m in [1/2, 1), so delta = M / 2^E
 * with M = m * 2^53 an integer in [2^52, 2^53) and E = 53 - e >= 1.
 * Stripping the common factors of two leaves num / 2^den_bits with
 * num < 2^53 (and num odd unless the fraction is integral).
 */
typedef struct {
  uint64_t num;
  u64 den_bits;
} lll_delta_t;

static lll_delta_t lll_delta_split(double delta)
{
  int e = 0;
  const double m = frexp(delta, &e);
  u64 M = (u64)(m * 9007199254740992.0); /* 2^53, exact */
  u64 den_bits = (u64)(53 - e);
  u64 s = (u64)__builtin_ctzll(M);
  if (s > den_bits) s = den_bits;

  lll_delta_t d;
  d.num = M >> s;
  d.den_bits = den_bits - s;
  return d;
}

/*
 * Lovasz test for the current k (1-based, lambda = lambda[k][k-1]):
 *
 *   swap iff d[k+1]*d[k-1] < delta*d[k]^2 - lambda^2.
 *
 * With delta = num / 2^den_bits this is evaluated as the exact integer
 * inequality
 *
 *   2^den_bits * d[k+1]*d[k-1] < num * d[k]^2 - 2^den_bits * lambda^2,
 *
 * so the reduction result is independent of any floating-point rounding.
 */
static rabin_err_t lll_lovasz(rmat_t* lam, const rz_t* dv, u64 k,
                              const lll_delta_t* dl, const rz_t* dden,
                              bool* swap, lll_scr_t* s)
{
  // lhs = d[k+1]*d[k-1]  (dv[k], dv[k-2])
  rabin_err_t err = rz_mul(&s->u, &dv[k], &dv[k - 2]);
  if (err != RABIN_SUCCESS) return err;
  // M = d[k]^2 (dv[k-1]), N = lambda^2
  if ((err = rz_sqr(&s->v, &dv[k - 1])) != RABIN_SUCCESS) return err;
  if ((err = rz_sqr(&s->w, RMAT_GET(lam, k - 1, k - 2))) != RABIN_SUCCESS)
    return err;

  // M <- num * M  (num < 2^53 always fits a positive i64)
  if ((err = rz_mul_i64(&s->v, &s->v, (i64)dl->num)) != RABIN_SUCCESS)
    return err;

  // lhs and N pick up the factor 2^den_bits
  if (dl->den_bits != 0) {
    if (dl->den_bits < 63) {
      const i64 den = (i64)1 << dl->den_bits;
      if ((err = rz_mul_i64(&s->w, &s->w, den)) != RABIN_SUCCESS) return err;
      if ((err = rz_mul_i64(&s->u, &s->u, den)) != RABIN_SUCCESS) return err;
    } else {
      if ((err = rz_mul(&s->w, &s->w, dden)) != RABIN_SUCCESS) return err;
      if ((err = rz_mul(&s->u, &s->u, dden)) != RABIN_SUCCESS) return err;
    }
  }

  if ((err = rz_sub(&s->v, &s->v, &s->w)) != RABIN_SUCCESS) return err;
  *swap = (rz_cmp(&s->u, &s->v) < 0);
  return RABIN_SUCCESS;
}

/* r = (column a of V) . (column b of V) */
static rabin_err_t lll_dot(rz_t* r, const rmat_t* V, u64 a, u64 b, rz_t* t)
{
  rabin_err_t err = rz_set_u64(r, 0);
  if (err != RABIN_SUCCESS) return err;

  for (u64 i = 0; i < V->rows; i++) {
    if ((err = rz_mul(t, RMAT_GET(V, i, a), RMAT_GET(V, i, b))) !=
        RABIN_SUCCESS)
      return err;
    if ((err = rz_add(r, r, t)) != RABIN_SUCCESS) return err;
  }
  return RABIN_SUCCESS;
}

/* ---------------------- LLL (Cohen 2.6.7) ------------------------- */

/* RED1: size-reduce vector k against vector l (0-based, l < k). */
static rabin_err_t lll_red1(rmat_t* V, rmat_t* H, rmat_t* lam, const rz_t* dv,
                            u64 k, u64 l, lll_scr_t* s)
{
  rz_t* lk = RMAT_GET(lam, k, l);

  // stop if |2 * lambda[k][l]| <= d[l+1]
  rabin_err_t err = rz_mul(&s->u, &s->two, lk);
  if (err != RABIN_SUCCESS) return err;
  if (rz_cmp_abs(&s->u, &dv[l + 1]) <= 0) return RABIN_SUCCESS;

  // q = nearest(lambda[k][l] / d[l+1]) =
  //     floor((2*lambda[k][l] + d[l+1]) / (2*d[l+1]))
  err = rz_add(&s->v, &s->u, &dv[l + 1]);
  if (err == RABIN_SUCCESS) err = rz_mul(&s->u, &s->two, &dv[l + 1]);
  if (err == RABIN_SUCCESS) err = rz_div_euclid(&s->w, &s->v, &s->u);
  if (err != RABIN_SUCCESS) return err;
  rz_t* q = &s->w;

  // V_k -= q * V_l
  for (u64 i = 0; i < V->rows; i++) {
    if ((err = rz_mul(&s->u, q, RMAT_GET(V, i, l))) != RABIN_SUCCESS)
      return err;
    if ((err = rz_sub(RMAT_GET(V, i, k), RMAT_GET(V, i, k), &s->u)) !=
        RABIN_SUCCESS)
      return err;
  }
  // H[:, k] -= q * H[:, l]
  for (u64 i = 0; i < H->rows; i++) {
    if ((err = rz_mul(&s->u, q, RMAT_GET(H, i, l))) != RABIN_SUCCESS)
      return err;
    if ((err = rz_sub(RMAT_GET(H, i, k), RMAT_GET(H, i, k), &s->u)) !=
        RABIN_SUCCESS)
      return err;
  }

  // lambda[k][l] -= q * d[l+1];  lambda[k][i] -= q * lambda[l][i]  (i < l)
  if ((err = rz_mul(&s->u, q, &dv[l + 1])) != RABIN_SUCCESS) return err;
  if ((err = rz_sub(lk, lk, &s->u)) != RABIN_SUCCESS) return err;
  for (u64 i = 0; i < l; i++) {
    if ((err = rz_mul(&s->u, q, RMAT_GET(lam, l, i))) != RABIN_SUCCESS)
      return err;
    if ((err = rz_sub(RMAT_GET(lam, k, i), RMAT_GET(lam, k, i), &s->u)) !=
        RABIN_SUCCESS)
      return err;
  }
  return RABIN_SUCCESS;
}

/* SWAP1: swap vectors k-1 and k (0-based, k >= 1). */
static rabin_err_t lll_swap1(rmat_t* V, rmat_t* H, rmat_t* lam, rz_t* dv, u64 k,
                             u64 kmax, lll_scr_t* s)
{
  rabin_err_t err;

  // swap V columns k and k-1
  for (u64 i = 0; i < V->rows; i++) {
    if ((err = rz_copy(&s->u, RMAT_GET(V, i, k))) != RABIN_SUCCESS) return err;
    if ((err = rz_copy(RMAT_GET(V, i, k), RMAT_GET(V, i, k - 1))) !=
        RABIN_SUCCESS)
      return err;
    if ((err = rz_copy(RMAT_GET(V, i, k - 1), &s->u)) != RABIN_SUCCESS)
      return err;
  }
  // swap H columns k and k-1
  for (u64 i = 0; i < H->rows; i++) {
    if ((err = rz_copy(&s->u, RMAT_GET(H, i, k))) != RABIN_SUCCESS) return err;
    if ((err = rz_copy(RMAT_GET(H, i, k), RMAT_GET(H, i, k - 1))) !=
        RABIN_SUCCESS)
      return err;
    if ((err = rz_copy(RMAT_GET(H, i, k - 1), &s->u)) != RABIN_SUCCESS)
      return err;
  }
  // swap lambda[k][j] <-> lambda[k-1][j]  (j < k-1)
  for (u64 j = 0; j + 1 < k; j++) {
    if ((err = rz_copy(&s->u, RMAT_GET(lam, k, j))) != RABIN_SUCCESS)
      return err;
    if ((err = rz_copy(RMAT_GET(lam, k, j), RMAT_GET(lam, k - 1, j))) !=
        RABIN_SUCCESS)
      return err;
    if ((err = rz_copy(RMAT_GET(lam, k - 1, j), &s->u)) != RABIN_SUCCESS)
      return err;
  }

  // A <- lambda[k][k-1]  (= lambda_{s+1,s} with s = k, 1-based)
  rz_t* A = RMAT_GET(lam, k, k - 1);

  // for i > k (up to kmax), 1-based I = i+1 > s+1:
  //   t <- lambda[i][k]   (old lambda_{I,s+1})
  //   lambda[i][k]   <- (d[k+1]*lambda[i][k-1] - A*t) / d[k] (lambda'_{I,s+1})
  //   lambda[i][k-1] <- (A*lambda[i][k-1] + d[k-1]*t) / d[k]   (lambda'_{I,s})
  for (u64 i = k + 1; i <= kmax; i++) {
    if ((err = rz_copy(&s->w, RMAT_GET(lam, i, k))) != RABIN_SUCCESS)
      return err;
    rz_t* t = &s->w;
    if ((err = rz_mul(&s->u, &dv[k + 1], RMAT_GET(lam, i, k - 1))) !=
        RABIN_SUCCESS)
      return err;
    if ((err = rz_mul(&s->v, A, t)) != RABIN_SUCCESS) return err;
    if ((err = rz_sub(&s->u, &s->u, &s->v)) != RABIN_SUCCESS) return err;
    if ((err = rz_div_exact(RMAT_GET(lam, i, k), &s->u, &dv[k])) !=
        RABIN_SUCCESS)
      return err;
    if ((err = rz_mul(&s->v, A, RMAT_GET(lam, i, k - 1))) != RABIN_SUCCESS)
      return err;
    if ((err = rz_mul(&s->u, &dv[k - 1], t)) != RABIN_SUCCESS) return err;
    if ((err = rz_add(&s->v, &s->v, &s->u)) != RABIN_SUCCESS) return err;
    if ((err = rz_div_exact(RMAT_GET(lam, i, k - 1), &s->v, &dv[k])) !=
        RABIN_SUCCESS)
      return err;
  }

  // d[k] <- (d[k+1]*d[k-1] + A^2) / d[k]   (new d_{s+1}, 1-based)
  err = rz_mul(&s->u, &dv[k + 1], &dv[k - 1]);
  if (err == RABIN_SUCCESS) err = rz_sqr(&s->v, A);
  if (err == RABIN_SUCCESS) err = rz_add(&s->u, &s->u, &s->v);
  if (err == RABIN_SUCCESS) err = rz_div_exact(&s->v, &s->u, &dv[k]);
  if (err != RABIN_SUCCESS) return err;
  return rz_copy(&dv[k], &s->v);
}

/*
 * Cohen 2.6.7 main loop, with 1-based k and kmax exactly as in Cohen
 * (0-based array indices are k-1 etc.). Requires full rank; returns
 * RABIN_ERR_INVALID_ARG when a d[k] becomes zero.
 */
static rabin_err_t rmat_lll_run(rmat_t* V, rmat_t* H, rmat_t* lam, rz_t* dv,
                                u64 n, double delta)
{
  lll_scr_t s;
  if (lll_scr_init(&s) != RABIN_SUCCESS) return RABIN_ERR_OUT_OF_MEMORY;

  rz_t dden;
  if (rz_init(&dden) != RABIN_SUCCESS) {
    lll_scr_clear(&s);
    return RABIN_ERR_OUT_OF_MEMORY;
  }

  // delta = num / 2^den_bits (exact); dden = 2^den_bits only needed when
  // the denominator no longer fits an i64 multiplier
  lll_delta_t dl = lll_delta_split(delta);
  rabin_err_t err = RABIN_SUCCESS;
  if (dl.den_bits >= 63) {
    if (rz_set_u64(&dden, 1) == RABIN_SUCCESS)
      err = rz_lshift(&dden, &dden, (int)dl.den_bits);
    if (err != RABIN_SUCCESS) goto cleanup;
  }

  TRY(rz_set_u64(&dv[0], 1));

  u64 k = 1;
  u64 kmax = 0;

step2:
  // [Initialize GS data for vector k]
  if (k > kmax) {
    kmax = k;
    for (u64 j = 1; j <= k; j++) {
      TRY(lll_dot(&s.u, V, k - 1, j - 1, &s.v));
      // U <- (d[i]*U - lambda[k][i]*lambda[j][i]) / d[i-1]   (1-based i < j,
      // d[1] = 1, d[i] stored in dv[i-1] (0-based))
      for (u64 i = 1; i < j; i++) {
        TRY(rz_mul(&s.w, &dv[i], &s.u));
        TRY(rz_mul(&s.v, RMAT_GET(lam, k - 1, i - 1),
                   RMAT_GET(lam, j - 1, i - 1)));
        TRY(rz_sub(&s.w, &s.w, &s.v));
        TRY(rz_div_exact(&s.u, &s.w, &dv[i - 1]));
      }
      if (j < k) {
        TRY(rz_copy(RMAT_GET(lam, k - 1, j - 1), &s.u));
      } else {
        TRY(rz_copy(&dv[k], &s.u));
        if (rz_is_zero(&dv[k])) {
          err = RABIN_ERR_INVALID_ARG;
          goto cleanup;
        }
      }
    }
  }

  // [LLL test: size-reduce k against k-1, check the Lovasz condition]
  if (k > 1) {
    TRY(lll_red1(V, H, lam, dv, k - 1, k - 2, &s));

    // Lovasz condition (1-based k, lambda = lambda[k][k-1], d[i] in dv[i-1]),
    // evaluated exactly for every double delta (see lll_lovasz):
    //   swap iff d[k+1]*d[k-1] < delta*d[k]^2 - lambda^2
    // (i.e. ||b*_k||^2 < (delta - mu^2) * ||b*_{k-1}||^2 with mu =
    // lambda/d[k]).
    bool swap = false;
    TRY(lll_lovasz(lam, dv, k, &dl, &dden, &swap, &s));

    if (swap) {
      TRY(lll_swap1(V, H, lam, dv, k - 1, kmax - 1, &s));
      k--;
      // the swapped-in vector now sits at k and must be re-examined
      // (size-reduce + Lovasz) before advancing: restart at step2.
      goto step2;
    }
  }

  // descending: RED1(k,l) updates lambda[k][i] for i < l, so l must decrease
  for (u64 l = k - 1; l >= 1; l--) {
    TRY(lll_red1(V, H, lam, dv, k - 1, l - 1, &s));
  }
  k++;
  if (k <= n) goto step2;

  err = RABIN_SUCCESS;
cleanup:
  rz_clear(&dden);
  lll_scr_clear(&s);
  return err;
}

rabin_err_t rmat_lll(rmat_t* B, u64 n, double delta, rmat_t* H)
{
  if (B == NULL || H == NULL) return RABIN_ERR_NULL_PTR;
  if (!(delta > 0.0 && delta <= 1.0)) return RABIN_ERR_INVALID_ARG;
  if (n != B->cols) return RABIN_ERR_INVALID_ARG;

  // nothing to reduce; H is re-initialized to the n x n identity
  rabin_err_t err = rmat_init(H, n, n);
  if (err != RABIN_SUCCESS) return err;
  err = rmat_id(H, n);
  if (err != RABIN_SUCCESS) return err;
  if (n <= 1) return RABIN_SUCCESS;

  rmat_t V = {0}, lam = {0};
  rz_t* dv = NULL;

  if (rmat_init(&V, B->rows, n) != RABIN_SUCCESS)
    return RABIN_ERR_OUT_OF_MEMORY;
  if (rmat_init(&lam, n, n) != RABIN_SUCCESS) {
    rmat_clear(&V);
    return RABIN_ERR_OUT_OF_MEMORY;
  }
  dv = malloc((n + 1) * sizeof(rz_t));
  if (dv == NULL) {
    rmat_clear(&V);
    rmat_clear(&lam);
    return RABIN_ERR_OUT_OF_MEMORY;
  }
  for (u64 i = 0; i <= n; i++) {
    if (rz_init(&dv[i]) != RABIN_SUCCESS) {
      for (u64 j = 0; j < i; j++) rz_clear(&dv[j]);
      free(dv);
      rmat_clear(&V);
      rmat_clear(&lam);
      return RABIN_ERR_OUT_OF_MEMORY;
    }
  }

  // H is already the identity (set above); on failure B is left untouched
  err = rmat_copy(&V, B);
  if (err == RABIN_SUCCESS) err = rmat_lll_run(&V, H, &lam, dv, n, delta);
  if (err == RABIN_SUCCESS) err = rmat_copy(B, &V);

  for (u64 i = 0; i <= n; i++) rz_clear(&dv[i]);
  free(dv);
  rmat_clear(&V);
  rmat_clear(&lam);
  return err;
}

/* ---------------------- MLLL (Cohen 2.6.8) ------------------------ */

/*
 * HNF via the exact column reduction (Cohen's Algorithm 2.4.5, same
 * loop as rmat_hermite_gcd()), additionally tracking the unimodular
 * column transformation.
 *
 * The reduction drives (n - p) columns to zero and drops them from the
 * left (p = rank of A), so on success:
 *
 *   W = A_work[:, k_final..n-1]   (m x p, HNF of the column lattice)
 *   A * U = [ 0_{m x (n-p)} | W ]
 *
 * Hence the first n - p columns of U form a Z-basis of the relation
 * space { x in Z^n : A x = 0 }, and the last p columns of U are
 * integer lifts of the HNF basis. W and U must point to empty
 * (zero-initialized) rmat_t objects; on success they hold the result.
 */
static rabin_err_t hnf_with_transform(rmat_t* W, rmat_t* U, const rmat_t* A)
{
  if (W == NULL || U == NULL || A == NULL) return RABIN_ERR_NULL_PTR;

  u64 m = A->rows;
  u64 n = A->cols;
  if (m == 0 || n == 0) return RABIN_ERR_INVALID_ARG;

  u64 i = m - 1;
  u64 j = n - 1;
  u64 k = n - 1;
  u64 l = 0;
  if (m > n) {
    l = m - n;
  }

  bool u_live = false;
  rmat_t A_work = {0};
  rz_t u, v, d, temp, temp2, b, q_k, q_j;
  rz_t* B = NULL;
  rz_t* C = NULL;
  rabin_err_t err = RABIN_SUCCESS;

  TRY(rmat_init(&A_work, m, n));
  TRY(rmat_copy(&A_work, A));

  TRY(rmat_init(U, n, n));
  u_live = true;
  TRY(rmat_id(U, n));

  rz_init_multi(&u, &v, &d, &temp, &temp2, &b, &q_k, &q_j, NULL);

  B = malloc(m * sizeof(rz_t));
  if (B == NULL) {
    err = RABIN_ERR_OUT_OF_MEMORY;
    goto cleanup;
  }
  for (u64 x = 0; x < m; x++) {
    rz_init(&B[x]);
  }

  C = malloc(n * sizeof(rz_t));
  if (C == NULL) {
    err = RABIN_ERR_OUT_OF_MEMORY;
    goto cleanup;
  }
  for (u64 x = 0; x < n; x++) {
    rz_init(&C[x]);
  }

// 2. [Check zero]
step2:
  if (j == 0) {
    goto step4;
  }
  j--;
  if (rz_is_zero(RMAT_GET(&A_work, i, j))) {
    goto step2;
  }

  // 3. [Euclidean step] (column ops applied to A_work and U)
  TRY(hnf_euclid_step(&A_work, U, B, C, i, j, k, &u, &v, &d, &q_k, &q_j, &temp,
                      &temp2));
  goto step2;

  // 4. [Final reduction]
step4:
  TRY(rz_copy(&b, RMAT_GET(&A_work, i, k)));
  if (b.is_neg) {  // zero is never negative, no extra rz_is_zero test needed
    for (u64 x = 0; x < m; x++) {
      TRY(rz_neg(RMAT_GET(&A_work, x, k), RMAT_GET(&A_work, x, k)));
    }
    for (u64 x = 0; x < n; x++) {
      TRY(rz_neg(RMAT_GET(U, x, k), RMAT_GET(U, x, k)));
    }
    TRY(rz_neg(&b, &b));
  }

  if (rz_is_zero(&b)) {
    k++;
    goto step5;
  }

  for (u64 j_0 = k + 1; j_0 < n; j_0++) {
    TRY(rz_div_euclid(&temp2, RMAT_GET(&A_work, i, j_0), &b));
    for (u64 x = 0; x < m; x++) {
      TRY(rz_mul(&temp, &temp2, RMAT_GET(&A_work, x, k)));
      TRY(rz_sub(RMAT_GET(&A_work, x, j_0), RMAT_GET(&A_work, x, j_0), &temp));
    }
    for (u64 x = 0; x < n; x++) {
      TRY(rz_mul(&temp, &temp2, RMAT_GET(U, x, k)));
      TRY(rz_sub(RMAT_GET(U, x, j_0), RMAT_GET(U, x, j_0), &temp));
    }
  }

step5:
  // 5. [Finished]
  if (i == l) {
    TRY(rmat_init(W, m, n - k));
    for (u64 x = 0; x < n - k; x++) {
      for (u64 y = 0; y < m; y++) {
        TRY(rz_copy(RMAT_GET(W, y, x), RMAT_GET(&A_work, y, x + k)));
      }
    }
    u_live = false;  // keep U for the caller
    goto cleanup;
  }

  i--;
  k--;
  j = k;
  goto step2;

cleanup:
  rz_clear_multi(&u, &v, &d, &temp, &temp2, &b, &q_k, &q_j, NULL);
  if (C != NULL) {
    for (u64 x = 0; x < n; x++) {
      rz_clear(&C[x]);
    }
    free(C);
  }
  if (B != NULL) {
    for (u64 x = 0; x < m; x++) {
      rz_clear(&B[x]);
    }
    free(B);
  }
  if (u_live) rmat_clear(U);
  rmat_clear(&A_work);
  return err;
}

rabin_err_t rmat_mlll(rmat_t* B, u64 n, double delta, rmat_t* H, u64* rank)
{
  if (B == NULL || H == NULL || rank == NULL) return RABIN_ERR_NULL_PTR;
  if (!(delta > 0.0 && delta <= 1.0)) return RABIN_ERR_INVALID_ARG;
  if (n != B->cols) return RABIN_ERR_INVALID_ARG;

  // the HNF machinery (and Cohen 2.6.8) is stated for d <= n: the
  // column reduction only inspects the bottom n rows
  if (B->rows > n) return RABIN_ERR_INVALID_ARG;

  *rank = 0;

  // no coordinates (or no vectors): every column is the zero vector
  if (n == 0 || B->rows == 0) {
    rabin_err_t e = rmat_init(H, n, n);
    if (e != RABIN_SUCCESS) return e;
    return rmat_id(H, n);
  }

  rmat_t A = {0};  // original generator matrix (kept for B_orig * H checks)
  rmat_t W = {0};  // HNF of the column lattice (m x p)
  rmat_t U = {0};  // unimodular: A * U = [ 0 | W ]
  rmat_t Q = {0};  // LLL transformation of W (p x p)
  rmat_t A_work = {0};
  rmat_t H_lift = {0};  // U[:, n-p..n-1] (n x p)
  rmat_t H_tail = {0};  // H_lift * Q (n x p)
  rmat_t B_new = {0};   // assembled result; swapped into B after the check

  rabin_err_t err = rmat_init(&A, B->rows, n);
  if (err != RABIN_SUCCESS) return err;
  err = rmat_copy(&A, B);
  if (err != RABIN_SUCCESS) goto cleanup;

  // 1. [Rank-revealing column reduction]
  err = hnf_with_transform(&W, &U, &A);
  if (err != RABIN_SUCCESS) goto cleanup;

  u64 p = W.cols;
  *rank = p;

  // 2. [LLL-reduce the lattice basis]
  if (p > 0) {
    TRY(rmat_lll(&W, p, delta, &Q));
  }

  // 3. [Assemble B = [ 0 | W ] and H = [ U_rel | U_lift * Q ]]
  err = rmat_init(H, n, n);
  if (err != RABIN_SUCCESS) goto cleanup;
  err = rmat_id(H, n);
  if (err != RABIN_SUCCESS) goto cleanup;
  for (u64 x = 0; x < n - p; x++) {
    for (u64 y = 0; y < n; y++) {
      TRY(rmat_set(H, RMAT_GET(&U, y, x), y, x));
    }
  }
  if (p > 0) {
    TRY(rmat_init(&H_lift, n, p));
    for (u64 x = 0; x < p; x++) {
      for (u64 y = 0; y < n; y++) {
        TRY(rmat_set(&H_lift, RMAT_GET(&U, y, n - p + x), y, x));
      }
    }
    TRY(rmat_init(&H_tail, n, p));
    TRY(rmat_mul(&H_tail, &H_lift, &Q));
    for (u64 x = 0; x < p; x++) {
      for (u64 y = 0; y < n; y++) {
        TRY(rmat_set(H, RMAT_GET(&H_tail, y, x), y, n - p + x));
      }
    }
  }
  // B_new = [ 0_{m x (n-p)} | W ] (the caller's B stays untouched until
  // the final commit, so a failure never corrupts it)
  TRY(rmat_init(&B_new, B->rows, n));
  for (u64 x = 0; x < p; x++) {
    for (u64 y = 0; y < B->rows; y++) {
      TRY(rmat_set(&B_new, RMAT_GET(&W, y, x), y, n - p + x));
    }
  }

  // 4. [Verify B_new == A_orig * H before committing]
  TRY(rmat_init(&A_work, B->rows, n));
  TRY(rmat_mul(&A_work, &A, H));
  if (!rmat_equal(&A_work, &B_new)) {
    err = RABIN_ERR_INVALID_ARG;
    goto cleanup;
  }

  // 5. [Commit: B now holds the reduced basis]
  rmat_swap(B, &B_new);

  err = RABIN_SUCCESS;
cleanup:
  rmat_clear(&B_new);
  rmat_clear(&A_work);
  rmat_clear(&H_tail);
  rmat_clear(&H_lift);
  rmat_clear(&Q);
  rmat_clear(&U);
  rmat_clear(&W);
  rmat_clear(&A);
  return err;
}

/* ------------- Kernel and image of a matrix (Cohen 2.7) ------------- */

/*
 * Extract columns [c0, c0 + len) of src into a fresh matrix dst.
 */
static rabin_err_t rmat_extract_cols(rmat_t* dst, const rmat_t* src, u64 c0,
                                     u64 len)
{
  if (dst == NULL || src == NULL) return RABIN_ERR_NULL_PTR;
  if (c0 + len > src->cols) return RABIN_ERR_INVALID_ARG;

  rabin_err_t err = rmat_init(dst, src->rows, len);
  if (err != RABIN_SUCCESS) return err;
  for (u64 x = 0; x < src->rows && err == RABIN_SUCCESS; x++) {
    for (u64 y = 0; y < len; y++) {
      err = rz_copy(RMAT_GET(dst, x, y), RMAT_GET(src, x, c0 + y));
    }
  }
  return err;
}

/*
 * Kernel-only LLL (Cohen 2.7.2): integer Gram-Schmidt with flags iota[]
 * marking the vectors that are nonzero in the Gram-Schmidt sense (iota = 1
 * iff B_k != 0). When u == 0 the convention d_k <- d_{k-1} keeps every d_k
 * nonzero, so no division ever hits zero and dependent vectors are simply
 * pushed to the front by SWAPK until iota[0..r-1] == 0 exactly delimits the
 * kernel: the first r columns of H then form a Z-basis of ker(A).
 */

/*
 * Exact integer Gram-Schmidt for one column (1-based row i of lam): fills
 * lambda[i][j] (j < i), d_i and iota[i-1] from the vectors in V and the GS
 * data of the rows below i. Mirrors the step-2 loop of rmat_lll_run.
 */
static rabin_err_t lll_gsrow(rmat_t* V, rmat_t* lam, rz_t* dv, bool* iota,
                             u64 i, lll_scr_t* s)
{
  rabin_err_t err;
  for (u64 j = 1; j <= i; j++) {
    if (!iota[j - 1] && j < i) {
      if ((err = rz_set_u64(RMAT_GET(lam, i - 1, j - 1), 0)) != RABIN_SUCCESS)
        return err;
      continue;
    }
    if ((err = lll_dot(&s->u, V, i - 1, j - 1, &s->v)) != RABIN_SUCCESS)
      return err;
    for (u64 x = 1; x < j; x++) {
      if (!iota[x - 1]) continue;
      if ((err = rz_mul(&s->w, &dv[x], &s->u)) != RABIN_SUCCESS) return err;
      if ((err = rz_mul(&s->v, RMAT_GET(lam, i - 1, x - 1),
                        RMAT_GET(lam, j - 1, x - 1))) != RABIN_SUCCESS)
        return err;
      if ((err = rz_sub(&s->w, &s->w, &s->v)) != RABIN_SUCCESS) return err;
      if ((err = rz_div_exact(&s->u, &s->w, &dv[x - 1])) != RABIN_SUCCESS)
        return err;
    }
    if (j < i) {
      if ((err = rz_copy(RMAT_GET(lam, i - 1, j - 1), &s->u)) != RABIN_SUCCESS)
        return err;
    } else {
      if (!rz_is_zero(&s->u)) {
        if ((err = rz_copy(&dv[i], &s->u)) != RABIN_SUCCESS) return err;
        iota[i - 1] = true;
      } else {
        if ((err = rz_copy(&dv[i], &dv[i - 1])) != RABIN_SUCCESS) return err;
        iota[i - 1] = false;
      }
    }
  }
  return RABIN_SUCCESS;
}

/*
 * SWAPK (Cohen 2.7.2): exchange the vectors at (1-based) positions k-1 and k.
 * On entry iota[k-2] is true and iota[k-1] is false (the vector at position
 * k is Gram-Schmidt-zero). Instead of updating lambda/d incrementally, the
 * GS data of rows k-1..kmax is recomputed from the (already swapped) vectors
 * with lll_gsrow(); this is exact by construction and keeps iota in sync.
 */
static rabin_err_t lll_swapk(rmat_t* V, rmat_t* H, rmat_t* lam, rz_t* dv,
                             bool* iota, u64 k, u64 kmax, lll_scr_t* s)
{
  rabin_err_t err;

  // swap the columns k-1 and k of V and H
  for (u64 i = 0; i < V->rows; i++) {
    if ((err = rz_copy(&s->u, RMAT_GET(V, i, k - 1))) != RABIN_SUCCESS)
      return err;
    if ((err = rz_copy(RMAT_GET(V, i, k - 1), RMAT_GET(V, i, k - 2))) !=
        RABIN_SUCCESS)
      return err;
    if ((err = rz_copy(RMAT_GET(V, i, k - 2), &s->u)) != RABIN_SUCCESS)
      return err;
  }
  for (u64 i = 0; i < H->rows; i++) {
    if ((err = rz_copy(&s->u, RMAT_GET(H, i, k - 1))) != RABIN_SUCCESS)
      return err;
    if ((err = rz_copy(RMAT_GET(H, i, k - 1), RMAT_GET(H, i, k - 2))) !=
        RABIN_SUCCESS)
      return err;
    if ((err = rz_copy(RMAT_GET(H, i, k - 2), &s->u)) != RABIN_SUCCESS)
      return err;
  }

  // recompute the GS data of all rows at or below the swap (ascending, so
  // each row sees the freshly updated data of the rows above it)
  for (u64 i = k - 1; i <= kmax; i++) {
    if ((err = lll_gsrow(V, lam, dv, iota, i, s)) != RABIN_SUCCESS) return err;
  }
  return RABIN_SUCCESS;
}

static rabin_err_t rmat_kernel_lll_run(rmat_t* V, rmat_t* H, rmat_t* lam,
                                       rz_t* dv, bool* iota, u64 n)
{
  lll_scr_t s;
  if (lll_scr_init(&s) != RABIN_SUCCESS) return RABIN_ERR_OUT_OF_MEMORY;

  rabin_err_t err = rz_set_u64(&dv[0], 1);  // d_0 = 1
  if (err != RABIN_SUCCESS) {
    lll_scr_clear(&s);
    return err;
  }

  u64 k = 1;
  u64 kmax = 0;

step2:
  // [Initialize Gram-Schmidt data for vector k]
  if (k > kmax) {
    kmax = k;
    TRY(lll_gsrow(V, lam, dv, iota, k, &s));
  }

step3:
  // [Test]
  if (k >= 2) {
    if (iota[k - 2]) {
      TRY(lll_red1(V, H, lam, dv, k - 1, k - 2, &s));
      if (!iota[k - 1]) {
        TRY(lll_swapk(V, H, lam, dv, iota, k, kmax, &s));
        if (k > 2) k--;  // k <- max(2, k-1)
        goto step3;
      }
    }
    // descending: RED(k,l) updates lambda[k][i] for i < l
    for (u64 l = k - 2; l >= 1; l--) {
      if (!iota[l - 1]) continue;
      TRY(lll_red1(V, H, lam, dv, k - 1, l - 1, &s));
    }
  }
  k++;
  if (k <= n) goto step2;

  err = RABIN_SUCCESS;
cleanup:
  lll_scr_clear(&s);
  return err;
}

/*
 * Nearest integer of the rational x (ties go up): with den > 0 the invariant
 * of rq_t,
 *
 *   round(num/den) = floor((2*num + den) / (2*den)),
 *
 * evaluated with the exact Euclidean quotient. t1, t2 are scratch.
 */
static rabin_err_t rz_round_nearest(rz_t* q, const rq_t* x, rz_t* t1, rz_t* t2)
{
  rabin_err_t err = rz_lshift(t1, &x->num, 1);  // 2*num
  if (err == RABIN_SUCCESS) err = rz_add(t1, t1, &x->den);
  if (err == RABIN_SUCCESS) err = rz_lshift(t2, &x->den, 1);  // 2*den > 0
  if (err == RABIN_SUCCESS) err = rz_div_euclid(q, t1, t2);
  return err;
}

/*
 * Invert the r x r rational matrix G by Gauss-Jordan elimination, returning
 * the inverse as a freshly allocated row-major array of r*r rationals
 * (freed by the caller with rq_clear + free). Fails with
 * RABIN_ERR_INVALID_ARG when G is singular.
 */
static rabin_err_t rq_mat_inverse(rq_t** Dout, const rmat_t* G, u64 r)
{
  if (Dout == NULL || G == NULL) return RABIN_ERR_NULL_PTR;
  *Dout = NULL;
  if (G->rows != r || G->cols != r) return RABIN_ERR_INVALID_ARG;
  if (r == 0) return RABIN_SUCCESS;

  rabin_err_t err = RABIN_SUCCESS;
  rq_t* A = calloc(r * 2 * r, sizeof(rq_t));
  if (A == NULL) return RABIN_ERR_OUT_OF_MEMORY;
  for (u64 x = 0; x < r * 2 * r; x++) {
    rq_init(&A[x]);
  }

  rq_t pv, f, t;
  rq_init_multi(&pv, &f, &t, NULL);

  rq_t* D = NULL;

  for (u64 x = 0; x < r && err == RABIN_SUCCESS; x++) {
    for (u64 y = 0; y < 2 * r; y++) {
      if (y < r) {
        err = rq_set(&A[x * 2 * r + y], RMAT_GET(G, x, y));
      } else {
        err = rq_set_u64(&A[x * 2 * r + y], y - r == x ? 1 : 0);
      }
    }
  }

  for (u64 c = 0; c < r && err == RABIN_SUCCESS; c++) {
    u64 p = c;
    while (p < r && rz_is_zero(&A[p * 2 * r + c].num)) p++;
    if (p == r) {
      err = RABIN_ERR_INVALID_ARG;  // singular
      break;
    }
    if (p != c) {
      for (u64 y = 0; y < 2 * r; y++) {
        err = rq_copy(&t, &A[p * 2 * r + y]);
        if (err == RABIN_SUCCESS)
          err = rq_copy(&A[p * 2 * r + y], &A[c * 2 * r + y]);
        if (err == RABIN_SUCCESS) err = rq_copy(&A[c * 2 * r + y], &t);
        if (err != RABIN_SUCCESS) break;
      }
    }
    if (err != RABIN_SUCCESS) break;

    // normalize row c so that the pivot becomes 1
    err = rq_copy(&pv, &A[c * 2 * r + c]);
    for (u64 y = c + 1; err == RABIN_SUCCESS && y < 2 * r; y++) {
      err = rq_div(&A[c * 2 * r + y], &A[c * 2 * r + y], &pv);
    }
    if (err == RABIN_SUCCESS) err = rq_set_u64(&A[c * 2 * r + c], 1);

    // clear the other rows
    for (u64 x = 0; x < r && err == RABIN_SUCCESS; x++) {
      if (x == c) continue;
      if (rz_is_zero(&A[x * 2 * r + c].num)) continue;
      err = rq_copy(&f, &A[x * 2 * r + c]);
      if (err == RABIN_SUCCESS) err = rq_set_u64(&A[x * 2 * r + c], 0);
      for (u64 y = c + 1; err == RABIN_SUCCESS && y < 2 * r; y++) {
        if (rz_is_zero(&A[c * 2 * r + y].num)) continue;
        err = rq_mul(&t, &f, &A[c * 2 * r + y]);
        if (err == RABIN_SUCCESS)
          err = rq_sub(&A[x * 2 * r + y], &A[x * 2 * r + y], &t);
      }
    }
  }
  if (err != RABIN_SUCCESS) goto cleanup;

  D = malloc(r * r * sizeof(rq_t));
  if (D == NULL) {
    err = RABIN_ERR_OUT_OF_MEMORY;
    goto cleanup;
  }
  for (u64 x = 0; x < r * r; x++) {
    rq_init(&D[x]);
  }
  for (u64 x = 0; x < r * r; x++) {
    if (rq_copy(&D[x], &A[(x / r) * 2 * r + r + (x % r)]) != RABIN_SUCCESS) {
      err = RABIN_ERR_OUT_OF_MEMORY;
      goto cleanup;
    }
  }
  *Dout = D;
  D = NULL;  // ownership transferred to the caller

cleanup:
  if (D != NULL) {
    for (u64 x = 0; x < r * r; x++) {
      rq_clear(&D[x]);
    }
    free(D);
  }
  rq_clear_multi(&pv, &f, &t, NULL);
  for (u64 x = 0; x < r * 2 * r; x++) {
    rq_clear(&A[x]);
  }
  free(A);
  return err;
}

/*
 * Steps 3-5 of Cohen 2.7.1: size-reduce every column of H with index
 * >= r against the (fixed) kernel basis formed by the first r columns.
 * H is updated in place.
 */
static rabin_err_t rmat_reduce_against_kernel(rmat_t* H, u64 n, u64 r)
{
  if (r == 0) return RABIN_SUCCESS;

  rabin_err_t err = RABIN_SUCCESS;
  rmat_t G = {0};
  rz_t* mr = NULL;
  rq_t* vq = NULL;
  rq_t* D = NULL;
  u64 mr_done = 0;
  rq_t acc, tmpq;
  rz_t t1, t2;
  bool temps = false;

  if ((err = rmat_init(&G, r, r)) != RABIN_SUCCESS) return err;
  mr = calloc(r, sizeof(rz_t));
  vq = calloc(r, sizeof(rq_t));
  if (mr == NULL || vq == NULL) {
    err = RABIN_ERR_OUT_OF_MEMORY;
    goto cleanup;
  }
  for (u64 j = 0; j < r; j++) {
    if (rz_init(&mr[j]) != RABIN_SUCCESS) goto cleanup;
    mr_done++;
  }
  for (u64 j = 0; j < r; j++) {
    if (rq_init(&vq[j]) != RABIN_SUCCESS) goto cleanup;
  }
  TRY(rq_init_multi(&acc, &tmpq, NULL));
  TRY(rz_init_multi(&t1, &t2, NULL));
  temps = true;

  // e = Gram matrix of the kernel basis columns
  for (u64 j = 0; j < r && err == RABIN_SUCCESS; j++) {
    for (u64 k = 0; k < r; k++) {
      err = lll_dot(RMAT_GET(&G, j, k), H, j, k, &t1);
    }
  }
  if (err == RABIN_SUCCESS) err = rq_mat_inverse(&D, &G, r);
  if (err != RABIN_SUCCESS) goto cleanup;

  for (u64 i = r; i < n && err == RABIN_SUCCESS; i++) {
    // v_j = H_i . H_j
    for (u64 j = 0; j < r; j++) {
      err = lll_dot(&t1, H, i, j, &t2);
      if (err != RABIN_SUCCESS) break;
      err = rq_set(&vq[j], &t1);
    }
    if (err != RABIN_SUCCESS) break;

    // m_j = nearest integer of x_j = sum_k D[j][k] * v_k
    for (u64 j = 0; j < r; j++) {
      err = rq_set_u64(&acc, 0);
      for (u64 k = 0; err == RABIN_SUCCESS && k < r; k++) {
        err = rq_mul(&tmpq, &D[j * r + k], &vq[k]);
        if (err == RABIN_SUCCESS) err = rq_add(&acc, &acc, &tmpq);
      }
      if (err != RABIN_SUCCESS) break;
      err = rz_round_nearest(&mr[j], &acc, &t1, &t2);
    }
    if (err != RABIN_SUCCESS) break;

    // H_i <- H_i - sum_j m_j * H_j
    for (u64 j = 0; j < r; j++) {
      if (rz_is_zero(&mr[j])) continue;
      for (u64 x = 0; x < n; x++) {
        err = rz_mul(&t1, &mr[j], RMAT_GET(H, x, j));
        if (err == RABIN_SUCCESS)
          err = rz_sub(RMAT_GET(H, x, i), RMAT_GET(H, x, i), &t1);
        if (err != RABIN_SUCCESS) break;
      }
      if (err != RABIN_SUCCESS) break;
    }
  }

cleanup:
  if (temps) {
    rq_clear_multi(&acc, &tmpq, NULL);
    rz_clear_multi(&t1, &t2, NULL);
  }
  if (D != NULL) {
    for (u64 x = 0; x < r * r; x++) rq_clear(&D[x]);
    free(D);
  }
  for (u64 j = 0; j < r && vq != NULL; j++) rq_clear(&vq[j]);
  free(vq);
  for (u64 j = 0; j < mr_done; j++) rz_clear(&mr[j]);
  free(mr);
  rmat_clear(&G);
  return err;
}

rabin_err_t rmat_kernel(const rmat_t* A, rmat_t* K)
{
  if (A == NULL || K == NULL) return RABIN_ERR_NULL_PTR;

  const u64 n = A->cols;
  if (n == 0) return rmat_init(K, 0, 0);

  rmat_t W = {0}, U = {0}, Uk = {0};
  rmat_t Kraw = {0};

  // no coordinates: every vector of Z^n is in the kernel
  if (A->rows == 0) {
    rabin_err_t err = rmat_init(K, n, n);
    if (err != RABIN_SUCCESS) return err;
    return rmat_id(K, n);
  }

  rabin_err_t err = hnf_with_transform(&W, &U, A);
  if (err != RABIN_SUCCESS) goto cleanup;

  const u64 r = n - W.cols;
  if (r == 0) {
    err = rmat_init(K, n, 0);
    goto cleanup;
  }

  err = rmat_extract_cols(&Kraw, &U, 0, r);
  if (err != RABIN_SUCCESS) goto cleanup;
  if (r > 1) {
    err = rmat_lll(&Kraw, r, 0.75, &Uk);
    if (err != RABIN_SUCCESS) goto cleanup;
  }
  rmat_swap(K, &Kraw);

cleanup:
  rmat_clear(&Kraw);
  rmat_clear(&Uk);
  rmat_clear(&U);
  rmat_clear(&W);
  return err;
}

rabin_err_t rmat_kernel_lll(const rmat_t* A, rmat_t* K)
{
  if (A == NULL || K == NULL) return RABIN_ERR_NULL_PTR;

  const u64 n = A->cols;
  if (n == 0) return rmat_init(K, 0, 0);

  rmat_t V = {0}, H = {0}, lam = {0}, Kraw = {0}, Uk = {0};
  rz_t* dv = NULL;
  bool* iota = NULL;

  // no coordinates: every vector of Z^n is in the kernel
  if (A->rows == 0) {
    rabin_err_t err = rmat_init(K, n, n);
    if (err != RABIN_SUCCESS) return err;
    return rmat_id(K, n);
  }

  rabin_err_t err = rmat_init(&V, A->rows, n);
  if (err != RABIN_SUCCESS) return err;
  err = rmat_copy(&V, A);
  if (err != RABIN_SUCCESS) goto cleanup;
  TRY(rmat_init(&H, n, n));
  TRY(rmat_id(&H, n));
  TRY(rmat_init(&lam, n, n));

  u64 dv_done = 0;
  dv = malloc((n + 1) * sizeof(rz_t));
  iota = calloc(n, sizeof(bool));
  if (dv == NULL || iota == NULL) {
    err = RABIN_ERR_OUT_OF_MEMORY;
    goto cleanup;
  }
  for (u64 x = 0; x <= n; x++) {
    if (rz_init(&dv[x]) != RABIN_SUCCESS) goto cleanup;
    dv_done++;
  }

  err = rmat_kernel_lll_run(&V, &H, &lam, dv, iota, n);
  if (err != RABIN_SUCCESS) goto cleanup;

  // r+1 is the first index with a nonzero flag: the first r columns of H
  // form a Z-basis of the kernel (they are independent)
  u64 r = 0;
  while (r < n && !iota[r]) r++;

  if (r == 0) {
    err = rmat_init(K, n, 0);
    goto cleanup;
  }
  err = rmat_extract_cols(&Kraw, &H, 0, r);
  if (err != RABIN_SUCCESS) goto cleanup;
  if (r > 1) {
    err = rmat_lll(&Kraw, r, 0.75, &Uk);
    if (err != RABIN_SUCCESS) goto cleanup;
  }
  rmat_swap(K, &Kraw);

cleanup:
  if (dv != NULL) {
    for (u64 x = 0; x < dv_done; x++) rz_clear(&dv[x]);
    free(dv);
  }
  free(iota);
  rmat_clear(&Kraw);
  rmat_clear(&Uk);
  rmat_clear(&lam);
  rmat_clear(&H);
  rmat_clear(&V);
  return err;
}

rabin_err_t rmat_kernel_image(const rmat_t* A, rmat_t* H, u64* p)
{
  if (A == NULL || H == NULL || p == NULL) return RABIN_ERR_NULL_PTR;
  *p = 0;

  const u64 n = A->cols;
  if (n == 0) return rmat_init(H, 0, 0);

  rmat_t W = {0}, U = {0}, Q = {0}, Uk = {0}, Kraw = {0}, Usub = {0};
  rmat_t Tail = {0}, Htmp = {0};
  rz_t det, one;
  bool det_live = false;

  // no coordinates: the whole Z^n is the kernel, the image is trivial
  if (A->rows == 0) {
    rabin_err_t err = rmat_init(H, n, n);
    if (err != RABIN_SUCCESS) return err;
    return rmat_id(H, n);
  }

  rabin_err_t err = hnf_with_transform(&W, &U, A);
  if (err != RABIN_SUCCESS) goto cleanup;
  const u64 rank = W.cols;
  const u64 r = n - rank;

  // LLL-reduce the HNF image basis W (in place), keeping the transform Q
  if (rank >= 1 && (err = rmat_lll(&W, rank, 0.75, &Q)) != RABIN_SUCCESS)
    goto cleanup;

  // LLL-reduce the kernel basis (columns 1..r of U, independent)
  if (r > 0) {
    err = rmat_extract_cols(&Kraw, &U, 0, r);
    if (err != RABIN_SUCCESS) goto cleanup;
    if (r > 1) err = rmat_lll(&Kraw, r, 0.75, &Uk);
    if (err != RABIN_SUCCESS) goto cleanup;
  }

  // image part: columns r..n-1 of U times Q
  if (rank >= 1) {
    err = rmat_extract_cols(&Usub, &U, r, rank);
    if (err == RABIN_SUCCESS) err = rmat_init(&Tail, n, rank);
    if (err == RABIN_SUCCESS) err = rmat_mul(&Tail, &Usub, &Q);
    if (err != RABIN_SUCCESS) goto cleanup;
  }

  // assemble H = [ K | U_tail * Q ] and commit only when every check passed
  TRY(rmat_init(&Htmp, n, n));
  for (u64 y = 0; y < n && err == RABIN_SUCCESS; y++) {
    for (u64 x = 0; x < r; x++) {
      err = rz_copy(RMAT_GET(&Htmp, y, x), RMAT_GET(&Kraw, y, x));
    }
    for (u64 x = 0; x < rank; x++) {
      if (err == RABIN_SUCCESS)
        err = rz_copy(RMAT_GET(&Htmp, y, r + x), RMAT_GET(&Tail, y, x));
    }
  }
  if (err != RABIN_SUCCESS) goto cleanup;

  // H must be unimodular
  TRY(rz_init_multi(&det, &one, NULL));
  det_live = true;
  TRY(rz_set_u64(&one, 1));
  TRY(rmat_det(&det, &Htmp));
  if (rz_cmp_abs(&det, &one) != 0) {
    err = RABIN_ERR_INVALID_ARG;
    goto cleanup;
  }

  // steps 3-5: keep the image columns small against the kernel basis
  TRY(rmat_reduce_against_kernel(&Htmp, n, r));

  rmat_swap(H, &Htmp);
  *p = rank;

cleanup:
  if (det_live) {
    rz_clear(&det);
    rz_clear(&one);
  }
  rmat_clear(&Htmp);
  rmat_clear(&Tail);
  rmat_clear(&Usub);
  rmat_clear(&Kraw);
  rmat_clear(&Uk);
  rmat_clear(&Q);
  rmat_clear(&U);
  rmat_clear(&W);
  return err;
}

/* ----------------- Dependence relations (Cohen 2.7.4) -------------- */

/* round scale*z to an integer, rejecting values outside the int64 range */
static rabin_err_t rz_round_scaled(rz_t* q, double z, double scale)
{
  if (!isfinite(z) || !(scale > 0.0) || !isfinite(scale) ||
      fabs(z) > 4.0e18 / scale)
    return RABIN_ERR_INVALID_ARG;
  return rz_set_i64(q, llround(scale * z));
}

/* best-effort double conversion of an integer, used for error residuals */
static double rz_as_double(const rz_t* z)
{
  double d = 0.0;
  for (size_t i = z->size; i > 0; i--) {
    d = d * 18446744073709551616.0 + (double)z->limbs[i - 1];
  }
  return z->is_neg ? -d : d;
}

rabin_err_t rmat_linear_dependence(const double* re, const double* im, u64 n,
                                   double scale, rmat_t* REL, double* err)
{
  if (re == NULL || REL == NULL) return RABIN_ERR_NULL_PTR;
  if (n < 2) return RABIN_ERR_INVALID_ARG;
  if (!(scale > 0.0) || !isfinite(scale)) return RABIN_ERR_INVALID_ARG;

  const bool complex_in = (im != NULL);
  for (u64 i = 0; i < n; i++) {
    if (!isfinite(re[i]) || (complex_in && !isfinite(im[i])))
      return RABIN_ERR_INVALID_ARG;
  }

  const u64 m = n + (complex_in ? 2 : 1);
  rmat_t M = {0};
  rabin_err_t e = rmat_init(&M, m, n);
  if (e != RABIN_SUCCESS) return e;

  // symmetric construction: rows 1..n hold the identity, the last rows the
  // scaled real (and imaginary) parts
  for (u64 i = 0; i < n && e == RABIN_SUCCESS; i++) {
    for (u64 j = 0; j < n; j++) {
      e = rz_set_u64(RMAT_GET(&M, i, j), i == j ? 1 : 0);
    }
  }
  for (u64 j = 0; j < n && e == RABIN_SUCCESS; j++) {
    e = rz_round_scaled(RMAT_GET(&M, n, j), re[j], scale);
  }
  if (e == RABIN_SUCCESS && complex_in) {
    for (u64 j = 0; j < n; j++) {
      e = rz_round_scaled(RMAT_GET(&M, n + 1, j), im[j], scale);
      if (e != RABIN_SUCCESS) break;
    }
  }

  // rmat_lll reduces M in place and returns the transformation matrix
  rmat_t Q = {0};
  if (e == RABIN_SUCCESS) e = rmat_lll(&M, n, 0.75, &Q);
  if (e == RABIN_SUCCESS) {
    e = rmat_init(REL, m, n);
    if (e == RABIN_SUCCESS) e = rmat_copy(REL, &M);
  }
  rmat_clear(&Q);
  if (e != RABIN_SUCCESS) {
    rmat_clear(&M);
    return e;
  }

  if (err != NULL) {
    // residual of every candidate relation: |sum_j REL[j][c] * z_j|
    for (u64 c = 0; c < n; c++) {
      long double sr = 0.0L, si = 0.0L;
      for (u64 j = 0; j < n; j++) {
        const double cj = rz_as_double(RMAT_GET(REL, j, c));
        sr += (long double)cj * (long double)re[j];
        if (complex_in) si += (long double)cj * (long double)im[j];
      }
      err[c] = (double)sqrtl(sr * sr + si * si);
    }
  }

  rmat_clear(&M);
  return RABIN_SUCCESS;
}

rabin_err_t rmat_algebraic_relation(double re, double im, u64 degree,
                                    double scale, rmat_t* poly, double* err)
{
  if (poly == NULL) return RABIN_ERR_NULL_PTR;
  if (degree < 1 || degree > 64) return RABIN_ERR_INVALID_ARG;
  if (!isfinite(re) || !isfinite(im) || !(scale > 0.0) || !isfinite(scale))
    return RABIN_ERR_INVALID_ARG;

  const u64 n = degree + 1;
  const bool complex_in = (im != 0.0);

  double* zre = malloc(n * sizeof(double));
  double* zim = complex_in ? malloc(n * sizeof(double)) : NULL;
  if (zre == NULL || (complex_in && zim == NULL)) {
    free(zre);
    free(zim);
    return RABIN_ERR_OUT_OF_MEMORY;
  }

  // z_i = a^i (1-based i, powers of the number)
  zre[0] = 1.0;
  if (complex_in) zim[0] = 0.0;
  for (u64 i = 1; i < n; i++) {
    if (complex_in) {
      const double pr = zre[i - 1], pi = zim[i - 1];
      zre[i] = pr * re - pi * im;
      zim[i] = pr * im + pi * re;
    } else {
      zre[i] = zre[i - 1] * re;
    }
    if (!isfinite(zre[i]) || (complex_in && !isfinite(zim[i]))) {
      free(zre);
      free(zim);
      return RABIN_ERR_INVALID_ARG;
    }
  }

  const rabin_err_t e = rmat_linear_dependence(zre, zim, n, scale, poly, err);
  free(zre);
  free(zim);
  return e;
}
