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

#include <omp.h>
#include <stdio.h>
#include <string.h>

#include "../../include/rns.h"
#include "../../include/rns_primes.h"
#include "../../include/rpol.h"
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

  rvec_t tmp = {0};
  rabin_err_t err = rvec_init(&tmp, A->cols);
  if (err != RABIN_SUCCESS) return err;

  rz_t dot;
  rz_init(&dot);

  for (u64 i = 0; i < A->rows; i++) {
    if ((err = rmat_get_row(&tmp, A, i)) != RABIN_SUCCESS) goto out;
    if ((err = rvec_dot(&dot, v, &tmp)) != RABIN_SUCCESS) goto out;
    if ((err = rvec_set(r, &dot, i)) != RABIN_SUCCESS) goto out;
  }

  err = RABIN_SUCCESS;
out:
  rz_clear(&dot);
  rvec_clear(&tmp);
  return err;
}

rabin_err_t rmat_vm(rvec_t* r, const rmat_t* A, rvec_t* v)
{
  if (r == NULL || A == NULL || v == NULL) return RABIN_ERR_NULL_PTR;
  if (r->size != A->cols || v->size != A->rows) return RABIN_ERR_MATRIX_DIM;

  rvec_t tmp = {0};
  rabin_err_t err = rvec_init(&tmp, A->rows);
  if (err != RABIN_SUCCESS) return err;

  rz_t dot;
  rz_init(&dot);

  for (u64 i = 0; i < A->cols; i++) {
    if ((err = rmat_get_col(&tmp, A, i)) != RABIN_SUCCESS) goto out;
    if ((err = rvec_dot(&dot, v, &tmp)) != RABIN_SUCCESS) goto out;
    if ((err = rvec_set(r, &dot, i)) != RABIN_SUCCESS) goto out;
  }

  err = RABIN_SUCCESS;
out:
  rz_clear(&dot);
  rvec_clear(&tmp);
  return err;
}

rabin_err_t rmat_hadamard(rz_t* r, const rmat_t* A)
{
  if (r == NULL || A == NULL) return RABIN_ERR_NULL_PTR;

  rabin_err_t err = rz_set_u64(r, 1);
  if (err != RABIN_SUCCESS) return err;

  rvec_t v = {0};
  if ((err = rvec_init(&v, A->rows)) != RABIN_SUCCESS) return err;
  rz_t tmp;
  rz_init(&tmp);

  for (u64 i = 0; i < A->cols; i++) {
    if ((err = rmat_get_col(&v, A, i)) != RABIN_SUCCESS) goto out;

    if ((err = rvec_norm(&tmp, &v)) != RABIN_SUCCESS) goto out;
    // add 1 since norm uses isqrt
    if ((err = rz_add_u64(&tmp, &tmp, 1)) != RABIN_SUCCESS) goto out;
    if ((err = rz_mul(r, r, &tmp)) != RABIN_SUCCESS) goto out;
  }

  err = RABIN_SUCCESS;
out:
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
  rmat_t B_T = {0};
  rabin_err_t err = rmat_init(&B_T, B->cols, B->rows);
  if (err != RABIN_SUCCESS) return err;

  for (u64 i = 0; i < B->rows; i++) {
    for (u64 j = 0; j < B->cols; j++) {
      if ((err = rz_copy(&B_T.data[j * B_T.cols + i],
                         &B->data[i * B->cols + j])) != RABIN_SUCCESS)
        goto out;
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

  err = RABIN_SUCCESS;
out:
  rmat_clear(&B_T);
  return err;
}

rabin_err_t rmat_det(rz_t* d, const rmat_t* A)
{
  if (d == NULL || A == NULL) return RABIN_ERR_NULL_PTR;
  if (A->cols != A->rows) return RABIN_ERR_MATRIX_DIM;

  rmat_t T = {0};
  rabin_err_t err = rmat_init(&T, A->cols, A->rows);
  if (err != RABIN_SUCCESS) return err;
  if ((err = rmat_copy(&T, A)) != RABIN_SUCCESS) goto out;

  u64 k = rns_estimate_det(&T);

  rns_ctx_t ctx;
  if ((err = rns_ctx_init(&ctx, RNS_PRIMES, k)) != RABIN_SUCCESS) goto out;

  if ((err = rmat_det_rns(d, &T, &ctx)) != RABIN_SUCCESS) goto out;

  rns_ctx_clear(&ctx);
  err = RABIN_SUCCESS;
out:
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
  rmat_t T = {0};
  rabin_err_t err = rmat_init(&T, n, n);
  if (err != RABIN_SUCCESS) return err;
  if ((err = rmat_copy(&T, A)) != RABIN_SUCCESS) goto out;

  rz_t p0;
  rz_init(&p0);
  if ((err = rz_set_u64(&p0, 1)) != RABIN_SUCCESS) goto out;

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
    err = rz_set_u64(det, 0);
  } else {
    if ((err = rz_copy(det, RMAT_GET(&T, n - 1, n - 1))) != RABIN_SUCCESS)
      goto out;
    if (sign < 0) {
      if ((err = rz_neg(det, det)) != RABIN_SUCCESS) goto out;
    }
  }

out:
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

  rmat_t C = {0}, T = {0};
  rabin_err_t err = rmat_init(&C, n, n);
  if (err != RABIN_SUCCESS) return err;
  if ((err = rmat_init(&T, n, n)) != RABIN_SUCCESS) goto fail_C;
  if ((err = rmat_id(&C, n)) != RABIN_SUCCESS) goto fail_T;

  rz_t temp, trace;
  rz_init_multi(&temp, &trace, NULL);
  rz_t* a = malloc((n + 1) * sizeof(rz_t));
  if (a == NULL) goto fail_T;
  for (u64 i = 0; i <= n; i++) {
    rz_init(&a[i]);
  }
  if ((err = rz_set_u64(&a[0], 1)) != RABIN_SUCCESS) goto fail_a;

  // 2. [Finished?]
  for (u64 i = 1; i <= n; i++) {
    // 3. [Compute next a_i and C]
    if ((err = rmat_mul(&T, A, &C)) != RABIN_SUCCESS) goto fail_a;

    if ((err = rmat_trace(&trace, &T)) != RABIN_SUCCESS) goto fail_a;

    if ((err = rz_set_u64(&temp, i)) != RABIN_SUCCESS) goto fail_a;
    rz_divmod_u64(&a[i], &trace, i);

    if ((err = rz_neg(&a[i], &a[i])) != RABIN_SUCCESS) goto fail_a;

    if (J && i == n) {
      if ((err = rmat_copy(J, &C)) != RABIN_SUCCESS) goto fail_a;
      if ((n - 1) & 1) {
        if ((err = rmat_neg(J)) != RABIN_SUCCESS) goto fail_a;
      }
    }
    if ((err = rmat_copy(&C, &T)) != RABIN_SUCCESS) goto fail_a;
    for (u64 j = 0; j < n; j++) {
      if ((err = rz_add(&temp, RMAT_GET(&C, j, j), &a[i])) != RABIN_SUCCESS)
        goto fail_a;
      if ((err = rmat_set(&C, &temp, j, j)) != RABIN_SUCCESS) goto fail_a;
    }
  }

  if ((err = rpol_alloc(p, n + 1)) != RABIN_SUCCESS) goto fail_a;

  for (u64 i = 0; i <= n; i++) {
    if ((err = rz_copy(&p->coeff[n - i], &a[i])) != RABIN_SUCCESS) goto fail_a;
    rz_clear(&a[i]);
  }

  p->deg = n;
  err = RABIN_SUCCESS;

fail_a:
  free(a);
fail_T:
  rmat_clear(&T);
fail_C:
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

  rz_t a, b;
  rz_init_multi(&a, &b, NULL);
  rabin_err_t err = rz_set_u64(&a, 1);
  if (err != RABIN_SUCCESS) goto out;
  if ((err = rz_set_u64(&b, 0)) != RABIN_SUCCESS) goto out;

  for (u64 i = 0; i < n; i++) {
    for (u64 j = 0; j < n; j++) {
      if (i == j) {
        if ((err = rmat_set(I, &a, i, j)) != RABIN_SUCCESS) goto out;
      } else {
        if ((err = rmat_set(I, &b, i, j)) != RABIN_SUCCESS) goto out;
      }
    }
  }

  err = RABIN_SUCCESS;
out:
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

  rmat_t A_work = {0};
  rabin_err_t err = rmat_init(&A_work, m, n);
  if (err != RABIN_SUCCESS) return err;
  if ((err = rmat_copy(&A_work, A)) != RABIN_SUCCESS) goto cleanup;

  rz_t b, q, rem, temp, one, t2;
  rz_init_multi(&b, &q, &rem, &temp, &one, &t2, NULL);
  if ((err = rz_set_u64(&one, 1)) != RABIN_SUCCESS) goto cleanup;

// 2. [Row finished?]
step2:
  if (rmat_row_prefix_is_zero(&A_work, i, k)) {
    if (RMAT_GET(&A_work, i, k)->is_neg) {
      if ((err = rmat_neg_col(&A_work, k)) != RABIN_SUCCESS) goto cleanup;
    }
    goto step5;
  }

  // 3. [Choose non-zero entry]
  u64 j_0 = rmat_min_abs_index(&A_work, i, k);
  if ((i64)j_0 < k) {
    if ((err = rmat_swap_col(&A_work, k, j_0)) != RABIN_SUCCESS) goto cleanup;
  }
  if (RMAT_GET(&A_work, i, k)->is_neg) {
    if ((err = rmat_neg_col(&A_work, k)) != RABIN_SUCCESS) goto cleanup;
  }

  if ((err = rz_copy(&b, RMAT_GET(&A_work, i, k))) != RABIN_SUCCESS)
    goto cleanup;

  // 4. [Reduce]
  for (u64 j = 0; j < (u64)k; j++) {
    if ((err = rz_div_euclid(&q, RMAT_GET(&A_work, i, j), &b)) != RABIN_SUCCESS)
      goto cleanup;
    for (u64 x = 0; x < m; x++) {
      if ((err = rz_mul(&temp, &q, RMAT_GET(&A_work, x, k))) != RABIN_SUCCESS)
        goto cleanup;
      if ((err = rz_sub(&t2, RMAT_GET(&A_work, x, j), &temp)) != RABIN_SUCCESS)
        goto cleanup;
      if ((err = rmat_set(&A_work, &t2, x, j)) != RABIN_SUCCESS) goto cleanup;
    }
  }
  goto step2;

// 5. [Final reductions]
step5:
  if ((err = rz_copy(&b, RMAT_GET(&A_work, i, k))) != RABIN_SUCCESS)
    goto cleanup;
  if (rz_is_zero(&b)) {
    k++;
    goto step6;
  }

  for (u64 j = k + 1; j < n; j++) {
    if ((err = rz_div_euclid(&q, RMAT_GET(&A_work, i, j), &b)) != RABIN_SUCCESS)
      goto cleanup;
    for (u64 x = 0; x < m; x++) {
      if ((err = rz_mul(&temp, &q, RMAT_GET(&A_work, x, k))) != RABIN_SUCCESS)
        goto cleanup;
      if ((err = rz_sub(&t2, RMAT_GET(&A_work, x, j), &temp)) != RABIN_SUCCESS)
        goto cleanup;
      if ((err = rmat_set(&A_work, &t2, x, j)) != RABIN_SUCCESS) goto cleanup;
    }
  }

// 6. [Finished?]
step6:
  if (i != (i64)l) {
    i--;
    k--;
    goto step2;
  }

  if ((err = rmat_init(W, m, n - k)) != RABIN_SUCCESS) goto cleanup;

  // copy result
  for (u64 j = 0; j < n - k; j++) {
    for (u64 x = 0; x < m; x++) {
      if ((err = rz_copy(RMAT_GET(W, x, j), RMAT_GET(&A_work, x, j + k))) !=
          RABIN_SUCCESS)
        goto cleanup;
    }
  }

  err = RABIN_SUCCESS;
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
    if ((err = rz_set_u64(u, 0)) != RABIN_SUCCESS) goto cleanup;
    if (rz_is_zero(b)) {
      if ((err = rz_set_u64(v, 0)) != RABIN_SUCCESS) goto cleanup;
    } else if (b->is_neg) {
      if ((err = rz_set_i64(v, -1)) != RABIN_SUCCESS) goto cleanup;
    } else {
      if ((err = rz_set_u64(v, 1)) != RABIN_SUCCESS) goto cleanup;
    }
    goto cleanup;
  }

  if (rz_is_zero(b)) {
    // u*sign(a) = 1 is forced; v = 0
    if (a->is_neg) {
      if ((err = rz_set_i64(u, -1)) != RABIN_SUCCESS) goto cleanup;
    } else {
      if ((err = rz_set_u64(u, 1)) != RABIN_SUCCESS) goto cleanup;
    }
    if ((err = rz_set_u64(v, 0)) != RABIN_SUCCESS) goto cleanup;
    goto cleanup;
  }

  if ((err = rz_copy(&p, a)) != RABIN_SUCCESS) goto cleanup;
  p.is_neg = false;

  // a | b  =>  d = |a|  =>  v = 0, u = sign(a)  (the book's essential
  // condition; it also keeps the Euclidean step from stalling)
  if (rz_cmp(&p, d) == 0) {
    if (a->is_neg) {
      if ((err = rz_set_i64(u, -1)) != RABIN_SUCCESS) goto cleanup;
    } else {
      if ((err = rz_set_u64(u, 1)) != RABIN_SUCCESS) goto cleanup;
    }
    if ((err = rz_set_u64(v, 0)) != RABIN_SUCCESS) goto cleanup;
    goto cleanup;
  }

  // centered half-period: v = v0 - t*(a/d) with |a/d| >= 2 and
  // v in (-|a/d|/2, |a/d|/2]; recompute u = (d - v*b)/a (exact)
  if ((err = rz_div_exact(&p, a, d)) != RABIN_SUCCESS) goto cleanup;
  if ((err = rz_copy(&pmag, &p)) != RABIN_SUCCESS) goto cleanup;
  pmag.is_neg = false;
  if ((err = rz_set_u64(&two, 2)) != RABIN_SUCCESS) goto cleanup;
  if ((err = rz_div(&half, &pmag, &two)) != RABIN_SUCCESS) goto cleanup;

  // r = v0 mod |p| in [0, |p|)
  if ((err = rz_divmod(&quot, &r, v, &pmag)) != RABIN_SUCCESS) goto cleanup;
  if (r.is_neg) {
    if ((err = rz_add(&r, &r, &pmag)) != RABIN_SUCCESS) goto cleanup;
  }
  if (rz_cmp(&r, &half) > 0) {
    if ((err = rz_sub(v, &r, &pmag)) != RABIN_SUCCESS) goto cleanup;
  } else {
    if ((err = rz_copy(v, &r)) != RABIN_SUCCESS) goto cleanup;
  }

  if ((err = rz_mul(&t, v, b)) != RABIN_SUCCESS) goto cleanup;
  if ((err = rz_sub(&t, d, &t)) != RABIN_SUCCESS) goto cleanup;
  if ((err = rz_div_exact(u, &t, a)) != RABIN_SUCCESS) goto cleanup;

cleanup:
  rz_clear_multi(&p, &pmag, &half, &quot, &r, &t, &two, NULL);
  return err;
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

  rmat_t A_work = {0};
  rabin_err_t err = rmat_init(&A_work, m, n);
  if (err != RABIN_SUCCESS) return err;
  if ((err = rmat_copy(&A_work, A)) != RABIN_SUCCESS) goto cleanup;

  rz_t R, u, v, d, q, t1, t2, t3, half;
  rz_init_multi(&R, &u, &v, &d, &q, &t1, &t2, &t3, &half, NULL);
  if ((err = rz_copy(&R, D)) != RABIN_SUCCESS) goto cleanup;
  R.is_neg = false;
  if ((err = rz_rshift(&half, &R, 1)) != RABIN_SUCCESS) goto cleanup;

  // reduce the input entries into (-R/2, R/2]; all column operations are
  // taken modulo R, so only the classes modulo R matter
  for (u64 x = 0; x < m * n; x++) {
    if ((err = rmat_hnf_centered_mod(&A_work.data[x], &A_work.data[x], &R,
                                     &half)) != RABIN_SUCCESS)
      goto cleanup;
  }

  rz_t* B;
  if (m == 0) {
    B = NULL;
  } else {
    B = malloc(m * sizeof(rz_t));
    if (B == NULL) {
      err = RABIN_ERR_OUT_OF_MEMORY;
      goto cleanup;
    }
  }
  for (u64 x = 0; x < m; x++) {
    rz_init(&B[x]);
  }

  rmat_t W_local = {0};
  if ((err = rmat_init(&W_local, m, m)) != RABIN_SUCCESS) goto cleanup_B;

  for (i64 i = (i64)m - 1; i >= 0; i--) {
    i64 k = (i64)n - (i64)m + i;

    // 2. [Check zero] / 3. [Euclidean step]
    for (i64 j = k - 1; j >= 0; j--) {
      if (rz_is_zero(RMAT_GET(&A_work, i, j))) {
        continue;
      }

      if ((err = rz_gcd_extended_lehmer(&u, &v, &d, RMAT_GET(&A_work, i, k),
                                        RMAT_GET(&A_work, i, j))) !=
          RABIN_SUCCESS)
        goto cleanup_B;
      if ((err = hnf_bezout_minimal(&u, &v, &d, RMAT_GET(&A_work, i, k),
                                    RMAT_GET(&A_work, i, j))) != RABIN_SUCCESS)
        goto cleanup_B;

      // B = u*A_k + v*A_j, reduced into (-R/2, R/2]
      for (u64 x = 0; x < m; x++) {
        if ((err = rz_mul(&B[x], &u, RMAT_GET(&A_work, x, k))) != RABIN_SUCCESS)
          goto cleanup_B;
        if ((err = rz_mul(&t2, &v, RMAT_GET(&A_work, x, j))) != RABIN_SUCCESS)
          goto cleanup_B;
        if ((err = rz_add(&B[x], &B[x], &t2)) != RABIN_SUCCESS) goto cleanup_B;
        if ((err = rmat_hnf_centered_mod(&B[x], &B[x], &R, &half)) !=
            RABIN_SUCCESS)
          goto cleanup_B;
      }

      // A_j = (a_{i,k}/d)*A_j - (a_{i,j}/d)*A_k, reduced into (-R/2, R/2]
      if ((err = rz_div_exact(&t2, RMAT_GET(&A_work, i, k), &d)) !=
          RABIN_SUCCESS)
        goto cleanup_B;
      if ((err = rz_div_exact(&q, RMAT_GET(&A_work, i, j), &d)) !=
          RABIN_SUCCESS)
        goto cleanup_B;
      for (u64 x = 0; x < m; x++) {
        if ((err = rz_mul(&t3, &t2, RMAT_GET(&A_work, x, j))) != RABIN_SUCCESS)
          goto cleanup_B;
        if ((err = rz_mul(&t1, &q, RMAT_GET(&A_work, x, k))) != RABIN_SUCCESS)
          goto cleanup_B;
        if ((err = rz_sub(&t1, &t3, &t1)) != RABIN_SUCCESS) goto cleanup_B;
        if ((err = rmat_hnf_centered_mod(RMAT_GET(&A_work, x, j), &t1, &R,
                                         &half)) != RABIN_SUCCESS)
          goto cleanup_B;
        if ((err = rz_copy(RMAT_GET(&A_work, x, k), &B[x])) != RABIN_SUCCESS)
          goto cleanup_B;
      }

      // a_{i,k} = u*a_{i,k} + v*a_{i,j} = d exactly
      if ((err = rz_copy(RMAT_GET(&A_work, i, k), &d)) != RABIN_SUCCESS)
        goto cleanup_B;
    }

    // 4. [Next row]
    // u*a_{i,k} + v*R = d = gcd(a_{i,k}, R)
    if ((err = rz_gcd_extended_lehmer(&u, &v, &d, RMAT_GET(&A_work, i, k),
                                      &R)) != RABIN_SUCCESS)
      goto cleanup_B;
    if ((err = hnf_bezout_minimal(&u, &v, &d, RMAT_GET(&A_work, i, k), &R)) !=
        RABIN_SUCCESS)
      goto cleanup_B;

    // W_i = u*A_k mod R, taken in [0, R)
    for (u64 x = 0; x < m; x++) {
      if ((err = rz_mul(&t1, &u, RMAT_GET(&A_work, x, k))) != RABIN_SUCCESS)
        goto cleanup_B;
      if ((err = rz_mod_pos(RMAT_GET(&W_local, x, i), &t1, &R)) !=
          RABIN_SUCCESS)
        goto cleanup_B;
    }

    // if d = R (i.e. R | a_{i,k}) the diagonal entry would be 0; set it to d
    if (rz_cmp(&d, &R) == 0) {
      if ((err = rz_copy(RMAT_GET(&W_local, i, i), &d)) != RABIN_SUCCESS)
        goto cleanup_B;
    }

    // final reductions: W_j -= floor(W_{i,j}/W_{i,i}) * W_i for j > i
    for (u64 j = (u64)i + 1; j < m; j++) {
      if (rz_is_zero(RMAT_GET(&W_local, i, j))) {
        continue;
      }
      if ((err = rz_div_euclid(&q, RMAT_GET(&W_local, i, j),
                               RMAT_GET(&W_local, i, i))) != RABIN_SUCCESS)
        goto cleanup_B;
      for (u64 x = 0; x <= (u64)i; x++) {
        if ((err = rz_mul(&t1, &q, RMAT_GET(&W_local, x, i))) != RABIN_SUCCESS)
          goto cleanup_B;
        if ((err = rz_sub(RMAT_GET(&W_local, x, j), RMAT_GET(&W_local, x, j),
                          &t1)) != RABIN_SUCCESS)
          goto cleanup_B;
      }
    }

    R.is_neg = false;
    if ((err = rz_div_exact(&R, &R, &d)) != RABIN_SUCCESS) goto cleanup_B;
    if ((err = rz_rshift(&half, &R, 1)) != RABIN_SUCCESS) goto cleanup_B;

    if (i > 0) {
      // working modulo R, a_{i-1,k-1} may have reduced to zero; replace it
      // by any nonzero multiple of R
      if (rz_is_zero(RMAT_GET(&A_work, i - 1, k - 1))) {
        if ((err = rz_copy(RMAT_GET(&A_work, i - 1, k - 1), &R)) !=
            RABIN_SUCCESS)
          goto cleanup_B;
      }
    }
  }

  if ((err = rmat_swap(W, &W_local)) != RABIN_SUCCESS) goto cleanup_B;
  rmat_clear(&W_local);
  err = RABIN_SUCCESS;

cleanup_B:
  for (u64 x = 0; x < m; x++) {
    if (B != NULL) rz_clear(&B[x]);
  }
  free(B);
cleanup:
  rz_clear_multi(&R, &u, &v, &d, &q, &t1, &t2, &t3, &half, NULL);
  rmat_clear(&A_work);
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

  rmat_t A_work = {0};
  rabin_err_t err = rmat_init(&A_work, m, n);
  if (err != RABIN_SUCCESS) return err;
  if ((err = rmat_copy(&A_work, A)) != RABIN_SUCCESS) goto cleanup;

  rz_t* B;
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
  if ((err = rz_set_u64(&one, 1)) != RABIN_SUCCESS) goto cleanup_B;

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
  if ((err = rz_gcd_extended_lehmer(&u, &v, &d, RMAT_GET(&A_work, i, k),
                                    RMAT_GET(&A_work, i, j))) != RABIN_SUCCESS)
    goto cleanup_B;
  if ((err = hnf_bezout_minimal(&u, &v, &d, RMAT_GET(&A_work, i, k),
                                RMAT_GET(&A_work, i, j))) != RABIN_SUCCESS)
    goto cleanup_B;

  for (u64 x = 0; x < m; x++) {
    if ((err = rz_mul(&B[x], &u, RMAT_GET(&A_work, x, k))) != RABIN_SUCCESS)
      goto cleanup_B;
    if ((err = rz_mul(&temp, &v, RMAT_GET(&A_work, x, j))) != RABIN_SUCCESS)
      goto cleanup_B;
    if ((err = rz_add(&B[x], &B[x], &temp)) != RABIN_SUCCESS) goto cleanup_B;
  }

  if ((err = rz_div_euclid(&q_k, RMAT_GET(&A_work, i, k), &d)) != RABIN_SUCCESS)
    goto cleanup_B;
  if ((err = rz_div_euclid(&q_j, RMAT_GET(&A_work, i, j), &d)) != RABIN_SUCCESS)
    goto cleanup_B;

  for (u64 x = 0; x < m; x++) {
    if ((err = rz_mul(&temp, &q_k, RMAT_GET(&A_work, x, j))) != RABIN_SUCCESS)
      goto cleanup_B;
    if ((err = rz_mul(&temp2, &q_j, RMAT_GET(&A_work, x, k))) != RABIN_SUCCESS)
      goto cleanup_B;
    if ((err = rz_sub(RMAT_GET(&A_work, x, j), &temp, &temp2)) != RABIN_SUCCESS)
      goto cleanup_B;
    if ((err = rz_copy(RMAT_GET(&A_work, x, k), &B[x])) != RABIN_SUCCESS)
      goto cleanup_B;
  }
  goto step2;

  // 4. [Final reduction]
step4:
  if ((err = rz_copy(&b, RMAT_GET(&A_work, i, k))) != RABIN_SUCCESS)
    goto cleanup_B;
  if (b.is_neg && !rz_is_zero(&b)) {
    for (u64 x = 0; x < m; x++) {
      if ((err = rz_neg(RMAT_GET(&A_work, x, k), RMAT_GET(&A_work, x, k))) !=
          RABIN_SUCCESS)
        goto cleanup_B;
    }
    if ((err = rz_neg(&b, &b)) != RABIN_SUCCESS) goto cleanup_B;
  }

  if (rz_is_zero(&b)) {
    k++;
    goto step5;
  }

  for (u64 j_0 = k + 1; j_0 < n; j_0++) {
    if ((err = rz_div_euclid(&temp2, RMAT_GET(&A_work, i, j_0), &b)) !=
        RABIN_SUCCESS)
      goto cleanup_B;
    for (u64 x = 0; x < m; x++) {
      if ((err = rz_mul(&temp, &temp2, RMAT_GET(&A_work, x, k))) !=
          RABIN_SUCCESS)
        goto cleanup_B;
      if ((err = rz_sub(RMAT_GET(&A_work, x, j_0), RMAT_GET(&A_work, x, j_0),
                        &temp)) != RABIN_SUCCESS)
        goto cleanup_B;
    }
  }

step5:
  // 5. [Finished]
  if (i == l) {
    if ((err = rmat_init(W, m, n - k)) != RABIN_SUCCESS) goto cleanup_B;
    for (u64 x = 0; x < n - k; x++) {
      for (u64 y = 0; y < m; y++) {
        if ((err = rz_copy(RMAT_GET(W, y, x), RMAT_GET(&A_work, y, x + k))) !=
            RABIN_SUCCESS)
          goto cleanup_B;
      }
    }
    err = RABIN_SUCCESS;
    goto cleanup_B;
  }

  i--;
  k--;
  j = k;
  goto step2;

cleanup_B:
  for (u64 x = 0; x < m; x++) {
    rz_clear(&B[x]);
  }
  free(B);
cleanup:
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
  rz_t R;
  rz_init(&R);
  rabin_err_t err = rmat_det(&R, A);
  if (err != RABIN_SUCCESS) return err;
  R.is_neg = false;

  if (n == 1) {
    if ((err = rmat_set(S, &R, 0, 0)) != RABIN_SUCCESS) goto fail_R1;
    err = RABIN_SUCCESS;
  fail_R1:
    rz_clear(&R);
    return err;
  }

  rz_t u, v, d, q_i, q_j, t1, t2, b;
  rz_init_multi(&u, &v, &d, &q_i, &q_j, &t1, &t2, &b, NULL);

  rmat_t A_work = {0};
  if ((err = rmat_init(&A_work, n, n)) != RABIN_SUCCESS) goto fail_R;
  if ((err = rmat_copy(&A_work, A)) != RABIN_SUCCESS) goto fail_AW;

  rz_t* B = malloc(n * sizeof(rz_t));
  if (B == NULL) {
    err = RABIN_ERR_OUT_OF_MEMORY;
    goto fail_AW;
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
    goto fail_B;

  // use Remark to find u,v minimal
  rz_t a_abs;
  rz_init(&a_abs);
  bool a_ii_live = false;
  rz_t a_ii_abs;
  rz_init(&a_ii_abs);
  if ((err = rz_copy(&a_abs, RMAT_GET(&A_work, i - 1, j - 1))) != RABIN_SUCCESS)
    goto fail_temps;
  a_abs.is_neg = false;
  if ((err = rz_copy(&a_ii_abs, RMAT_GET(&A_work, i - 1, i - 1))) !=
      RABIN_SUCCESS)
    goto fail_temps;
  a_ii_live = true;
  a_ii_abs.is_neg = false;

  // if d == |a_{i,i}|, force u = sign(a_{i,i}), v = 0
  if (rz_cmp(&d, &a_ii_abs) == 0) {
    if ((err = rz_set_u64(&u, 1)) != RABIN_SUCCESS) goto fail_temps;
    if ((err = rz_set_u64(&v, 0)) != RABIN_SUCCESS) goto fail_temps;
    if (RMAT_GET(&A_work, i - 1, i - 1)->is_neg)
      if ((err = rz_set_i64(&u, -1)) != RABIN_SUCCESS) goto fail_temps;
  }
  // if d == |a_{i,j}|, force u = 0, v = sign(a_{i,j})
  else if (rz_cmp(&d, &a_abs) == 0) {
    if ((err = rz_set_u64(&u, 0)) != RABIN_SUCCESS) goto fail_temps;
    if ((err = rz_set_u64(&v, 1)) != RABIN_SUCCESS) goto fail_temps;
    if (RMAT_GET(&A_work, i - 1, j - 1)->is_neg)
      if ((err = rz_set_i64(&v, -1)) != RABIN_SUCCESS) goto fail_temps;
  }

  rz_clear(&a_abs);
  rz_clear(&a_ii_abs);
  a_ii_live = false;

  // B = uA_i + vA_j
  for (u64 x = 0; x < (u64)n; x++) {
    if ((err = rz_mul(&B[x], &u, RMAT_GET(&A_work, x, i - 1))) != RABIN_SUCCESS)
      goto fail_B;
    if ((err = rz_mul(&t1, &v, RMAT_GET(&A_work, x, j - 1))) != RABIN_SUCCESS)
      goto fail_B;
    if ((err = rz_add(&B[x], &B[x], &t1)) != RABIN_SUCCESS) goto fail_B;
  }

  // A_j = ((a_ii / d)A_j - (a_ij / d)A_i) mod R
  if ((err = rz_div(&q_i, RMAT_GET(&A_work, i - 1, i - 1), &d)) !=
      RABIN_SUCCESS)
    goto fail_B;
  if ((err = rz_div(&q_j, RMAT_GET(&A_work, i - 1, j - 1), &d)) !=
      RABIN_SUCCESS)
    goto fail_B;
  for (u64 x = 0; x < (u64)n; x++) {
    if ((err = rz_mul(&t1, &q_i, RMAT_GET(&A_work, x, j - 1))) != RABIN_SUCCESS)
      goto fail_B;
    if ((err = rz_mul(&t2, &q_j, RMAT_GET(&A_work, x, i - 1))) != RABIN_SUCCESS)
      goto fail_B;
    if ((err = rz_sub(&t1, &t1, &t2)) != RABIN_SUCCESS) goto fail_B;
    if ((err = rz_mod(&t1, &t1, &R)) != RABIN_SUCCESS) goto fail_B;
    if ((err = rz_copy(RMAT_GET(&A_work, x, j - 1), &t1)) != RABIN_SUCCESS)
      goto fail_B;
  }

  // A_i = B mod R
  for (u64 x = 0; x < (u64)n; x++) {
    if ((err = rz_mod(RMAT_GET(&A_work, x, i - 1), &B[x], &R)) != RABIN_SUCCESS)
      goto fail_B;
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
    goto fail_B;

  // use Remark to find u,v minimal safely
  rz_t a_abs_col;
  rz_init(&a_abs_col);
  bool a_ii_col_live = false;
  rz_t a_ii_abs_col;
  rz_init(&a_ii_abs_col);
  if ((err = rz_copy(&a_abs_col, RMAT_GET(&A_work, j - 1, i - 1))) !=
      RABIN_SUCCESS)
    goto fail_temps_col;
  a_abs_col.is_neg = false;
  if ((err = rz_copy(&a_ii_abs_col, RMAT_GET(&A_work, i - 1, i - 1))) !=
      RABIN_SUCCESS)
    goto fail_temps_col;
  a_ii_col_live = true;
  a_ii_abs_col.is_neg = false;

  // if d == |a_{i,i}|, force u = sign(a_{i,i}), v = 0
  if (rz_cmp(&d, &a_ii_abs_col) == 0) {
    if ((err = rz_set_u64(&u, 1)) != RABIN_SUCCESS) goto fail_temps_col;
    if ((err = rz_set_u64(&v, 0)) != RABIN_SUCCESS) goto fail_temps_col;
    if (RMAT_GET(&A_work, i - 1, i - 1)->is_neg)
      if ((err = rz_set_i64(&u, -1)) != RABIN_SUCCESS) goto fail_temps_col;
  }
  // if d == |a_{j,i}|, force u = 0, v = sign(a_{j,i})
  else if (rz_cmp(&d, &a_abs_col) == 0) {
    if ((err = rz_set_u64(&u, 0)) != RABIN_SUCCESS) goto fail_temps_col;
    if ((err = rz_set_u64(&v, 1)) != RABIN_SUCCESS) goto fail_temps_col;
    if (RMAT_GET(&A_work, j - 1, i - 1)->is_neg)
      if ((err = rz_set_i64(&v, -1)) != RABIN_SUCCESS) goto fail_temps_col;
  }

  rz_clear(&a_abs_col);
  rz_clear(&a_ii_abs_col);
  a_ii_col_live = false;

  // B = uA'_i + vA'_j
  for (u64 x = 0; x < (u64)n; x++) {
    if ((err = rz_mul(&B[x], &u, RMAT_GET(&A_work, i - 1, x))) != RABIN_SUCCESS)
      goto fail_B;
    if ((err = rz_mul(&t1, &v, RMAT_GET(&A_work, j - 1, x))) != RABIN_SUCCESS)
      goto fail_B;
    if ((err = rz_add(&B[x], &B[x], &t1)) != RABIN_SUCCESS) goto fail_B;
  }

  // A'_j = ((a_ii / d)A'_j - (a_ij / d)A'_i) mod R
  if ((err = rz_div(&q_i, RMAT_GET(&A_work, i - 1, i - 1), &d)) !=
      RABIN_SUCCESS)
    goto fail_B;
  if ((err = rz_div(&q_j, RMAT_GET(&A_work, j - 1, i - 1), &d)) !=
      RABIN_SUCCESS)
    goto fail_B;
  for (u64 x = 0; x < (u64)n; x++) {
    if ((err = rz_mul(&t1, &q_i, RMAT_GET(&A_work, j - 1, x))) != RABIN_SUCCESS)
      goto fail_B;
    if ((err = rz_mul(&t2, &q_j, RMAT_GET(&A_work, i - 1, x))) != RABIN_SUCCESS)
      goto fail_B;
    if ((err = rz_sub(&t1, &t1, &t2)) != RABIN_SUCCESS) goto fail_B;
    if ((err = rz_mod(&t1, &t1, &R)) != RABIN_SUCCESS) goto fail_B;
    if ((err = rz_copy(RMAT_GET(&A_work, j - 1, x), &t1)) != RABIN_SUCCESS)
      goto fail_B;
  }

  // A'_i = B mod R
  for (u64 x = 0; x < (u64)n; x++) {
    if ((err = rz_mod(RMAT_GET(&A_work, i - 1, x), &B[x], &R)) != RABIN_SUCCESS)
      goto fail_B;
  }

  c++;
  goto step6;

// 8. [Repeat stage i?]
step8:
  if (c > 0) {
    goto step2;
  }

  // 9. [Check the rest of the matrix]
  if ((err = rz_copy(&b, RMAT_GET(&A_work, i - 1, i - 1))) != RABIN_SUCCESS)
    goto fail_B;
  for (i64 k_0 = 0; k_0 < i - 1; k_0++) {
    for (i64 l_0 = 0; l_0 < i - 1; l_0++) {
      bool not_divisible = false;
      if (rz_is_zero(&b)) {
        if (!rz_is_zero(RMAT_GET(&A_work, k_0, l_0))) {
          not_divisible = true;
        }
      } else {
        if ((err = rz_mod(&t1, RMAT_GET(&A_work, k_0, l_0), &b)) !=
            RABIN_SUCCESS)
          goto fail_B;
        if (!rz_is_zero(&t1)) {
          not_divisible = true;
        }
      }

      if (not_divisible) {
        for (i64 x = 0; x < n; x++) {
          if ((err = rz_add(RMAT_GET(&A_work, i - 1, x),
                            RMAT_GET(&A_work, i - 1, x),
                            RMAT_GET(&A_work, k_0, x))) != RABIN_SUCCESS)
            goto fail_B;
        }
        goto step2;
      }
    }
  }

  // 10. [Next stage]
  if ((err = rz_gcd_lehmer(&t1, RMAT_GET(&A_work, i - 1, i - 1), &R)) !=
      RABIN_SUCCESS)
    goto fail_B;
  if ((err = rmat_set(&A_work, &t1, i - 1, i - 1)) != RABIN_SUCCESS)
    goto fail_B;
  if ((err = rz_div(&R, &R, &t1)) != RABIN_SUCCESS) goto fail_B;

  if (i == 2) {
    if ((err = rz_gcd_lehmer(&t2, RMAT_GET(&A_work, 0, 0), &R)) !=
        RABIN_SUCCESS)
      goto fail_B;
    if ((err = rmat_set(&A_work, &t2, 0, 0)) != RABIN_SUCCESS) goto fail_B;

    // order diagonal entries in proper order
    rz_t zero;
    rz_init(&zero);
    for (i64 r = 0; r < n; r++) {
      for (i64 c = 0; c < n; c++) {
        if (r == c) {
          if ((err = rmat_set(S, RMAT_GET(&A_work, n - 1 - r, n - 1 - r), r,
                              r)) != RABIN_SUCCESS)
            goto fail_zero;
        } else {
          if ((err = rmat_set(S, &zero, r, c)) != RABIN_SUCCESS) goto fail_zero;
        }
      }
    }
    rz_clear(&zero);
    err = RABIN_SUCCESS;
    goto fail_B;
  }

  i--;
  goto step2;

fail_zero:
fail_temps_col:
  rz_clear(&a_abs_col);
  if (a_ii_col_live) rz_clear(&a_ii_abs_col);
fail_temps:
  rz_clear(&a_abs);
  if (a_ii_live) rz_clear(&a_ii_abs);
fail_B:
  for (i64 x = 0; x < n; x++) {
    rz_clear(&B[x]);
  }
  free(B);
fail_AW:
  rmat_clear(&A_work);
fail_R:
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
