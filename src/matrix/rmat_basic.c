/*
 * rmat_basic.c
 *
 * Constructors and elementwise building blocks for rmat_t: in-place
 * constants, diagonal construction, row/column swaps, stacking and
 * integer matrix powers.
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

rabin_err_t rmat_zero(rmat_t* A)
{
  if (A == NULL) return RABIN_ERR_NULL_PTR;

  for (u64 i = 0; i < A->rows * A->cols; i++) {
    rabin_err_t err = rz_set_u64(&A->data[i], 0);
    if (err != RABIN_SUCCESS) return err;
  }
  return RABIN_SUCCESS;
}

rabin_err_t rmat_ones(rmat_t* A)
{
  if (A == NULL) return RABIN_ERR_NULL_PTR;

  for (u64 i = 0; i < A->rows * A->cols; i++) {
    rabin_err_t err = rz_set_u64(&A->data[i], 1);
    if (err != RABIN_SUCCESS) return err;
  }
  return RABIN_SUCCESS;
}

rabin_err_t rmat_diag(rmat_t* D, const rvec_t* v)
{
  if (D == NULL || v == NULL) return RABIN_ERR_NULL_PTR;
  if (D->rows != v->size || D->cols != v->size) return RABIN_ERR_MATRIX_DIM;

  rabin_err_t err = RABIN_SUCCESS;

  for (u64 i = 0; i < D->rows; i++) {
    for (u64 j = 0; j < D->cols; j++) {
      if (i == j) {
        TRY(rz_copy(RMAT_GET(D, i, j), &v->data[i]));
      } else {
        TRY(rz_set_u64(RMAT_GET(D, i, j), 0));
      }
    }
  }

cleanup:
  return err;
}

rabin_err_t rmat_swap_rows(rmat_t* A, u64 i, u64 j)
{
  if (A == NULL) return RABIN_ERR_NULL_PTR;
  if (i >= A->rows || j >= A->rows) return RABIN_ERR_INVALID_ARG;
  if (i == j) return RABIN_SUCCESS;

  rabin_err_t err = RABIN_SUCCESS;

  for (u64 c = 0; c < A->cols; c++) {
    TRY(rz_swap(RMAT_GET(A, i, c), RMAT_GET(A, j, c)));
  }

cleanup:
  return err;
}

rabin_err_t rmat_swap_cols(rmat_t* A, u64 i, u64 j)
{
  if (A == NULL) return RABIN_ERR_NULL_PTR;
  if (i >= A->cols || j >= A->cols) return RABIN_ERR_INVALID_ARG;
  if (i == j) return RABIN_SUCCESS;

  rabin_err_t err = RABIN_SUCCESS;

  for (u64 r = 0; r < A->rows; r++) {
    TRY(rz_swap(RMAT_GET(A, r, i), RMAT_GET(A, r, j)));
  }

cleanup:
  return err;
}

rabin_err_t rmat_hstack(rmat_t* R, const rmat_t* A, const rmat_t* B)
{
  if (R == NULL || A == NULL || B == NULL) return RABIN_ERR_NULL_PTR;
  if (A->rows != B->rows) return RABIN_ERR_MATRIX_DIM;
  if (A->cols > SIZE_MAX - B->cols) return RABIN_ERR_OVERFLOW;
  if (R->rows != A->rows || R->cols != A->cols + B->cols) {
    return RABIN_ERR_MATRIX_DIM;
  }

  rabin_err_t err = RABIN_SUCCESS;

  for (u64 i = 0; i < A->rows; i++) {
    for (u64 j = 0; j < A->cols; j++) {
      TRY(rz_copy(RMAT_GET(R, i, j), RMAT_GET(A, i, j)));
    }
    for (u64 j = 0; j < B->cols; j++) {
      TRY(rz_copy(RMAT_GET(R, i, A->cols + j), RMAT_GET(B, i, j)));
    }
  }

cleanup:
  return err;
}

rabin_err_t rmat_vstack(rmat_t* R, const rmat_t* A, const rmat_t* B)
{
  if (R == NULL || A == NULL || B == NULL) return RABIN_ERR_NULL_PTR;
  if (A->cols != B->cols) return RABIN_ERR_MATRIX_DIM;
  if (A->rows > SIZE_MAX - B->rows) return RABIN_ERR_OVERFLOW;
  if (R->rows != A->rows + B->rows || R->cols != A->cols) {
    return RABIN_ERR_MATRIX_DIM;
  }

  rabin_err_t err = RABIN_SUCCESS;

  for (u64 i = 0; i < A->rows; i++) {
    for (u64 j = 0; j < A->cols; j++) {
      TRY(rz_copy(RMAT_GET(R, i, j), RMAT_GET(A, i, j)));
    }
  }
  for (u64 i = 0; i < B->rows; i++) {
    for (u64 j = 0; j < A->cols; j++) {
      TRY(rz_copy(RMAT_GET(R, A->rows + i, j), RMAT_GET(B, i, j)));
    }
  }

cleanup:
  return err;
}

rabin_err_t rmat_pow(rmat_t* R, const rmat_t* A, u64 e)
{
  if (R == NULL || A == NULL) return RABIN_ERR_NULL_PTR;
  if (!rmat_is_square(A)) return RABIN_ERR_MATRIX_DIM;
  if (R->rows != A->rows || R->cols != A->cols) return RABIN_ERR_MATRIX_DIM;

  const u64 n = A->rows;
  if (n == 0) return RABIN_SUCCESS;

  if (e == 0) return rmat_id(R, n);

  rabin_err_t err = RABIN_SUCCESS;
  rmat_t acc = {0};
  rmat_t base = {0};
  rmat_t tmp = {0};

  TRY(rmat_init(&acc, n, n));
  TRY(rmat_init(&base, n, n));
  TRY(rmat_init(&tmp, n, n));
  TRY(rmat_id(&acc, n));
  TRY(rmat_copy(&base, A));

  while (e > 0) {
    if (e & 1) {
      TRY(rmat_mul(&tmp, &acc, &base));
      rmat_swap(&acc, &tmp);
    }
    e >>= 1;
    if (e > 0) {
      TRY(rmat_mul(&tmp, &base, &base));
      rmat_swap(&base, &tmp);
    }
  }

  TRY(rmat_copy(R, &acc));

cleanup:
  rmat_clear(&acc);
  rmat_clear(&base);
  rmat_clear(&tmp);
  return err;
}
