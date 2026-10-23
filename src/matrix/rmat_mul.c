/*
 * rmat_mul.c
 *
 * Multiplication for rmat_t: schoolbook matrix product, matrix-vector and
 * vector-matrix products, and the Hadamard product.
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

#include <omp.h>
#include <string.h>

#include "../../include/rmat.h"
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
    TRY(rz_set_u64(&dot, 0));
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
    TRY(rz_set_u64(&dot, 0));
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
