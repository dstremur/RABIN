/*
 * rmatq.c
 *
 * Rational matrices (common denominator): lifecycle, element access, arithmetic
 * and printing.
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

#include "../../include/rmatq.h"

#include <stdio.h>
#include <string.h>
rabin_err_t rmatq_init(rmatq_t* A, u64 m, u64 n)
{
  if (A == NULL) return RABIN_ERR_NULL_PTR;

  rabin_err_t err = rmat_init(&A->M, m, n);
  if (err != RABIN_SUCCESS) return err;

  rz_init(&A->den);
  return rz_set_u64(&A->den, 1);
}

rabin_err_t rmatq_clear(rmatq_t* A)
{
  if (A == NULL) return RABIN_ERR_NULL_PTR;

  rmat_clear(&A->M);
  rz_clear(&A->den);
  return RABIN_SUCCESS;
}

rabin_err_t rmatq_copy(rmatq_t* R, const rmatq_t* A)
{
  if (R == NULL || A == NULL) return RABIN_ERR_NULL_PTR;

  rabin_err_t err = rmat_copy(&R->M, &A->M);
  if (err != RABIN_SUCCESS) return err;

  return rz_copy(&R->den, &A->den);
}

rabin_err_t rmatq_copyZ(rmatq_t* R, rmat_t* A)
{
  if (R == NULL || A == NULL) return RABIN_ERR_NULL_PTR;

  rabin_err_t err = rmat_copy(&R->M, A);
  if (err != RABIN_SUCCESS) return err;

  return rz_set_u64(&R->den, 1);
}

rabin_err_t rmatq_add(rmatq_t* R, const rmatq_t* A, const rmatq_t* B)
{
  if (R == NULL || A == NULL || B == NULL) return RABIN_ERR_NULL_PTR;
  if (A->M.rows != B->M.rows || A->M.cols != B->M.cols)
    return RABIN_ERR_MATRIX_DIM;
  if (R->M.rows != A->M.rows || R->M.cols != A->M.cols)
    return RABIN_ERR_MATRIX_DIM;

  rmat_t TA = {0}, TB = {0};
  rz_t l, sA, sB;
  rz_init_multi(&l, &sA, &sB, NULL);
  rabin_err_t err = RABIN_SUCCESS;

  TRY(rmat_init(&TA, A->M.rows, A->M.cols));
  TRY(rmat_init(&TB, A->M.rows, A->M.cols));

  TRY(rz_lcm(&l, &A->den, &B->den));

  // find scaling factor
  TRY(rz_div_exact(&sA, &l, &A->den));
  TRY(rz_div_exact(&sB, &l, &B->den));

  // scale the matrices
  TRY(rmat_scalar(&TA, &A->M, &sA));
  TRY(rmat_scalar(&TB, &B->M, &sB));

  TRY(rmat_add(&R->M, &TA, &TB));
  TRY(rz_copy(&R->den, &l));

cleanup:
  rmat_clear(&TB);
  rmat_clear(&TA);
  rz_clear_multi(&l, &sA, &sB, NULL);
  return err;
}

rabin_err_t rmatq_get(rq_t* a, const rmatq_t* A, u64 r, u64 c)
{
  if (a == NULL || A == NULL) return RABIN_ERR_NULL_PTR;
  if (r >= A->M.rows || c >= A->M.cols) return RABIN_ERR_INVALID_ARG;

  rabin_err_t err = rz_copy(&a->num, RMAT_GET(&A->M, r, c));
  if (err != RABIN_SUCCESS) return err;

  err = rz_copy(&a->den, &A->den);
  if (err != RABIN_SUCCESS) return err;

  return rq_normalize(a);
}

rabin_err_t rmatq_set(rmatq_t* A, const rq_t* a, u64 r, u64 c)
{
  if (A == NULL || a == NULL) return RABIN_ERR_NULL_PTR;
  if (r >= A->M.rows || c >= A->M.cols) return RABIN_ERR_INVALID_ARG;
  if (rz_is_zero(&a->den)) return RABIN_ERR_INVALID_ARG;

  rz_t l, s, t;
  rz_init_multi(&l, &s, &t, NULL);
  rabin_err_t err = RABIN_SUCCESS;

  // Grow the common denominator so the new entry fits without rounding.
  TRY(rz_lcm(&l, &A->den, &a->den));

  if (rz_cmp(&l, &A->den) != 0) {
    TRY(rz_div_exact(&s, &l, &A->den));

    for (u64 k = 0; k < A->M.rows * A->M.cols; k++) {
      TRY(rz_mul(&t, &A->M.data[k], &s));
      TRY(rz_copy(&A->M.data[k], &t));
    }

    TRY(rz_copy(&A->den, &l));
  }

  // Store a = a.num / a.den as (a.num * den / a.den) / den.
  TRY(rz_div_exact(&s, &A->den, &a->den));
  TRY(rz_mul(&t, &a->num, &s));
  TRY(rz_copy(RMAT_GET(&A->M, r, c), &t));

  TRY(rmatq_normalize(A));

cleanup:
  rz_clear_multi(&l, &s, &t, NULL);
  return err;
}

rabin_err_t rmatq_set_zero(rmatq_t* A, u64 r, u64 c)
{
  if (A == NULL) return RABIN_ERR_NULL_PTR;
  if (r >= A->M.rows || c >= A->M.cols) return RABIN_ERR_INVALID_ARG;

  return rz_set_u64(RMAT_GET(&A->M, r, c), 0);
}

rabin_err_t rmatq_print(const rmatq_t* A)
{
  if (A == NULL) return RABIN_ERR_NULL_PTR;

  rz_t temp;
  rz_init(&temp);

  for (u64 i = 0; i < A->M.rows; i++) {
    for (u64 j = 0; j < A->M.cols; j++) {
      rmat_get(&temp, &A->M, i, j);
      rz_print(&temp);
      printf("/");
      rz_print(&A->den);
      printf(" ");
    }
    printf("\n");
  }

  rz_clear(&temp);
  return RABIN_SUCCESS;
}

rabin_err_t rmatq_sub(rmatq_t* R, const rmatq_t* A, const rmatq_t* B)
{
  if (R == NULL || A == NULL || B == NULL) return RABIN_ERR_NULL_PTR;
  if (A->M.rows != B->M.rows || A->M.cols != B->M.cols)
    return RABIN_ERR_MATRIX_DIM;
  if (R->M.rows != A->M.rows || R->M.cols != A->M.cols)
    return RABIN_ERR_MATRIX_DIM;

  rmat_t TA = {0}, TB = {0};
  rz_t l, sA, sB;
  rz_init_multi(&l, &sA, &sB, NULL);
  rabin_err_t err = RABIN_SUCCESS;

  TRY(rmat_init(&TA, A->M.rows, A->M.cols));
  TRY(rmat_init(&TB, A->M.rows, A->M.cols));

  TRY(rz_lcm(&l, &A->den, &B->den));

  // find scaling factor
  TRY(rz_div_exact(&sA, &l, &A->den));
  TRY(rz_div_exact(&sB, &l, &B->den));

  // scale the matrices
  TRY(rmat_scalar(&TA, &A->M, &sA));
  TRY(rmat_scalar(&TB, &B->M, &sB));

  TRY(rmat_sub(&R->M, &TA, &TB));
  TRY(rz_copy(&R->den, &l));

cleanup:
  rmat_clear(&TB);
  rmat_clear(&TA);
  rz_clear_multi(&l, &sA, &sB, NULL);
  return err;
}

rabin_err_t rmatq_normalize(rmatq_t* M)
{
  if (M == NULL) return RABIN_ERR_NULL_PTR;

  rz_t d;
  rz_init(&d);

  // Compute GCD of m->den and all entries in m->num
  rabin_err_t err = RABIN_SUCCESS;
  TRY(rmat_gcd_all(&d, &M->M));
  TRY(rz_gcd(&d, &d, &M->den));

  // If GCD > 1, divide both numerator matrix entries and denominator by g
  if (!rz_is_eq_i64(&d, 1)) {
    TRY(rmat_div_exact_scalar(&M->M, &M->M, &d));
    TRY(rz_div_exact(&M->den, &M->den, &d));
  }

  err = RABIN_SUCCESS;
cleanup:
  rz_clear(&d);
  return err;
}

rabin_err_t rmatq_mul(rmatq_t* R, const rmatq_t* A, const rmatq_t* B)
{
  if (R == NULL || A == NULL || B == NULL) return RABIN_ERR_NULL_PTR;
  if (A->M.cols != B->M.rows) return RABIN_ERR_MATRIX_DIM;
  if (R->M.rows != A->M.rows || R->M.cols != B->M.cols)
    return RABIN_ERR_MATRIX_DIM;

  rabin_err_t err = RABIN_SUCCESS;
  rz_t t;
  rz_init(&t);

  for (u64 i = 0; i < R->M.rows; i++) {
    for (u64 j = 0; j < R->M.cols; j++) {
      TRY(rz_set_u64(RMAT_GET(&R->M, i, j), 0));
      for (u64 k = 0; k < A->M.cols; k++) {
        TRY(rz_mul(&t, RMAT_GET(&A->M, i, k), RMAT_GET(&B->M, k, j)));
        TRY(rz_add(RMAT_GET(&R->M, i, j), RMAT_GET(&R->M, i, j), &t));
      }
    }
  }

  TRY(rz_mul(&R->den, &A->den, &B->den));
  TRY(rmatq_normalize(R));

cleanup:
  rz_clear(&t);
  return err;
}
