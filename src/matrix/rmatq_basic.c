/*
 * rmatq_basic.c
 *
 * Constructors and elementwise building blocks for rmatq_t: in-place
 * constants, diagonal construction, row/column swaps, stacking, scalar
 * multiplication, powers and rational matrix predicates.
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
#include "../../include/rvec.h"

/* Cross-multiplied entry compare: A[r][c] == B[r][c] without touching
 * the common denominators of either matrix. */
static bool rmatq_entry_eq(const rmatq_t* A, const rmatq_t* B, u64 r, u64 c)
{
  rz_t lhs = {0};
  rz_t rhs = {0};
  bool eq = false;
  if (rz_init(&lhs) != RABIN_SUCCESS) return false;
  if (rz_init(&rhs) != RABIN_SUCCESS) {
    rz_clear(&lhs);
    return false;
  }
  if (rz_mul(&lhs, RMAT_GET(&A->M, r, c), &B->den) != RABIN_SUCCESS) {
    goto done;
  }
  if (rz_mul(&rhs, RMAT_GET(&B->M, r, c), &A->den) != RABIN_SUCCESS) {
    goto done;
  }
  eq = (rz_cmp(&lhs, &rhs) == 0);
done:
  rz_clear(&lhs);
  rz_clear(&rhs);
  return eq;
}

rabin_err_t rmatq_zero(rmatq_t* A)
{
  if (A == NULL) return RABIN_ERR_NULL_PTR;

  rabin_err_t err = RABIN_SUCCESS;
  for (u64 i = 0; i < A->M.rows && err == RABIN_SUCCESS; i++) {
    for (u64 j = 0; j < A->M.cols && err == RABIN_SUCCESS; j++) {
      err = rz_set_u64(RMAT_GET(&A->M, i, j), 0);
    }
  }
  if (err == RABIN_SUCCESS) err = rz_set_u64(&A->den, 1);
  return err;
}

rabin_err_t rmatq_ones(rmatq_t* A)
{
  if (A == NULL) return RABIN_ERR_NULL_PTR;

  rabin_err_t err = RABIN_SUCCESS;
  for (u64 i = 0; i < A->M.rows && err == RABIN_SUCCESS; i++) {
    for (u64 j = 0; j < A->M.cols && err == RABIN_SUCCESS; j++) {
      err = rz_set_u64(RMAT_GET(&A->M, i, j), 1);
    }
  }
  if (err == RABIN_SUCCESS) err = rz_set_u64(&A->den, 1);
  return err;
}

rabin_err_t rmatq_id(rmatq_t* I, u64 n)
{
  if (I == NULL) return RABIN_ERR_NULL_PTR;
  if (I->M.rows != n || I->M.cols != n) return RABIN_ERR_MATRIX_DIM;

  rabin_err_t err = RABIN_SUCCESS;
  for (u64 i = 0; i < n && err == RABIN_SUCCESS; i++) {
    for (u64 j = 0; j < n && err == RABIN_SUCCESS; j++) {
      err = rz_set_u64(RMAT_GET(&I->M, i, j), i == j ? 1 : 0);
    }
  }
  if (err == RABIN_SUCCESS) err = rz_set_u64(&I->den, 1);
  return err;
}

rabin_err_t rmatq_diag(rmatq_t* D, const rvec_t* v)
{
  if (D == NULL || v == NULL) return RABIN_ERR_NULL_PTR;
  if (D->M.rows != v->size || D->M.cols != v->size) {
    return RABIN_ERR_MATRIX_DIM;
  }

  rabin_err_t err = RABIN_SUCCESS;
  for (u64 i = 0; i < D->M.rows && err == RABIN_SUCCESS; i++) {
    for (u64 j = 0; j < D->M.cols && err == RABIN_SUCCESS; j++) {
      if (i == j) {
        err = rz_copy(RMAT_GET(&D->M, i, j), &v->data[i]);
      } else {
        err = rz_set_u64(RMAT_GET(&D->M, i, j), 0);
      }
    }
  }
  if (err == RABIN_SUCCESS) err = rz_set_u64(&D->den, 1);
  return err;
}

rabin_err_t rmatq_swap_rows(rmatq_t* A, u64 i, u64 j)
{
  if (A == NULL) return RABIN_ERR_NULL_PTR;
  return rmat_swap_rows(&A->M, i, j);
}

rabin_err_t rmatq_swap_cols(rmatq_t* A, u64 i, u64 j)
{
  if (A == NULL) return RABIN_ERR_NULL_PTR;
  return rmat_swap_cols(&A->M, i, j);
}

rabin_err_t rmatq_transpose(rmatq_t* T, const rmatq_t* A)
{
  if (T == NULL || A == NULL) return RABIN_ERR_NULL_PTR;
  if (T->M.rows != A->M.cols || T->M.cols != A->M.rows) {
    return RABIN_ERR_MATRIX_DIM;
  }

  rabin_err_t err = RABIN_SUCCESS;
  rmat_t tmp = {0};

  if (T == A) {
    TRY(rmat_init(&tmp, A->M.rows, A->M.cols));
    TRY(rmat_copy(&tmp, &A->M));
    TRY(rmat_transpose(&T->M, &tmp));
  } else {
    TRY(rmat_transpose(&T->M, &A->M));
  }
  TRY(rz_copy(&T->den, &A->den));

cleanup:
  rmat_clear(&tmp);
  return err;
}

rabin_err_t rmatq_hstack(rmatq_t* R, const rmatq_t* A, const rmatq_t* B)
{
  if (R == NULL || A == NULL || B == NULL) return RABIN_ERR_NULL_PTR;
  if (A->M.rows != B->M.rows) return RABIN_ERR_MATRIX_DIM;
  if (R->M.rows != A->M.rows) return RABIN_ERR_MATRIX_DIM;
  if (A->M.cols > SIZE_MAX - B->M.cols) return RABIN_ERR_OVERFLOW;
  if (R->M.cols != A->M.cols + B->M.cols) return RABIN_ERR_MATRIX_DIM;

  rabin_err_t err = RABIN_SUCCESS;
  rq_t v = {0};
  TRY(rq_init(&v));
  for (u64 i = 0; i < A->M.rows; i++) {
    for (u64 j = 0; j < A->M.cols; j++) {
      TRY(rmatq_get(&v, A, i, j));
      TRY(rmatq_set(R, &v, i, j));
    }
    for (u64 j = 0; j < B->M.cols; j++) {
      TRY(rmatq_get(&v, B, i, j));
      TRY(rmatq_set(R, &v, i, A->M.cols + j));
    }
  }
cleanup:
  rq_clear(&v);
  return err;
}

rabin_err_t rmatq_vstack(rmatq_t* R, const rmatq_t* A, const rmatq_t* B)
{
  if (R == NULL || A == NULL || B == NULL) return RABIN_ERR_NULL_PTR;
  if (A->M.cols != B->M.cols) return RABIN_ERR_MATRIX_DIM;
  if (A->M.rows > SIZE_MAX - B->M.rows) return RABIN_ERR_OVERFLOW;
  if (R->M.rows != A->M.rows + B->M.rows) return RABIN_ERR_MATRIX_DIM;
  if (R->M.cols != A->M.cols) return RABIN_ERR_MATRIX_DIM;

  rabin_err_t err = RABIN_SUCCESS;
  rq_t v = {0};
  TRY(rq_init(&v));
  for (u64 i = 0; i < A->M.rows; i++) {
    for (u64 j = 0; j < A->M.cols; j++) {
      TRY(rmatq_get(&v, A, i, j));
      TRY(rmatq_set(R, &v, i, j));
    }
  }
  for (u64 i = 0; i < B->M.rows; i++) {
    for (u64 j = 0; j < B->M.cols; j++) {
      TRY(rmatq_get(&v, B, i, j));
      TRY(rmatq_set(R, &v, A->M.rows + i, j));
    }
  }
cleanup:
  rq_clear(&v);
  return err;
}

rabin_err_t rmatq_scalar(rmatq_t* R, const rmatq_t* A, const rq_t* a)
{
  if (R == NULL || A == NULL || a == NULL) return RABIN_ERR_NULL_PTR;
  if (R->M.rows != A->M.rows || R->M.cols != A->M.cols) {
    return RABIN_ERR_MATRIX_DIM;
  }

  rabin_err_t err = RABIN_SUCCESS;
  rq_t v = {0};
  rq_t w = {0};
  TRY(rq_init(&v));
  TRY(rq_init(&w));
  for (u64 i = 0; i < A->M.rows; i++) {
    for (u64 j = 0; j < A->M.cols; j++) {
      TRY(rmatq_get(&v, A, i, j));
      TRY(rq_mul(&w, &v, a));
      TRY(rmatq_set(R, &w, i, j));
    }
  }
cleanup:
  rq_clear(&v);
  rq_clear(&w);
  return err;
}

rabin_err_t rmatq_pow(rmatq_t* R, const rmatq_t* A, u64 e)
{
  if (R == NULL || A == NULL) return RABIN_ERR_NULL_PTR;
  if (!rmat_is_square(&A->M)) return RABIN_ERR_MATRIX_DIM;
  if (R->M.rows != A->M.rows || R->M.cols != A->M.cols) {
    return RABIN_ERR_MATRIX_DIM;
  }

  const u64 n = A->M.rows;
  rabin_err_t err = RABIN_SUCCESS;
  rmatq_t base = {0};
  rmatq_t tmp = {0};

  TRY(rmatq_init(&base, n, n));
  TRY(rmatq_copy(&base, A));
  TRY(rmatq_id(R, n));
  if (n == 0) goto cleanup;
  TRY(rmatq_init(&tmp, n, n));

  u64 k = e;
  while (k > 0) {
    if (k & 1) {
      TRY(rmatq_mul(&tmp, R, &base));
      TRY(rmatq_copy(R, &tmp));
    }
    k >>= 1;
    if (k > 0) {
      TRY(rmatq_mul(&tmp, &base, &base));
      TRY(rmatq_copy(&base, &tmp));
    }
  }

cleanup:
  rmatq_clear(&base);
  rmatq_clear(&tmp);
  return err;
}

bool rmatq_equal(const rmatq_t* A, const rmatq_t* B)
{
  if (A == NULL || B == NULL) return false;
  if (A->M.rows != B->M.rows || A->M.cols != B->M.cols) return false;
  for (u64 i = 0; i < A->M.rows; i++) {
    for (u64 j = 0; j < A->M.cols; j++) {
      if (!rmatq_entry_eq(A, B, i, j)) return false;
    }
  }
  return true;
}

bool rmatq_is_square(const rmatq_t* A)
{
  if (A == NULL) return false;
  return rmat_is_square(&A->M);
}

bool rmatq_is_zero(const rmatq_t* A)
{
  if (A == NULL) return false;
  for (u64 i = 0; i < A->M.rows; i++) {
    for (u64 j = 0; j < A->M.cols; j++) {
      if (!rz_is_zero(RMAT_GET(&A->M, i, j))) return false;
    }
  }
  return true;
}

bool rmatq_is_identity(const rmatq_t* A)
{
  if (A == NULL || !rmat_is_square(&A->M)) return false;
  if (rz_is_zero(&A->den)) return false;
  for (u64 i = 0; i < A->M.rows; i++) {
    for (u64 j = 0; j < A->M.cols; j++) {
      if (i == j) {
        if (rz_cmp(RMAT_GET(&A->M, i, j), &A->den) != 0) return false;
      } else if (!rz_is_zero(RMAT_GET(&A->M, i, j))) {
        return false;
      }
    }
  }
  return true;
}

bool rmatq_is_symmetric(const rmatq_t* A)
{
  if (A == NULL || !rmat_is_square(&A->M)) return false;
  for (u64 i = 0; i < A->M.rows; i++) {
    for (u64 j = i + 1; j < A->M.cols; j++) {
      if (rz_cmp(RMAT_GET(&A->M, i, j), RMAT_GET(&A->M, j, i)) != 0) {
        return false;
      }
    }
  }
  return true;
}

rabin_err_t rmatq_trace(rq_t* t, const rmatq_t* A)
{
  if (t == NULL || A == NULL) return RABIN_ERR_NULL_PTR;
  if (!rmat_is_square(&A->M)) return RABIN_ERR_MATRIX_DIM;

  rabin_err_t err = RABIN_SUCCESS;
  rq_t acc = {0};
  rq_t v = {0};
  TRY(rq_init(&acc));
  TRY(rq_init(&v));
  TRY(rq_set_u64(&acc, 0));
  for (u64 i = 0; i < A->M.rows; i++) {
    TRY(rmatq_get(&v, A, i, i));
    TRY(rq_add(&acc, &acc, &v));
  }
  TRY(rq_copy(t, &acc));

cleanup:
  rq_clear(&acc);
  rq_clear(&v);
  return err;
}

rabin_err_t rmatq_det(rq_t* d, const rmatq_t* A)
{
  if (d == NULL || A == NULL) return RABIN_ERR_NULL_PTR;
  if (!rmat_is_square(&A->M)) return RABIN_ERR_MATRIX_DIM;

  const u64 n = A->M.rows;
  rabin_err_t err = RABIN_SUCCESS;
  rz_t num = {0};
  rz_t pown = {0};
  rz_t e = {0};
  rq_t rn = {0};
  rq_t rd = {0};

  TRY(rz_init(&num));
  TRY(rz_init(&pown));
  TRY(rz_init(&e));
  TRY(rq_init(&rn));
  TRY(rq_init(&rd));

  TRY(rmat_det_bareiss(&num, &A->M));
  if (n == 0) {
    TRY(rq_set_u64(d, 1));
    goto cleanup;
  }
  TRY(rz_set_u64(&e, n));
  TRY(rz_pow(&pown, &A->den, &e));
  TRY(rq_set(&rn, &num));
  TRY(rq_set(&rd, &pown));
  TRY(rq_div(d, &rn, &rd));

cleanup:
  rz_clear(&num);
  rz_clear(&pown);
  rz_clear(&e);
  rq_clear(&rn);
  rq_clear(&rd);
  return err;
}
