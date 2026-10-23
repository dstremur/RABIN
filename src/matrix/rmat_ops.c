/*
 * rmat_ops.c
 *
 * Component-wise and scalar operations over rmat_t matrices: scalar mul/div,
 * add/sub, trace, swap, gcd, identity, transpose and equality.
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
#include "../../include/rzmath.h"
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

rabin_err_t rmat_transpose(rmat_t* T, const rmat_t* A)
{
  if (T == NULL || A == NULL) return RABIN_ERR_NULL_PTR;
  if (T->rows != A->cols || T->cols != A->rows) return RABIN_ERR_MATRIX_DIM;
  if (A->rows > SIZE_MAX / (A->cols + 1)) return RABIN_ERR_OVERFLOW;

  rabin_err_t err = RABIN_SUCCESS;

#pragma omp for collapse(2)
  for (u64 i = 0; i < A->rows; i++) {
    for (u64 j = 0; j < A->cols; j++) {
      err = rz_copy(RMAT_GET(T, j, i), RMAT_GET(A, i, j));
    }
  }

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
