/*
 * rmatq_elim.c
 *
 * Exact linear algebra over Q for rmatq_t: rank, reduced row echelon
 * form, LU decomposition with row pivoting and linear system solving.
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

rabin_err_t rmatq_rank(u64* rank, const rmatq_t* A)
{
  if (rank == NULL || A == NULL) return RABIN_ERR_NULL_PTR;
  /* Scaling every entry by the common denominator is invertible over Q,
   * so the numerator matrix has the same rank. */
  return rmat_rank(rank, &A->M);
}

rabin_err_t rmatq_rref(rmatq_t* R, u64* rank)
{
  if (R == NULL) return RABIN_ERR_NULL_PTR;

  const u64 m = R->M.rows;
  const u64 n = R->M.cols;
  rabin_err_t err = RABIN_SUCCESS;
  rq_t piv = {0};
  rq_t f = {0};
  rq_t u = {0};
  rq_t t = {0};
  rq_t w = {0};
  rq_t z = {0};

  TRY(rq_init(&piv));
  TRY(rq_init(&f));
  TRY(rq_init(&u));
  TRY(rq_init(&t));
  TRY(rq_init(&w));
  TRY(rq_init(&z));
  TRY(rq_set_u64(&z, 0));

  u64 r = 0;
  for (u64 col = 0; col < n && r < m; col++) {
    u64 pr = 0;
    bool found = false;
    for (u64 i = r; i < m; i++) {
      TRY(rmatq_get(&u, R, i, col));
      if (!rq_is_zero(&u)) {
        pr = i;
        found = true;
        break;
      }
    }
    if (!found) continue;

    if (pr != r) {
      TRY(rmatq_swap_rows(R, pr, r));
    }

    TRY(rmatq_get(&piv, R, r, col));
    for (u64 j = col; j < n; j++) {
      TRY(rmatq_get(&u, R, r, j));
      TRY(rq_div(&w, &u, &piv));
      TRY(rmatq_set(R, &w, r, j));
    }

    for (u64 i = 0; i < m; i++) {
      if (i == r) continue;
      TRY(rmatq_get(&f, R, i, col));
      if (rq_is_zero(&f)) continue;
      for (u64 j = col; j < n; j++) {
        TRY(rmatq_get(&u, R, i, j));
        TRY(rmatq_get(&t, R, r, j));
        TRY(rq_mul(&w, &f, &t));
        TRY(rq_sub(&u, &u, &w));
        TRY(rmatq_set(R, &u, i, j));
      }
      TRY(rmatq_set(R, &z, i, col));
    }
    r++;
  }

  if (rank != NULL) *rank = r;

cleanup:
  rq_clear(&piv);
  rq_clear(&f);
  rq_clear(&u);
  rq_clear(&t);
  rq_clear(&w);
  rq_clear(&z);
  return err;
}

rabin_err_t rmatq_lu(rmatq_t* L, rmatq_t* U, u64* P, const rmatq_t* A)
{
  if (L == NULL || U == NULL || P == NULL || A == NULL) {
    return RABIN_ERR_NULL_PTR;
  }
  if (!rmat_is_square(&A->M)) return RABIN_ERR_MATRIX_DIM;
  if (L->M.rows != A->M.rows || L->M.cols != A->M.cols) {
    return RABIN_ERR_MATRIX_DIM;
  }
  if (U->M.rows != A->M.rows || U->M.cols != A->M.cols) {
    return RABIN_ERR_MATRIX_DIM;
  }

  const u64 n = A->M.rows;
  rabin_err_t err = RABIN_SUCCESS;
  u64* perm = NULL;
  rmatq_t W = {0};
  rq_t piv = {0};
  rq_t f = {0};
  rq_t m_ik = {0};
  rq_t u = {0};
  rq_t t = {0};
  rq_t w = {0};
  rq_t z = {0};

  perm = calloc(n, sizeof(u64));
  if (perm == NULL) return RABIN_ERR_OUT_OF_MEMORY;
  for (u64 i = 0; i < n; i++) perm[i] = i;

  TRY(rmatq_init(&W, n, n));
  TRY(rmatq_copy(&W, A));
  TRY(rmatq_id(L, n));
  TRY(rq_init(&piv));
  TRY(rq_init(&f));
  TRY(rq_init(&m_ik));
  TRY(rq_init(&u));
  TRY(rq_init(&t));
  TRY(rq_init(&w));
  TRY(rq_init(&z));
  TRY(rq_set_u64(&z, 0));

  for (u64 k = 0; k < n; k++) {
    u64 pr = 0;
    bool found = false;
    for (u64 i = k; i < n; i++) {
      TRY(rmatq_get(&u, &W, i, k));
      if (!rq_is_zero(&u)) {
        pr = i;
        found = true;
        break;
      }
    }
    if (!found) continue;

    if (pr != k) {
      TRY(rmatq_swap_rows(&W, pr, k));
      TRY(rmatq_swap_rows(L, pr, k));
      const u64 tmp = perm[pr];
      perm[pr] = perm[k];
      perm[k] = tmp;
    }

    TRY(rmatq_get(&piv, &W, k, k));
    for (u64 i = k + 1; i < n; i++) {
      TRY(rmatq_get(&f, &W, i, k));
      if (rq_is_zero(&f)) continue;
      TRY(rq_div(&m_ik, &f, &piv));
      TRY(rmatq_set(L, &m_ik, i, k));
      for (u64 j = k; j < n; j++) {
        TRY(rmatq_get(&u, &W, i, j));
        TRY(rmatq_get(&t, &W, k, j));
        TRY(rq_mul(&w, &m_ik, &t));
        TRY(rq_sub(&u, &u, &w));
        TRY(rmatq_set(&W, &u, i, j));
      }
      TRY(rmatq_set(&W, &z, i, k));
    }
  }

  for (u64 i = 0; i < n; i++) {
    P[i] = perm[i];
    for (u64 j = 0; j < n; j++) {
      TRY(rmatq_get(&u, &W, i, j));
      TRY(rmatq_set(U, &u, i, j));
    }
  }

cleanup:
  free(perm);
  rmatq_clear(&W);
  rq_clear(&piv);
  rq_clear(&f);
  rq_clear(&m_ik);
  rq_clear(&u);
  rq_clear(&t);
  rq_clear(&w);
  rq_clear(&z);
  return err;
}

rabin_err_t rmatq_solve(rmatq_t* X, const rmatq_t* A, const rmatq_t* B)
{
  if (X == NULL || A == NULL || B == NULL) return RABIN_ERR_NULL_PTR;
  if (!rmat_is_square(&A->M)) return RABIN_ERR_MATRIX_DIM;
  if (B->M.rows != A->M.rows) return RABIN_ERR_MATRIX_DIM;
  if (X->M.rows != A->M.rows || X->M.cols != B->M.cols) {
    return RABIN_ERR_MATRIX_DIM;
  }

  const u64 n = A->M.rows;
  const u64 k = B->M.cols;
  if (n == 0 || k == 0) return RABIN_SUCCESS;

  rabin_err_t err = RABIN_SUCCESS;
  rmatq_t C = {0};
  rq_t v = {0};

  TRY(rmatq_init(&C, n, n + k));
  TRY(rmatq_hstack(&C, A, B));

  u64 rank_aug = 0;
  TRY(rmatq_rref(&C, &rank_aug));

  u64 rank_a = 0;
  TRY(rmatq_rank(&rank_a, A));

  if (rank_aug > rank_a) {
    err = RABIN_ERR_INCONSISTENT;
    goto cleanup;
  }
  if (rank_a < n) {
    err = RABIN_ERR_SINGULAR;
    goto cleanup;
  }

  TRY(rq_init(&v));
  for (u64 i = 0; i < n; i++) {
    for (u64 j = 0; j < k; j++) {
      TRY(rmatq_get(&v, &C, i, n + j));
      TRY(rmatq_set(X, &v, i, j));
    }
  }

cleanup:
  rmatq_clear(&C);
  rq_clear(&v);
  return err;
}
