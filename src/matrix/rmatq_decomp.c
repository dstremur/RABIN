/*
 * rmatq_decomp.c
 *
 * Rational matrix decompositions: exact Cholesky factorization and inversion.
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

#include <string.h>

#include "../../include/rmatq.h"
/**
 * Exact rational square root: sets \f$r = \sqrt a\f$ when both \f$a.num\f$
 * and \f$a.den\f$ are perfect squares, and fails otherwise. Internal to
 * rmatq_cholesky(), which never calls it with a negative or zero value.
 */
static rabin_err_t rq_sqrt_exact(rq_t* r, const rq_t* a)
{
  if (r == NULL || a == NULL) return RABIN_ERR_NULL_PTR;
  if (rz_is_zero(&a->den) || a->num.is_neg) return RABIN_ERR_INVALID_ARG;

  rz_t sn, sd;
  rz_init_multi(&sn, &sd, NULL);

  rabin_err_t err = RABIN_SUCCESS;

  if (!rz_is_square(&sn, &a->num) || !rz_is_square(&sd, &a->den)) {
    err = RABIN_ERR_INVALID_ARG;
    goto cleanup;
  }

  TRY(rz_copy(&r->num, &sn));
  TRY(rz_copy(&r->den, &sd));

cleanup:
  rz_clear_multi(&sn, &sd, NULL);
  return err;
}

rabin_err_t rmatq_cholesky(rmatq_t* C, const rmatq_t* A)
{
  if (C == NULL || A == NULL) return RABIN_ERR_NULL_PTR;
  if (A->M.rows != A->M.cols) return RABIN_ERR_MATRIX_DIM;
  if (C->M.rows != A->M.rows || C->M.cols != A->M.cols)
    return RABIN_ERR_MATRIX_DIM;

  const u64 n = A->M.rows;
  if (n > SIZE_MAX / n) return RABIN_ERR_OVERFLOW;
  if (n == 0) return rz_set_u64(&C->den, 1);

  // Scratch upper triangle Q of rationals, row-major n x n.
  rq_t* Q = malloc(n * n * sizeof(rq_t));
  if (Q == NULL) return RABIN_ERR_OUT_OF_MEMORY;

  rq_t tmp;
  rz_t D, s, t;
  rq_init(&tmp);
  rz_init_multi(&D, &s, &t, NULL);

  rabin_err_t err = RABIN_SUCCESS;
  u64 inited = 0;

  for (u64 k = 0; k < n * n; k++) {
    TRY(rq_init(&Q[k]));
    inited = k + 1;
  }

  // 1. [Initialize]
  for (u64 i = 0; i < n; i++) {
    for (u64 j = i; j < n; j++) {
      TRY(rmatq_get(&Q[i * n + j], A, i, j));
    }
  }

  // 2. [Loop on i]
  for (u64 i = 0; i < n; i++) {
    rq_t* pivot = &Q[i * n + i];

    if (rz_is_zero(&pivot->num)) {
      err = RABIN_ERR_DIV_BY_ZERO;
      goto cleanup;
    }
    if (pivot->num.is_neg) {
      err = RABIN_ERR_INVALID_ARG;
      goto cleanup;
    }

    for (u64 j = i + 1; j < n; j++) {
      // q_{j,i} <- q_{i,j}, then q_{i,j} <- q_{i,j} / q_{i,i}
      TRY(rq_copy(&Q[j * n + i], &Q[i * n + j]));
      TRY(rq_div(&Q[i * n + j], &Q[i * n + j], pivot));
    }

    // 3. [Main loop]
    for (u64 k = i + 1; k < n; k++) {
      if (rz_is_zero(&Q[k * n + i].num)) continue;

      for (u64 l = k; l < n; l++) {
        if (rz_is_zero(&Q[i * n + l].num)) continue;

        // q_{k,l} <- q_{k,l} - q_{k,i} * q_{i,l}
        TRY(rq_mul(&tmp, &Q[k * n + i], &Q[i * n + l]));
        TRY(rq_sub(&Q[k * n + l], &Q[k * n + l], &tmp));
      }
    }
  }

  // 4. [Find matrix R]
  for (u64 i = 0; i < n; i++) {
    TRY(rq_sqrt_exact(&Q[i * n + i], &Q[i * n + i]));

    for (u64 j = i + 1; j < n; j++) {
      TRY(rq_mul(&Q[i * n + j], &Q[i * n + j], &Q[i * n + i]));
    }
  }

  // Flush R into C on a single common denominator; C is untouched on error.
  TRY(rz_set_u64(&D, 1));
  for (u64 i = 0; i < n; i++) {
    for (u64 j = i; j < n; j++) {
      TRY(rz_lcm(&D, &D, &Q[i * n + j].den));
    }
  }

  for (u64 i = 0; i < n; i++) {
    for (u64 j = 0; j < n; j++) {
      if (j < i) {
        TRY(rz_set_u64(RMAT_GET(&C->M, i, j), 0));
        continue;
      }

      TRY(rz_div_exact(&s, &D, &Q[i * n + j].den));
      TRY(rz_mul(&t, &Q[i * n + j].num, &s));
      TRY(rz_copy(RMAT_GET(&C->M, i, j), &t));
    }
  }

  TRY(rz_copy(&C->den, &D));
  TRY(rmatq_normalize(C));

cleanup:
  if (Q != NULL) {
    for (u64 k = 0; k < inited; k++) rq_clear(&Q[k]);
    free(Q);
  }
  rq_clear(&tmp);
  rz_clear_multi(&D, &s, &t, NULL);
  return err;
}

rabin_err_t rmatq_inv(rmatq_t* B, const rmatq_t* A)
{
  if (B == NULL || A == NULL) return RABIN_ERR_NULL_PTR;
  if (A->M.rows != A->M.cols) return RABIN_ERR_MATRIX_DIM;
  if (B->M.rows != A->M.rows || B->M.cols != A->M.cols)
    return RABIN_ERR_MATRIX_DIM;

  const u64 n = A->M.rows;
  if (n != 0 && 2 * n > SIZE_MAX / n) return RABIN_ERR_OVERFLOW;
  if (n == 0) return rz_set_u64(&B->den, 1);

  rmatq_t P = {0};
  rz_t D = {0}, s = {0}, t = {0};

  // Augmented system [A | I] over the rationals.
  rq_t* W = malloc(n * 2 * n * sizeof(rq_t));
  if (W == NULL) return RABIN_ERR_OUT_OF_MEMORY;

  rq_t tmp, fac;
  rq_init_multi(&tmp, &fac, NULL);

  rabin_err_t err = RABIN_SUCCESS;
  u64 inited = 0;
  for (u64 k = 0; k < 2 * n * n; k++) {
    TRY(rq_init(&W[k]));
    inited = k + 1;
  }

  for (u64 i = 0; i < n; i++) {
    for (u64 j = 0; j < n; j++) {
      TRY(rmatq_get(&W[i * 2 * n + j], A, i, j));
    }
    for (u64 j = 0; j < n; j++) {
      TRY(rq_set_u64(&W[i * 2 * n + n + j], i == j ? 1 : 0));
    }
  }

  for (u64 k = 0; k < n; k++) {
    u64 p = k;
    while (p < n && rz_is_zero(&W[p * 2 * n + k].num)) p++;
    if (p == n) {
      err = RABIN_ERR_DIV_BY_ZERO;
      goto cleanup;
    }
    if (p != k) {
      for (u64 j = 0; j < 2 * n; j++) {
        TRY(rq_copy(&tmp, &W[k * 2 * n + j]));
        TRY(rq_copy(&W[k * 2 * n + j], &W[p * 2 * n + j]));
        TRY(rq_copy(&W[p * 2 * n + j], &tmp));
      }
    }

    TRY(rq_copy(&fac, &W[k * 2 * n + k]));
    for (u64 j = 0; j < 2 * n; j++) {
      TRY(rq_div(&W[k * 2 * n + j], &W[k * 2 * n + j], &fac));
    }

    for (u64 i = 0; i < n; i++) {
      if (i == k) continue;
      TRY(rq_copy(&fac, &W[i * 2 * n + k]));
      if (rq_is_zero(&fac)) continue;
      for (u64 j = 0; j < 2 * n; j++) {
        if (rz_is_zero(&W[k * 2 * n + j].num)) continue;
        TRY(rq_mul(&tmp, &fac, &W[k * 2 * n + j]));
        TRY(rq_sub(&W[i * 2 * n + j], &W[i * 2 * n + j], &tmp));
      }
    }
  }

  // Flush the right block into B on one common denominator.
  TRY(rz_set_u64(&D, 1));
  for (u64 i = 0; i < n; i++) {
    for (u64 j = 0; j < n; j++) {
      TRY(rz_lcm(&D, &D, &W[i * 2 * n + n + j].den));
    }
  }
  for (u64 i = 0; i < n; i++) {
    for (u64 j = 0; j < n; j++) {
      const rq_t* e = &W[i * 2 * n + n + j];
      TRY(rz_div_exact(&s, &D, &e->den));
      TRY(rz_mul(&t, &e->num, &s));
      TRY(rz_copy(RMAT_GET(&B->M, i, j), &t));
    }
  }
  TRY(rz_copy(&B->den, &D));
  TRY(rmatq_normalize(B));

  // Verify B * A = I before trusting the result.
  TRY(rmatq_init(&P, n, n));
  TRY(rmatq_mul(&P, B, A));
  for (u64 i = 0; i < n && err == RABIN_SUCCESS; i++) {
    for (u64 j = 0; j < n; j++) {
      const int eq = rz_cmp(RMAT_GET(&P.M, i, j), &P.den);
      if ((i == j && eq != 0) ||
          (i != j && !rz_is_zero(RMAT_GET(&P.M, i, j)))) {
        err = RABIN_ERR_INVALID_ARG;
        break;
      }
    }
  }

cleanup:
  if (W != NULL) {
    for (u64 k = 0; k < inited; k++) rq_clear(&W[k]);
    free(W);
  }
  rmatq_clear(&P);
  rz_clear_multi(&D, &s, &t, NULL);
  rq_clear_multi(&tmp, &fac, NULL);
  return err;
}
