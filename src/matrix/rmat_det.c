/*
 * rmat_det.c
 *
 * Determinants and adjugate for rmat_t: RNS-based and Bareiss determinants,
 * characteristic polynomial and unimodular inversion.
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
#include "../../include/rmatq.h"
#include "../../include/rns.h"
#include "../../include/rns_primes.h"
#include "../../include/rpol.h"
#include "../../include/rzmath.h"
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

rabin_err_t rmat_inv_unimodular(rmat_t* B, const rmat_t* A)
{
  if (B == NULL || A == NULL) return RABIN_ERR_NULL_PTR;
  if (!rmat_is_square(A)) return RABIN_ERR_MATRIX_DIM;
  if (B->rows != A->rows || B->cols != A->cols) return RABIN_ERR_MATRIX_DIM;

  const u64 n = A->rows;
  if (n != 0 && n > SIZE_MAX / 2 / n) return RABIN_ERR_OVERFLOW;
  if (n == 0) return RABIN_SUCCESS;

  rz_t det;
  rz_init(&det);
  rmat_t W = {0};
  rmat_t P = {0};
  rabin_err_t err = RABIN_SUCCESS;

  TRY(rmat_det(&det, A));
  if (!rz_is_eq_i64(&det, 1) && !rz_is_eq_i64(&det, -1)) {
    err = RABIN_ERR_INVALID_ARG;
    goto cleanup;
  }

  // Augmented system [A | I]; fraction-free Gauss-Jordan keeps every
  // intermediate value integral. Row i of the right block ends equal to
  // W[i][i] * (A^-1)[i][:], so each row is divided by its pivot at the end.
  TRY(rmat_init(&W, n, 2 * n));
  for (u64 i = 0; i < n; i++) {
    for (u64 j = 0; j < n; j++) {
      TRY(rz_copy(RMAT_GET(&W, i, j), RMAT_GET(A, i, j)));
      TRY(rz_set_u64(RMAT_GET(&W, i, n + j), i == j ? 1 : 0));
    }
  }

  rz_t prev, piv, fac;
  rz_init_multi(&prev, &piv, &fac, NULL);
  TRY(rz_set_u64(&prev, 1));

  bool fail = false;
  for (u64 k = 0; k < n && !fail; k++) {
    u64 p = k;
    while (p < n && rz_is_zero(RMAT_GET(&W, p, k))) p++;
    if (p == n) {
      err = RABIN_ERR_INVALID_ARG;
      goto cleanup2;
    }
    if (p != k) {
      for (u64 j = 0; j < 2 * n; j++) {
        TRY(rz_swap(RMAT_GET(&W, k, j), RMAT_GET(&W, p, j)));
      }
    }

    TRY(rz_copy(&piv, RMAT_GET(&W, k, k)));
    for (u64 i = 0; i < n; i++) {
      if (i == k) continue;
      TRY(rz_copy(&fac, RMAT_GET(&W, i, k)));

      rz_t t1, t2;
      rz_init_multi(&t1, &t2, NULL);
      for (u64 j = 0; j < 2 * n; j++) {
        TRY(rz_mul(&t1, &piv, RMAT_GET(&W, i, j)));
        TRY(rz_mul(&t2, &fac, RMAT_GET(&W, k, j)));
        TRY(rz_sub(RMAT_GET(&W, i, j), &t1, &t2));
        TRY(rz_div_exact(RMAT_GET(&W, i, j), RMAT_GET(&W, i, j), &prev));
      }
      rz_clear_multi(&t1, &t2, NULL);
      TRY(rz_set_u64(RMAT_GET(&W, i, k), 0));
    }
    TRY(rz_copy(&prev, &piv));
  }

  // Normalize each row by its pivot; a unimodular A ends with pivots ±1.
  TRY(rmat_init(&P, n, n));
  for (u64 i = 0; i < n; i++) {
    if (!rz_is_eq_i64(RMAT_GET(&W, i, i), 1) &&
        !rz_is_eq_i64(RMAT_GET(&W, i, i), -1)) {
      err = RABIN_ERR_INVALID_ARG;
      goto cleanup2;
    }
    const bool flip = rz_is_eq_i64(RMAT_GET(&W, i, i), -1);
    for (u64 j = 0; j < n; j++) {
      TRY(rz_copy(RMAT_GET(&P, i, j), RMAT_GET(&W, i, n + j)));
      if (flip) TRY(rz_neg(RMAT_GET(&P, i, j), RMAT_GET(&P, i, j)));
    }
  }

  // Guarantee correctness even against elimination bugs: A * P = I.
  rmat_t V = {0};
  TRY(rmat_init(&V, n, n));
  TRY(rmat_mul(&V, A, &P));
  if (!rmat_is_identity(&V)) {
    err = RABIN_ERR_INVALID_ARG;
    goto cleanup2;
  }

  TRY(rmat_copy(B, &P));

cleanup2:
  rz_clear_multi(&prev, &piv, &fac, NULL);
  rmat_clear(&P);
  rmat_clear(&V);
  rmat_clear(&W);
cleanup:
  rz_clear(&det);
  return err;
}
