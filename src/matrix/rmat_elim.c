/*
 * rmat_elim.c
 *
 * Elimination core and derived functions for rmat_t: rank, reduced row
 * echelon form, exact integer LU with row pivoting, rational linear
 * systems, cofactor/adjugate and singularity tests.
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

#include "../../include/rmat.h"
#include "../../include/rmatq.h"
#include "../../include/rq.h"

/*
 * Division-free column elimination on W (m rows, W->cols columns; the
 * last `extra = W->cols - ncols` columns are carried along unchanged as
 * side blocks). For each pivot column a nonzero entry of smallest
 * absolute value is swapped to the pivot row and the entries below are
 * reduced with integer row operations only (Euclid descent), so every
 * value stays integral and the pivot ends as the gcd of the column
 * section. perm (optional, size rows) tracks the row permutation; the
 * pivot row count is stored in rank when non-NULL.
 */
static rabin_err_t rmat_elim_core(rmat_t* W, u64* perm, u64 ncols, u64* rank)
{
  if (ncols == 0 || W->rows == 0) {
    if (rank != NULL) *rank = 0;
    return RABIN_SUCCESS;
  }

  rabin_err_t err = RABIN_SUCCESS;
  rz_t q, t, best;
  rz_init_multi(&q, &t, &best, NULL);

  u64 r = 0;
  for (u64 col = 0; col < ncols && r < W->rows; col++) {
    for (;;) {
      // pick the smallest nonzero entry of the column below row r
      u64 p = W->rows;
      bool have = false;
      for (u64 i = r; i < W->rows; i++) {
        const rz_t* e = RMAT_GET(W, i, col);
        if (rz_is_zero(e)) continue;
        if (!have || rz_cmp_abs(&best, e) > 0) {
          TRY(rz_copy(&best, e));
          p = i;
          have = true;
        }
      }
      if (!have) break;  // empty column: advance
      if (p != r) {
        for (u64 j = 0; j < W->cols; j++) {
          TRY(rz_swap(RMAT_GET(W, r, j), RMAT_GET(W, p, j)));
        }
        if (perm != NULL) {
          const u64 s = perm[r];
          perm[r] = perm[p];
          perm[p] = s;
        }
      }

      bool reduced = false;
      for (u64 i = r + 1; i < W->rows; i++) {
        if (rz_is_zero(RMAT_GET(W, i, col))) continue;
        TRY(rz_divmod(&q, RMAT_GET(W, i, col), RMAT_GET(W, i, col),
                      RMAT_GET(W, r, col)));
        for (u64 j = col + 1; j < W->cols; j++) {
          TRY(rz_mul(&t, &q, RMAT_GET(W, r, j)));
          TRY(rz_sub(RMAT_GET(W, i, j), RMAT_GET(W, i, j), &t));
        }
        reduced = true;
      }
      if (!reduced) {
        r++;
        break;
      }
    }
  }

  if (rank != NULL) *rank = r;

cleanup:
  rz_clear_multi(&q, &t, &best, NULL);
  return err;
}

rabin_err_t rmat_rank(u64* rank, const rmat_t* A)
{
  if (rank == NULL || A == NULL) return RABIN_ERR_NULL_PTR;

  rabin_err_t err = RABIN_SUCCESS;
  rmat_t W = {0};

  TRY(rmat_init(&W, A->rows, A->cols));
  for (u64 i = 0; i < A->rows; i++) {
    for (u64 j = 0; j < A->cols; j++) {
      TRY(rz_copy(RMAT_GET(&W, i, j), RMAT_GET(A, i, j)));
    }
  }
  TRY(rmat_elim_core(&W, NULL, A->cols, rank));

cleanup:
  rmat_clear(&W);
  return err;
}

rabin_err_t rmat_rref(rmatq_t* R, u64* rank, const rmat_t* A)
{
  if (R == NULL || A == NULL) return RABIN_ERR_NULL_PTR;
  if (R->M.rows != A->rows || R->M.cols != A->cols) {
    return RABIN_ERR_MATRIX_DIM;
  }

  const u64 m = A->rows;
  const u64 n = A->cols;
  if (m == 0 || n == 0) {
    if (rank != NULL) *rank = 0;
    return RABIN_SUCCESS;
  }
  if (m > SIZE_MAX / n) return RABIN_ERR_OVERFLOW;

  rabin_err_t err = RABIN_SUCCESS;
  rq_t* g = calloc(m * n, sizeof(rq_t));
  if (g == NULL) return RABIN_ERR_OUT_OF_MEMORY;

  rq_t piv, s, factor, t;
  rq_init_multi(&piv, &s, &factor, &t, NULL);
  rq_t zero;
  rq_init(&zero);
  TRY(rq_set_u64(&zero, 0));
  rq_t one;
  rq_init(&one);
  TRY(rq_set_u64(&one, 1));

  for (u64 i = 0; i < m * n; i++) TRY(rq_init(&g[i]));
  for (u64 i = 0; i < m; i++) {
    for (u64 j = 0; j < n; j++) {
      TRY(rq_set(&g[i * n + j], RMAT_GET(A, i, j)));
    }
  }

  u64 r = 0;
  for (u64 col = 0; col < n && r < m; col++) {
    u64 p = m;
    for (u64 i = r; i < m; i++) {
      if (!rq_is_zero(&g[i * n + col])) {
        p = i;
        break;
      }
    }
    if (p == m) continue;

    if (p != r) {
      for (u64 j = 0; j < n; j++) {
        TRY(rq_copy(&t, &g[r * n + j]));
        TRY(rq_copy(&g[r * n + j], &g[p * n + j]));
        TRY(rq_copy(&g[p * n + j], &t));
      }
    }

    TRY(rq_copy(&piv, &g[r * n + col]));
    for (u64 i = 0; i < m; i++) {
      if (i == r) continue;
      if (rq_is_zero(&g[i * n + col])) continue;
      TRY(rq_div(&factor, &g[i * n + col], &piv));
      for (u64 j = col + 1; j < n; j++) {
        TRY(rq_mul(&t, &factor, &g[r * n + j]));
        TRY(rq_sub(&g[i * n + j], &g[i * n + j], &t));
      }
      TRY(rq_copy(&g[i * n + col], &zero));
    }
    for (u64 j = col; j < n; j++) {
      TRY(rq_div(&g[r * n + j], &g[r * n + j], &piv));
    }
    TRY(rq_copy(&g[r * n + col], &one));
    r++;
  }

  for (u64 i = 0; i < m; i++) {
    for (u64 j = 0; j < n; j++) {
      TRY(rmatq_set(R, &g[i * n + j], i, j));
    }
  }
  if (rank != NULL) *rank = r;
  TRY(rmatq_normalize(R));

cleanup:
  if (g != NULL) {
    for (u64 i = 0; i < m * n; i++) rq_clear(&g[i]);
    free(g);
  }
  rq_clear_multi(&piv, &s, &factor, &t, NULL);
  rq_clear(&zero);
  rq_clear(&one);
  return err;
}

rabin_err_t rmat_lu(rmat_t* L, rmat_t* U, u64* P, const rmat_t* A)
{
  if (L == NULL || U == NULL || P == NULL || A == NULL) {
    return RABIN_ERR_NULL_PTR;
  }
  if (!rmat_is_square(A)) return RABIN_ERR_MATRIX_DIM;
  if (L->rows != A->rows || L->cols != A->cols) return RABIN_ERR_MATRIX_DIM;
  if (U->rows != A->rows || U->cols != A->cols) return RABIN_ERR_MATRIX_DIM;

  const u64 n = A->rows;

  rabin_err_t err = RABIN_SUCCESS;
  u64* perm = NULL;
  rmat_t W = {0};  // augmented [A | I]
  rmat_t S = {0};  // E * Perm^-1
  rmat_t Ei = {0};

  if (n == 0) return RABIN_SUCCESS;
  if (n > SIZE_MAX / 2 / n) return RABIN_ERR_OVERFLOW;

  perm = calloc(n, sizeof(u64));
  if (perm == NULL) return RABIN_ERR_OUT_OF_MEMORY;
  for (u64 i = 0; i < n; i++) perm[i] = i;

  TRY(rmat_init(&W, n, 2 * n));
  for (u64 i = 0; i < n; i++) {
    for (u64 j = 0; j < n; j++) {
      TRY(rz_copy(RMAT_GET(&W, i, j), RMAT_GET(A, i, j)));
      TRY(rz_set_u64(RMAT_GET(&W, i, n + j), i == j ? 1 : 0));
    }
  }
  TRY(rmat_elim_core(&W, perm, n, NULL));

  // U is the left block. The right block E satisfies E * A = U.
  // With the permutation matrix
  // Perm[i][perm[i]] = 1 the product S = E * Perm^T is lower triangular
  // with diagonal entries d_i = +/-1, S = D * S~ with S~ unit lower
  // triangular, and S * (Perm * A) = U. Hence
  //     Perm * A = S~^-1 * (D * U) = L * U
  // is an exact unit lower / upper factorization over the integers.
  TRY(rmat_init(&S, n, n));
  TRY(rmat_init(&Ei, n, n));
  for (u64 i = 0; i < n; i++) {
    for (u64 j = 0; j < n; j++) {
      TRY(rz_copy(RMAT_GET(&S, i, j), RMAT_GET(&W, i, n + perm[j])));
    }
  }

  // Ei <- S~^-1 by forward substitution (S~[i][k] = S[i][k] * d_i).
  // S is always invertible (E unimodular), so d_i = +/-1.
  for (u64 i = 0; i < n; i++) {
    TRY(rz_set_u64(RMAT_GET(&Ei, i, i), 1));
    for (u64 j = 0; j < i; j++) {
      rz_t acc;
      rz_init(&acc);
      for (u64 k = j; k < i; k++) {
        rz_t t;
        rz_init(&t);
        TRY(rz_mul(&t, RMAT_GET(&S, i, k), RMAT_GET(&S, i, i)));
        rz_t u;
        rz_init(&u);
        TRY(rz_mul(&u, &t, RMAT_GET(&Ei, k, j)));
        TRY(rz_sub(&acc, &acc, &u));
        rz_clear(&t);
        rz_clear(&u);
      }
      TRY(rz_copy(RMAT_GET(&Ei, i, j), &acc));
      rz_clear(&acc);
    }
  }

  TRY(rmat_copy(L, &Ei));
  for (u64 i = 0; i < n; i++) {
    P[i] = perm[i];
    const bool flip = !rz_is_one(RMAT_GET(&S, i, i));
    for (u64 j = 0; j < n; j++) {
      if (flip) {
        TRY(rz_neg(RMAT_GET(U, i, j), RMAT_GET(&W, i, j)));
      } else {
        TRY(rz_copy(RMAT_GET(U, i, j), RMAT_GET(&W, i, j)));
      }
    }
  }

cleanup:
  free(perm);
  rmat_clear(&W);
  rmat_clear(&S);
  rmat_clear(&Ei);
  return err;
}

rabin_err_t rmat_solve(rmatq_t* X, const rmat_t* A, const rmat_t* B)
{
  if (X == NULL || A == NULL || B == NULL) return RABIN_ERR_NULL_PTR;
  if (!rmat_is_square(A)) return RABIN_ERR_MATRIX_DIM;
  if (B->rows != A->rows) return RABIN_ERR_MATRIX_DIM;
  if (X->M.rows != A->rows || X->M.cols != B->cols) {
    return RABIN_ERR_MATRIX_DIM;
  }

  const u64 n = A->rows;
  const u64 k = B->cols;
  if (n == 0 || k == 0) return RABIN_SUCCESS;
  if (n > SIZE_MAX / (n + k + 1)) return RABIN_ERR_OVERFLOW;

  rabin_err_t err = RABIN_SUCCESS;
  rmat_t C = {0};
  rmatq_t R = {0};

  TRY(rmat_init(&C, n, n + k));
  TRY(rmat_hstack(&C, A, B));
  TRY(rmatq_init(&R, n, n + k));

  u64 rank_aug = 0;
  TRY(rmat_rref(&R, &rank_aug, &C));

  u64 rank_a = 0;
  TRY(rmat_rank(&rank_a, A));

  if (rank_aug > rank_a) {
    err = RABIN_ERR_INCONSISTENT;
    goto cleanup;
  }
  if (rank_a < n) {
    err = RABIN_ERR_SINGULAR;
    goto cleanup;
  }

  rq_t v;
  rq_init(&v);
  for (u64 i = 0; i < n; i++) {
    for (u64 j = 0; j < k; j++) {
      TRY(rmatq_get(&v, &R, i, n + j));
      TRY(rmatq_set(X, &v, i, j));
    }
  }
  rq_clear(&v);

cleanup:
  rmat_clear(&C);
  rmatq_clear(&R);
  return err;
}

rabin_err_t rmat_cofactor(rmat_t* C, const rmat_t* A)
{
  if (C == NULL || A == NULL) return RABIN_ERR_NULL_PTR;
  if (!rmat_is_square(A)) return RABIN_ERR_MATRIX_DIM;
  if (C->rows != A->rows || C->cols != A->cols) return RABIN_ERR_MATRIX_DIM;

  const u64 n = A->rows;
  if (n == 0) return RABIN_SUCCESS;

  rabin_err_t err = RABIN_SUCCESS;
  rmat_t M = {0};
  rz_t d;
  rz_init(&d);

  if (n > 1) TRY(rmat_init(&M, n - 1, n - 1));

  for (u64 i = 0; i < n; i++) {
    for (u64 j = 0; j < n; j++) {
      if (n == 1) {
        TRY(rz_set_u64(RMAT_GET(C, i, j), 1));
        continue;
      }
      for (u64 r = 0; r < n - 1; r++) {
        for (u64 c = 0; c < n - 1; c++) {
          TRY(rz_copy(RMAT_GET(&M, r, c),
                      RMAT_GET(A, r < i ? r : r + 1, c < j ? c : c + 1)));
        }
      }
      TRY(rmat_det_bareiss(&d, &M));
      if ((i + j) & 1) TRY(rz_neg(&d, &d));
      TRY(rz_copy(RMAT_GET(C, i, j), &d));
    }
  }

cleanup:
  rmat_clear(&M);
  rz_clear(&d);
  return err;
}

rabin_err_t rmat_adj(rmat_t* ADJ, const rmat_t* A)
{
  if (ADJ == NULL || A == NULL) return RABIN_ERR_NULL_PTR;
  if (!rmat_is_square(A)) return RABIN_ERR_MATRIX_DIM;
  if (ADJ->rows != A->cols || ADJ->cols != A->rows) {
    return RABIN_ERR_MATRIX_DIM;
  }

  rabin_err_t err = RABIN_SUCCESS;
  rmat_t C = {0};

  TRY(rmat_init(&C, A->rows, A->cols));
  TRY(rmat_cofactor(&C, A));
  TRY(rmat_transpose(ADJ, &C));

cleanup:
  rmat_clear(&C);
  return err;
}

bool rmat_is_singular(const rmat_t* A)
{
  if (A == NULL) return false;
  if (!rmat_is_square(A)) return false;

  bool singular = true;
  rz_t d;
  rz_init(&d);
  if (rmat_det(&d, A) == RABIN_SUCCESS) singular = rz_is_zero(&d);
  rz_clear(&d);
  return singular;
}
