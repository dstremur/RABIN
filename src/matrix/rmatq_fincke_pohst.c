/*
 * rmatq_fincke_pohst.c
 *
 * Fincke-Pohst enumeration for shortest vectors of rational quadratic forms.
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

#include <math.h>
#include <string.h>

#include "../../include/rmatq.h"
#include "../../include/rzmath.h"
/* floor of the square root of a non-negative rational. */
static rabin_err_t rq_sqrt_floor(rz_t* m, const rq_t* r)
{
  if (rz_is_zero(&r->num)) return rz_set_u64(m, 0);
  if (r->num.is_neg) return rz_set_u64(m, 0);

  rz_t t;
  rz_init(&t);
  rabin_err_t err = RABIN_SUCCESS;
  TRY(rz_mul(&t, &r->num, &r->den));
  TRY(rz_isqrt(m, &t));
  TRY(rz_div(m, m, &r->den));
cleanup:
  rz_clear(&t);
  return err;
}

/*
 * m = floor(sqrt(r) + c) for r >= 0 and rational c. With Zf = floor(sqrt r)
 * and uf = floor(c), floor(sqrt r + c) lies in {Zf+uf, Zf+uf+1}, so a single
 * comparison decides: the candidate cand = Zf+uf+1 is correct exactly when
 * cand <= sqrt(r) + c. The difference cand - c = Zf + 1 - (c - uf) is always
 * positive, so that is equivalent to (cand - c)^2 <= r.
 */
static rabin_err_t rq_sqrt_shift_floor(rz_t* m, const rq_t* r, const rq_t* c)
{
  rq_t t, sq;
  rq_init_multi(&t, &sq, NULL);
  rz_t one, cand, uf, prod;
  rz_init_multi(&one, &cand, &uf, &prod, NULL);

  rabin_err_t err = RABIN_SUCCESS;
  TRY(rz_set_u64(&one, 1));
  TRY(rq_sqrt_floor(m, r));

  // uf = floor(c): rz_div truncates toward zero, so a negative non-exact
  // quotient needs one correction.
  TRY(rz_copy(&uf, &c->num));
  TRY(rz_div(&uf, &uf, &c->den));
  if (c->num.is_neg) {
    TRY(rz_mul(&prod, &uf, &c->den));
    if (rz_cmp(&prod, &c->num) != 0) {
      TRY(rz_sub(&uf, &uf, &one));
    }
  }

  TRY(rz_add(m, m, &uf));       // m = Zf + uf
  TRY(rz_add(&cand, m, &one));  // cand = Zf + uf + 1

  // accept cand when cand - c <= sqrt(r), i.e. (cand - c)^2 <= r
  TRY(rq_set(&t, &cand));
  TRY(rq_sub(&t, &t, c));
  TRY(rq_sqr(&sq, &t));
  if (rq_cmp(&sq, r) <= 0) TRY(rz_copy(m, &cand));

cleanup:
  rq_clear_multi(&t, &sq, NULL);
  rz_clear_multi(&one, &cand, &uf, &prod, NULL);
  return err;
}

typedef struct {
  u64 n;
  rq_t* q;        // Cholesky coefficients, row-major n x n
  rz_t* xv;       // current integer coordinates
  const rq_t* C;  // enumeration threshold
  rq_t* best;     // best Q' value found
  rz_t* ybest;    // coordinates of the best vector
  bool have;
} fp_enum_ctx_t;

/*
 * Algorithm 2.7.5 (short vectors): at level i enumerate every x_i with
 * q_ii (x_i + shift_i)^2 <= T, where shift_i = sum_{j>i} q_ij x_j. The
 * coordinates are visited from n-1 down to 0 so that shift_i only depends
 * on coordinates already fixed by the enclosing levels. T is the budget
 * left by the enclosing levels, and T - q_ii (x_i + shift_i)^2 is passed
 * on; a negative budget can never be repaid (all terms are >= 0), so such
 * a node is pruned.
 */
static rabin_err_t fp_enum(fp_enum_ctx_t* f, u64 i, const rq_t* T)
{
  if (T->num.is_neg && !rz_is_zero(&T->num)) return RABIN_SUCCESS;

  const u64 n = f->n;
  const rq_t* qii = &f->q[i * n + i];

  rq_t budget, shift, term, val, rem, tmp;
  rq_init_multi(&budget, &shift, &term, &val, &rem, &tmp, NULL);
  rz_t lo, hi, one;
  rz_init_multi(&lo, &hi, &one, NULL);

  rabin_err_t err = RABIN_SUCCESS;
  TRY(rz_set_u64(&one, 1));
  TRY(rq_set_u64(&shift, 0));

  TRY(rq_div(&budget, T, qii));

  for (u64 j = i + 1; j < n; j++) {
    if (rz_is_zero(&f->q[i * n + j].num)) continue;
    TRY(rq_set(&tmp, &f->xv[j]));
    TRY(rq_mul(&tmp, &tmp, &f->q[i * n + j]));
    TRY(rq_add(&shift, &shift, &tmp));
  }

  // bounds: ceil(-sqrt(budget) - shift) = -floor(sqrt(budget) + shift)
  //         <= x_i <= floor(sqrt(budget) - shift)
  TRY(rq_sqrt_shift_floor(&hi, &budget, &shift));
  TRY(rz_neg(&lo, &hi));
  TRY(rq_copy(&tmp, &shift));
  TRY(rq_neg(&tmp));
  TRY(rq_sqrt_shift_floor(&hi, &budget, &tmp));

  for (; rz_cmp(&lo, &hi) <= 0; rz_add(&lo, &lo, &one)) {
    TRY(rz_copy(&f->xv[i], &lo));

    TRY(rq_set(&term, &lo));
    TRY(rq_add(&term, &term, &shift));
    TRY(rq_sqr(&val, &term));
    TRY(rq_mul(&val, &val, qii));

    if (i == 0) {
      // canonical representative of {y, -y}: skip when the first
      // non-zero coordinate is negative, skip y = 0
      u64 k = 0;
      while (k < n && rz_is_zero(&f->xv[k])) k++;
      if (k == n) continue;
      if (f->xv[k].is_neg) continue;

      // Q' = C - rem; val <= T by construction of the bounds, so Q' <= C
      TRY(rq_sub(&rem, T, &val));
      TRY(rq_sub(&term, f->C, &rem));
      if (!f->have || rq_cmp(&term, f->best) < 0) {
        TRY(rq_copy(f->best, &term));
        for (u64 k = 0; k < n; k++) {
          TRY(rz_copy(&f->ybest[k], &f->xv[k]));
        }
        f->have = true;
      }
    } else {
      TRY(rq_sub(&rem, T, &val));
      TRY(fp_enum(f, i - 1, &rem));
    }
  }

cleanup:
  rq_clear_multi(&budget, &shift, &term, &val, &rem, &tmp, NULL);
  rz_clear_multi(&lo, &hi, &one, NULL);
  return err;
}

/*
 * num = sum_k U[k][ci] * AU[k][cj], where AU = A U shares the denominator
 * AU->den; the numerator of entry (ci, cj) of U^t A U follows by dividing
 * by AU->den. ci and cj index the columns of U.
 */
static rabin_err_t fp_gram_entry(rz_t* num, const rmat_t* U, const rmatq_t* AU,
                                 u64 n, u64 ci, u64 cj)
{
  rz_t t;
  rz_init(&t);
  rabin_err_t err = RABIN_SUCCESS;
  TRY(rz_set_u64(num, 0));
  for (u64 k = 0; k < n; k++) {
    TRY(rz_mul(&t, RMAT_GET(U, k, ci), RMAT_GET(&AU->M, k, cj)));
    TRY(rz_add(num, num, &t));
  }
cleanup:
  rz_clear(&t);
  return err;
}

/* norms[i] = (U^t A U)_{i,i} = Q(u_i) = ||s_i||^2, AU = A U. */
static rabin_err_t fp_column_norms(rq_t* norms, const rmat_t* U,
                                   const rmatq_t* AU, u64 n)
{
  rabin_err_t err = RABIN_SUCCESS;
  for (u64 i = 0; i < n; i++) {
    TRY(fp_gram_entry(&norms[i].num, U, AU, n, i, i));
    TRY(rz_copy(&norms[i].den, &AU->den));
    TRY(rq_normalize(&norms[i]));
  }
  return err;
cleanup:
  return err;
}

/* Stable selection sort of ord[] on increasing norms[ord[]]. */
static void fp_order_columns(u64* ord, const rq_t* norms, u64 n)
{
  for (u64 i = 0; i < n; i++) ord[i] = i;
  for (u64 i = 0; i + 1 < n; i++) {
    u64 best = i;
    for (u64 j = i + 1; j < n; j++) {
      if (rq_cmp(&norms[ord[j]], &norms[ord[best]]) < 0) best = j;
    }
    if (best != i) {
      u64 tmp = ord[i];
      ord[i] = ord[best];
      ord[best] = tmp;
    }
  }
}

/* Gp = (U P)^t A (U P) = P^t (U^t A U) P from AU = A U. */
static rabin_err_t fp_permuted_gram(rmatq_t* Gp, const rmat_t* U,
                                    const rmatq_t* AU, const u64* ord, u64 n)
{
  rz_t num;
  rz_init(&num);
  rabin_err_t err = RABIN_SUCCESS;
  for (u64 i = 0; i < n; i++) {
    for (u64 j = 0; j < n; j++) {
      TRY(fp_gram_entry(&num, U, AU, n, ord[i], ord[j]));
      TRY(rz_copy(RMAT_GET(&Gp->M, i, j), &num));
    }
  }
  TRY(rz_copy(&Gp->den, &AU->den));
  TRY(rmatq_normalize(Gp));
cleanup:
  rz_clear(&num);
  return err;
}

/*
 * Cholesky coefficients of Gp (steps 1-3 of Algorithm 2.7.6, no square
 * roots): they express Q'(y) = y^t Gp y as sum_i q_ii (y_i +
 * sum_{j>i} q_ij y_j)^2, where the strict upper triangle holds the
 * multipliers G[i][j] / pivot and the strict lower triangle the
 * eliminated values. q must hold n*n initialized rationals. A pivot that
 * is not strictly positive means Gp is not positive definite.
 */
static rabin_err_t fp_cholesky_coeffs(rq_t* q, const rmatq_t* Gp, u64 n)
{
  rq_t t;
  rq_init(&t);
  rabin_err_t err = RABIN_SUCCESS;
  for (u64 i = 0; i < n; i++) {
    for (u64 j = 0; j < n; j++) {
      TRY(rmatq_get(&q[i * n + j], Gp, i, j));
    }
  }
  for (u64 i = 0; i < n; i++) {
    rq_t* pivot = &q[i * n + i];
    if (rz_is_zero(&pivot->num)) {
      err = RABIN_ERR_DIV_BY_ZERO;
      goto cleanup;
    }
    if (pivot->num.is_neg) {
      err = RABIN_ERR_INVALID_ARG;
      goto cleanup;
    }
    for (u64 j = i + 1; j < n; j++) {
      TRY(rq_copy(&q[j * n + i], &q[i * n + j]));
      TRY(rq_div(&q[i * n + j], &q[i * n + j], pivot));
    }
    for (u64 k = i + 1; k < n; k++) {
      if (rz_is_zero(&q[k * n + i].num)) continue;
      for (u64 l = k; l < n; l++) {
        if (rz_is_zero(&q[i * n + l].num)) continue;
        TRY(rq_mul(&t, &q[k * n + i], &q[i * n + l]));
        TRY(rq_sub(&q[k * n + l], &q[k * n + l], &t));
      }
    }
  }
cleanup:
  rq_clear(&t);
  return err;
}

rabin_err_t rmatq_fincke_pohst(rvec_t* x, rq_t* norm, const rmatq_t* A)
{
  if (x == NULL || norm == NULL || A == NULL) return RABIN_ERR_NULL_PTR;
  if (A->M.rows != A->M.cols) return RABIN_ERR_MATRIX_DIM;

  const u64 n = A->M.rows;
  if (n == 0) return RABIN_ERR_INVALID_ARG;
  if (x->size != n) return RABIN_ERR_MATRIX_DIM;

  rmatq_t Ainv = {0}, Uq = {0}, AU = {0}, Gp = {0};
  rmat_t B = {0}, H = {0}, Ht = {0}, U = {0}, Ur = {0};
  rvec_t yv = {0};
  rq_t* q = NULL;
  rq_t* norms = NULL;
  rz_t* xv = NULL;
  rz_t* ybest = NULL;
  u64* ord = NULL;
  rq_t best, C;
  u64 norms_inited = 0, q_inited = 0, xv_inited = 0, yb_inited = 0;

  rq_init_multi(&best, &C, NULL);

  rabin_err_t err = RABIN_SUCCESS;

  // 1. [LLL reduction of the dual lattice] rows of A^-1 generate the
  // dual of (Z^n, Q), the same lattice as the rows of Cholesky's R^-1,
  // without needing any square root.
  TRY(rmatq_init(&Ainv, n, n));
  TRY(rmatq_inv(&Ainv, A));
  TRY(rmat_init(&B, n, n));
  TRY(rmat_transpose(&B, &Ainv.M));
  TRY(rmat_lll(&B, n, 0.75, &H));

  // 2. U = (H^t)^-1 is unimodular, so its columns are a Z-basis of Z^n;
  // S^-1 = U^-1 R^-1 and S = R U.
  TRY(rmat_init(&Ht, n, n));
  TRY(rmat_transpose(&Ht, &H));
  TRY(rmat_init(&U, n, n));
  TRY(rmat_inv_unimodular(&U, &Ht));

  // 3. Order the columns u_i of U on increasing norm Q(u_i) = ||s_i||^2
  // and build the permuted Gram matrix Gp = (U P)^t A (U P). AU = A U
  // already holds a common denominator, so AU->den is shared by all
  // numerators of U^t A U.
  TRY(rmatq_init(&Uq, n, n));
  TRY(rmatq_copyZ(&Uq, &U));
  TRY(rmatq_init(&AU, n, n));
  TRY(rmatq_mul(&AU, A, &Uq));

  norms = malloc(n * sizeof(rq_t));
  ord = malloc(n * sizeof(u64));
  if (norms == NULL || ord == NULL) {
    err = RABIN_ERR_OUT_OF_MEMORY;
    goto cleanup;
  }
  for (u64 i = 0; i < n; i++) {
    TRY(rq_init(&norms[i]));
    norms_inited = i + 1;
  }
  TRY(fp_column_norms(norms, &U, &AU, n));
  fp_order_columns(ord, norms, n);
  TRY(rq_copy(&C, &norms[ord[0]]));

  TRY(rmatq_init(&Gp, n, n));
  TRY(fp_permuted_gram(&Gp, &U, &AU, ord, n));

  // 4. Cholesky coefficients of Gp (steps 1-3 of Algorithm 2.7.6, no
  // square roots).
  q = malloc(n * n * sizeof(rq_t));
  if (q == NULL) {
    err = RABIN_ERR_OUT_OF_MEMORY;
    goto cleanup;
  }
  for (u64 k = 0; k < n * n; k++) {
    TRY(rq_init(&q[k]));
    q_inited = k + 1;
  }
  TRY(fp_cholesky_coeffs(q, &Gp, n));

  // 5. [Enumeration] bound C = min ||s_i||^2; the minimum is attained
  // since the s_i are genuine lattice vectors.
  xv = malloc(n * sizeof(rz_t));
  ybest = malloc(n * sizeof(rz_t));
  if (xv == NULL || ybest == NULL) {
    err = RABIN_ERR_OUT_OF_MEMORY;
    goto cleanup;
  }
  for (u64 k = 0; k < n; k++) {
    TRY(rz_init(&xv[k]));
    xv_inited = k + 1;
  }
  for (u64 k = 0; k < n; k++) {
    TRY(rz_init(&ybest[k]));
    yb_inited = k + 1;
  }
  TRY(rq_copy(&best, &C));  // any solution <= C; C itself bounds the best

  fp_enum_ctx_t f = {n, q, xv, &C, &best, ybest, false};
  TRY(fp_enum(&f, n - 1, &C));
  if (!f.have) {
    err = RABIN_ERR_INVALID_ARG;
    goto cleanup;
  }

  // 6. Back to the original coordinates: x = U P y.
  TRY(rmat_init(&Ur, n, n));
  for (u64 j = 0; j < n; j++) {
    for (u64 i = 0; i < n; i++) {
      TRY(rz_copy(RMAT_GET(&Ur, i, j), RMAT_GET(&U, i, ord[j])));
    }
  }
  TRY(rvec_init(&yv, n));
  for (u64 k = 0; k < n; k++) {
    TRY(rz_copy(&yv.data[k], &ybest[k]));
  }
  TRY(rmat_mv(x, &Ur, &yv));
  TRY(rq_copy(norm, &best));

cleanup:
  rmatq_clear(&Ainv);
  rmatq_clear(&Uq);
  rmatq_clear(&AU);
  rmatq_clear(&Gp);
  rmat_clear(&B);
  rmat_clear(&H);
  rmat_clear(&Ht);
  rmat_clear(&U);
  rmat_clear(&Ur);
  rvec_clear(&yv);
  if (q != NULL) {
    for (u64 k = 0; k < q_inited; k++) rq_clear(&q[k]);
    free(q);
  }
  if (norms != NULL) {
    for (u64 k = 0; k < norms_inited; k++) rq_clear(&norms[k]);
    free(norms);
  }
  if (xv != NULL) {
    for (u64 k = 0; k < xv_inited; k++) rz_clear(&xv[k]);
    free(xv);
  }
  if (ybest != NULL) {
    for (u64 k = 0; k < yb_inited; k++) rz_clear(&ybest[k]);
    free(ybest);
  }
  free(ord);
  rq_clear_multi(&best, &C, NULL);
  return err;
}
