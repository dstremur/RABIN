/*
 * rmat_lll.c
 *
 * LLL basis reduction (Cohen 2.6.7) and modular LLL (Cohen 2.6.8) for rmat_t
 * lattices.
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

#include "../../include/rmat.h"
#include "../../include/rzmath.h"
#include "rmat_internal.h"
/* ------------------------------------------------------------------ */
/* LLL basis reduction                                                 */
/* ------------------------------------------------------------------ */
/*
 * rmat_lll()  implements Cohen's Algorithm 2.6.7 (all-integer LLL,
 *             full rank only).
 * rmat_mlll() handles zero / dependent columns as a two-phase rebuild:
 *             (1) a rank-revealing HNF column reduction (Cohen 2.4.5)
 *             tracking the unimodular transform U with
 *                 B_orig * U = [ 0_{m x (n-p)} | W ],
 *             so the first n-p columns of U are a Z-basis of the
 *             relation space and W (m x p) a Z-basis of the column
 *             lattice;
 *             (2) rmat_lll() on W, W <- W * Q, and
 *                 B <- [ 0 | W ],  H <- [ U_rel | U_lift * Q ].
 *
 * rmat_lll() maintains the current generator vectors V (initially a
 * copy of B) together with the unimodular transform H such that
 * V = B_orig * H holds at all times; on success B is overwritten with
 * V. It stores the scaled Gram-Schmidt data (1-based indices as in
 * Cohen; 0-based here)
 *
 *     lambda[i][j] = d[j+1] * mu[i+1][j+1]     (0 <= j < i < n)
 *     d[0] = 1,  d[i+1] = ||b_{i+1}||^2 * d[i]
 *
 * so that all quantities are integers. In rmat_lll(), RED(k, l)
 * size-reduces vector k against vector l and SWAP exchanges vectors
 * k-1 and k (0-based).
 */

rabin_err_t lll_scr_init(lll_scr_t* s)
{
  if (rz_init_multi(&s->u, &s->v, &s->w, &s->one, &s->two, NULL) !=
      RABIN_SUCCESS)
    return RABIN_ERR_OUT_OF_MEMORY;

  rabin_err_t err = rz_set_u64(&s->one, 1);
  if (err == RABIN_SUCCESS) err = rz_set_u64(&s->two, 2);
  return err;
}

void lll_scr_clear(lll_scr_t* s)
{
  rz_clear_multi(&s->u, &s->v, &s->w, &s->one, &s->two, NULL);
}

/*
 * delta as an exact binary rational. Any double delta in (0, 1] is one:
 * frexp() gives delta = m * 2^e with m in [1/2, 1), so delta = M / 2^E
 * with M = m * 2^53 an integer in [2^52, 2^53) and E = 53 - e >= 1.
 * Stripping the common factors of two leaves num / 2^den_bits with
 * num < 2^53 (and num odd unless the fraction is integral).
 */
typedef struct {
  uint64_t num;
  u64 den_bits;
} lll_delta_t;

static lll_delta_t lll_delta_split(double delta)
{
  int e = 0;
  const double m = frexp(delta, &e);
  u64 M = (u64)(m * 9007199254740992.0); /* 2^53, exact */
  u64 den_bits = (u64)(53 - e);
  u64 s = (u64)__builtin_ctzll(M);
  if (s > den_bits) s = den_bits;

  lll_delta_t d;
  d.num = M >> s;
  d.den_bits = den_bits - s;
  return d;
}

/*
 * Lovasz test for the current k (1-based, lambda = lambda[k][k-1]):
 *
 *   swap iff d[k+1]*d[k-1] < delta*d[k]^2 - lambda^2.
 *
 * With delta = num / 2^den_bits this is evaluated as the exact integer
 * inequality
 *
 *   2^den_bits * d[k+1]*d[k-1] < num * d[k]^2 - 2^den_bits * lambda^2,
 *
 * so the reduction result is independent of any floating-point rounding.
 */
static rabin_err_t lll_lovasz(rmat_t* lam, const rz_t* dv, u64 k,
                              const lll_delta_t* dl, const rz_t* dden,
                              bool* swap, lll_scr_t* s)
{
  // lhs = d[k+1]*d[k-1]  (dv[k], dv[k-2])
  rabin_err_t err = rz_mul(&s->u, &dv[k], &dv[k - 2]);
  if (err != RABIN_SUCCESS) return err;
  // M = d[k]^2 (dv[k-1]), N = lambda^2
  if ((err = rz_sqr(&s->v, &dv[k - 1])) != RABIN_SUCCESS) return err;
  if ((err = rz_sqr(&s->w, RMAT_GET(lam, k - 1, k - 2))) != RABIN_SUCCESS)
    return err;

  // M <- num * M  (num < 2^53 always fits a positive i64)
  if ((err = rz_mul_i64(&s->v, &s->v, (i64)dl->num)) != RABIN_SUCCESS)
    return err;

  // lhs and N pick up the factor 2^den_bits
  if (dl->den_bits != 0) {
    if (dl->den_bits < 63) {
      const i64 den = (i64)1 << dl->den_bits;
      if ((err = rz_mul_i64(&s->w, &s->w, den)) != RABIN_SUCCESS) return err;
      if ((err = rz_mul_i64(&s->u, &s->u, den)) != RABIN_SUCCESS) return err;
    } else {
      if ((err = rz_mul(&s->w, &s->w, dden)) != RABIN_SUCCESS) return err;
      if ((err = rz_mul(&s->u, &s->u, dden)) != RABIN_SUCCESS) return err;
    }
  }

  if ((err = rz_sub(&s->v, &s->v, &s->w)) != RABIN_SUCCESS) return err;
  *swap = (rz_cmp(&s->u, &s->v) < 0);
  return RABIN_SUCCESS;
}

/* r = (column a of V) . (column b of V) */
rabin_err_t lll_dot(rz_t* r, const rmat_t* V, u64 a, u64 b, rz_t* t)
{
  rabin_err_t err = rz_set_u64(r, 0);
  if (err != RABIN_SUCCESS) return err;

  for (u64 i = 0; i < V->rows; i++) {
    if ((err = rz_mul(t, RMAT_GET(V, i, a), RMAT_GET(V, i, b))) !=
        RABIN_SUCCESS)
      return err;
    if ((err = rz_add(r, r, t)) != RABIN_SUCCESS) return err;
  }
  return RABIN_SUCCESS;
}

/* ---------------------- LLL (Cohen 2.6.7) ------------------------- */

/* RED1: size-reduce vector k against vector l (0-based, l < k). */
rabin_err_t lll_red1(rmat_t* V, rmat_t* H, rmat_t* lam, const rz_t* dv, u64 k,
                     u64 l, lll_scr_t* s)
{
  rz_t* lk = RMAT_GET(lam, k, l);

  // stop if |2 * lambda[k][l]| <= d[l+1]
  rabin_err_t err = rz_mul(&s->u, &s->two, lk);
  if (err != RABIN_SUCCESS) return err;
  if (rz_cmp_abs(&s->u, &dv[l + 1]) <= 0) return RABIN_SUCCESS;

  // q = nearest(lambda[k][l] / d[l+1]) =
  //     floor((2*lambda[k][l] + d[l+1]) / (2*d[l+1]))
  err = rz_add(&s->v, &s->u, &dv[l + 1]);
  if (err == RABIN_SUCCESS) err = rz_mul(&s->u, &s->two, &dv[l + 1]);
  if (err == RABIN_SUCCESS) err = rz_div_euclid(&s->w, &s->v, &s->u);
  if (err != RABIN_SUCCESS) return err;
  rz_t* q = &s->w;

  // V_k -= q * V_l
  for (u64 i = 0; i < V->rows; i++) {
    if ((err = rz_mul(&s->u, q, RMAT_GET(V, i, l))) != RABIN_SUCCESS)
      return err;
    if ((err = rz_sub(RMAT_GET(V, i, k), RMAT_GET(V, i, k), &s->u)) !=
        RABIN_SUCCESS)
      return err;
  }
  // H[:, k] -= q * H[:, l]
  for (u64 i = 0; i < H->rows; i++) {
    if ((err = rz_mul(&s->u, q, RMAT_GET(H, i, l))) != RABIN_SUCCESS)
      return err;
    if ((err = rz_sub(RMAT_GET(H, i, k), RMAT_GET(H, i, k), &s->u)) !=
        RABIN_SUCCESS)
      return err;
  }

  // lambda[k][l] -= q * d[l+1];  lambda[k][i] -= q * lambda[l][i]  (i < l)
  if ((err = rz_mul(&s->u, q, &dv[l + 1])) != RABIN_SUCCESS) return err;
  if ((err = rz_sub(lk, lk, &s->u)) != RABIN_SUCCESS) return err;
  for (u64 i = 0; i < l; i++) {
    if ((err = rz_mul(&s->u, q, RMAT_GET(lam, l, i))) != RABIN_SUCCESS)
      return err;
    if ((err = rz_sub(RMAT_GET(lam, k, i), RMAT_GET(lam, k, i), &s->u)) !=
        RABIN_SUCCESS)
      return err;
  }
  return RABIN_SUCCESS;
}

/* SWAP1: swap vectors k-1 and k (0-based, k >= 1). */
static rabin_err_t lll_swap1(rmat_t* V, rmat_t* H, rmat_t* lam, rz_t* dv, u64 k,
                             u64 kmax, lll_scr_t* s)
{
  rabin_err_t err;

  // swap V columns k and k-1
  for (u64 i = 0; i < V->rows; i++) {
    if ((err = rz_copy(&s->u, RMAT_GET(V, i, k))) != RABIN_SUCCESS) return err;
    if ((err = rz_copy(RMAT_GET(V, i, k), RMAT_GET(V, i, k - 1))) !=
        RABIN_SUCCESS)
      return err;
    if ((err = rz_copy(RMAT_GET(V, i, k - 1), &s->u)) != RABIN_SUCCESS)
      return err;
  }
  // swap H columns k and k-1
  for (u64 i = 0; i < H->rows; i++) {
    if ((err = rz_copy(&s->u, RMAT_GET(H, i, k))) != RABIN_SUCCESS) return err;
    if ((err = rz_copy(RMAT_GET(H, i, k), RMAT_GET(H, i, k - 1))) !=
        RABIN_SUCCESS)
      return err;
    if ((err = rz_copy(RMAT_GET(H, i, k - 1), &s->u)) != RABIN_SUCCESS)
      return err;
  }
  // swap lambda[k][j] <-> lambda[k-1][j]  (j < k-1)
  for (u64 j = 0; j + 1 < k; j++) {
    if ((err = rz_copy(&s->u, RMAT_GET(lam, k, j))) != RABIN_SUCCESS)
      return err;
    if ((err = rz_copy(RMAT_GET(lam, k, j), RMAT_GET(lam, k - 1, j))) !=
        RABIN_SUCCESS)
      return err;
    if ((err = rz_copy(RMAT_GET(lam, k - 1, j), &s->u)) != RABIN_SUCCESS)
      return err;
  }

  // A <- lambda[k][k-1]  (= lambda_{s+1,s} with s = k, 1-based)
  rz_t* A = RMAT_GET(lam, k, k - 1);

  // for i > k (up to kmax), 1-based I = i+1 > s+1:
  //   t <- lambda[i][k]   (old lambda_{I,s+1})
  //   lambda[i][k]   <- (d[k+1]*lambda[i][k-1] - A*t) / d[k] (lambda'_{I,s+1})
  //   lambda[i][k-1] <- (A*lambda[i][k-1] + d[k-1]*t) / d[k]   (lambda'_{I,s})
  for (u64 i = k + 1; i <= kmax; i++) {
    if ((err = rz_copy(&s->w, RMAT_GET(lam, i, k))) != RABIN_SUCCESS)
      return err;
    rz_t* t = &s->w;
    if ((err = rz_mul(&s->u, &dv[k + 1], RMAT_GET(lam, i, k - 1))) !=
        RABIN_SUCCESS)
      return err;
    if ((err = rz_mul(&s->v, A, t)) != RABIN_SUCCESS) return err;
    if ((err = rz_sub(&s->u, &s->u, &s->v)) != RABIN_SUCCESS) return err;
    if ((err = rz_div_exact(RMAT_GET(lam, i, k), &s->u, &dv[k])) !=
        RABIN_SUCCESS)
      return err;
    if ((err = rz_mul(&s->v, A, RMAT_GET(lam, i, k - 1))) != RABIN_SUCCESS)
      return err;
    if ((err = rz_mul(&s->u, &dv[k - 1], t)) != RABIN_SUCCESS) return err;
    if ((err = rz_add(&s->v, &s->v, &s->u)) != RABIN_SUCCESS) return err;
    if ((err = rz_div_exact(RMAT_GET(lam, i, k - 1), &s->v, &dv[k])) !=
        RABIN_SUCCESS)
      return err;
  }

  // d[k] <- (d[k+1]*d[k-1] + A^2) / d[k]   (new d_{s+1}, 1-based)
  err = rz_mul(&s->u, &dv[k + 1], &dv[k - 1]);
  if (err == RABIN_SUCCESS) err = rz_sqr(&s->v, A);
  if (err == RABIN_SUCCESS) err = rz_add(&s->u, &s->u, &s->v);
  if (err == RABIN_SUCCESS) err = rz_div_exact(&s->v, &s->u, &dv[k]);
  if (err != RABIN_SUCCESS) return err;
  return rz_copy(&dv[k], &s->v);
}

/*
 * Cohen 2.6.7 main loop, with 1-based k and kmax exactly as in Cohen
 * (0-based array indices are k-1 etc.). Requires full rank; returns
 * RABIN_ERR_INVALID_ARG when a d[k] becomes zero.
 */
static rabin_err_t rmat_lll_run(rmat_t* V, rmat_t* H, rmat_t* lam, rz_t* dv,
                                u64 n, double delta)
{
  lll_scr_t s;
  if (lll_scr_init(&s) != RABIN_SUCCESS) return RABIN_ERR_OUT_OF_MEMORY;

  rz_t dden;
  if (rz_init(&dden) != RABIN_SUCCESS) {
    lll_scr_clear(&s);
    return RABIN_ERR_OUT_OF_MEMORY;
  }

  // delta = num / 2^den_bits (exact); dden = 2^den_bits only needed when
  // the denominator no longer fits an i64 multiplier
  lll_delta_t dl = lll_delta_split(delta);
  rabin_err_t err = RABIN_SUCCESS;
  if (dl.den_bits >= 63) {
    if (rz_set_u64(&dden, 1) == RABIN_SUCCESS)
      err = rz_lshift(&dden, &dden, (int)dl.den_bits);
    if (err != RABIN_SUCCESS) goto cleanup;
  }

  TRY(rz_set_u64(&dv[0], 1));

  u64 k = 1;
  u64 kmax = 0;

step2:
  // [Initialize GS data for vector k]
  if (k > kmax) {
    kmax = k;
    for (u64 j = 1; j <= k; j++) {
      TRY(lll_dot(&s.u, V, k - 1, j - 1, &s.v));
      // U <- (d[i]*U - lambda[k][i]*lambda[j][i]) / d[i-1]   (1-based i < j,
      // d[1] = 1, d[i] stored in dv[i-1] (0-based))
      for (u64 i = 1; i < j; i++) {
        TRY(rz_mul(&s.w, &dv[i], &s.u));
        TRY(rz_mul(&s.v, RMAT_GET(lam, k - 1, i - 1),
                   RMAT_GET(lam, j - 1, i - 1)));
        TRY(rz_sub(&s.w, &s.w, &s.v));
        TRY(rz_div_exact(&s.u, &s.w, &dv[i - 1]));
      }
      if (j < k) {
        TRY(rz_copy(RMAT_GET(lam, k - 1, j - 1), &s.u));
      } else {
        TRY(rz_copy(&dv[k], &s.u));
        if (rz_is_zero(&dv[k])) {
          err = RABIN_ERR_INVALID_ARG;
          goto cleanup;
        }
      }
    }
  }

  // [LLL test: size-reduce k against k-1, check the Lovasz condition]
  if (k > 1) {
    TRY(lll_red1(V, H, lam, dv, k - 1, k - 2, &s));

    // Lovasz condition (1-based k, lambda = lambda[k][k-1], d[i] in dv[i-1]),
    // evaluated exactly for every double delta (see lll_lovasz):
    //   swap iff d[k+1]*d[k-1] < delta*d[k]^2 - lambda^2
    // (i.e. ||b*_k||^2 < (delta - mu^2) * ||b*_{k-1}||^2 with mu =
    // lambda/d[k]).
    bool swap = false;
    TRY(lll_lovasz(lam, dv, k, &dl, &dden, &swap, &s));

    if (swap) {
      TRY(lll_swap1(V, H, lam, dv, k - 1, kmax - 1, &s));
      k--;
      // the swapped-in vector now sits at k and must be re-examined
      // (size-reduce + Lovasz) before advancing: restart at step2.
      goto step2;
    }
  }

  // descending: RED1(k,l) updates lambda[k][i] for i < l, so l must decrease
  for (u64 l = k - 1; l >= 1; l--) {
    TRY(lll_red1(V, H, lam, dv, k - 1, l - 1, &s));
  }
  k++;
  if (k <= n) goto step2;

  err = RABIN_SUCCESS;
cleanup:
  rz_clear(&dden);
  lll_scr_clear(&s);
  return err;
}

rabin_err_t rmat_lll(rmat_t* B, u64 n, double delta, rmat_t* H)
{
  if (B == NULL || H == NULL) return RABIN_ERR_NULL_PTR;
  if (!(delta > 0.0 && delta <= 1.0)) return RABIN_ERR_INVALID_ARG;
  if (n != B->cols) return RABIN_ERR_INVALID_ARG;

  // nothing to reduce; H is re-initialized to the n x n identity
  rabin_err_t err = rmat_init(H, n, n);
  if (err != RABIN_SUCCESS) return err;
  err = rmat_id(H, n);
  if (err != RABIN_SUCCESS) {
    rmat_clear(H);
    *H = (rmat_t){0};
    return err;
  }
  if (n <= 1) return RABIN_SUCCESS;

  rmat_t V = {0}, lam = {0};
  rz_t* dv = NULL;

  if (rmat_init(&V, B->rows, n) != RABIN_SUCCESS) {
    rmat_clear(H);
    *H = (rmat_t){0};
    return RABIN_ERR_OUT_OF_MEMORY;
  }
  if (rmat_init(&lam, n, n) != RABIN_SUCCESS) {
    rmat_clear(&V);
    rmat_clear(H);
    *H = (rmat_t){0};
    return RABIN_ERR_OUT_OF_MEMORY;
  }
  dv = malloc((n + 1) * sizeof(rz_t));
  if (dv == NULL) {
    rmat_clear(&V);
    rmat_clear(&lam);
    rmat_clear(H);
    *H = (rmat_t){0};
    return RABIN_ERR_OUT_OF_MEMORY;
  }
  for (u64 i = 0; i <= n; i++) {
    if (rz_init(&dv[i]) != RABIN_SUCCESS) {
      for (u64 j = 0; j < i; j++) rz_clear(&dv[j]);
      free(dv);
      rmat_clear(&V);
      rmat_clear(&lam);
      rmat_clear(H);
      *H = (rmat_t){0};
      return RABIN_ERR_OUT_OF_MEMORY;
    }
  }

  // H is already the identity (set above); on failure B is left untouched
  // and H is reset to an empty matrix so callers own no allocations
  err = rmat_copy(&V, B);
  if (err == RABIN_SUCCESS) err = rmat_lll_run(&V, H, &lam, dv, n, delta);
  if (err == RABIN_SUCCESS) err = rmat_copy(B, &V);

  for (u64 i = 0; i <= n; i++) rz_clear(&dv[i]);
  free(dv);
  rmat_clear(&V);
  rmat_clear(&lam);
  if (err != RABIN_SUCCESS) {
    rmat_clear(H);
    *H = (rmat_t){0};
  }
  return err;
}

/* ---------------------- MLLL (Cohen 2.6.8) ------------------------ */

/*
 * HNF via the exact column reduction (Cohen's Algorithm 2.4.5, same
 * loop as rmat_hermite_gcd()), additionally tracking the unimodular
 * column transformation.
 *
 * The reduction drives (n - p) columns to zero and drops them from the
 * left (p = rank of A), so on success:
 *
 *   W = A_work[:, k_final..n-1]   (m x p, HNF of the column lattice)
 *   A * U = [ 0_{m x (n-p)} | W ]
 *
 * Hence the first n - p columns of U form a Z-basis of the relation
 * space { x in Z^n : A x = 0 }, and the last p columns of U are
 * integer lifts of the HNF basis. W and U must point to empty
 * (zero-initialized) rmat_t objects; on success they hold the result.
 */
rabin_err_t hnf_with_transform(rmat_t* W, rmat_t* U, const rmat_t* A)
{
  if (W == NULL || U == NULL || A == NULL) return RABIN_ERR_NULL_PTR;

  u64 m = A->rows;
  u64 n = A->cols;
  if (m == 0 || n == 0) return RABIN_ERR_INVALID_ARG;

  u64 i = m - 1;
  u64 j = n - 1;
  u64 k = n - 1;
  u64 l = 0;
  if (m > n) {
    l = m - n;
  }

  bool u_live = false;
  rmat_t A_work = {0};
  rz_t u, v, d, temp, temp2, b, q_k, q_j;
  rz_t* B = NULL;
  rz_t* C = NULL;
  rabin_err_t err = RABIN_SUCCESS;

  TRY(rmat_init(&A_work, m, n));
  TRY(rmat_copy(&A_work, A));

  TRY(rmat_init(U, n, n));
  u_live = true;
  TRY(rmat_id(U, n));

  rz_init_multi(&u, &v, &d, &temp, &temp2, &b, &q_k, &q_j, NULL);

  B = malloc(m * sizeof(rz_t));
  if (B == NULL) {
    err = RABIN_ERR_OUT_OF_MEMORY;
    goto cleanup;
  }
  for (u64 x = 0; x < m; x++) {
    rz_init(&B[x]);
  }

  C = malloc(n * sizeof(rz_t));
  if (C == NULL) {
    err = RABIN_ERR_OUT_OF_MEMORY;
    goto cleanup;
  }
  for (u64 x = 0; x < n; x++) {
    rz_init(&C[x]);
  }

// 2. [Check zero]
step2:
  if (j == 0) {
    goto step4;
  }
  j--;
  if (rz_is_zero(RMAT_GET(&A_work, i, j))) {
    goto step2;
  }

  // 3. [Euclidean step] (column ops applied to A_work and U)
  TRY(hnf_euclid_step(&A_work, U, B, C, i, j, k, &u, &v, &d, &q_k, &q_j, &temp,
                      &temp2));
  goto step2;

  // 4. [Final reduction]
step4:
  TRY(rz_copy(&b, RMAT_GET(&A_work, i, k)));
  if (b.is_neg) {  // zero is never negative, no extra rz_is_zero test needed
    for (u64 x = 0; x < m; x++) {
      TRY(rz_neg(RMAT_GET(&A_work, x, k), RMAT_GET(&A_work, x, k)));
    }
    for (u64 x = 0; x < n; x++) {
      TRY(rz_neg(RMAT_GET(U, x, k), RMAT_GET(U, x, k)));
    }
    TRY(rz_neg(&b, &b));
  }

  if (rz_is_zero(&b)) {
    k++;
    goto step5;
  }

  for (u64 j_0 = k + 1; j_0 < n; j_0++) {
    TRY(rz_div_euclid(&temp2, RMAT_GET(&A_work, i, j_0), &b));
    for (u64 x = 0; x < m; x++) {
      TRY(rz_mul(&temp, &temp2, RMAT_GET(&A_work, x, k)));
      TRY(rz_sub(RMAT_GET(&A_work, x, j_0), RMAT_GET(&A_work, x, j_0), &temp));
    }
    for (u64 x = 0; x < n; x++) {
      TRY(rz_mul(&temp, &temp2, RMAT_GET(U, x, k)));
      TRY(rz_sub(RMAT_GET(U, x, j_0), RMAT_GET(U, x, j_0), &temp));
    }
  }

step5:
  // 5. [Finished]
  if (i == l) {
    TRY(rmat_init(W, m, n - k));
    for (u64 x = 0; x < n - k; x++) {
      for (u64 y = 0; y < m; y++) {
        TRY(rz_copy(RMAT_GET(W, y, x), RMAT_GET(&A_work, y, x + k)));
      }
    }
    u_live = false;  // keep U for the caller
    goto cleanup;
  }

  i--;
  k--;
  j = k;
  goto step2;

cleanup:
  rz_clear_multi(&u, &v, &d, &temp, &temp2, &b, &q_k, &q_j, NULL);
  if (C != NULL) {
    for (u64 x = 0; x < n; x++) {
      rz_clear(&C[x]);
    }
    free(C);
  }
  if (B != NULL) {
    for (u64 x = 0; x < m; x++) {
      rz_clear(&B[x]);
    }
    free(B);
  }
  if (u_live) rmat_clear(U);
  rmat_clear(&A_work);
  return err;
}

rabin_err_t rmat_mlll(rmat_t* B, u64 n, double delta, rmat_t* H, u64* rank)
{
  if (B == NULL || H == NULL || rank == NULL) return RABIN_ERR_NULL_PTR;
  if (!(delta > 0.0 && delta <= 1.0)) return RABIN_ERR_INVALID_ARG;
  if (n != B->cols) return RABIN_ERR_INVALID_ARG;

  // the HNF machinery (and Cohen 2.6.8) is stated for d <= n: the
  // column reduction only inspects the bottom n rows
  if (B->rows > n) return RABIN_ERR_INVALID_ARG;

  *rank = 0;

  // no coordinates (or no vectors): every column is the zero vector
  if (n == 0 || B->rows == 0) {
    rabin_err_t e = rmat_init(H, n, n);
    if (e != RABIN_SUCCESS) return e;
    return rmat_id(H, n);
  }

  rmat_t A = {0};  // original generator matrix (kept for B_orig * H checks)
  rmat_t W = {0};  // HNF of the column lattice (m x p)
  rmat_t U = {0};  // unimodular: A * U = [ 0 | W ]
  rmat_t Q = {0};  // LLL transformation of W (p x p)
  rmat_t A_work = {0};
  rmat_t H_lift = {0};  // U[:, n-p..n-1] (n x p)
  rmat_t H_tail = {0};  // H_lift * Q (n x p)
  rmat_t B_new = {0};   // assembled result; swapped into B after the check

  rabin_err_t err = rmat_init(&A, B->rows, n);
  if (err != RABIN_SUCCESS) return err;
  err = rmat_copy(&A, B);
  if (err != RABIN_SUCCESS) goto cleanup;

  // 1. [Rank-revealing column reduction]
  err = hnf_with_transform(&W, &U, &A);
  if (err != RABIN_SUCCESS) goto cleanup;

  u64 p = W.cols;
  *rank = p;

  // 2. [LLL-reduce the lattice basis]
  if (p > 0) {
    TRY(rmat_lll(&W, p, delta, &Q));
  }

  // 3. [Assemble B = [ 0 | W ] and H = [ U_rel | U_lift * Q ]]
  err = rmat_init(H, n, n);
  if (err != RABIN_SUCCESS) goto cleanup;
  err = rmat_id(H, n);
  if (err != RABIN_SUCCESS) goto cleanup;
  for (u64 x = 0; x < n - p; x++) {
    for (u64 y = 0; y < n; y++) {
      TRY(rmat_set(H, RMAT_GET(&U, y, x), y, x));
    }
  }
  if (p > 0) {
    TRY(rmat_init(&H_lift, n, p));
    for (u64 x = 0; x < p; x++) {
      for (u64 y = 0; y < n; y++) {
        TRY(rmat_set(&H_lift, RMAT_GET(&U, y, n - p + x), y, x));
      }
    }
    TRY(rmat_init(&H_tail, n, p));
    TRY(rmat_mul(&H_tail, &H_lift, &Q));
    for (u64 x = 0; x < p; x++) {
      for (u64 y = 0; y < n; y++) {
        TRY(rmat_set(H, RMAT_GET(&H_tail, y, x), y, n - p + x));
      }
    }
  }
  // B_new = [ 0_{m x (n-p)} | W ] (the caller's B stays untouched until
  // the final commit, so a failure never corrupts it)
  TRY(rmat_init(&B_new, B->rows, n));
  for (u64 x = 0; x < p; x++) {
    for (u64 y = 0; y < B->rows; y++) {
      TRY(rmat_set(&B_new, RMAT_GET(&W, y, x), y, n - p + x));
    }
  }

  // 4. [Verify B_new == A_orig * H before committing]
  TRY(rmat_init(&A_work, B->rows, n));
  TRY(rmat_mul(&A_work, &A, H));
  if (!rmat_equal(&A_work, &B_new)) {
    err = RABIN_ERR_INVALID_ARG;
    goto cleanup;
  }

  // 5. [Commit: B now holds the reduced basis]
  rmat_swap(B, &B_new);

  err = RABIN_SUCCESS;
cleanup:
  rmat_clear(&B_new);
  rmat_clear(&A_work);
  rmat_clear(&H_tail);
  rmat_clear(&H_lift);
  rmat_clear(&Q);
  rmat_clear(&U);
  rmat_clear(&W);
  rmat_clear(&A);
  return err;
}
