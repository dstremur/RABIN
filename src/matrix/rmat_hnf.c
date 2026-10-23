/*
 * rmat_hnf.c
 *
 * Hermite normal form for rmat_t: exact and modular HNF with optional
 * transformation tracking, plus HNF verification predicates.
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
static rabin_err_t rmat_neg_col(rmat_t* A, u64 j)
{
  if (A == NULL) return RABIN_ERR_NULL_PTR;
  if (j >= A->cols) return RABIN_ERR_INVALID_ARG;

  for (u64 i = 0; i < A->rows; i++) {
    rabin_err_t err = rz_neg(RMAT_GET(A, i, j), RMAT_GET(A, i, j));
    if (err != RABIN_SUCCESS) return err;
  }
  return RABIN_SUCCESS;
}

static bool rmat_row_prefix_is_zero(const rmat_t* A, u64 i, u64 k)
{
  for (u64 j = 0; j < k; j++) {
    if (!rz_is_zero(RMAT_GET(A, i, j))) {
      return false;
    }
  }
  return true;
}

static u64 rmat_min_abs_index(const rmat_t* A, int i, int k)
{
  int j, j0 = -1;
  for (j = 0; j <= k; j++) {
    if (rz_is_zero(RMAT_GET(A, i, j))) continue;

    if (j0 == -1 || rz_cmp_abs(RMAT_GET(A, i, j), RMAT_GET(A, i, j0)) < 0)
      j0 = j;
  }

  return j0;
}

static rabin_err_t rmat_swap_col(rmat_t* A, u64 k1, u64 k2)
{
  if (A == NULL) return RABIN_ERR_NULL_PTR;
  if (k1 >= A->cols || k2 >= A->cols) return RABIN_ERR_INVALID_ARG;

  for (u64 r = 0; r < A->rows; r++) {
    rabin_err_t err = rz_swap(RMAT_GET(A, r, k1), RMAT_GET(A, r, k2));
    if (err != RABIN_SUCCESS) return err;
  }
  return RABIN_SUCCESS;
}

rabin_err_t rmat_hermite(rmat_t* W, const rmat_t* A)
{
  if (W == NULL || A == NULL) return RABIN_ERR_NULL_PTR;

  // 1. [Initialize]
  u64 m = A->rows;
  u64 n = A->cols;
  if (m == 0 || n == 0) return RABIN_SUCCESS;

  i64 i = m - 1;
  i64 k = n - 1;
  u64 l = 0;
  if (m > n) {
    l = m - n;
  }

  rabin_err_t err = RABIN_SUCCESS;
  rmat_t A_work = {0};
  TRY(rmat_init(&A_work, m, n));
  TRY(rmat_copy(&A_work, A));

  rz_t b, q, rem, temp, one, t2;
  rz_init_multi(&b, &q, &rem, &temp, &one, &t2, NULL);
  TRY(rz_set_u64(&one, 1));

// 2. [Row finished?]
step2:
  if (rmat_row_prefix_is_zero(&A_work, i, k)) {
    if (RMAT_GET(&A_work, i, k)->is_neg) {
      TRY(rmat_neg_col(&A_work, k));
    }
    goto step5;
  }

  // 3. [Choose non-zero entry]
  u64 j_0 = rmat_min_abs_index(&A_work, i, k);
  if ((i64)j_0 < k) {
    TRY(rmat_swap_col(&A_work, k, j_0));
  }
  if (RMAT_GET(&A_work, i, k)->is_neg) {
    TRY(rmat_neg_col(&A_work, k));
  }

  TRY(rz_copy(&b, RMAT_GET(&A_work, i, k)));

  // 4. [Reduce]
  for (u64 j = 0; j < (u64)k; j++) {
    TRY(rz_div_euclid(&q, RMAT_GET(&A_work, i, j), &b));
    for (u64 x = 0; x < m; x++) {
      TRY(rz_mul(&temp, &q, RMAT_GET(&A_work, x, k)));
      TRY(rz_sub(&t2, RMAT_GET(&A_work, x, j), &temp));
      TRY(rmat_set(&A_work, &t2, x, j));
    }
  }
  goto step2;

// 5. [Final reductions]
step5:
  TRY(rz_copy(&b, RMAT_GET(&A_work, i, k)));
  if (rz_is_zero(&b)) {
    k++;
    goto step6;
  }

  for (u64 j = k + 1; j < n; j++) {
    TRY(rz_div_euclid(&q, RMAT_GET(&A_work, i, j), &b));
    for (u64 x = 0; x < m; x++) {
      TRY(rz_mul(&temp, &q, RMAT_GET(&A_work, x, k)));
      TRY(rz_sub(&t2, RMAT_GET(&A_work, x, j), &temp));
      TRY(rmat_set(&A_work, &t2, x, j));
    }
  }

// 6. [Finished?]
step6:
  if (i != (i64)l) {
    i--;
    k--;
    goto step2;
  }

  TRY(rmat_init(W, m, n - k));

  // copy result
  for (u64 j = 0; j < n - k; j++) {
    for (u64 x = 0; x < m; x++) {
      TRY(rz_copy(RMAT_GET(W, x, j), RMAT_GET(&A_work, x, j + k)));
    }
  }

cleanup:
  rz_clear_multi(&b, &q, &rem, &temp, &one, &t2, NULL);
  rmat_clear(&A_work);
  return err;
}

static rabin_err_t rmat_hnf_centered_mod(rz_t* r, const rz_t* x, const rz_t* R,
                                         const rz_t* half)
{
  if (r == NULL || x == NULL || R == NULL || half == NULL)
    return RABIN_ERR_NULL_PTR;

  if (rz_cmp_abs(x, half) <= 0) {
    return rz_copy(r, x);
  }

  rabin_err_t err = rz_mod_pos(r, x, R);
  if (err != RABIN_SUCCESS) return err;
  if (rz_cmp(r, half) > 0) {
    if ((err = rz_sub(r, r, R)) != RABIN_SUCCESS) return err;
  }
  return RABIN_SUCCESS;
}

// Cohen's "Important Remark" after Algorithm 2.4.5: the Euclidean steps of
// the HNF algorithms need the Bezout pair (u, v) of u*a + v*b = d with
// small coefficients. All solutions are (u0 + t*(b/d), v0 - t*(a/d)), so
// this rewrites the pair in place to the unique one with v in the centered
// half-period (-|a/d|/2, |a/d|/2] (which minimizes |v|); when a | b
// (d = |a|) the book's essential condition v = 0, u = sign(a) applies
// instead. The identity u*a + v*b = d is preserved in all cases.
static rabin_err_t hnf_bezout_minimal(rz_t* u, rz_t* v, const rz_t* d,
                                      const rz_t* a, const rz_t* b)
{
  rz_t p, pmag, half, quot, r, t, two;
  rz_init_multi(&p, &pmag, &half, &quot, &r, &t, &two, NULL);
  rabin_err_t err = RABIN_SUCCESS;

  if (rz_is_zero(a)) {
    // v*sign(b) = d/|b| = 1 is forced; u is free, take u = 0
    TRY(rz_set_u64(u, 0));
    if (rz_is_zero(b)) {
      TRY(rz_set_u64(v, 0));
    } else if (b->is_neg) {
      TRY(rz_set_i64(v, -1));
    } else {
      TRY(rz_set_u64(v, 1));
    }
    goto cleanup;
  }

  if (rz_is_zero(b)) {
    // u*sign(a) = 1 is forced; v = 0
    if (a->is_neg) {
      TRY(rz_set_i64(u, -1));
    } else {
      TRY(rz_set_u64(u, 1));
    }
    TRY(rz_set_u64(v, 0));
    goto cleanup;
  }

  TRY(rz_copy(&p, a));
  p.is_neg = false;

  // a | b  =>  d = |a|  =>  v = 0, u = sign(a)  (the book's essential
  // condition; it also keeps the Euclidean step from stalling)
  if (rz_cmp(&p, d) == 0) {
    if (a->is_neg) {
      TRY(rz_set_i64(u, -1));
    } else {
      TRY(rz_set_u64(u, 1));
    }
    TRY(rz_set_u64(v, 0));
    goto cleanup;
  }

  // centered half-period: v = v0 - t*(a/d) with |a/d| >= 2 and
  // v in (-|a/d|/2, |a/d|/2]; recompute u = (d - v*b)/a (exact)
  TRY(rz_div_exact(&p, a, d));
  TRY(rz_copy(&pmag, &p));
  pmag.is_neg = false;
  TRY(rz_set_u64(&two, 2));
  TRY(rz_div(&half, &pmag, &two));

  // r = v0 mod |p| in [0, |p|)
  TRY(rz_divmod(&quot, &r, v, &pmag));
  if (r.is_neg) {
    TRY(rz_add(&r, &r, &pmag));
  }
  if (rz_cmp(&r, &half) > 0) {
    TRY(rz_sub(v, &r, &pmag));
  } else {
    TRY(rz_copy(v, &r));
  }

  TRY(rz_mul(&t, v, b));
  TRY(rz_sub(&t, d, &t));
  TRY(rz_div_exact(u, &t, a));

cleanup:
  rz_clear_multi(&p, &pmag, &half, &quot, &r, &t, &two, NULL);
  return err;
}

/*
 * One [Euclidean step] of Cohen 2.4.5 (shared by rmat_hermite_gcd and
 * hnf_with_transform): with pivot row i and columns j < k, compute the
 * minimal Bezout pair (u, v, d) of (a_{i,k}, a_{i,j}) and apply the
 * column operation
 *
 *   col_k <- u*col_k + v*col_j                       (places d at row i)
 *   col_j <- (a_{i,k}/d)*col_j - (a_{i,j}/d)*col_k  (places 0 at row i)
 *
 * buf must hold A_work->rows initialized scratch rz_t values. If
 * U != NULL (with scratch vector C of U->rows values, both non-NULL or
 * both NULL) the identical operations are applied to columns k and j of
 * U, maintaining the invariant orig * U = A_work.
 */
rabin_err_t hnf_euclid_step(rmat_t* A_work, rmat_t* U, rz_t* buf, rz_t* C,
                            u64 i, u64 j, u64 k, rz_t* u, rz_t* v, rz_t* d,
                            rz_t* q_k, rz_t* q_j, rz_t* temp, rz_t* temp2)
{
  if (A_work == NULL || buf == NULL || u == NULL || v == NULL || d == NULL ||
      q_k == NULL || q_j == NULL || temp == NULL || temp2 == NULL ||
      (U == NULL) != (C == NULL))
    return RABIN_ERR_NULL_PTR;

  rabin_err_t err = rz_gcd_extended_lehmer(u, v, d, RMAT_GET(A_work, i, k),
                                           RMAT_GET(A_work, i, j));
  if (err != RABIN_SUCCESS) return err;
  if ((err = hnf_bezout_minimal(u, v, d, RMAT_GET(A_work, i, k),
                                RMAT_GET(A_work, i, j))) != RABIN_SUCCESS)
    return err;

  // buf[x] = u*A[x][k] + v*A[x][j]
  for (u64 x = 0; x < A_work->rows; x++) {
    if ((err = rz_mul(&buf[x], u, RMAT_GET(A_work, x, k))) != RABIN_SUCCESS)
      return err;
    if ((err = rz_mul(temp, v, RMAT_GET(A_work, x, j))) != RABIN_SUCCESS)
      return err;
    if ((err = rz_add(&buf[x], &buf[x], temp)) != RABIN_SUCCESS) return err;
  }

  if ((err = rz_div_euclid(q_k, RMAT_GET(A_work, i, k), d)) != RABIN_SUCCESS)
    return err;
  if ((err = rz_div_euclid(q_j, RMAT_GET(A_work, i, j), d)) != RABIN_SUCCESS)
    return err;

  // col_j <- q_k*col_j - q_j*col_k,  col_k <- buf
  for (u64 x = 0; x < A_work->rows; x++) {
    if ((err = rz_mul(temp, q_k, RMAT_GET(A_work, x, j))) != RABIN_SUCCESS)
      return err;
    if ((err = rz_mul(temp2, q_j, RMAT_GET(A_work, x, k))) != RABIN_SUCCESS)
      return err;
    if ((err = rz_sub(RMAT_GET(A_work, x, j), temp, temp2)) != RABIN_SUCCESS)
      return err;
    if ((err = rz_copy(RMAT_GET(A_work, x, k), &buf[x])) != RABIN_SUCCESS)
      return err;
  }

  if (U != NULL) {
    for (u64 x = 0; x < U->rows; x++) {
      if ((err = rz_mul(&C[x], u, RMAT_GET(U, x, k))) != RABIN_SUCCESS)
        return err;
      if ((err = rz_mul(temp, v, RMAT_GET(U, x, j))) != RABIN_SUCCESS)
        return err;
      if ((err = rz_add(&C[x], &C[x], temp)) != RABIN_SUCCESS) return err;
    }
    for (u64 x = 0; x < U->rows; x++) {
      if ((err = rz_mul(temp, q_k, RMAT_GET(U, x, j))) != RABIN_SUCCESS)
        return err;
      if ((err = rz_mul(temp2, q_j, RMAT_GET(U, x, k))) != RABIN_SUCCESS)
        return err;
      if ((err = rz_sub(RMAT_GET(U, x, j), temp, temp2)) != RABIN_SUCCESS)
        return err;
      if ((err = rz_copy(RMAT_GET(U, x, k), &C[x])) != RABIN_SUCCESS)
        return err;
    }
  }
  return RABIN_SUCCESS;
}

rabin_err_t rmat_hermite_mod_d(rmat_t* W, const rmat_t* A, const rz_t* D)
{
  if (W == NULL || A == NULL || D == NULL) return RABIN_ERR_NULL_PTR;
  if (rz_is_zero(D)) return RABIN_ERR_DIV_BY_ZERO;

  // 1. [Initialize]
  u64 m = A->rows;
  u64 n = A->cols;
  if (m == 0 || n == 0) return RABIN_SUCCESS;

  if (m > n) {
    return rmat_hermite_gcd(W, A);
  }

  rabin_err_t err = RABIN_SUCCESS;
  rmat_t A_work = {0}, W_local = {0};
  rz_t R, u, v, d, q, t1, t2, t3, half;
  rz_init_multi(&R, &u, &v, &d, &q, &t1, &t2, &t3, &half, NULL);
  rz_t* B = NULL;

  TRY(rmat_init(&A_work, m, n));
  TRY(rmat_copy(&A_work, A));
  TRY(rz_copy(&R, D));
  R.is_neg = false;
  TRY(rz_rshift(&half, &R, 1));

  // reduce the input entries into (-R/2, R/2]; all column operations are
  // taken modulo R, so only the classes modulo R matter
  for (u64 x = 0; x < m * n; x++) {
    TRY(rmat_hnf_centered_mod(&A_work.data[x], &A_work.data[x], &R, &half));
  }

  if (m > 0) {
    B = malloc(m * sizeof(rz_t));
    if (B == NULL) {
      err = RABIN_ERR_OUT_OF_MEMORY;
      goto cleanup;
    }
    for (u64 x = 0; x < m; x++) {
      rz_init(&B[x]);
    }
  }

  TRY(rmat_init(&W_local, m, m));

  for (i64 i = (i64)m - 1; i >= 0; i--) {
    i64 k = (i64)n - (i64)m + i;

    // 2. [Check zero] / 3. [Euclidean step]
    for (i64 j = k - 1; j >= 0; j--) {
      if (rz_is_zero(RMAT_GET(&A_work, i, j))) {
        continue;
      }

      TRY(rz_gcd_extended_lehmer(&u, &v, &d, RMAT_GET(&A_work, i, k),
                                 RMAT_GET(&A_work, i, j)));
      TRY(hnf_bezout_minimal(&u, &v, &d, RMAT_GET(&A_work, i, k),
                             RMAT_GET(&A_work, i, j)));

      // B = u*A_k + v*A_j, reduced into (-R/2, R/2]
      for (u64 x = 0; x < m; x++) {
        TRY(rz_mul(&B[x], &u, RMAT_GET(&A_work, x, k)));
        TRY(rz_mul(&t2, &v, RMAT_GET(&A_work, x, j)));
        TRY(rz_add(&B[x], &B[x], &t2));
        TRY(rmat_hnf_centered_mod(&B[x], &B[x], &R, &half));
      }

      // A_j = (a_{i,k}/d)*A_j - (a_{i,j}/d)*A_k, reduced into (-R/2, R/2]
      TRY(rz_div_exact(&t2, RMAT_GET(&A_work, i, k), &d));
      TRY(rz_div_exact(&q, RMAT_GET(&A_work, i, j), &d));
      for (u64 x = 0; x < m; x++) {
        TRY(rz_mul(&t3, &t2, RMAT_GET(&A_work, x, j)));
        TRY(rz_mul(&t1, &q, RMAT_GET(&A_work, x, k)));
        TRY(rz_sub(&t1, &t3, &t1));
        TRY(rmat_hnf_centered_mod(RMAT_GET(&A_work, x, j), &t1, &R, &half));
        TRY(rz_copy(RMAT_GET(&A_work, x, k), &B[x]));
      }

      // a_{i,k} = u*a_{i,k} + v*a_{i,j} = d exactly
      TRY(rz_copy(RMAT_GET(&A_work, i, k), &d));
    }

    // 4. [Next row]
    // u*a_{i,k} + v*R = d = gcd(a_{i,k}, R)
    TRY(rz_gcd_extended_lehmer(&u, &v, &d, RMAT_GET(&A_work, i, k), &R));
    TRY(hnf_bezout_minimal(&u, &v, &d, RMAT_GET(&A_work, i, k), &R));

    // W_i = u*A_k mod R, taken in [0, R)
    for (u64 x = 0; x < m; x++) {
      TRY(rz_mul(&t1, &u, RMAT_GET(&A_work, x, k)));
      TRY(rz_mod_pos(RMAT_GET(&W_local, x, i), &t1, &R));
    }

    // if d = R (i.e. R | a_{i,k}) the diagonal entry would be 0; set it to d
    if (rz_cmp(&d, &R) == 0) {
      TRY(rz_copy(RMAT_GET(&W_local, i, i), &d));
    }

    // final reductions: W_j -= floor(W_{i,j}/W_{i,i}) * W_i for j > i
    for (u64 j = (u64)i + 1; j < m; j++) {
      if (rz_is_zero(RMAT_GET(&W_local, i, j))) {
        continue;
      }
      TRY(rz_div_euclid(&q, RMAT_GET(&W_local, i, j),
                        RMAT_GET(&W_local, i, i)));
      for (u64 x = 0; x <= (u64)i; x++) {
        TRY(rz_mul(&t1, &q, RMAT_GET(&W_local, x, i)));
        TRY(rz_sub(RMAT_GET(&W_local, x, j), RMAT_GET(&W_local, x, j), &t1));
      }
    }

    R.is_neg = false;
    TRY(rz_div_exact(&R, &R, &d));
    TRY(rz_rshift(&half, &R, 1));

    if (i > 0) {
      // working modulo R, a_{i-1,k-1} may have reduced to zero; replace it
      // by any nonzero multiple of R
      if (rz_is_zero(RMAT_GET(&A_work, i - 1, k - 1))) {
        TRY(rz_copy(RMAT_GET(&A_work, i - 1, k - 1), &R));
      }
    }
  }

  TRY(rmat_swap(W, &W_local));
  rmat_clear(&W_local);

cleanup:
  if (B != NULL) {
    for (u64 x = 0; x < m; x++) {
      rz_clear(&B[x]);
    }
    free(B);
  }
  rz_clear_multi(&R, &u, &v, &d, &q, &t1, &t2, &t3, &half, NULL);
  rmat_clear(&A_work);
  rmat_clear(&W_local);
  return err;
}

rabin_err_t rmat_hermite_gcd(rmat_t* W, const rmat_t* A)
{
  if (W == NULL || A == NULL) return RABIN_ERR_NULL_PTR;

  // 1. [Initialize]
  u64 m = A->rows;
  u64 n = A->cols;
  if (m == 0 || n == 0) return RABIN_SUCCESS;

  u64 i = m - 1;
  u64 j = n - 1;
  u64 k = n - 1;
  u64 l = 0;
  if (m > n) {
    l = m - n;
  }

  rabin_err_t err = RABIN_SUCCESS;
  rmat_t A_work = {0};
  rz_t* B = NULL;

  TRY(rmat_init(&A_work, m, n));
  TRY(rmat_copy(&A_work, A));

  B = malloc(m * sizeof(rz_t));
  if (B == NULL) {
    err = RABIN_ERR_OUT_OF_MEMORY;
    goto cleanup;
  }
  for (u64 x = 0; x < m; x++) {
    rz_init(&B[x]);
  }

  rz_t u, v, d, temp, temp2, temp3, one, b, q_k, q_j;
  rz_init_multi(&u, &v, &d, &temp, &temp2, &temp3, &one, &b, &q_k, &q_j, NULL);
  TRY(rz_set_u64(&one, 1));

// 2. [Check zero]
step2:
  if (j == 0) {
    goto step4;
  }
  j--;
  if (rz_is_zero(RMAT_GET(&A_work, i, j))) {
    goto step2;
  }

  // 3. [Euclidean step]
  TRY(hnf_euclid_step(&A_work, NULL, B, NULL, i, j, k, &u, &v, &d, &q_k, &q_j,
                      &temp, &temp2));
  goto step2;

  // 4. [Final reduction]
step4:
  TRY(rz_copy(&b, RMAT_GET(&A_work, i, k)));
  if (b.is_neg) {  // zero is never negative, no extra rz_is_zero test needed
    for (u64 x = 0; x < m; x++) {
      TRY(rz_neg(RMAT_GET(&A_work, x, k), RMAT_GET(&A_work, x, k)));
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
    goto cleanup;
  }

  i--;
  k--;
  j = k;
  goto step2;

cleanup:
  if (B != NULL) {
    for (u64 x = 0; x < m; x++) {
      rz_clear(&B[x]);
    }
    free(B);
  }
  rmat_clear(&A_work);
  rz_clear_multi(&u, &v, &d, &temp, &temp2, &temp3, &one, &b, &q_k, &q_j, NULL);
  return err;
}

bool rmat_hnf_check_structure(const rmat_t* H)
{
  if (H == NULL) return false;

  u64 r = H->rows;
  u64 c = H->cols;

  // once a zero row is hit, every later row must be zero as well
  bool zero_row_seen = false;

  for (u64 i = 0; i < r; i++) {
    rz_t* row = &H->data[i * c];

    // upper triangular: H[i][j] == 0 for all j < i (capped at c for
    // rows below the last column)
    u64 left = RZ_MIN(i, c);
    for (u64 j = 0; j < left; j++) {
      if (!rz_is_zero(&row[j])) {
        return false;
      }
    }

    // zero rows only in a trailing block: from here on everything must
    // be zero
    if (zero_row_seen) {
      for (u64 j = 0; j < c; j++) {
        if (!rz_is_zero(&row[j])) {
          return false;
        }
      }
      continue;
    }

    // a nonzero row needs a diagonal pivot, so rows at or beyond the
    // last column can only be zero rows
    if (i >= c) {
      for (u64 j = 0; j < c; j++) {
        if (!rz_is_zero(&row[j])) {
          return false;
        }
      }
      zero_row_seen = true;
      continue;
    }

    rz_t* pivot = &row[i];

    if (rz_is_zero(pivot)) {
      // a zero pivot may only start the trailing zero block: everything
      // right of the pivot must also be zero
      for (u64 j = i + 1; j < c; j++) {
        if (!rz_is_zero(&row[j])) {
          return false;
        }
      }
      zero_row_seen = true;
      continue;
    }

    // positive pivot
    if (pivot->is_neg) {
      return false;
    }

    // reduced entries: the entries right of the pivot, H[i][j] for
    // j > i, must lie in [0, pivot) — this is the reduction invariant
    // established by the column-operation HNF algorithms in this
    // library (each row is reduced modulo its own pivot)
    for (u64 j = i + 1; j < c; j++) {
      if (row[j].is_neg || rz_cmp(&row[j], pivot) >= 0) {
        return false;
      }
    }
  }

  return true;
}

bool rmat_hnf_check_transformation(const rmat_t* A, const rmat_t* H,
                                   const rmat_t* U)
{
  // column-operation convention: H = A * U; reject shape mismatches
  // before allocating the product
  if (A == NULL || H == NULL || U == NULL) return false;
  if (A->rows != H->rows || A->cols != U->rows || U->cols != H->cols) {
    return false;
  }

  rmat_t P = {0};
  if (rmat_init(&P, H->rows, H->cols) != RABIN_SUCCESS) return false;

  if (rmat_mul(&P, A, U) != RABIN_SUCCESS) {
    rmat_clear(&P);
    return false;
  }

  bool ok = rmat_equal(&P, H);

  rmat_clear(&P);
  return ok;
}

bool rmat_hnf_check_unimodular(const rmat_t* U)
{
  if (U == NULL) return false;
  if (!rmat_is_square(U)) {
    return false;
  }
  // the empty matrix acts as the identity transformation
  if (U->rows == 0) {
    return true;
  }

  rz_t det;
  rz_init(&det);

  if (rmat_det(&det, U) != RABIN_SUCCESS) {
    rz_clear(&det);
    return false;
  }

  // |det| == 1: clear the sign flag and use the O(1) one-check
  det.is_neg = false;
  bool ok = rz_is_one(&det);

  rz_clear(&det);
  return ok;
}

bool rmat_hnf_verify(const rmat_t* A, const rmat_t* H, const rmat_t* U)
{
  if (A == NULL || H == NULL || U == NULL) return false;

  // order matters: the cheap structure check runs first, the
  // multiplication and the determinant only if it passes
  return rmat_hnf_check_structure(H) &&
         rmat_hnf_check_transformation(A, H, U) && rmat_hnf_check_unimodular(U);
}
