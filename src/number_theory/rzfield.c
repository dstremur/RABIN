/*
 * rzfield.c
 *
 * Arithmetic in Z_m and in the polynomial ring Z_m[x]/(q).
 *
 * This file implements a higher layer on top of the rz_t and polynomial
 * types:
 *
 *   - field_ctx_t: arithmetic modulo m. When m is odd and > 1 the elements
 *     live in the Montgomery domain and the operations use REDC; for even
 *     m the elements live in plain form [0, m) and the operations use
 *     plain multiply + reduce. The representation is an internal detail.
 *   - rpol_divmod / rpol_xgcd: polynomial long division and extended
 *     GCD over the coefficient ring.
 *   - poly_ring_t: the ring of polynomials modulo q with coefficients in
 *     Z_m.
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

#include "../../include/rzfield.h"

#include <assert.h>
#include <stdlib.h>
#include <string.h>

/*===========================================================================
 *  context management
 *=========================================================================*/

rabin_err_t field_ctx_init(field_ctx_t* ctx, const rz_t* m)
{
  if (ctx == NULL || m == NULL) return RABIN_ERR_NULL_PTR;
  if (m->is_neg || rz_is_zero(m) || rz_is_eq_i64(m, 1))
    return RABIN_ERR_INVALID_ARG;

  ctx->q = malloc(sizeof(rz_t));
  if (ctx->q == NULL) return RABIN_ERR_OUT_OF_MEMORY;
  rz_init(ctx->q);

  rabin_err_t err = rz_copy(ctx->q, m);
  if (err != RABIN_SUCCESS) goto out_free;

  ctx->mont = !rz_is_even(m);

  if (ctx->mont) {
    if ((err = rz_mont_ctx_init(&ctx->mctx, m)) != RABIN_SUCCESS) goto out_free;
  }

  err = RABIN_SUCCESS;
  return err;
out_free:
  rz_clear(ctx->q);
  free(ctx->q);
  ctx->q = NULL;
  return err;
}

rabin_err_t field_ctx_clear(field_ctx_t* ctx)
{
  if (ctx == NULL) return RABIN_ERR_NULL_PTR;

  if (ctx->mont) {
    rz_mont_ctx_clear(&ctx->mctx);
  }
  rz_clear(ctx->q);
  free(ctx->q);
  ctx->q = NULL;
  ctx->mont = false;
  return RABIN_SUCCESS;
}

/*===========================================================================
 *  Z_m element operations
 *=========================================================================*/

rabin_err_t field_in(rz_t* r, const rz_t* a, field_ctx_t* ctx)
{
  if (r == NULL || a == NULL || ctx == NULL) return RABIN_ERR_NULL_PTR;

  if (ctx->mont) {
    return rz_mont_in(r, a, &ctx->mctx);
  }

  rabin_err_t err = rz_mod(r, a, ctx->q);
  if (err != RABIN_SUCCESS) return err;
  if (r->is_neg) {
    err = rz_add(r, r, ctx->q);
  }
  return err;
}

rabin_err_t field_out(rz_t* r, const rz_t* a, field_ctx_t* ctx)
{
  if (r == NULL || a == NULL || ctx == NULL) return RABIN_ERR_NULL_PTR;

  if (ctx->mont) {
    return rz_mont_out(r, a, &ctx->mctx);
  }

  return rz_copy(r, a);
}

rabin_err_t field_set_u64(rz_t* r, u64 v, field_ctx_t* ctx)
{
  if (r == NULL || ctx == NULL) return RABIN_ERR_NULL_PTR;

  rz_t tmp;
  rz_init(&tmp);
  rabin_err_t err = rz_set_u64(&tmp, v);
  if (err == RABIN_SUCCESS) {
    err = field_in(r, &tmp, ctx);
  }
  rz_clear(&tmp);
  return err;
}

rabin_err_t field_add(rz_t* r, const rz_t* a, const rz_t* b, field_ctx_t* ctx)
{
  if (r == NULL || a == NULL || b == NULL || ctx == NULL)
    return RABIN_ERR_NULL_PTR;

  rabin_err_t err = rz_add(r, a, b);
  if (err != RABIN_SUCCESS) return err;
  if (rz_cmp(r, ctx->q) >= 0) {
    err = rz_sub_abs(r, r, ctx->q);
  }
  return err;
}

rabin_err_t field_sub(rz_t* r, const rz_t* a, const rz_t* b, field_ctx_t* ctx)
{
  if (r == NULL || a == NULL || b == NULL || ctx == NULL)
    return RABIN_ERR_NULL_PTR;

  rabin_err_t err = rz_sub(r, a, b);
  if (err != RABIN_SUCCESS) return err;
  if (r->is_neg) {
    err = rz_add(r, r, ctx->q);
  }
  return err;
}

rabin_err_t field_neg(rz_t* r, const rz_t* a, field_ctx_t* ctx)
{
  if (r == NULL || a == NULL || ctx == NULL) return RABIN_ERR_NULL_PTR;

  rabin_err_t err = rz_sub(r, ctx->q, a);
  if (err != RABIN_SUCCESS) return err;
  if (rz_cmp(r, ctx->q) >= 0) {
    err = rz_sub_abs(r, r, ctx->q);
  }
  return err;
}

rabin_err_t field_mul(rz_t* r, const rz_t* a, const rz_t* b, field_ctx_t* ctx)
{
  if (r == NULL || a == NULL || b == NULL || ctx == NULL)
    return RABIN_ERR_NULL_PTR;

  if (ctx->mont) {
    return rz_mont_mul(r, a, b, &ctx->mctx);
  }

  rz_t t;
  rz_init(&t);
  rabin_err_t err = rz_mul(&t, a, b);
  if (err == RABIN_SUCCESS) {
    err = rz_mod(r, &t, ctx->q);
  }
  rz_clear(&t);
  return err;
}

/**
 * @brief Modular inverse via the classical extended Euclidean algorithm.
 *
 * Let n = m->size, measured in 64-bit limbs.
 *
 * Computes res = a^(-1) mod m (in [0, m)) using the iterative extended
 * Euclidean algorithm, which is valid for ANY modulus (odd or even) as
 * long as gcd(a, m) == 1. Unlike the binary extended GCD
 * (rz_mod_inverse), it never divides a tracked coefficient by 2, so it
 * is correct for even moduli where 2 is not invertible.
 *
 * Complexity:
 *   Time: O(n^3) - O(n) division steps, each an O(n^2) rz_t division
 *   Auxiliary memory: O(n) limbs
 *   Output memory: O(n) limbs
 *
 * @param[out] res Receives a^(-1) mod m if it exists.
 * @param[in]  a   Value to invert (nonnegative).
 * @param[in]  m   Modulus (nonzero).
 *
 * @return true  If a is a unit (gcd(a, m) == 1); res is set.
 * @return false If a is zero, m is zero or one, or a and m are not
 *               coprime; res is left unchanged.
 */
static rabin_err_t u64_mod_inverse_euclid(rz_t* res, const rz_t* a,
                                          const rz_t* m)
{
  if (res == NULL || a == NULL || m == NULL) return RABIN_ERR_NULL_PTR;
  if (rz_is_zero(a) || rz_is_zero(m) || rz_is_eq_i64(m, 1))
    return RABIN_ERR_INVALID_ARG;

  // (old_r, r) = (a, m); (old_x, x) = (1, 0)
  // Invariant: old_x * a + old_y * m = old_r (old_y implicit)
  rz_t old_r, r, old_x, x, q, tmp;
  rz_init_multi(&old_r, &r, &old_x, &x, &q, &tmp, NULL);

  rabin_err_t err = rz_copy(&old_r, a);
  if (err == RABIN_SUCCESS) err = rz_copy(&r, m);
  if (err == RABIN_SUCCESS) err = rz_set_u64(&old_x, 1);
  if (err == RABIN_SUCCESS) err = rz_set_u64(&x, 0);
  if (err != RABIN_SUCCESS) goto out;

  while (!rz_is_zero(&r)) {
    if ((err = rz_div(&q, &old_r, &r)) != RABIN_SUCCESS) goto out;

    // (old_r, r) = (r, old_r - q*r)
    if ((err = rz_mul(&tmp, &q, &r)) != RABIN_SUCCESS) goto out;
    if ((err = rz_sub(&tmp, &old_r, &tmp)) != RABIN_SUCCESS) goto out;
    if ((err = rz_copy(&old_r, &r)) != RABIN_SUCCESS) goto out;
    if ((err = rz_copy(&r, &tmp)) != RABIN_SUCCESS) goto out;

    // (old_x, x) = (x, old_x - q*x)
    if ((err = rz_mul(&tmp, &q, &x)) != RABIN_SUCCESS) goto out;
    if ((err = rz_sub(&tmp, &old_x, &tmp)) != RABIN_SUCCESS) goto out;
    if ((err = rz_copy(&old_x, &x)) != RABIN_SUCCESS) goto out;
    if ((err = rz_copy(&x, &tmp)) != RABIN_SUCCESS) goto out;
  }

  // old_r = gcd(a, m); old_x is the Bezout coefficient of a
  if (!rz_is_eq_i64(&old_r, 1)) {
    err = RABIN_ERR_INVALID_ARG;
    goto out;
  }

  if ((err = rz_mod(res, &old_x, m)) != RABIN_SUCCESS) goto out;
  if (res->is_neg) {
    if ((err = rz_add(res, res, m)) != RABIN_SUCCESS) goto out;
  }

  err = RABIN_SUCCESS;
out:
  rz_clear_multi(&old_r, &r, &old_x, &x, &q, &tmp, NULL);
  return err;
}

rabin_err_t field_inv(rz_t* r, const rz_t* a, field_ctx_t* ctx)
{
  if (r == NULL || a == NULL || ctx == NULL) return RABIN_ERR_NULL_PTR;
  if (rz_is_zero(a)) return RABIN_ERR_INVALID_ARG;

  if (!ctx->mont) {
    return u64_mod_inverse_euclid(r, a, ctx->q);
  }

  rz_t plain, inv;
  rz_init_multi(&plain, &inv, NULL);

  rabin_err_t err = rz_mont_out(&plain, a, &ctx->mctx);
  if (err == RABIN_SUCCESS) {
    if ((err = rz_mod_inverse(&inv, &plain, ctx->q)) == RABIN_SUCCESS) {
      err = rz_mont_in(r, &inv, &ctx->mctx);
    }
  }

  rz_clear_multi(&plain, &inv, NULL);
  return err;
}

rabin_err_t field_div(rz_t* r, const rz_t* a, const rz_t* b, field_ctx_t* ctx)
{
  if (r == NULL || a == NULL || b == NULL || ctx == NULL)
    return RABIN_ERR_NULL_PTR;

  rz_t inv;
  rz_init(&inv);

  rabin_err_t err = field_inv(&inv, b, ctx);
  if (err == RABIN_SUCCESS) {
    err = field_mul(r, a, &inv, ctx);
  }

  rz_clear(&inv);
  return err;
}

rabin_err_t field_pow(rz_t* r, const rz_t* a, const rz_t* e, field_ctx_t* ctx)
{
  if (r == NULL || a == NULL || e == NULL || ctx == NULL)
    return RABIN_ERR_NULL_PTR;

  if (ctx->mont) {
    return rz_mont_exp(r, a, e, &ctx->mctx);
  }

  return rz_mod_exp(r, a, e, ctx->q);
}

bool field_is_zero(const rz_t* a) { return a != NULL && rz_is_zero(a); }

bool field_equal(const rz_t* a, const rz_t* b)
{
  return a != NULL && b != NULL && rz_cmp(a, b) == 0;
}

/*===========================================================================
 *  polynomial helpers over the coefficient ring (no reduction mod q)
 *=========================================================================*/

/**
 * @brief Test whether a polynomial is the zero polynomial.
 *
 * Complexity:
 *   Time: O(d) where d = p->deg
 *   Auxiliary memory: O(1)
 *   Output memory: O(1)
 *
 * @param[in] p Polynomial.
 *
 * @return true If all coefficients up to the degree are zero.
 */
static bool poly_is_zero(const rpol_t* p)
{
  for (u64 i = 0; i <= p->deg; i++) {
    if (!rz_is_zero(&p->coeff[i])) {
      return false;
    }
  }
  return true;
}

/**
 * @brief Deep-copy src into dst (dst grows as needed, deg is set).
 *
 * Complexity:
 *   Time: O(d * n) where d = src->deg and n the coefficient size
 *   Auxiliary memory: O(1)
 *   Output memory: O(d) bignums
 *
 * @param[out] dst Destination polynomial.
 * @param[in]  src Source polynomial.
 */
static rabin_err_t poly_assign(rpol_t* dst, const rpol_t* src)
{
  if (dst == NULL || src == NULL) return RABIN_ERR_NULL_PTR;

  rabin_err_t err = rpol_alloc(dst, src->deg);
  if (err != RABIN_SUCCESS) return err;
  for (u64 i = 0; i <= src->deg; i++) {
    if ((err = rz_copy(&dst->coeff[i], &src->coeff[i])) != RABIN_SUCCESS)
      return err;
  }
  dst->deg = src->deg;
  return RABIN_SUCCESS;
}

/**
 * @brief r = a + b in Z_m[x] (coefficients reduced mod m, no q-reduction).
 *
 * Complexity:
 *   Time: O(d * n) where d = max(a->deg, b->deg)
 *   Auxiliary memory: O(1)
 *   Output memory: O(d) bignums
 *
 * @param[out] r Result polynomial.
 * @param[in]  a First polynomial.
 * @param[in]  b Second polynomial.
 * @param[in]  fctx Coefficient ring.
 */
static rabin_err_t poly_fadd(rpol_t* r, const rpol_t* a, const rpol_t* b,
                             field_ctx_t* fctx)
{
  if (r == NULL || a == NULL || b == NULL || fctx == NULL)
    return RABIN_ERR_NULL_PTR;

  u64 max = RZ_MAX(a->deg, b->deg);
  rabin_err_t err = rpol_alloc(r, max);
  if (err != RABIN_SUCCESS) return err;

  for (u64 i = 0; i <= max; i++) {
    if (i <= a->deg && i <= b->deg) {
      if ((err = field_add(&r->coeff[i], &a->coeff[i], &b->coeff[i], fctx)) !=
          RABIN_SUCCESS)
        return err;
    } else if (i <= a->deg) {
      if ((err = rz_copy(&r->coeff[i], &a->coeff[i])) != RABIN_SUCCESS)
        return err;
    } else {
      if ((err = rz_copy(&r->coeff[i], &b->coeff[i])) != RABIN_SUCCESS)
        return err;
    }
  }

  r->deg = max;
  rpol_trim(r);
  return RABIN_SUCCESS;
}

/**
 * @brief r = a - b in Z_m[x] (coefficients reduced mod m, no q-reduction).
 *
 * Complexity:
 *   Time: O(d * n) where d = max(a->deg, b->deg)
 *   Auxiliary memory: O(1)
 *   Output memory: O(d) bignums
 *
 * @param[out] r Result polynomial.
 * @param[in]  a First polynomial.
 * @param[in]  b Second polynomial.
 * @param[in]  fctx Coefficient ring.
 */
static rabin_err_t poly_fsub(rpol_t* r, const rpol_t* a, const rpol_t* b,
                             field_ctx_t* fctx)
{
  if (r == NULL || a == NULL || b == NULL || fctx == NULL)
    return RABIN_ERR_NULL_PTR;

  u64 max = RZ_MAX(a->deg, b->deg);
  rabin_err_t err = rpol_alloc(r, max);
  if (err != RABIN_SUCCESS) return err;

  for (u64 i = 0; i <= max; i++) {
    if (i <= a->deg && i <= b->deg) {
      if ((err = field_sub(&r->coeff[i], &a->coeff[i], &b->coeff[i], fctx)) !=
          RABIN_SUCCESS)
        return err;
    } else if (i <= a->deg) {
      if ((err = rz_copy(&r->coeff[i], &a->coeff[i])) != RABIN_SUCCESS)
        return err;
    } else {
      if ((err = field_neg(&r->coeff[i], &b->coeff[i], fctx)) != RABIN_SUCCESS)
        return err;
    }
  }

  r->deg = max;
  rpol_trim(r);
  return RABIN_SUCCESS;
}

/**
 * @brief r = a * b in Z_m[x] (coefficients reduced mod m, no q-reduction).
 *
 * Schoolbook convolution with the coefficient ring operations.
 *
 * Complexity:
 *   Time: O((a->deg + 1) * (b->deg + 1) * n^2)
 *   Auxiliary memory: O(1)
 *   Output memory: O(a->deg + b->deg) bignums
 *
 * @param[out] r Result polynomial.
 * @param[in]  a First polynomial.
 * @param[in]  b Second polynomial.
 * @param[in]  fctx Coefficient ring.
 */
static rabin_err_t poly_fmul(rpol_t* r, const rpol_t* a, const rpol_t* b,
                             field_ctx_t* fctx)
{
  if (r == NULL || a == NULL || b == NULL || fctx == NULL)
    return RABIN_ERR_NULL_PTR;
  if (a->deg > UINT64_MAX - b->deg) return RABIN_ERR_OVERFLOW;

  u64 r_deg = a->deg + b->deg;
  rabin_err_t err = rpol_alloc(r, r_deg);
  if (err != RABIN_SUCCESS) return err;

  for (u64 i = 0; i <= r_deg; i++) {
    if ((err = rz_set_u64(&r->coeff[i], 0)) != RABIN_SUCCESS) return err;
  }

  rz_t t;
  rz_init(&t);

  for (u64 i = 0; i <= a->deg; i++) {
    for (u64 j = 0; j <= b->deg; j++) {
      if ((err = field_mul(&t, &a->coeff[i], &b->coeff[j], fctx)) !=
          RABIN_SUCCESS)
        goto out;
      if ((err = field_add(&r->coeff[i + j], &r->coeff[i + j], &t, fctx)) !=
          RABIN_SUCCESS)
        goto out;
    }
  }

  r->deg = r_deg;
  rpol_trim(r);
  err = RABIN_SUCCESS;
out:
  rz_clear(&t);
  return err;
}

/*===========================================================================
 *  polynomial division and extended GCD
 *=========================================================================*/

rabin_err_t rpol_divmod(rpol_t* q, rpol_t* r, const rpol_t* a, const rpol_t* b,
                        field_ctx_t* fctx)
{
  if (a == NULL || b == NULL || fctx == NULL) return RABIN_ERR_NULL_PTR;

  // write q = 0, r = a on the "not a valid divisor" errors
  if (poly_is_zero(b)) {
    if (r != NULL) {
      rabin_err_t err = poly_assign(r, a);
      if (err != RABIN_SUCCESS) return err;
    }
    if (q != NULL) {
      rabin_err_t err = rpol_alloc(q, 0);
      if (err == RABIN_SUCCESS) err = rz_set_u64(&q->coeff[0], 0);
      if (err == RABIN_SUCCESS) q->deg = 0;
      if (err != RABIN_SUCCESS) return err;
    }
    return RABIN_ERR_DIV_BY_ZERO;
  }

  rz_t inv_lead;
  rz_init(&inv_lead);
  rabin_err_t err = field_inv(&inv_lead, &b->coeff[b->deg], fctx);
  if (err != RABIN_SUCCESS) {
    rz_clear(&inv_lead);
    if (r != NULL) {
      err = poly_assign(r, a);
      if (err != RABIN_SUCCESS) return err;
    }
    if (q != NULL) {
      err = rpol_alloc(q, 0);
      if (err == RABIN_SUCCESS) err = rz_set_u64(&q->coeff[0], 0);
      if (err == RABIN_SUCCESS) q->deg = 0;
    }
    return err;
  }

  rpol_t rem = {0}, quo = {0};
  if ((err = rpol_init(&rem)) != RABIN_SUCCESS) goto out_inv;
  if ((err = rpol_init(&quo)) != RABIN_SUCCESS) goto out_rem;
  if ((err = poly_assign(&rem, a)) != RABIN_SUCCESS) goto out_rem_q;

  if ((err = rpol_alloc(&quo, a->deg)) != RABIN_SUCCESS) goto out_rem_q;
  for (u64 i = 0; i <= a->deg; i++) {
    if ((err = rz_set_u64(&quo.coeff[i], 0)) != RABIN_SUCCESS) goto out_rem_q;
  }
  quo.deg = 0;

  rz_t factor, t;
  rz_init_multi(&factor, &t, NULL);

  while (rem.deg >= b->deg && !poly_is_zero(&rem)) {
    u64 shift = rem.deg - b->deg;

    // factor = r[deg r] / b[deg b]
    if ((err = field_mul(&factor, &rem.coeff[rem.deg], &inv_lead, fctx)) !=
        RABIN_SUCCESS)
      goto out_rt;

    if ((err = field_add(&quo.coeff[shift], &quo.coeff[shift], &factor,
                         fctx)) != RABIN_SUCCESS)
      goto out_rt;
    if (shift > quo.deg) {
      quo.deg = shift;
    }

    // r -= factor * x^shift * b
    for (u64 i = 0; i <= b->deg; i++) {
      if ((err = field_mul(&t, &factor, &b->coeff[i], fctx)) != RABIN_SUCCESS)
        goto out_rt;
      if ((err = field_sub(&rem.coeff[i + shift], &rem.coeff[i + shift], &t,
                           fctx)) != RABIN_SUCCESS)
        goto out_rt;
    }

    rpol_trim(&rem);
  }

  rpol_trim(&quo);

  if (q != NULL) {
    if ((err = poly_assign(q, &quo)) != RABIN_SUCCESS) goto out_rt;
  }
  if (r != NULL) {
    if ((err = poly_assign(r, &rem)) != RABIN_SUCCESS) goto out_rt;
  }

  rz_clear_multi(&factor, &t, NULL);
out_rt:
  rpol_clear(&quo);
out_rem_q:
  rpol_clear(&rem);
out_rem:
  rpol_clear(&rem);
out_inv:
  rz_clear(&inv_lead);
  return err;
}

rabin_err_t rpol_xgcd(rpol_t* g, rpol_t* x, rpol_t* y, const rpol_t* a,
                      const rpol_t* b, field_ctx_t* fctx)
{
  if (a == NULL || b == NULL || fctx == NULL) return RABIN_ERR_NULL_PTR;

  rpol_t r0 = {0}, r1 = {0}, r2 = {0}, s0 = {0}, s1 = {0}, s2 = {0}, t0 = {0},
         t1 = {0}, t2 = {0}, qq = {0}, tmp = {0};
  rabin_err_t err = rpol_init(&r0);
  if (err == RABIN_SUCCESS) err = rpol_init(&r1);
  if (err == RABIN_SUCCESS) err = rpol_init(&r2);
  if (err == RABIN_SUCCESS) err = rpol_init(&s0);
  if (err == RABIN_SUCCESS) err = rpol_init(&s1);
  if (err == RABIN_SUCCESS) err = rpol_init(&s2);
  if (err == RABIN_SUCCESS) err = rpol_init(&t0);
  if (err == RABIN_SUCCESS) err = rpol_init(&t1);
  if (err == RABIN_SUCCESS) err = rpol_init(&t2);
  if (err == RABIN_SUCCESS) err = rpol_init(&qq);
  if (err == RABIN_SUCCESS) err = rpol_init(&tmp);
  if (err != RABIN_SUCCESS) goto out_r0;

  if ((err = poly_assign(&r0, a)) != RABIN_SUCCESS) goto out_all;
  if ((err = poly_assign(&r1, b)) != RABIN_SUCCESS) goto out_all;

  if ((err = rpol_alloc(&s0, 0)) != RABIN_SUCCESS) goto out_all;
  if ((err = rz_set_u64(&s0.coeff[0], 1)) != RABIN_SUCCESS) goto out_all;
  s0.deg = 0;
  if ((err = rpol_alloc(&t0, 0)) != RABIN_SUCCESS) goto out_all;
  if ((err = rz_set_u64(&t0.coeff[0], 0)) != RABIN_SUCCESS) goto out_all;
  t0.deg = 0;
  if ((err = rpol_alloc(&s1, 0)) != RABIN_SUCCESS) goto out_all;
  if ((err = rz_set_u64(&s1.coeff[0], 0)) != RABIN_SUCCESS) goto out_all;
  s1.deg = 0;
  if ((err = rpol_alloc(&t1, 0)) != RABIN_SUCCESS) goto out_all;
  if ((err = rz_set_u64(&t1.coeff[0], 1)) != RABIN_SUCCESS) goto out_all;
  t1.deg = 0;

  while (!poly_is_zero(&r1)) {
    if ((err = rpol_divmod(&qq, &r2, &r0, &r1, fctx)) != RABIN_SUCCESS)
      goto out_all;

    if ((err = poly_assign(&r0, &r1)) != RABIN_SUCCESS) goto out_all;
    if ((err = poly_assign(&r1, &r2)) != RABIN_SUCCESS) goto out_all;

    // s2 = s0 - qq * s1
    if ((err = poly_fmul(&tmp, &qq, &s1, fctx)) != RABIN_SUCCESS) goto out_all;
    if ((err = poly_fsub(&s2, &s0, &tmp, fctx)) != RABIN_SUCCESS) goto out_all;
    if ((err = poly_assign(&s0, &s1)) != RABIN_SUCCESS) goto out_all;
    if ((err = poly_assign(&s1, &s2)) != RABIN_SUCCESS) goto out_all;

    // t2 = t0 - qq * t1
    if ((err = poly_fmul(&tmp, &qq, &t1, fctx)) != RABIN_SUCCESS) goto out_all;
    if ((err = poly_fsub(&t2, &t0, &tmp, fctx)) != RABIN_SUCCESS) goto out_all;
    if ((err = poly_assign(&t0, &t1)) != RABIN_SUCCESS) goto out_all;
    if ((err = poly_assign(&t1, &t2)) != RABIN_SUCCESS) goto out_all;
  }

  if (g != NULL) {
    if ((err = poly_assign(g, &r0)) != RABIN_SUCCESS) goto out_all;
  }
  if (x != NULL) {
    if ((err = poly_assign(x, &s0)) != RABIN_SUCCESS) goto out_all;
  }
  if (y != NULL) {
    if ((err = poly_assign(y, &t0)) != RABIN_SUCCESS) goto out_all;
  }

  err = RABIN_SUCCESS;
out_all:
  rpol_clear(&tmp);
  rpol_clear(&qq);
  rpol_clear(&t2);
  rpol_clear(&t1);
  rpol_clear(&t0);
  rpol_clear(&s2);
  rpol_clear(&s1);
  rpol_clear(&s0);
  rpol_clear(&r2);
  rpol_clear(&r1);
  rpol_clear(&r0);
  return err;
out_r0:
  rpol_clear(&r0);
  return err;
}

/*===========================================================================
 *  ring Z_m[x]/(q) operations
 *=========================================================================*/

rabin_err_t poly_ring_init(poly_ring_t* ctx, const rpol_t* q, field_ctx_t* fctx)
{
  if (ctx == NULL || q == NULL || fctx == NULL) return RABIN_ERR_NULL_PTR;

  ctx->q = malloc(sizeof(rpol_t));
  if (ctx->q == NULL) return RABIN_ERR_OUT_OF_MEMORY;
  *ctx->q = (rpol_t){0};
  rabin_err_t err = rpol_init(ctx->q);
  if (err == RABIN_SUCCESS) {
    err = poly_assign(ctx->q, q);
  }
  if (err != RABIN_SUCCESS) {
    rpol_clear(ctx->q);
    free(ctx->q);
    ctx->q = NULL;
    ctx->fctx = NULL;
    return err;
  }
  ctx->fctx = fctx;
  return RABIN_SUCCESS;
}

rabin_err_t poly_ring_clear(poly_ring_t* ctx)
{
  if (ctx == NULL) return RABIN_ERR_NULL_PTR;

  rpol_clear(ctx->q);
  free(ctx->q);
  ctx->q = NULL;
  ctx->fctx = NULL;
  return RABIN_SUCCESS;
}

/**
 * @brief Reduce a polynomial modulo the ring modulus (in place into r).
 *
 * Thin wrapper around rpol_divmod() keeping the remainder.
 *
 * Complexity:
 *   Time: O((a->deg - q->deg + 1) * q->deg * n^2)
 *   Auxiliary memory: O(a->deg) bignums
 *   Output memory: O(q->deg) bignums
 *
 * @param[out] r Remainder a mod q.
 * @param[in]  a Polynomial to reduce.
 * @param[in]  ctx Ring context.
 */
static rabin_err_t poly_reduce_mod_q(rpol_t* r, const rpol_t* a,
                                     poly_ring_t* ctx)
{
  return rpol_divmod(NULL, r, a, ctx->q, ctx->fctx);
}

rabin_err_t poly_ring_add(rpol_t* r, const rpol_t* a, const rpol_t* b,
                          poly_ring_t* ctx)
{
  if (r == NULL || a == NULL || b == NULL || ctx == NULL)
    return RABIN_ERR_NULL_PTR;

  rpol_t t = {0};
  rabin_err_t err = rpol_init(&t);
  if (err == RABIN_SUCCESS) err = poly_fadd(&t, a, b, ctx->fctx);
  if (err == RABIN_SUCCESS) err = poly_reduce_mod_q(r, &t, ctx);
  rpol_clear(&t);
  return err;
}

rabin_err_t poly_ring_sub(rpol_t* r, const rpol_t* a, const rpol_t* b,
                          poly_ring_t* ctx)
{
  if (r == NULL || a == NULL || b == NULL || ctx == NULL)
    return RABIN_ERR_NULL_PTR;

  rpol_t t = {0};
  rabin_err_t err = rpol_init(&t);
  if (err == RABIN_SUCCESS) err = poly_fsub(&t, a, b, ctx->fctx);
  if (err == RABIN_SUCCESS) err = poly_reduce_mod_q(r, &t, ctx);
  rpol_clear(&t);
  return err;
}

rabin_err_t poly_ring_neg(rpol_t* r, const rpol_t* a, poly_ring_t* ctx)
{
  if (r == NULL || a == NULL || ctx == NULL) return RABIN_ERR_NULL_PTR;

  rpol_t zero = {0}, t = {0};
  rabin_err_t err = rpol_init(&zero);
  if (err == RABIN_SUCCESS) err = rpol_init(&t);
  if (err != RABIN_SUCCESS) goto out;

  if ((err = rpol_alloc(&zero, 0)) != RABIN_SUCCESS) goto out;
  if ((err = rz_set_u64(&zero.coeff[0], 0)) != RABIN_SUCCESS) goto out;
  zero.deg = 0;

  if ((err = poly_fsub(&t, &zero, a, ctx->fctx)) != RABIN_SUCCESS) goto out;
  err = poly_reduce_mod_q(r, &t, ctx);

out:
  rpol_clear(&zero);
  rpol_clear(&t);
  return err;
}

rabin_err_t poly_ring_mul(rpol_t* r, const rpol_t* a, const rpol_t* b,
                          poly_ring_t* ctx)
{
  if (r == NULL || a == NULL || b == NULL || ctx == NULL)
    return RABIN_ERR_NULL_PTR;

  rpol_t t = {0};
  rabin_err_t err = rpol_init(&t);
  if (err == RABIN_SUCCESS) err = poly_fmul(&t, a, b, ctx->fctx);
  if (err == RABIN_SUCCESS) err = poly_reduce_mod_q(r, &t, ctx);
  rpol_clear(&t);
  return err;
}

rabin_err_t poly_ring_inv(rpol_t* r, const rpol_t* a, poly_ring_t* ctx)
{
  if (r == NULL || a == NULL || ctx == NULL) return RABIN_ERR_NULL_PTR;

  rpol_t g = {0}, x = {0}, y = {0};
  rabin_err_t err = rpol_init(&g);
  if (err == RABIN_SUCCESS) err = rpol_init(&x);
  if (err == RABIN_SUCCESS) err = rpol_init(&y);
  if (err != RABIN_SUCCESS) goto out;

  err = rpol_xgcd(&g, &x, &y, a, ctx->q, ctx->fctx);
  if (err != RABIN_SUCCESS) goto out;

  if (g.deg != 0 || rz_is_zero(&g.coeff[0])) {
    err = RABIN_ERR_INVALID_ARG;
    goto out;
  }

  rz_t g_inv;
  rz_init(&g_inv);
  if ((err = field_inv(&g_inv, &g.coeff[0], ctx->fctx)) != RABIN_SUCCESS) {
    rz_clear(&g_inv);
    goto out;
  }

  rpol_t scaled = {0};
  if ((err = rpol_init(&scaled)) != RABIN_SUCCESS) {
    rz_clear(&g_inv);
    goto out;
  }
  if ((err = rpol_alloc(&scaled, x.deg)) == RABIN_SUCCESS) {
    for (u64 i = 0; i <= x.deg; i++) {
      if ((err = field_mul(&scaled.coeff[i], &x.coeff[i], &g_inv, ctx->fctx)) !=
          RABIN_SUCCESS)
        break;
    }
    if (err == RABIN_SUCCESS) {
      scaled.deg = x.deg;
      rpol_trim(&scaled);
      err = poly_reduce_mod_q(r, &scaled, ctx);
    }
  }
  rpol_clear(&scaled);
  rz_clear(&g_inv);

out:
  rpol_clear(&y);
  rpol_clear(&x);
  rpol_clear(&g);
  return err;
}

rabin_err_t poly_ring_div(rpol_t* r, const rpol_t* a, const rpol_t* b,
                          poly_ring_t* ctx)
{
  if (r == NULL || a == NULL || b == NULL || ctx == NULL)
    return RABIN_ERR_NULL_PTR;

  rpol_t b_inv = {0};
  rabin_err_t err = rpol_init(&b_inv);
  if (err == RABIN_SUCCESS) err = poly_ring_inv(&b_inv, b, ctx);
  if (err == RABIN_SUCCESS) err = poly_ring_mul(r, a, &b_inv, ctx);
  rpol_clear(&b_inv);
  return err;
}

rabin_err_t poly_ring_pow(rpol_t* r, const rpol_t* a, const rz_t* e,
                          poly_ring_t* ctx)
{
  if (r == NULL || a == NULL || e == NULL || ctx == NULL)
    return RABIN_ERR_NULL_PTR;

  rpol_t result = {0}, base = {0}, tmp = {0};
  rabin_err_t err = rpol_init(&result);
  if (err == RABIN_SUCCESS) err = rpol_init(&base);
  if (err == RABIN_SUCCESS) err = rpol_init(&tmp);
  if (err != RABIN_SUCCESS) goto out;

  // result = 1
  if ((err = rpol_alloc(&result, 0)) != RABIN_SUCCESS) goto out;
  if ((err = rz_set_u64(&result.coeff[0], 1)) != RABIN_SUCCESS) goto out;
  result.deg = 0;

  if ((err = poly_assign(&base, a)) != RABIN_SUCCESS) goto out;

  rz_t exp;
  rz_init(&exp);
  if ((err = rz_copy(&exp, e)) != RABIN_SUCCESS) {
    rz_clear(&exp);
    goto out;
  }

  while (!rz_is_zero(&exp)) {
    if (!rz_is_even(&exp)) {
      if ((err = poly_ring_mul(&tmp, &result, &base, ctx)) != RABIN_SUCCESS)
        goto out_exp;
      if ((err = poly_assign(&result, &tmp)) != RABIN_SUCCESS) goto out_exp;
    }
    if ((err = poly_ring_mul(&tmp, &base, &base, ctx)) != RABIN_SUCCESS)
      goto out_exp;
    if ((err = poly_assign(&base, &tmp)) != RABIN_SUCCESS) goto out_exp;
    if ((err = rz_rshift1(&exp)) != RABIN_SUCCESS) goto out_exp;
  }

  if ((err = poly_assign(r, &result)) != RABIN_SUCCESS) goto out_exp;

  err = RABIN_SUCCESS;
out_exp:
  rz_clear(&exp);
out:
  rpol_clear(&tmp);
  rpol_clear(&base);
  rpol_clear(&result);
  return err;
}

bool poly_ring_is_zero(const rpol_t* a) { return a != NULL && poly_is_zero(a); }

bool poly_ring_equal(const rpol_t* a, const rpol_t* b)
{
  if (a == NULL || b == NULL) return false;
  u64 max = RZ_MAX(a->deg, b->deg);
  for (u64 i = 0; i <= max; i++) {
    bool ai_zero = (i > a->deg) || rz_is_zero(&a->coeff[i]);
    bool bi_zero = (i > b->deg) || rz_is_zero(&b->coeff[i]);
    if (ai_zero != bi_zero) {
      return false;
    }
    if (!ai_zero && rz_cmp(&a->coeff[i], &b->coeff[i]) != 0) {
      return false;
    }
  }
  return true;
}
