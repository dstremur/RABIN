/*
 * bigfield.c
 *
 * Arithmetic in Z_m and in the polynomial ring Z_m[x]/(q).
 *
 * This file implements a higher layer on top of the bignum and polynomial
 * types:
 *
 *   - field_ctx: arithmetic modulo m. When m is odd and > 1 the elements
 *     live in the Montgomery domain and the operations use REDC; for even
 *     m the elements live in plain form [0, m) and the operations use
 *     plain multiply + reduce. The representation is an internal detail.
 *   - bigpoly_divmod / bigpoly_xgcd: polynomial long division and extended
 *     GCD over the coefficient ring.
 *   - poly_ring: the ring of polynomials modulo q with coefficients in
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

#include "../../include/bigfield.h"

#include <assert.h>
#include <stdlib.h>
#include <string.h>

/*===========================================================================
 *  context management
 *=========================================================================*/

void field_ctx_init(field_ctx* ctx, const bignum* m)
{
  ctx->q = malloc(sizeof(bignum));
  bn_init(ctx->q);
  bn_copy(ctx->q, m);

  ctx->mont = !bn_is_even(m) && !bn_is_zero(m) && !bn_is_eq_i64(m, 1);

  if (ctx->mont) {
    bn_mont_ctx_init(&ctx->mctx, m);
  }
}

void field_ctx_free(field_ctx* ctx)
{
  if (ctx->mont) {
    bn_mont_ctx_free(&ctx->mctx);
  }
  bn_free(ctx->q);
  free(ctx->q);
  ctx->q = NULL;
  ctx->mont = false;
}

/*===========================================================================
 *  Z_m element operations
 *=========================================================================*/

void field_in(bignum* r, const bignum* a, field_ctx* ctx)
{
  if (ctx->mont) {
    bn_mont_in(r, a, &ctx->mctx);
    return;
  }

  bn_mod(r, a, ctx->q);
  if (r->is_neg) {
    bn_add(r, r, ctx->q);
  }
}

void field_out(bignum* r, const bignum* a, field_ctx* ctx)
{
  if (ctx->mont) {
    bn_mont_out(r, a, &ctx->mctx);
    return;
  }

  bn_copy(r, a);
}

void field_set_u64(bignum* r, u64 v, field_ctx* ctx)
{
  bignum tmp;
  bn_init(&tmp);
  bn_set_u64(&tmp, v);
  field_in(r, &tmp, ctx);
  bn_free(&tmp);
}

void field_add(bignum* r, const bignum* a, const bignum* b, field_ctx* ctx)
{
  bn_add(r, a, b);
  if (bn_cmp(r, ctx->q) >= 0) {
    bn_sub_abs(r, r, ctx->q);
  }
}

void field_sub(bignum* r, const bignum* a, const bignum* b, field_ctx* ctx)
{
  bn_sub(r, a, b);
  if (r->is_neg) {
    bn_add(r, r, ctx->q);
  }
}

void field_neg(bignum* r, const bignum* a, field_ctx* ctx)
{
  bn_sub(r, ctx->q, a);
  if (bn_cmp(r, ctx->q) >= 0) {
    bn_sub_abs(r, r, ctx->q);
  }
}

void field_mul(bignum* r, const bignum* a, const bignum* b, field_ctx* ctx)
{
  if (ctx->mont) {
    bn_mont_mul(r, a, b, &ctx->mctx);
    return;
  }

  bignum t;
  bn_init(&t);
  bn_mul(&t, a, b);
  bn_mod(r, &t, ctx->q);
  bn_free(&t);
}

/**
 * @brief Modular inverse via the classical extended Euclidean algorithm.
 *
 * Let n = m->size, measured in 64-bit limbs.
 *
 * Computes res = a^(-1) mod m (in [0, m)) using the iterative extended
 * Euclidean algorithm, which is valid for ANY modulus (odd or even) as
 * long as gcd(a, m) == 1. Unlike the binary extended GCD
 * (bn_mod_inverse), it never divides a tracked coefficient by 2, so it
 * is correct for even moduli where 2 is not invertible.
 *
 * Complexity:
 *   Time: O(n^3) - O(n) division steps, each an O(n^2) bignum division
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
static bool mod_inverse_euclid(bignum* res, const bignum* a, const bignum* m)
{
  if (bn_is_zero(a) || bn_is_zero(m) || bn_is_eq_i64(m, 1)) {
    return false;
  }

  // (old_r, r) = (a, m); (old_x, x) = (1, 0)
  // Invariant: old_x * a + old_y * m = old_r (old_y implicit)
  bignum old_r, r, old_x, x, q, tmp;
  bn_init_multi(&old_r, &r, &old_x, &x, &q, &tmp, NULL);

  bn_copy(&old_r, a);
  bn_copy(&r, m);
  bn_set_u64(&old_x, 1);
  bn_set_u64(&x, 0);

  while (!bn_is_zero(&r)) {
    bn_div(&q, &old_r, &r);

    // (old_r, r) = (r, old_r - q*r)
    bn_mul(&tmp, &q, &r);
    bn_sub(&tmp, &old_r, &tmp);
    bn_copy(&old_r, &r);
    bn_copy(&r, &tmp);

    // (old_x, x) = (x, old_x - q*x)
    bn_mul(&tmp, &q, &x);
    bn_sub(&tmp, &old_x, &tmp);
    bn_copy(&old_x, &x);
    bn_copy(&x, &tmp);
  }

  // old_r = gcd(a, m); old_x is the Bezout coefficient of a
  bool ok = bn_is_eq_i64(&old_r, 1);
  if (ok) {
    bn_mod(res, &old_x, m);
    if (res->is_neg) {
      bn_add(res, res, m);
    }
  }

  bn_free_multi(&old_r, &r, &old_x, &x, &q, &tmp, NULL);
  return ok;
}

bool field_inv(bignum* r, const bignum* a, field_ctx* ctx)
{
  if (bn_is_zero(a)) {
    return false;
  }

  if (!ctx->mont) {
    return mod_inverse_euclid(r, a, ctx->q);
  }

  bignum plain, inv;
  bn_init_multi(&plain, &inv, NULL);

  bn_mont_out(&plain, a, &ctx->mctx);
  bool ok = bn_mod_inverse(&inv, &plain, ctx->q);
  if (ok) {
    bn_mont_in(r, &inv, &ctx->mctx);
  }

  bn_free_multi(&plain, &inv, NULL);
  return ok;
}

bool field_div(bignum* r, const bignum* a, const bignum* b, field_ctx* ctx)
{
  bignum inv;
  bn_init(&inv);

  bool ok = field_inv(&inv, b, ctx);
  if (ok) {
    field_mul(r, a, &inv, ctx);
  }

  bn_free(&inv);
  return ok;
}

void field_pow(bignum* r, const bignum* a, const bignum* e, field_ctx* ctx)
{
  if (ctx->mont) {
    bn_mont_exp(r, a, e, &ctx->mctx);
    return;
  }

  bn_mod_exp(r, a, e, ctx->q);
}

bool field_is_zero(const bignum* a) { return bn_is_zero(a); }

bool field_equal(const bignum* a, const bignum* b) { return bn_cmp(a, b) == 0; }

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
static bool poly_is_zero(const bigpoly* p)
{
  for (u64 i = 0; i <= p->deg; i++) {
    if (!bn_is_zero(&p->coeff[i])) {
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
static void poly_assign(bigpoly* dst, const bigpoly* src)
{
  if (!bigpoly_alloc(dst, src->deg)) {
    return;
  }
  for (u64 i = 0; i <= src->deg; i++) {
    bn_copy(&dst->coeff[i], &src->coeff[i]);
  }
  dst->deg = src->deg;
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
static void poly_fadd(bigpoly* r, const bigpoly* a, const bigpoly* b,
                      field_ctx* fctx)
{
  u64 max = MAX(a->deg, b->deg);
  if (!bigpoly_alloc(r, max)) {
    return;
  }

  for (u64 i = 0; i <= max; i++) {
    if (i <= a->deg && i <= b->deg) {
      field_add(&r->coeff[i], &a->coeff[i], &b->coeff[i], fctx);
    } else if (i <= a->deg) {
      bn_copy(&r->coeff[i], &a->coeff[i]);
    } else {
      bn_copy(&r->coeff[i], &b->coeff[i]);
    }
  }

  r->deg = max;
  bigpoly_trim(r);
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
static void poly_fsub(bigpoly* r, const bigpoly* a, const bigpoly* b,
                      field_ctx* fctx)
{
  u64 max = MAX(a->deg, b->deg);
  if (!bigpoly_alloc(r, max)) {
    return;
  }

  for (u64 i = 0; i <= max; i++) {
    if (i <= a->deg && i <= b->deg) {
      field_sub(&r->coeff[i], &a->coeff[i], &b->coeff[i], fctx);
    } else if (i <= a->deg) {
      bn_copy(&r->coeff[i], &a->coeff[i]);
    } else {
      field_neg(&r->coeff[i], &b->coeff[i], fctx);
    }
  }

  r->deg = max;
  bigpoly_trim(r);
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
static void poly_fmul(bigpoly* r, const bigpoly* a, const bigpoly* b,
                      field_ctx* fctx)
{
  u64 r_deg = a->deg + b->deg;
  if (!bigpoly_alloc(r, r_deg)) {
    return;
  }

  for (u64 i = 0; i <= r_deg; i++) {
    bn_set_u64(&r->coeff[i], 0);
  }

  bignum t;
  bn_init(&t);

  for (u64 i = 0; i <= a->deg; i++) {
    for (u64 j = 0; j <= b->deg; j++) {
      field_mul(&t, &a->coeff[i], &b->coeff[j], fctx);
      field_add(&r->coeff[i + j], &r->coeff[i + j], &t, fctx);
    }
  }

  bn_free(&t);

  r->deg = r_deg;
  bigpoly_trim(r);
}

/*===========================================================================
 *  polynomial division and extended GCD
 *=========================================================================*/

bool bigpoly_divmod(bigpoly* q, bigpoly* r, const bigpoly* a, const bigpoly* b,
                    field_ctx* fctx)
{
  if (poly_is_zero(b)) {
    if (r) poly_assign(r, a);
    if (q) {
      bigpoly_alloc(q, 0);
      bn_set_u64(&q->coeff[0], 0);
      q->deg = 0;
    }
    return false;
  }

  bignum inv_lead;
  bn_init(&inv_lead);
  bool ok = field_inv(&inv_lead, &b->coeff[b->deg], fctx);
  if (!ok) {
    if (r) poly_assign(r, a);
    if (q) {
      bigpoly_alloc(q, 0);
      bn_set_u64(&q->coeff[0], 0);
      q->deg = 0;
    }
    bn_free(&inv_lead);
    return false;
  }

  bigpoly rem, quo;
  bigpoly_init(&rem);
  bigpoly_init(&quo);
  poly_assign(&rem, a);

  bigpoly_alloc(&quo, a->deg);
  for (u64 i = 0; i <= a->deg; i++) {
    bn_set_u64(&quo.coeff[i], 0);
  }
  quo.deg = 0;

  bignum factor, t;
  bn_init_multi(&factor, &t, NULL);

  while (rem.deg >= b->deg && !poly_is_zero(&rem)) {
    u64 shift = rem.deg - b->deg;

    // factor = r[deg r] / b[deg b]
    field_mul(&factor, &rem.coeff[rem.deg], &inv_lead, fctx);

    field_add(&quo.coeff[shift], &quo.coeff[shift], &factor, fctx);
    if (shift > quo.deg) {
      quo.deg = shift;
    }

    // r -= factor * x^shift * b
    for (u64 i = 0; i <= b->deg; i++) {
      field_mul(&t, &factor, &b->coeff[i], fctx);
      field_sub(&rem.coeff[i + shift], &rem.coeff[i + shift], &t, fctx);
    }

    bigpoly_trim(&rem);
  }

  bigpoly_trim(&quo);

  if (q) poly_assign(q, &quo);
  if (r) poly_assign(r, &rem);

  bigpoly_free(&quo);
  bigpoly_free(&rem);
  bn_free_multi(&factor, &t, NULL);
  bn_free(&inv_lead);

  return true;
}

bool bigpoly_xgcd(bigpoly* g, bigpoly* x, bigpoly* y, const bigpoly* a,
                  const bigpoly* b, field_ctx* fctx)
{
  bigpoly r0, r1, r2, s0, s1, s2, t0, t1, t2, qq, tmp;
  bigpoly_init(&r0);
  bigpoly_init(&r1);
  bigpoly_init(&r2);
  bigpoly_init(&s0);
  bigpoly_init(&s1);
  bigpoly_init(&s2);
  bigpoly_init(&t0);
  bigpoly_init(&t1);
  bigpoly_init(&t2);
  bigpoly_init(&qq);
  bigpoly_init(&tmp);

  poly_assign(&r0, a);
  poly_assign(&r1, b);

  bigpoly_alloc(&s0, 0);
  bn_set_u64(&s0.coeff[0], 1);
  s0.deg = 0;
  bigpoly_alloc(&t0, 0);
  bn_set_u64(&t0.coeff[0], 0);
  t0.deg = 0;
  bigpoly_alloc(&s1, 0);
  bn_set_u64(&s1.coeff[0], 0);
  s1.deg = 0;
  bigpoly_alloc(&t1, 0);
  bn_set_u64(&t1.coeff[0], 1);
  t1.deg = 0;

  bool ok = true;
  while (!poly_is_zero(&r1)) {
    if (!bigpoly_divmod(&qq, &r2, &r0, &r1, fctx)) {
      ok = false;
      break;
    }

    poly_assign(&r0, &r1);
    poly_assign(&r1, &r2);

    // s2 = s0 - qq * s1
    poly_fmul(&tmp, &qq, &s1, fctx);
    poly_fsub(&s2, &s0, &tmp, fctx);
    poly_assign(&s0, &s1);
    poly_assign(&s1, &s2);

    // t2 = t0 - qq * t1
    poly_fmul(&tmp, &qq, &t1, fctx);
    poly_fsub(&t2, &t0, &tmp, fctx);
    poly_assign(&t0, &t1);
    poly_assign(&t1, &t2);
  }

  if (ok) {
    if (g) poly_assign(g, &r0);
    if (x) poly_assign(x, &s0);
    if (y) poly_assign(y, &t0);
  }

  bigpoly_free(&r0);
  bigpoly_free(&r1);
  bigpoly_free(&r2);
  bigpoly_free(&s0);
  bigpoly_free(&s1);
  bigpoly_free(&s2);
  bigpoly_free(&t0);
  bigpoly_free(&t1);
  bigpoly_free(&t2);
  bigpoly_free(&qq);
  bigpoly_free(&tmp);

  return ok;
}

/*===========================================================================
 *  ring Z_m[x]/(q) operations
 *=========================================================================*/

void poly_ring_init(poly_ring* ctx, const bigpoly* q, field_ctx* fctx)
{
  ctx->q = malloc(sizeof(bigpoly));
  bigpoly_init(ctx->q);
  poly_assign(ctx->q, q);
  ctx->fctx = fctx;
}

void poly_ring_free(poly_ring* ctx)
{
  bigpoly_free(ctx->q);
  free(ctx->q);
  ctx->q = NULL;
  ctx->fctx = NULL;
}

/**
 * @brief Reduce a polynomial modulo the ring modulus (in place into r).
 *
 * Thin wrapper around bigpoly_divmod() keeping the remainder.
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
static void poly_reduce_mod_q(bigpoly* r, const bigpoly* a, poly_ring* ctx)
{
  bigpoly_divmod(NULL, r, a, ctx->q, ctx->fctx);
}

void poly_ring_add(bigpoly* r, const bigpoly* a, const bigpoly* b,
                   poly_ring* ctx)
{
  bigpoly t;
  bigpoly_init(&t);
  poly_fadd(&t, a, b, ctx->fctx);
  poly_reduce_mod_q(r, &t, ctx);
  bigpoly_free(&t);
}

void poly_ring_sub(bigpoly* r, const bigpoly* a, const bigpoly* b,
                   poly_ring* ctx)
{
  bigpoly t;
  bigpoly_init(&t);
  poly_fsub(&t, a, b, ctx->fctx);
  poly_reduce_mod_q(r, &t, ctx);
  bigpoly_free(&t);
}

void poly_ring_neg(bigpoly* r, const bigpoly* a, poly_ring* ctx)
{
  bigpoly zero, t;
  bigpoly_init(&zero);
  bigpoly_init(&t);

  bigpoly_alloc(&zero, 0);
  bn_set_u64(&zero.coeff[0], 0);
  zero.deg = 0;

  poly_fsub(&t, &zero, a, ctx->fctx);
  poly_reduce_mod_q(r, &t, ctx);

  bigpoly_free(&zero);
  bigpoly_free(&t);
}

void poly_ring_mul(bigpoly* r, const bigpoly* a, const bigpoly* b,
                   poly_ring* ctx)
{
  bigpoly t;
  bigpoly_init(&t);
  poly_fmul(&t, a, b, ctx->fctx);
  poly_reduce_mod_q(r, &t, ctx);
  bigpoly_free(&t);
}

bool poly_ring_inv(bigpoly* r, const bigpoly* a, poly_ring* ctx)
{
  bigpoly g, x, y;
  bigpoly_init(&g);
  bigpoly_init(&x);
  bigpoly_init(&y);

  bool ok = bigpoly_xgcd(&g, &x, &y, a, ctx->q, ctx->fctx);

  if (ok && g.deg == 0 && !bn_is_zero(&g.coeff[0])) {
    bignum g_inv;
    bn_init(&g_inv);
    ok = field_inv(&g_inv, &g.coeff[0], ctx->fctx);

    if (ok) {
      bigpoly scaled;
      bigpoly_init(&scaled);
      bigpoly_alloc(&scaled, x.deg);
      for (u64 i = 0; i <= x.deg; i++) {
        field_mul(&scaled.coeff[i], &x.coeff[i], &g_inv, ctx->fctx);
      }
      scaled.deg = x.deg;
      bigpoly_trim(&scaled);

      poly_reduce_mod_q(r, &scaled, ctx);
      bigpoly_free(&scaled);
    }

    bn_free(&g_inv);
  } else {
    ok = false;
  }

  bigpoly_free(&g);
  bigpoly_free(&x);
  bigpoly_free(&y);

  return ok;
}

bool poly_ring_div(bigpoly* r, const bigpoly* a, const bigpoly* b,
                   poly_ring* ctx)
{
  bigpoly b_inv;
  bigpoly_init(&b_inv);

  bool ok = poly_ring_inv(&b_inv, b, ctx);
  if (ok) {
    poly_ring_mul(r, a, &b_inv, ctx);
  }

  bigpoly_free(&b_inv);
  return ok;
}

void poly_ring_pow(bigpoly* r, const bigpoly* a, const bignum* e,
                   poly_ring* ctx)
{
  bigpoly result, base, tmp;
  bigpoly_init(&result);
  bigpoly_init(&base);
  bigpoly_init(&tmp);

  // result = 1
  bigpoly_alloc(&result, 0);
  bn_set_u64(&result.coeff[0], 1);
  result.deg = 0;

  poly_assign(&base, a);

  bignum exp;
  bn_init(&exp);
  bn_copy(&exp, e);

  while (!bn_is_zero(&exp)) {
    if (!bn_is_even(&exp)) {
      poly_ring_mul(&tmp, &result, &base, ctx);
      poly_assign(&result, &tmp);
    }
    poly_ring_mul(&tmp, &base, &base, ctx);
    poly_assign(&base, &tmp);
    bn_rshift1(&exp);
  }

  poly_assign(r, &result);

  bn_free(&exp);
  bigpoly_free(&result);
  bigpoly_free(&base);
  bigpoly_free(&tmp);
}

bool poly_ring_is_zero(const bigpoly* a) { return poly_is_zero(a); }

bool poly_ring_equal(const bigpoly* a, const bigpoly* b)
{
  u64 max = MAX(a->deg, b->deg);
  for (u64 i = 0; i <= max; i++) {
    bool ai_zero = (i > a->deg) || bn_is_zero(&a->coeff[i]);
    bool bi_zero = (i > b->deg) || bn_is_zero(&b->coeff[i]);
    if (ai_zero != bi_zero) {
      return false;
    }
    if (!ai_zero && bn_cmp(&a->coeff[i], &b->coeff[i]) != 0) {
      return false;
    }
  }
  return true;
}
