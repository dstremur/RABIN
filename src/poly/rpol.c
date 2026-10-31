/*
 * rpol.c
 *
 * Polynomials over bignums.
 *
 * This file implements polynomials with rz_t coefficients:
 * initialization, setting from coefficient arrays, copying,
 * comparison, growing the coefficient array, trimming, printing,
 * addition, subtraction, multiplication (schoolbook, rz_t NTT and
 * u64 NTT), pseudo-division and exact division, content, primitive
 * part, gcd, evaluation and derivation, fixed-width decomposition
 * and recomposition, and carry propagation.
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

#include "../../include/rpol.h"

#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../../include/rntt.h"
#include "../../include/rzlogic.h"
#include "../../include/rzmath.h"
#include "../../include/u64.h"
#include "../core/strparse_internal.h"

rabin_err_t rpol_init(rpol_t* p)
{
  if (p == NULL) return RABIN_ERR_NULL_PTR;

  // idempotent init: release any previous allocation
  if (p->coeff != NULL) {
    rpol_clear(p);
  }

  p->coeff = NULL;
  p->deg = 0;
  p->size = 0;
  return RABIN_SUCCESS;
}

rabin_err_t rpol_set(rpol_t* p, const rz_t* coeff, u64 deg)
{
  if (p == NULL || coeff == NULL) return RABIN_ERR_NULL_PTR;

  rabin_err_t err = rpol_alloc(p, deg);
  if (err != RABIN_SUCCESS) return err;

  for (u64 i = 0; i <= deg; i++) {
    if ((err = rz_copy(&p->coeff[i], &coeff[i])) != RABIN_SUCCESS) return err;
  }

  p->deg = deg;
  return RABIN_SUCCESS;
}

rabin_err_t rpol_set_i64(rpol_t* p, const i64* coeff, u64 deg)
{
  if (p == NULL || coeff == NULL) return RABIN_ERR_NULL_PTR;

  rabin_err_t err = rpol_alloc(p, deg);
  if (err != RABIN_SUCCESS) return err;

  for (u64 i = 0; i <= deg; i++) {
    if ((err = rz_set_i64(&p->coeff[i], coeff[i])) != RABIN_SUCCESS) return err;
  }

  p->deg = deg;
  return rpol_trim(p);
}

rabin_err_t rpol_copy(rpol_t* p, const rpol_t* q)
{
  if (p == NULL || q == NULL) return RABIN_ERR_NULL_PTR;

  if (p->size < q->deg + 1) return RABIN_ERR_OVERFLOW;

  for (u64 i = 0; i <= q->deg; i++) {
    rabin_err_t err = rz_copy(&p->coeff[i], &q->coeff[i]);
    if (err != RABIN_SUCCESS) return err;
  }
  return RABIN_SUCCESS;
}

bool rpol_equal(const rpol_t* a, const rpol_t* b)
{
  if (a == NULL || b == NULL) return false;

  if (a->deg != b->deg) return false;
  for (u64 i = 0; i <= a->deg; i++) {
    if (rz_cmp(&a->coeff[i], &b->coeff[i]) != 0) {
      return false;
    }
  }
  return true;
}

rabin_err_t rpol_clear(rpol_t* p)
{
  if (p == NULL) return RABIN_ERR_NULL_PTR;
  if (!p->coeff) return RABIN_SUCCESS;

  for (u64 i = 0; i < p->size; i++) {
    rz_clear(&p->coeff[i]);
  }

  free(p->coeff);
  p->coeff = NULL;
  p->deg = 0;
  p->size = 0;
  return RABIN_SUCCESS;
}

rabin_err_t rpol_alloc(rpol_t* p, u64 deg)
{
  if (p == NULL) return RABIN_ERR_NULL_PTR;

  // deg + 1 must not overflow
  if (deg == UINT64_MAX) return RABIN_ERR_OVERFLOW;

  u64 needed = deg + 1;
  if (p->size >= needed) return RABIN_SUCCESS;

  u64 new_cap = (p->size == 0) ? needed : (p->size * 2);
  if (p->size != 0 && new_cap < p->size) {  // doubling overflowed
    new_cap = needed;
    if (new_cap < p->size) return RABIN_ERR_OVERFLOW;
  }
  if (new_cap < needed) {
    new_cap = needed;
  }
  if (new_cap > SIZE_MAX / sizeof(rz_t)) return RABIN_ERR_OVERFLOW;

  rz_t* new = realloc(p->coeff, new_cap * sizeof(rz_t));
  if (!new) return RABIN_ERR_OUT_OF_MEMORY;

  p->coeff = new;

  for (u64 i = p->size; i < new_cap; i++) {
    rz_init(&p->coeff[i]);
  }

  p->size = new_cap;
  return RABIN_SUCCESS;
}

rabin_err_t rpol_trim(rpol_t* p)
{
  if (p == NULL) return RABIN_ERR_NULL_PTR;

  while (p->deg > 0 && rz_is_zero(&p->coeff[p->deg])) {
    p->deg--;
  }
  return RABIN_SUCCESS;
}

rabin_err_t rpol_print(const rpol_t* p)
{
  if (p == NULL || p->coeff == NULL || p->size == 0) {
    printf("0\n");
    return RABIN_SUCCESS;
  }

  rz_t mag;
  rz_init(&mag);
  rabin_err_t err;

  bool first = true;

  for (i64 i = (i64)p->deg; i >= 0; i--) {
    const rz_t* c = &p->coeff[i];
    if (rz_is_zero(c)) continue;

    const bool neg = c->is_neg;

    if (first) {
      if (neg) putchar('-');
    } else {
      putchar(' ');
      putchar(neg ? '-' : '+');
      putchar(' ');
    }

    TRY(rz_copy(&mag, c));
    if (neg) {
      TRY(rz_neg(&mag, &mag));
    }

    /* hide the coefficient 1 on every term except the constant one */
    if (i == 0 || !rz_is_eq_i64(&mag, 1)) rz_print(&mag);

    if (i >= 2)
      printf("x^%" PRId64, i);
    else if (i == 1)
      putchar('x');

    first = false;
  }

  if (first) putchar('0');

  printf("\n");
  err = RABIN_SUCCESS;
cleanup:
  rz_clear(&mag);
  return err;
}

rabin_err_t rpol_add(rpol_t* r, const rpol_t* p, const rpol_t* q)
{
  if (r == NULL || p == NULL || q == NULL) return RABIN_ERR_NULL_PTR;

  u64 min = RZ_MIN(p->deg, q->deg);
  u64 max = RZ_MAX(p->deg, q->deg);

  rabin_err_t err = rpol_alloc(r, max);
  if (err != RABIN_SUCCESS) return err;

  // add upto the minimum degree
  for (u64 i = 0; i <= min; i++) {
    if ((err = rz_add(&r->coeff[i], &p->coeff[i], &q->coeff[i])) !=
        RABIN_SUCCESS)
      return err;
  }

  const rpol_t* longer = (p->deg > q->deg) ? p : q;
  for (u64 i = min + 1; i <= max; i++) {
    if ((err = rz_copy(&r->coeff[i], &longer->coeff[i])) != RABIN_SUCCESS)
      return err;
  }

  r->deg = max;
  return rpol_trim(r);
}

rabin_err_t rpol_sub(rpol_t* r, const rpol_t* p, const rpol_t* q)
{
  if (r == NULL || p == NULL || q == NULL) return RABIN_ERR_NULL_PTR;

  u64 min = RZ_MIN(p->deg, q->deg);
  u64 max = RZ_MAX(p->deg, q->deg);

  rabin_err_t err = rpol_alloc(r, max);
  if (err != RABIN_SUCCESS) return err;

  // subtract upto the minimum degree
  for (u64 i = 0; i <= min; i++) {
    if ((err = rz_sub(&r->coeff[i], &p->coeff[i], &q->coeff[i])) !=
        RABIN_SUCCESS)
      return err;
  }

  // coefficients beyond min(p,q) belong to the longer operand only: p - q
  // keeps p's own sign there, and negates q's
  if (p->deg > q->deg) {
    for (u64 i = min + 1; i <= max; i++) {
      if ((err = rz_copy(&r->coeff[i], &p->coeff[i])) != RABIN_SUCCESS)
        return err;
    }
  } else if (q->deg > p->deg) {
    for (u64 i = min + 1; i <= max; i++) {
      if ((err = rz_neg(&r->coeff[i], &q->coeff[i])) != RABIN_SUCCESS)
        return err;
    }
  }

  r->deg = max;
  return rpol_trim(r);
}

rabin_err_t rpol_mul_ntt(rpol_t* r, const rpol_t* p, const rpol_t* q)
{
  if (r == NULL || p == NULL || q == NULL) return RABIN_ERR_NULL_PTR;
  if (p->deg > UINT64_MAX - q->deg - 1) return RABIN_ERR_OVERFLOW;

  rntt_ctx_t ctx;

  u64 required_len = p->deg + q->deg + 1;
  u64 ntt_size = 1;
  u64 k = 0;
  // pad to next power of 2
  while (ntt_size < required_len) {
    ntt_size <<= 1;
    k++;
  }

  rabin_err_t err = rpol_alloc(r, required_len);
  if (err != RABIN_SUCCESS) return err;

  if ((err = rntt_ctx_init_golden(&ctx, k)) != RABIN_SUCCESS) return err;

  rpol_t p_hat = {0}, q_hat = {0}, r_hat = {0};
  rpol_init(&p_hat);
  rpol_init(&q_hat);
  rpol_init(&r_hat);

  TRY(rntt_cyclic_forward(&p_hat, p, &ctx));
  TRY(rntt_cyclic_forward(&q_hat, q, &ctx));
  TRY(rpol_alloc(&r_hat, ntt_size));

  for (u64 i = 0; i < ntt_size; i++) {
    // MontMul(pR, qR) = (pR * qR * R^-1) mod q = (p*q)R mod q
    TRY(rz_mont_mul(&r_hat.coeff[i], &p_hat.coeff[i], &q_hat.coeff[i],
                    &ctx.mctx));
  }
  r_hat.deg = ntt_size - 1;

  rpol_t r_ntt = {0};
  rpol_init(&r_ntt);
  TRY(rntt_cyclic_inverse_mont_in(&r_ntt, &r_hat, &ctx));

  rpol_clear(r);
  rpol_init(r);
  TRY(rpol_alloc(r, required_len));

  for (u64 i = 0; i < required_len; i++) {
    TRY(rz_copy(&r->coeff[i], &r_ntt.coeff[i]));
  }

  r->deg = required_len - 1;

  TRY(rpol_trim(r));

  err = RABIN_SUCCESS;
cleanup:
  rpol_clear(&p_hat);
  rpol_clear(&q_hat);
  rpol_clear(&r_hat);
  rpol_clear(&r_ntt);
  rntt_ctx_clear(&ctx);
  return err;
}

rabin_err_t rpol_mul_ntt_u64(rpol_t* r, const rpol_t* p, const rpol_t* q)
{
  if (r == NULL || p == NULL || q == NULL) return RABIN_ERR_NULL_PTR;
  if (p->deg > UINT64_MAX - q->deg - 1) return RABIN_ERR_OVERFLOW;

  u64 required_len = p->deg + q->deg + 1;
  u64 ntt_size = 1;
  u64 k = 0;

  // Pad to next power of 2
  while (ntt_size < required_len) {
    ntt_size <<= 1;
    k++;
  }

  u64_ntt_ctx_t ctx;
  rabin_err_t err = u64_ntt_ctx_init_golden(&ctx, k);
  if (err != RABIN_SUCCESS) return err;

  // 1. One arena block for all six flat u64 arrays (zeroed, like calloc)
  if (ntt_size > SIZE_MAX / (6 * sizeof(u64))) return RABIN_ERR_OVERFLOW;
  u64* buf = rz_scratch_get(6 * ntt_size);
  if (buf == NULL) {
    return RABIN_ERR_OUT_OF_MEMORY;
  }
  memset(buf, 0, 6 * ntt_size * sizeof(u64));
  u64* p_arr = buf;
  u64* q_arr = buf + ntt_size;
  u64* p_hat = buf + 2 * ntt_size;
  u64* q_hat = buf + 3 * ntt_size;
  u64* r_hat = buf + 4 * ntt_size;
  u64* r_arr = buf + 5 * ntt_size;

  // 2. Extract u64 values from the rz_t polynomials
  for (u64 i = 0; i <= p->deg; i++) {
    p_arr[i] = p->coeff[i].limbs[0];
  }
  for (u64 i = 0; i <= q->deg; i++) {
    q_arr[i] = q->coeff[i].limbs[0];
  }

  // 3. Perform Forward NTTs
  u64_ntt_cyclic_forward(p_hat, p_arr, &ctx);
  u64_ntt_cyclic_forward(q_hat, q_arr, &ctx);

  // 4. Pointwise Multiplication
  for (u64 i = 0; i < ntt_size; i++) {
    r_hat[i] = u64_mont_mul(p_hat[i], q_hat[i], &ctx.mctx);
  }

  // 5. Perform Inverse NTT
  u64_ntt_cyclic_inverse_montgomery_in(r_arr, r_hat, &ctx);

  err = rpol_clear(r);
  if (err != RABIN_SUCCESS) goto cleanup;
  rpol_init(r);
  TRY(rpol_alloc(r, required_len));

  for (u64 i = 0; i < required_len; i++) {
    TRY(rz_set_u64(&r->coeff[i], r_arr[i]));
  }

  r->deg = required_len - 1;
  TRY(rpol_trim(r));

  err = RABIN_SUCCESS;
cleanup:
  // 7. Cleanup
  rz_scratch_release();
  u64_ntt_ctx_clear(&ctx);

  return err;
}

rabin_err_t rpol_mul_school(rpol_t* r, const rpol_t* p, const rpol_t* q)
{
  if (r == NULL || p == NULL || q == NULL) return RABIN_ERR_NULL_PTR;
  if (p->deg > UINT64_MAX - q->deg) return RABIN_ERR_OVERFLOW;

  u64 r_deg = p->deg + q->deg;

  rpol_t temp = {0};
  rpol_init(&temp);
  rabin_err_t err = rpol_alloc(&temp, r_deg);
  if (err != RABIN_SUCCESS) return err;

  temp.deg = r_deg;

  rz_t t;
  rz_init(&t);

  for (u64 i = 0; i <= p->deg; i++) {
    for (u64 j = 0; j <= q->deg; j++) {
      TRY(rz_mul(&t, &p->coeff[i], &q->coeff[j]));

      TRY(rz_add(&temp.coeff[i + j], &temp.coeff[i + j], &t));
    }
  }

  TRY(rpol_trim(&temp));

  TRY(rpol_alloc(r, temp.deg));
  for (u64 i = 0; i <= temp.deg; i++) {
    TRY(rz_copy(&r->coeff[i], &temp.coeff[i]));
  }
  r->deg = temp.deg;

  err = RABIN_SUCCESS;
cleanup:
  rz_clear(&t);
  rpol_clear(&temp);
  return err;
}

rabin_err_t rpol_mul(rpol_t* r, const rpol_t* p, const rpol_t* q)
{
  return rpol_mul_school(r, p, q);
}

rabin_err_t rz_decompose(rpol_t* r, const rz_t* n, u64 width)
{
  if (r == NULL || n == NULL) return RABIN_ERR_NULL_PTR;
  if (width == 0) return RABIN_ERR_INVALID_ARG;

  rz_t tmp, mask, digit;
  rz_init_multi(&tmp, &mask, &digit, NULL);

  rabin_err_t err = rz_copy(&tmp, n);
  if (err != RABIN_SUCCESS) goto cleanup;

  // mask = 0xFFFF..
  TRY(rz_lshift(&mask, &RZ_ONE, width));
  TRY(rz_sub(&mask, &mask, &RZ_ONE));

  u64 i = 0;
  while (!rz_is_zero(&tmp)) {
    TRY(rz_and(&digit, &tmp, &mask));

    TRY(rpol_alloc(r, i + 1));
    TRY(rz_copy(&r->coeff[i], &digit));

    TRY(rz_rshift(&tmp, &tmp, width));

    i++;
  }

  r->deg = (i > 0) ? i - 1 : 0;

  err = RABIN_SUCCESS;
cleanup:
  rz_clear_multi(&tmp, &mask, &digit, NULL);
  return err;
}

rabin_err_t rpol_carry_propagation(rpol_t* r, u64 bit_width)
{
  if (r == NULL) return RABIN_ERR_NULL_PTR;
  if (bit_width == 0) return RABIN_ERR_INVALID_ARG;

  rz_t carry, base, mask, total;
  rz_init_multi(&carry, &base, &mask, &total, NULL);

  rabin_err_t err = rz_lshift(&base, &RZ_ONE, bit_width);
  if (err != RABIN_SUCCESS) goto cleanup;
  TRY(rz_sub(&mask, &base, &RZ_ONE));
  TRY(rz_set_u64(&carry, 0));

  u64 i = 0;
  // Iterate through all coefficients plus any remaining carries
  while (i <= r->deg || !rz_is_zero(&carry)) {
    if (i > r->deg) {
      // Expand poly if carry exceeds current deg
      TRY(rpol_alloc(r, i + 1));
      r->deg = i;
    }

    // total = coeff[i] + carry
    TRY(rz_add(&total, &r->coeff[i], &carry));

    // carry = total >> bit_width
    TRY(rz_rshift(&carry, &total, bit_width));

    // coeff[i] = total & mask
    TRY(rz_and(&r->coeff[i], &total, &mask));

    i++;
  }

  err = RABIN_SUCCESS;
cleanup:
  rz_clear_multi(&carry, &base, &mask, &total, NULL);
  return err;
}

rabin_err_t rz_recompose(rz_t* n, const rpol_t* p, u64 bit_width)
{
  if (n == NULL || p == NULL) return RABIN_ERR_NULL_PTR;

  rabin_err_t err = rz_set_u64(n, 0);
  if (err != RABIN_SUCCESS) return err;

  rz_t term;
  rz_init(&term);

  for (u64 i = 0; i <= p->deg; i++) {
    // term = coeff[i] << (i * bit_width)
    TRY(rz_lshift(&term, &p->coeff[i], i * bit_width));
    // n += term
    TRY(rz_add(n, n, &term));
  }

  err = RABIN_SUCCESS;
cleanup:
  rz_clear(&term);
  return err;
}

/* ------------------------------------------------------------------------ */
/* Division, content and gcd over Z[x]                                      */
/* ------------------------------------------------------------------------ */

/*
 * True if p represents the zero polynomial. p may be freshly initialized
 * (size == 0) or hold untrimmed leading zeros.
 */
static bool rpol_is_zero_poly(const rpol_t* p)
{
  if (p->size == 0) return true;

  for (u64 i = 0; i <= p->deg; i++) {
    if (!rz_is_zero(&p->coeff[i])) return false;
  }
  return true;
}

/* Set p to the canonical zero polynomial (deg 0, single zero coefficient). */
static rabin_err_t rpol_set_zero_poly(rpol_t* p)
{
  rabin_err_t err = rpol_alloc(p, 0);
  if (err != RABIN_SUCCESS) return err;

  p->deg = 0;
  return rz_set_u64(&p->coeff[0], 0);
}

/*
 * Deep copy src into dst (allocating dst as needed) and set dst->deg.
 * An empty src (size == 0) becomes the canonical zero polynomial.
 */
static rabin_err_t rpol_copy_full(rpol_t* dst, const rpol_t* src)
{
  if (src->size == 0) return rpol_set_zero_poly(dst);

  rabin_err_t err = rpol_alloc(dst, src->deg);
  if (err != RABIN_SUCCESS) return err;

  err = rpol_copy(dst, src);
  if (err != RABIN_SUCCESS) return err;

  dst->deg = src->deg;
  return RABIN_SUCCESS;
}

/* Flip the sign of every coefficient, in place. */
static rabin_err_t rpol_negate(rpol_t* p)
{
  for (u64 i = 0; i <= p->deg; i++) {
    rabin_err_t err = rz_neg(&p->coeff[i], &p->coeff[i]);
    if (err != RABIN_SUCCESS) return err;
  }
  return RABIN_SUCCESS;
}

/*
 * Divide p by its content and normalize the leading coefficient to be
 * positive, in place. p must be trimmed, nonzero and have size >= deg + 1;
 * scratch is a zero-initialized bignum used as working space.
 */
static rabin_err_t rpol_make_primitive(rpol_t* p, rz_t* scratch)
{
  rz_set_u64(scratch, 0);

  for (u64 i = 0; i <= p->deg; i++) {
    rabin_err_t err = rz_gcd(scratch, scratch, &p->coeff[i]);
    if (err != RABIN_SUCCESS) return err;
    if (rz_is_eq_i64(scratch, 1)) break;
  }

  if (!rz_is_eq_i64(scratch, 1)) {
    for (u64 i = 0; i <= p->deg; i++) {
      rabin_err_t err = rz_div(&p->coeff[i], &p->coeff[i], scratch);
      if (err != RABIN_SUCCESS) return err;
    }
  }

  if (p->coeff[p->deg].is_neg) return rpol_negate(p);
  return RABIN_SUCCESS;
}

rabin_err_t rpol_pseudo_div(rpol_t* q, rpol_t* r, const rpol_t* a,
                            const rpol_t* b)
{
  if (q == NULL || r == NULL || a == NULL || b == NULL)
    return RABIN_ERR_NULL_PTR;

  rabin_err_t err;
  rpol_t bb = {0}, R = {0}, Q = {0};
  rz_t lc_b, lc_pow, expo, coef, prod;
  rpol_init(&bb);
  rpol_init(&R);
  rpol_init(&Q);
  rz_init_multi(&lc_b, &lc_pow, &expo, &coef, &prod, NULL);

  /* Local, trimmed copy of the divisor */
  err = rpol_copy_full(&bb, b);
  if (err != RABIN_SUCCESS) goto cleanup;
  err = rpol_trim(&bb);
  if (err != RABIN_SUCCESS) goto cleanup;

  if (rpol_is_zero_poly(&bb)) {
    err = RABIN_ERR_DIV_BY_ZERO;
    goto cleanup;
  }

  /* Local, trimmed copy of the dividend */
  err = rpol_copy_full(&R, a);
  if (err != RABIN_SUCCESS) goto cleanup;
  err = rpol_trim(&R);
  if (err != RABIN_SUCCESS) goto cleanup;

  /* a = 0 or deg a < deg b: q = 0, r = a, identity a = 0*b + a */
  if (rpol_is_zero_poly(&R) || R.deg < bb.deg) {
    rpol_clear(q);
    rpol_init(q);
    err = rpol_set_zero_poly(q);
    if (err != RABIN_SUCCESS) goto cleanup;

    rpol_clear(r);
    rpol_init(r);
    err = rpol_copy_full(r, &R);
    goto cleanup;
  }

  /* R := lc(b)^(deg a - deg b + 1) * a, so all later divisions are exact */
  u64 delta = R.deg - bb.deg;
  err = rz_copy(&lc_b, &bb.coeff[bb.deg]);
  if (err != RABIN_SUCCESS) goto cleanup;
  err = rz_set_u64(&expo, delta + 1);
  if (err != RABIN_SUCCESS) goto cleanup;
  err = rz_pow(&lc_pow, &lc_b, &expo);
  if (err != RABIN_SUCCESS) goto cleanup;

  for (u64 i = 0; i <= R.deg; i++) {
    err = rz_mul(&R.coeff[i], &R.coeff[i], &lc_pow);
    if (err != RABIN_SUCCESS) goto cleanup;
  }

  /* Q holds the quotient coefficients q_delta ... q_0 (zero-initialized) */
  err = rpol_alloc(&Q, delta);
  if (err != RABIN_SUCCESS) goto cleanup;
  Q.deg = delta;

  while (!rpol_is_zero_poly(&R) && R.deg >= bb.deg) {
    /* q_k := lc(R) / lc(b); exact by the pseudo-division invariant */
    u64 k = R.deg - bb.deg;
    err = rz_div(&coef, &R.coeff[R.deg], &lc_b);
    if (err != RABIN_SUCCESS) goto cleanup;
    err = rz_copy(&Q.coeff[k], &coef);
    if (err != RABIN_SUCCESS) goto cleanup;

    /* R := R - q_k * b * x^k (the leading term cancels exactly) */
    for (u64 j = 0; j <= bb.deg; j++) {
      if (rz_is_zero(&bb.coeff[j])) continue;

      err = rz_mul(&prod, &coef, &bb.coeff[j]);
      if (err != RABIN_SUCCESS) goto cleanup;
      err = rz_sub(&R.coeff[k + j], &R.coeff[k + j], &prod);
      if (err != RABIN_SUCCESS) goto cleanup;
    }

    err = rpol_trim(&R);
    if (err != RABIN_SUCCESS) goto cleanup;
  }

  err = rpol_trim(&Q);
  if (err != RABIN_SUCCESS) goto cleanup;

  rpol_clear(q);
  rpol_init(q);
  err = rpol_copy_full(q, &Q);
  if (err != RABIN_SUCCESS) goto cleanup;

  rpol_clear(r);
  rpol_init(r);
  err = rpol_copy_full(r, &R);

cleanup:
  rpol_clear(&bb);
  rpol_clear(&R);
  rpol_clear(&Q);
  rz_clear_multi(&lc_b, &lc_pow, &expo, &coef, &prod, NULL);
  return err;
}

rabin_err_t rpol_div_exact(rpol_t* q, const rpol_t* a, const rpol_t* b)
{
  if (q == NULL || a == NULL || b == NULL) return RABIN_ERR_NULL_PTR;

  rabin_err_t err;
  rpol_t bb = {0}, R = {0}, Q = {0};
  rz_t lc_b, coef, prod, rem;
  rpol_init(&bb);
  rpol_init(&R);
  rpol_init(&Q);
  rz_init_multi(&lc_b, &coef, &prod, &rem, NULL);

  /* Local, trimmed copy of the divisor */
  err = rpol_copy_full(&bb, b);
  if (err != RABIN_SUCCESS) goto cleanup;
  err = rpol_trim(&bb);
  if (err != RABIN_SUCCESS) goto cleanup;

  if (rpol_is_zero_poly(&bb)) {
    err = RABIN_ERR_DIV_BY_ZERO;
    goto cleanup;
  }

  /* Local, trimmed copy of the dividend */
  err = rpol_copy_full(&R, a);
  if (err != RABIN_SUCCESS) goto cleanup;
  err = rpol_trim(&R);
  if (err != RABIN_SUCCESS) goto cleanup;

  bool exact = true;

  if (rpol_is_zero_poly(&R)) {
    /* 0 = 0 * b is exact */
  } else if (R.deg < bb.deg) {
    exact = false;
  } else {
    u64 delta = R.deg - bb.deg;

    err = rz_copy(&lc_b, &bb.coeff[bb.deg]);
    if (err != RABIN_SUCCESS) goto cleanup;

    err = rpol_alloc(&Q, delta);
    if (err != RABIN_SUCCESS) goto cleanup;
    Q.deg = delta;

    while (!rpol_is_zero_poly(&R) && R.deg >= bb.deg) {
      /* q_k := lc(R) / lc(b), but only if lc(b) divides lc(R) */
      u64 k = R.deg - bb.deg;
      err = rz_mod(&rem, &R.coeff[R.deg], &lc_b);
      if (err != RABIN_SUCCESS) goto cleanup;
      if (!rz_is_zero(&rem)) {
        exact = false;
        break;
      }

      err = rz_div(&coef, &R.coeff[R.deg], &lc_b);
      if (err != RABIN_SUCCESS) goto cleanup;
      err = rz_copy(&Q.coeff[k], &coef);
      if (err != RABIN_SUCCESS) goto cleanup;

      /* R := R - q_k * b * x^k */
      for (u64 j = 0; j <= bb.deg; j++) {
        if (rz_is_zero(&bb.coeff[j])) continue;

        err = rz_mul(&prod, &coef, &bb.coeff[j]);
        if (err != RABIN_SUCCESS) goto cleanup;
        err = rz_sub(&R.coeff[k + j], &R.coeff[k + j], &prod);
        if (err != RABIN_SUCCESS) goto cleanup;
      }

      err = rpol_trim(&R);
      if (err != RABIN_SUCCESS) goto cleanup;
    }

    /* a = q*b holds over Z[x] only if the remainder vanished */
    if (exact && !rpol_is_zero_poly(&R)) exact = false;
  }

  if (!exact) {
    err = RABIN_ERR_INVALID_ARG;
    goto cleanup;
  }

  err = rpol_trim(&Q);
  if (err != RABIN_SUCCESS) goto cleanup;

  rpol_clear(q);
  rpol_init(q);
  err = rpol_copy_full(q, &Q);

cleanup:
  rpol_clear(&bb);
  rpol_clear(&R);
  rpol_clear(&Q);
  rz_clear_multi(&lc_b, &coef, &prod, &rem, NULL);
  return err;
}

rabin_err_t rpol_content(rz_t* g, const rpol_t* p)
{
  if (g == NULL || p == NULL) return RABIN_ERR_NULL_PTR;

  rabin_err_t err = rz_set_u64(g, 0);
  if (err != RABIN_SUCCESS) return err;

  for (u64 i = 0; i <= p->deg; i++) {
    err = rz_gcd(g, g, &p->coeff[i]);
    if (err != RABIN_SUCCESS) return err;
    if (rz_is_eq_i64(g, 1)) break;
  }
  return RABIN_SUCCESS;
}

rabin_err_t rpol_primitive_part(rpol_t* pp, const rpol_t* p)
{
  if (pp == NULL || p == NULL) return RABIN_ERR_NULL_PTR;

  rabin_err_t err;
  rpol_t T = {0};
  rz_t c;
  rpol_init(&T);
  rz_init(&c);

  err = rpol_copy_full(&T, p);
  if (err != RABIN_SUCCESS) goto cleanup;
  err = rpol_trim(&T);
  if (err != RABIN_SUCCESS) goto cleanup;

  if (rpol_is_zero_poly(&T)) {
    rpol_clear(pp);
    rpol_init(pp);
    err = rpol_set_zero_poly(pp);
    goto cleanup;
  }

  err = rpol_make_primitive(&T, &c);
  if (err != RABIN_SUCCESS) goto cleanup;

  rpol_clear(pp);
  rpol_init(pp);
  err = rpol_copy_full(pp, &T);

cleanup:
  rpol_clear(&T);
  rz_clear(&c);
  return err;
}

rabin_err_t rpol_gcd(rpol_t* g, const rpol_t* a, const rpol_t* b)
{
  if (g == NULL || a == NULL || b == NULL) return RABIN_ERR_NULL_PTR;

  rabin_err_t err;
  rpol_t aa = {0}, bb = {0}, Q = {0}, R = {0};
  rz_t ca, cb, ct, tmp;
  rpol_init(&aa);
  rpol_init(&bb);
  rpol_init(&Q);
  rpol_init(&R);
  rz_init_multi(&ca, &cb, &ct, &tmp, NULL);

  err = rpol_copy_full(&aa, a);
  if (err != RABIN_SUCCESS) goto cleanup;
  err = rpol_trim(&aa);
  if (err != RABIN_SUCCESS) goto cleanup;

  err = rpol_copy_full(&bb, b);
  if (err != RABIN_SUCCESS) goto cleanup;
  err = rpol_trim(&bb);
  if (err != RABIN_SUCCESS) goto cleanup;

  bool aa_zero = rpol_is_zero_poly(&aa);
  bool bb_zero = rpol_is_zero_poly(&bb);

  /* gcd(0, 0) = 0; gcd(0, b) = b with positive leading coefficient */
  if (aa_zero && bb_zero) {
    rpol_clear(g);
    rpol_init(g);
    err = rpol_set_zero_poly(g);
    goto cleanup;
  }

  if (aa_zero || bb_zero) {
    rpol_t* nz = aa_zero ? &bb : &aa;

    if (nz->coeff[nz->deg].is_neg) {
      err = rpol_negate(nz);
      if (err != RABIN_SUCCESS) goto cleanup;
    }

    rpol_clear(g);
    rpol_init(g);
    err = rpol_copy_full(g, nz);
    goto cleanup;
  }

  /* gcd(a, b) = gcd(cont a, cont b) * gcd(pp a, pp b), Gauss's lemma */
  err = rpol_content(&ca, &aa);
  if (err != RABIN_SUCCESS) goto cleanup;
  err = rpol_content(&cb, &bb);
  if (err != RABIN_SUCCESS) goto cleanup;
  err = rz_gcd(&ct, &ca, &cb);
  if (err != RABIN_SUCCESS) goto cleanup;

  err = rpol_make_primitive(&aa, &tmp);
  if (err != RABIN_SUCCESS) goto cleanup;
  err = rpol_make_primitive(&bb, &tmp);
  if (err != RABIN_SUCCESS) goto cleanup;

  if (aa.deg < bb.deg) {
    rpol_t swap = aa;
    aa = bb;
    bb = swap;
  }

  /* Reduced pseudo-remainder sequence: strip the content of every
   * remainder so coefficients cannot blow up unchecked. aa always holds
   * the last nonzero remainder, so the rotation must happen before the
   * zero-remainder exit. */
  while (!rpol_is_zero_poly(&bb)) {
    err = rpol_pseudo_div(&Q, &R, &aa, &bb);
    if (err != RABIN_SUCCESS) goto cleanup;

    err = rpol_copy_full(&aa, &bb);
    if (err != RABIN_SUCCESS) goto cleanup;

    if (rpol_is_zero_poly(&R)) break;

    err = rpol_make_primitive(&R, &tmp);
    if (err != RABIN_SUCCESS) goto cleanup;

    err = rpol_copy_full(&bb, &R);
    if (err != RABIN_SUCCESS) goto cleanup;
  }

  /* aa is the primitive gcd with positive leading coefficient */
  if (!rz_is_eq_i64(&ct, 1)) {
    for (u64 i = 0; i <= aa.deg; i++) {
      err = rz_mul(&aa.coeff[i], &aa.coeff[i], &ct);
      if (err != RABIN_SUCCESS) goto cleanup;
    }
  }

  rpol_clear(g);
  rpol_init(g);
  err = rpol_copy_full(g, &aa);

cleanup:
  rpol_clear(&aa);
  rpol_clear(&bb);
  rpol_clear(&Q);
  rpol_clear(&R);
  rz_clear_multi(&ca, &cb, &ct, &tmp, NULL);
  return err;
}

rabin_err_t rpol_eval(rz_t* r, const rpol_t* p, const rz_t* x)
{
  if (r == NULL || p == NULL || x == NULL) return RABIN_ERR_NULL_PTR;

  rabin_err_t err;
  rpol_t T = {0};
  rz_t acc;
  rpol_init(&T);
  rz_init(&acc);

  err = rpol_copy_full(&T, p);
  if (err != RABIN_SUCCESS) goto cleanup;
  err = rpol_trim(&T);
  if (err != RABIN_SUCCESS) goto cleanup;

  if (rpol_is_zero_poly(&T)) {
    err = rz_set_u64(r, 0);
    goto cleanup;
  }

  /* Horner's scheme: acc = (...((c_d*x + c_{d-1})*x + ...)*x + c_0) */
  err = rz_copy(&acc, &T.coeff[T.deg]);
  if (err != RABIN_SUCCESS) goto cleanup;

  for (u64 i = T.deg; i-- > 0;) {
    err = rz_mul(&acc, &acc, x);
    if (err != RABIN_SUCCESS) goto cleanup;
    err = rz_add(&acc, &acc, &T.coeff[i]);
    if (err != RABIN_SUCCESS) goto cleanup;
  }

  err = rz_copy(r, &acc);

cleanup:
  rpol_clear(&T);
  rz_clear(&acc);
  return err;
}

rabin_err_t rpol_derive(rpol_t* dp, const rpol_t* p)
{
  if (dp == NULL || p == NULL) return RABIN_ERR_NULL_PTR;

  rabin_err_t err;
  rpol_t T = {0}, D = {0};
  rz_t k;
  rpol_init(&T);
  rpol_init(&D);
  rz_init(&k);

  err = rpol_copy_full(&T, p);
  if (err != RABIN_SUCCESS) goto cleanup;
  err = rpol_trim(&T);
  if (err != RABIN_SUCCESS) goto cleanup;

  /* The derivative of a constant (or of 0) is 0 */
  if (rpol_is_zero_poly(&T) || T.deg == 0) {
    rpol_clear(dp);
    rpol_init(dp);
    err = rpol_set_zero_poly(dp);
    goto cleanup;
  }

  err = rpol_alloc(&D, T.deg - 1);
  if (err != RABIN_SUCCESS) goto cleanup;
  D.deg = T.deg - 1;

  for (u64 i = 1; i <= T.deg; i++) {
    err = rz_set_u64(&k, i);
    if (err != RABIN_SUCCESS) goto cleanup;
    err = rz_mul(&D.coeff[i - 1], &T.coeff[i], &k);
    if (err != RABIN_SUCCESS) goto cleanup;
  }

  err = rpol_trim(&D);
  if (err != RABIN_SUCCESS) goto cleanup;

  rpol_clear(dp);
  rpol_init(dp);
  err = rpol_copy_full(dp, &D);

cleanup:
  rpol_clear(&T);
  rpol_clear(&D);
  rz_clear(&k);
  return err;
}

size_t rpol_to_string_len(const rpol_t* p)
{
  if (p == NULL) return 0;

  const u64 count = (p->coeff == NULL || p->size == 0) ? 0 : p->deg + 1;
  size_t len = 3; /* '[' + ']' + NUL */
  for (u64 i = 0; i < count; i++) {
    char* ns = rz_to_string(&p->coeff[i]);
    if (ns == NULL) return 0;
    len += strlen(ns);
    free(ns);
    if (i + 1 < count) len++; /* comma */
  }
  return len;
}

rabin_err_t rpol_to_string(const rpol_t* p, char* out, size_t cap)
{
  if (p == NULL || out == NULL) return RABIN_ERR_NULL_PTR;

  size_t need = rpol_to_string_len(p);
  if (need == 0) return RABIN_ERR_OUT_OF_MEMORY;
  if (need > cap) return RABIN_ERR_INVALID_ARG;

  const u64 count = (p->coeff == NULL || p->size == 0) ? 0 : p->deg + 1;

  size_t w = 0;
  out[w++] = '[';
  for (u64 i = 0; i < count; i++) {
    if (i > 0) out[w++] = ',';
    char* ns = rz_to_string(&p->coeff[i]);
    if (ns == NULL) return RABIN_ERR_OUT_OF_MEMORY;
    const size_t l = strlen(ns);
    memcpy(out + w, ns, l);
    w += l;
    free(ns);
  }
  out[w++] = ']';
  out[w] = '\0';
  return RABIN_SUCCESS;
}

rabin_err_t rpol_init_val(rpol_t* p, const char* str)
{
  if (p == NULL || str == NULL) return RABIN_ERR_NULL_PTR;

  strpx_lexer lx;
  strpx_init(&lx, str);
  if (!strpx_eat(&lx, '[')) return RABIN_ERR_INVALID_ARG;

  /* Pass 1: validate the grammar and count coefficients. */
  u64 count = 0;
  rabin_err_t err = RABIN_SUCCESS;
  rz_t scratch;
  if (rz_init(&scratch) != RABIN_SUCCESS) return RABIN_ERR_OUT_OF_MEMORY;

  if (!strpx_eat(&lx, ']')) {
    for (;;) {
      if (strpx_int(&lx, &scratch) != RABIN_SUCCESS) {
        err = RABIN_ERR_INVALID_ARG;
        goto cleanup;
      }
      count++;
      if (strpx_eat(&lx, ',')) continue;
      break;
    }
    if (!strpx_eat(&lx, ']')) {
      err = RABIN_ERR_INVALID_ARG;
      goto cleanup;
    }
  }
  if (!strpx_done(&lx)) {
    err = RABIN_ERR_INVALID_ARG;
    goto cleanup;
  }

  if (count == 0) {
    rpol_init(p); /* empty state: the zero polynomial */
    goto cleanup;
  }

  /* Pass 2: parse into a temporary array, then commit and trim. */
  rz_t* tmp = malloc(count * sizeof(rz_t));
  if (tmp == NULL) {
    err = RABIN_ERR_OUT_OF_MEMORY;
    goto cleanup;
  }
  u64 ninit = 0;
  for (u64 i = 0; i < count; i++) {
    if (rz_init(&tmp[i]) != RABIN_SUCCESS) {
      err = RABIN_ERR_OUT_OF_MEMORY;
      for (u64 j = 0; j < ninit; j++) rz_clear(&tmp[j]);
      free(tmp);
      goto cleanup;
    }
    ninit++;
  }

  strpx_init(&lx, str);
  strpx_eat(&lx, '[');
  for (u64 i = 0; i < count; i++) {
    if (strpx_int(&lx, &tmp[i]) != RABIN_SUCCESS) {
      err = RABIN_ERR_INVALID_ARG;
      goto abort_commit;
    }
    if (i + 1 < count && !strpx_eat(&lx, ',')) {
      err = RABIN_ERR_INVALID_ARG;
      goto abort_commit;
    }
  }

  err = rpol_set(p, tmp, count - 1);
  if (err == RABIN_SUCCESS) err = rpol_trim(p);

abort_commit:
  for (u64 i = 0; i < ninit; i++) rz_clear(&tmp[i]);
  free(tmp);

cleanup:
  rz_clear(&scratch);
  return err;
}
