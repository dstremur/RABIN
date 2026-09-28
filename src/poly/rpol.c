/*
 * rpol.c
 *
 * Polynomials over bignums.
 *
 * This file implements polynomials with rz_t coefficients:
 * initialization, setting from coefficient arrays, copying,
 * comparison, growing the coefficient array, trimming, printing,
 * addition, subtraction, multiplication (schoolbook, rz_t NTT and
 * u64 NTT), fixed-width decomposition and recomposition, and carry
 * propagation.
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
#include <string.h>

#include "../../include/rntt.h"
#include "../../include/rzlogic.h"
#include "../../include/u64.h"

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

    if ((err = rz_copy(&mag, c)) != RABIN_SUCCESS) goto out;
    if (neg) {
      if ((err = rz_neg(&mag, &mag)) != RABIN_SUCCESS) goto out;
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
out:
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

  const rpol_t* longer = (p->deg > q->deg) ? p : q;
  for (u64 i = min + 1; i <= max; i++) {
    if ((err = rz_copy(&r->coeff[i], &longer->coeff[i])) != RABIN_SUCCESS)
      return err;
    r->coeff[i].is_neg = true;
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

  if ((err = rntt_cyclic_forward(&p_hat, p, &ctx)) != RABIN_SUCCESS) goto out;
  if ((err = rntt_cyclic_forward(&q_hat, q, &ctx)) != RABIN_SUCCESS) goto out;
  if ((err = rpol_alloc(&r_hat, ntt_size)) != RABIN_SUCCESS) goto out;

  for (u64 i = 0; i < ntt_size; i++) {
    // MontMul(pR, qR) = (pR * qR * R^-1) mod q = (p*q)R mod q
    if ((err = rz_mont_mul(&r_hat.coeff[i], &p_hat.coeff[i], &q_hat.coeff[i],
                           &ctx.mctx)) != RABIN_SUCCESS)
      goto out;
  }
  r_hat.deg = ntt_size - 1;

  rpol_t r_ntt = {0};
  rpol_init(&r_ntt);
  if ((err = rntt_cyclic_inverse_mont_in(&r_ntt, &r_hat, &ctx)) !=
      RABIN_SUCCESS)
    goto out;

  rpol_clear(r);
  rpol_init(r);
  if ((err = rpol_alloc(r, required_len)) != RABIN_SUCCESS) goto out;

  for (u64 i = 0; i < required_len; i++) {
    if ((err = rz_copy(&r->coeff[i], &r_ntt.coeff[i])) != RABIN_SUCCESS)
      goto out;
  }

  r->deg = required_len - 1;

  if ((err = rpol_trim(r)) != RABIN_SUCCESS) goto out;

  err = RABIN_SUCCESS;
out:
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

  u64_ntt_ctx_t* ctx = u64_ntt_ctx_golden_cached(k);
  if (!ctx) {
    return RABIN_ERR_OUT_OF_MEMORY;
  }

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
  u64_ntt_cyclic_forward(p_hat, p_arr, ctx);
  u64_ntt_cyclic_forward(q_hat, q_arr, ctx);

  // 4. Pointwise Multiplication
  for (u64 i = 0; i < ntt_size; i++) {
    r_hat[i] = u64_mont_mul(p_hat[i], q_hat[i], &ctx->mctx);
  }

  // 5. Perform Inverse NTT
  u64_ntt_cyclic_inverse_montgomery_in(r_arr, r_hat, ctx);

  rabin_err_t err = rpol_clear(r);
  if (err != RABIN_SUCCESS) goto out;
  rpol_init(r);
  if ((err = rpol_alloc(r, required_len)) != RABIN_SUCCESS) goto out;

  for (u64 i = 0; i < required_len; i++) {
    if ((err = rz_set_u64(&r->coeff[i], r_arr[i])) != RABIN_SUCCESS) goto out;
  }

  r->deg = required_len - 1;
  if ((err = rpol_trim(r)) != RABIN_SUCCESS) goto out;

  err = RABIN_SUCCESS;
out:
  // 7. Cleanup
  rz_scratch_release();
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
      if ((err = rz_mul(&t, &p->coeff[i], &q->coeff[j])) != RABIN_SUCCESS)
        goto out;

      if ((err = rz_add(&temp.coeff[i + j], &temp.coeff[i + j], &t)) !=
          RABIN_SUCCESS)
        goto out;
    }
  }

  if ((err = rpol_trim(&temp)) != RABIN_SUCCESS) goto out;

  if ((err = rpol_alloc(r, temp.deg)) != RABIN_SUCCESS) goto out;
  for (u64 i = 0; i <= temp.deg; i++) {
    if ((err = rz_copy(&r->coeff[i], &temp.coeff[i])) != RABIN_SUCCESS)
      goto out;
  }
  r->deg = temp.deg;

  err = RABIN_SUCCESS;
out:
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
  if (err != RABIN_SUCCESS) goto out;

  // mask = 0xFFFF..
  if ((err = rz_lshift(&mask, &RZ_ONE, width)) != RABIN_SUCCESS) goto out;
  if ((err = rz_sub(&mask, &mask, &RZ_ONE)) != RABIN_SUCCESS) goto out;

  u64 i = 0;
  while (!rz_is_zero(&tmp)) {
    if ((err = rz_and(&digit, &tmp, &mask)) != RABIN_SUCCESS) goto out;

    if ((err = rpol_alloc(r, i + 1)) != RABIN_SUCCESS) goto out;
    if ((err = rz_copy(&r->coeff[i], &digit)) != RABIN_SUCCESS) goto out;

    if ((err = rz_rshift(&tmp, &tmp, width)) != RABIN_SUCCESS) goto out;

    i++;
  }

  r->deg = (i > 0) ? i - 1 : 0;

  err = RABIN_SUCCESS;
out:
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
  if (err != RABIN_SUCCESS) goto out;
  if ((err = rz_sub(&mask, &base, &RZ_ONE)) != RABIN_SUCCESS) goto out;
  if ((err = rz_set_u64(&carry, 0)) != RABIN_SUCCESS) goto out;

  u64 i = 0;
  // Iterate through all coefficients plus any remaining carries
  while (i <= r->deg || !rz_is_zero(&carry)) {
    if (i > r->deg) {
      // Expand poly if carry exceeds current deg
      if ((err = rpol_alloc(r, i + 1)) != RABIN_SUCCESS) goto out;
      r->deg = i;
    }

    // total = coeff[i] + carry
    if ((err = rz_add(&total, &r->coeff[i], &carry)) != RABIN_SUCCESS) goto out;

    // carry = total >> bit_width
    if ((err = rz_rshift(&carry, &total, bit_width)) != RABIN_SUCCESS) goto out;

    // coeff[i] = total & mask
    if ((err = rz_and(&r->coeff[i], &total, &mask)) != RABIN_SUCCESS) goto out;

    i++;
  }

  err = RABIN_SUCCESS;
out:
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
    if ((err = rz_lshift(&term, &p->coeff[i], i * bit_width)) != RABIN_SUCCESS)
      goto out;
    // n += term
    if ((err = rz_add(n, n, &term)) != RABIN_SUCCESS) goto out;
  }

  err = RABIN_SUCCESS;
out:
  rz_clear(&term);
  return err;
}

rabin_err_t rpol_test()
{
  printf("--- BigPoly Test --- \n");

  rpol_t p = {0}, q = {0}, r = {0};
  rpol_init(&p);
  rpol_init(&q);
  rpol_init(&r);

  // Let p(x) = 2x^2 + 3x + 1
  i64 p_vals[] = {1, 3, 2};
  if (rpol_set_i64(&p, p_vals, 2) != RABIN_SUCCESS) return RABIN_ERR_OVERFLOW;

  // Let q(x) = 4x + 5
  i64 q_vals[] = {5, 4};
  if (rpol_set_i64(&q, q_vals, 1) != RABIN_SUCCESS) return RABIN_ERR_OVERFLOW;

  printf("\nPolynomial P(x):\n");
  rpol_print(&p);

  printf("\nPolynomial Q(x):\n");
  rpol_print(&q);

  // Test Addition: r = p + q = 2x^2 + 7x + 6
  printf("\n--- Test: Addition (P + Q) ---\n");
  rpol_add(&r, &p, &q);
  rpol_print(&r);

  // Test Multiplication: r = p * q = 8x^3 + 22x^2 + 19x + 5
  printf("\n--- Test: Multiplication (P * Q) ---\n");
  rpol_mul_school(&r, &p, &q);
  rpol_print(&r);

  rpol_clear(&r);
  rpol_init(&r);

  rpol_mul_ntt_u64(&r, &p, &q);
  rpol_print(&r);

  // Cleanup
  rpol_clear(&p);
  rpol_clear(&q);
  rpol_clear(&r);

  printf("\nTests complete.\n");
  return RABIN_SUCCESS;
}
