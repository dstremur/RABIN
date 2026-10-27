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

// Knuth's Pseudo-Division: l^(m - n + 1) * a = q * b + r
// rabin_err_t rpol_div(rpol_t* q_out, rpol_t* r_out, rpol_t* a, rpol_t* b)
// {
//   if (!a || !b || !q_out || !r_out) return RABIN_ERR_INVALID_ARG;

//   i64 m = a->deg;
//   i64 n = b->deg;

//   if (n < 0 || b->coeff[n] == 0)
//     return RABIN_ERR_INVALID_ARG;  // Division by zero
//   if (m < n) {
//     rpol_zero(q_out);
//     rpol_copy(r_out, a);
//     return RABIN_SUCCESS;
//   }

//   rz_t l = b->coeff[n];  // Leading coefficient of b

//   rpol_t Q, R;
//   rpol_init(&Q);
//   rpol_copy(&R, a);  // R starts as a copy of the dividend 'a'

//   // 1. [Initialize] k = m - n
//   i64 k = m - n;

//   // 2. [Loop] Repeat while k >= 0
//   while (k >= 0) {
//     // 3. [Multiply by lc(v)] Q = l * Q, R = l * R
//     rpol_scale(&Q, l);
//     rpol_scale(&R, l);

//     // 4. [Division step] If deg(R) == n + k, eliminate leading term
//     if (R.deg == n + k) {
//       rz_t c = R.coeff[R.deg];  // Leading coefficient of current remainder

//       // Add term c * x^k to Q
//       rpol_add_term(&Q, c, k);

//       // Subtract (c * x^k) * b from R
//       // This cancels out the leading coefficient cleanly without fractions
//       rpol_sub_scaled_poly(&R, b, c, k);
//     }

//     // 5. [Decrease k]
//     k--;
//   }

//   // Copy results back to user-provided pointers
//   rpol_copy(q_out, &Q);
//   rpol_copy(r_out, &R);

//   // Clean up local polynomials if needed
//   rpol_clear(&Q);

//   return RABIN_SUCCESS;
// }
