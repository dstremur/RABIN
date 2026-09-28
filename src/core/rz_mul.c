/*
 * rz_mul.c
 *
 * rz_t multiplication routines.
 *
 * This file implements signed rz_t multiplication with a schoolbook
 * and a Karatsuba path (dispatched by size), squaring, and an NTT-based
 * fast multiplication path for very large operands. The raw-limb
 * kernels live in bighelper.c.
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

#include <stddef.h>
#include <stdio.h>
#include <string.h>

#include "../../include/precomp.h"
#include "../../include/rabin.h"
#include "../../include/rntt.h"
#include "../../include/rzlimb.h"
#include "../../include/u64.h"
#include "rz_internal.h"

rabin_err_t rz_sqr(rz_t* r, const rz_t* a)
{
  if (r == NULL || a == NULL) return RABIN_ERR_NULL_PTR;

  if (a->size == 0) {
    r->size = 0;
    if (r->capacity > 0) r->limbs[0] = 0;
    return RABIN_SUCCESS;
  }

  if (r == a) {
    rz_t tmp;
    rz_init(&tmp);
    rabin_err_t err = rz_sqr(&tmp, a);
    if (err != RABIN_SUCCESS) {
      rz_clear(&tmp);
      return err;
    }
    err = rz_copy(r, &tmp);
    rz_clear(&tmp);
    return err;
  }

  rabin_err_t err = rz_alloc(r, 2 * a->size);
  if (err != RABIN_SUCCESS) return err;

  u64 scratch_size = 8 * a->size + 8;  // conservative heuristic
  u64* scratch = rz_scratch_get(scratch_size);
  if (scratch == NULL) return RABIN_ERR_OUT_OF_MEMORY;

  limbs_sqr_karatsuba(r->limbs, a->limbs, a->size, scratch);

  rz_scratch_release();

  r->size = 2 * a->size;
  r->is_neg = false;
  rz_trim(r);
  return RABIN_SUCCESS;
}

rabin_err_t rz_mul(rz_t* r, const rz_t* a, const rz_t* b)
{
  if (r == NULL || a == NULL || b == NULL) return RABIN_ERR_NULL_PTR;

  if (a->size == 0 || b->size == 0) {
    r->size = 0;
    if (r->capacity > 0) r->limbs[0] = 0;
    return RABIN_SUCCESS;
  }

  if (r == a || r == b) {
    rz_t tmp;
    rz_init(&tmp);
    rabin_err_t err = rz_mul(&tmp, a, b);
    if (err != RABIN_SUCCESS) {
      rz_clear(&tmp);
      return err;
    }
    err = rz_copy(r, &tmp);
    rz_clear(&tmp);
    return err;
  }

  if (a == b) {
    return rz_sqr(r, a);
  }

  u64 max_len = RZ_MAX(a->size, b->size);
  rabin_err_t err = rz_alloc(r, 2 * max_len);
  if (err != RABIN_SUCCESS) return err;

  // use karasuba if both inputs are larger then the limit
  if (a->size >= RZ_KARATSUBA_LIMIT && b->size >= RZ_KARATSUBA_LIMIT) {
    // One arena block: pad_a, pad_b, then the Karatsuba scratch
    uint64_t* buf = rz_scratch_get(10 * max_len);
    if (buf == NULL) return RABIN_ERR_OUT_OF_MEMORY;
    uint64_t* pad_a = buf;
    uint64_t* pad_b = buf + max_len;
    uint64_t* scratch = buf + 2 * max_len;

    // Normalize asymmetric arrays with zero-padding up front
    memset(pad_a, 0, max_len * sizeof(uint64_t));
    memset(pad_b, 0, max_len * sizeof(uint64_t));
    memcpy(pad_a, a->limbs, a->size * sizeof(uint64_t));
    memcpy(pad_b, b->limbs, b->size * sizeof(uint64_t));

    limbs_mul_karatsuba(r->limbs, pad_a, pad_b, max_len, scratch);

    rz_scratch_release();
  } else {
    limbs_mul_school(r->limbs, a->limbs, a->size, b->limbs, b->size);
  }

  r->size = a->size + b->size;
  r->is_neg = a->is_neg ^ b->is_neg;
  rz_trim(r);
  return RABIN_SUCCESS;
}

rabin_err_t rz_mul_school(rz_t* r, const rz_t* a, const rz_t* b)
{
  if (r == NULL || a == NULL || b == NULL) return RABIN_ERR_NULL_PTR;

  u64 max = a->size + b->size;
  rabin_err_t err = rz_alloc(r, max);
  if (err != RABIN_SUCCESS) return err;

  memset(r->limbs, 0, max * sizeof(u64));
  r->size = max;
  r->is_neg = a->is_neg ^ b->is_neg;

  for (u64 i = 0; i < a->size; i++) {
    if (a->limbs[i] == 0) continue;

    rz_mul_add_inner(&r->limbs[i], b->limbs, a->limbs[i], b->size);
  }

  rz_trim(r);
  return RABIN_SUCCESS;
}

rabin_err_t rz_mul_fast(rz_t* res, const rz_t* a, const rz_t* b)
{
  if (res == NULL || a == NULL || b == NULL) return RABIN_ERR_NULL_PTR;

  if (rz_is_zero(a) || rz_is_zero(b)) {
    return rz_set_u64(res, 0);
  }

  u64 bit_width = 16;
  rpol_t poly_a = {0}, poly_b = {0}, poly_res = {0};
  rpol_init(&poly_a);
  rpol_init(&poly_b);
  rpol_init(&poly_res);

  rabin_err_t err;

  // 1. Decompose integers into polynomials
  err = rz_decompose(&poly_a, a, bit_width);
  if (err != RABIN_SUCCESS) goto cleanup;
  err = rz_decompose(&poly_b, b, bit_width);
  if (err != RABIN_SUCCESS) goto cleanup;

  // 2. Perform Polynomial NTT Multiplication (Your existing function)
  // Note: rpol_mul_ntt internally uses rntt_ctx_init_golden
  err = rpol_mul_ntt(&poly_res, &poly_a, &poly_b);
  if (err != RABIN_SUCCESS) goto cleanup;

  // 3. Ripple carries
  err = rpol_carry_propagation(&poly_res, bit_width);
  if (err != RABIN_SUCCESS) goto cleanup;

  // 4. Convert back to rz_t
  err = rz_recompose(res, &poly_res, bit_width);

cleanup:
  rpol_clear(&poly_a);
  rpol_clear(&poly_b);
  rpol_clear(&poly_res);
  return err;
}

rabin_err_t rz_mul_i64(rz_t* r, const rz_t* a, const i64 c)
{
  if (r == NULL || a == NULL) return RABIN_ERR_NULL_PTR;

  u64 carry = 0;
  r->size = a->size;

  rabin_err_t err = rz_alloc(r, a->size + 1);
  if (err != RABIN_SUCCESS) return err;

  // In two's complement, 0 - (u64)c safely computes |c| for all negative c,
  // including INT64_MIN.
  u64 abs_c = (c < 0) ? (0 - (u64)c) : (u64)c;

  for (u64 i = 0; i < a->size; i++) {
    u128 prod = (u128)a->limbs[i] * abs_c + carry;
    r->limbs[i] = (u64)prod;
    carry = (u64)(prod >> 64);
  }

  if (carry) {
    r->limbs[r->size++] = carry;
  }

  bool neg = (c < 0);
  r->is_neg = a->is_neg ^ neg;
  return RABIN_SUCCESS;
}

bool rz_is_square(rz_t* q, const rz_t* a)
{
  if (a == NULL) return false;

  // 1. [Test 64]
  u64 t = a->limbs[0] & 63;

  if (q64[t] == 0) return false;

  u64 r = 0;
  for (int i = a->size - 1; i >= 0; --i)
    r = (r * 16 + (a->limbs[i] % 45045)) % 45045;

  // 2. [Test 63
  if (q63[r % 63] == 0) return false;
  // 3. [Test 65]
  if (q65[r % 65] == 0) return false;
  // 4. [Test 11]
  if (q11[r % 11] == 0) return false;

  rz_t q_1, q_2;
  rz_init_multi(&q_1, &q_2, NULL);

  rabin_err_t err = rz_isqrt(&q_1, a);
  if (err != RABIN_SUCCESS) {
    rz_clear_multi(&q_1, &q_2, NULL);
    return false;
  }
  err = rz_sqr(&q_2, &q_1);
  if (err != RABIN_SUCCESS) {
    rz_clear_multi(&q_1, &q_2, NULL);
    return false;
  }
  if (rz_cmp(&q_2, a) != 0) {
    rz_clear_multi(&q_1, &q_2, NULL);
    return false;
  }

  if (q) {
    rz_copy(q, &q_1);
  }

  rz_clear_multi(&q_1, &q_2, NULL);
  return true;
}

bool rz_is_prime_power(rz_t* p, const rz_t* n)
{
  if (n == NULL) return false;

  rz_t a, b, p_0, temp, r;
  rz_init_multi(&a, &b, &p_0, &temp, &r, NULL);
  if (rz_set_u64(&a, 1) != RABIN_SUCCESS) {
    rz_clear_multi(&a, &b, &p_0, &temp, &r, NULL);
    return false;
  }

  bool verdict = false;

  // 2. [Compute GCD]
  while (true) {
    if (rz_add_u64(&a, &a, 1) != RABIN_SUCCESS) break;
    if (rz_mod_exp(&b, &a, n, n) != RABIN_SUCCESS) break;

    rz_sub(&temp, &b, &a);
    temp.is_neg = false;

    if (rz_gcd(&p_0, &temp, n) != RABIN_SUCCESS) break;

    // 3. [Finished?]
    if (rz_is_eq_i64(&p_0, 1)) {
      verdict = false;
      break;
    }

    if (!rz_bpsw(&p_0)) {
      continue;
    }

    // 4. [Final Test]
    if (rz_copy(&temp, n) != RABIN_SUCCESS) break;
    // repeatadly divide n by p
    while (true) {
      if (rz_divmod(&b, &r, &temp, &p_0) != RABIN_SUCCESS) break;
      if (!rz_is_zero(&r)) {
        break;
      }
      if (rz_copy(&temp, &b) != RABIN_SUCCESS) break;
    }

    if (rz_is_eq_i64(&temp, 1)) {
      verdict = true;
      if (p) {
        rz_copy(p, &p_0);
      }
    } else {
      verdict = false;
    }

    break;
  }

  rz_clear_multi(&a, &b, &p_0, &temp, &r, NULL);
  return verdict;
}