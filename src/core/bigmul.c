/*
 * bigmul.c
 *
 * bignum multiplication routines.
 *
 * This file implements signed bignum multiplication with a schoolbook
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

#include "../../include/bighelper.h"
#include "../../include/bigntt.h"
#include "../../include/bignum.h"
#include "../../include/precomp.h"
#include "../../include/u64.h"

void bn_sqr(bignum* r, const bignum* a)
{
  if (a->size == 0) {
    r->size = 0;
    if (r->capacity > 0) r->limbs[0] = 0;
    return;
  }

  if (r == a) {
    bignum tmp;
    bn_init(&tmp);
    bn_sqr(&tmp, a);
    bn_copy(r, &tmp);
    bn_free(&tmp);
    return;
  }

  bn_alloc(r, 2 * a->size);

  u64 scratch_size = 8 * a->size + 8;  // conservative heuristic
  u64* scratch = bn_scratch_get(scratch_size);

  limbs_sqr_karatsuba(r->limbs, a->limbs, a->size, scratch);

  bn_scratch_release();

  r->size = 2 * a->size;
  r->is_neg = false;
  bn_trim(r);
}

void bn_mul(bignum* r, const bignum* a, const bignum* b)
{
  if (a->size == 0 || b->size == 0) {
    r->size = 0;
    if (r->capacity > 0) r->limbs[0] = 0;
    return;
  }

  if (r == a || r == b) {
    bignum tmp;
    bn_init(&tmp);
    bn_mul(&tmp, a, b);
    bn_copy(r, &tmp);
    bn_free(&tmp);
    return;
  }

  if (a == b) {
    bn_sqr(r, a);
    return;
  }

  u64 max_len = MAX(a->size, b->size);
  bn_alloc(r, 2 * max_len);

  // use karasuba if both inputs are larger then the limit
  if (a->size >= KARATSUBA_LIMIT && b->size >= KARATSUBA_LIMIT) {
    // One arena block: pad_a, pad_b, then the Karatsuba scratch
    uint64_t* buf = bn_scratch_get(10 * max_len);
    uint64_t* pad_a = buf;
    uint64_t* pad_b = buf + max_len;
    uint64_t* scratch = buf + 2 * max_len;

    // Normalize asymmetric arrays with zero-padding up front
    memset(pad_a, 0, max_len * sizeof(uint64_t));
    memset(pad_b, 0, max_len * sizeof(uint64_t));
    memcpy(pad_a, a->limbs, a->size * sizeof(uint64_t));
    memcpy(pad_b, b->limbs, b->size * sizeof(uint64_t));

    limbs_mul_karatsuba(r->limbs, pad_a, pad_b, max_len, scratch);

    bn_scratch_release();
  } else {
    limbs_mul_school(r->limbs, a->limbs, a->size, b->limbs, b->size);
  }

  r->size = a->size + b->size;
  r->is_neg = a->is_neg ^ b->is_neg;
  bn_trim(r);
}

void bn_mul_school(bignum* r, const bignum* a, const bignum* b)
{
  u64 max = a->size + b->size;
  bn_alloc(r, max);

  memset(r->limbs, 0, max * sizeof(u64));
  r->size = max;
  r->is_neg = a->is_neg ^ b->is_neg;

  for (u64 i = 0; i < a->size; i++) {
    if (a->limbs[i] == 0) continue;

    bn_mul_add_inner(&r->limbs[i], b->limbs, a->limbs[i], b->size);
  }

  bn_trim(r);
}

void bn_mul_fast(bignum* res, const bignum* a, const bignum* b)
{
  if (bn_is_zero(a) || bn_is_zero(b)) {
    bn_set_u64(res, 0);
    return;
  }

  u64 bit_width = 16;
  bigpoly poly_a, poly_b, poly_res;
  bigpoly_init(&poly_a);
  bigpoly_init(&poly_b);
  bigpoly_init(&poly_res);

  // 1. Decompose integers into polynomials
  bn_decompose(&poly_a, a, bit_width);
  bn_decompose(&poly_b, b, bit_width);

  // 2. Perform Polynomial NTT Multiplication (Your existing function)
  // Note: bigpoly_mul_ntt internally uses bigntt_ctx_init_golden
  bigpoly_mul_ntt(&poly_res, &poly_a, &poly_b);

  // 3. Ripple carries
  poly_carry_propagation(&poly_res, bit_width);

  // 4. Convert back to bignum
  bn_recompose(res, &poly_res, bit_width);

  // Cleanup
  bigpoly_free(&poly_a);
  bigpoly_free(&poly_b);
  bigpoly_free(&poly_res);
}

void bn_mul_i64(bignum* r, const bignum* a, const i64 c)
{
  u64 carry = 0;
  r->size = a->size;

  bn_alloc(r, a->size + 1);

  u64 abs_c = (c < 0) ? (u64)(-c) : (u64)c;

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
}

bool bn_is_square(bignum* q, const bignum* a)
{
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

  bignum q_1, q_2;
  bn_init_multi(&q_1, &q_2, NULL);
  bn_isqrt(&q_1, a);
  bn_sqr(&q_2, &q_1);
  if (bn_cmp(&q_2, a) != 0) {
    bn_free_multi(&q_1, &q_2, NULL);
    return false;
  }

  if (q) {
    bn_copy(q, &q_1);
  }

  bn_free_multi(&q_1, &q_2, NULL);
  return true;
}

bool bn_is_prime_power(bignum* p, const bignum* n)
{
  bignum a, b, p_0, temp, r;
  bn_init_multi(&a, &b, &p_0, &temp, &r, NULL);
  bn_set_u64(&a, 1);

  bool verdict = false;

  // 2. [Compute GCD]
  while (true) {
    bn_add_u64(&a, &a, 1);
    bn_mod_exp(&b, &a, n, n);

    bn_sub(&temp, &b, &a);
    temp.is_neg = false;

    bn_gcd(&p_0, &temp, n);

    // 3. [Finished?]
    if (bn_is_eq_i64(&p_0, 1)) {
      verdict = false;
      break;
    }

    if (!bn_bpsw(&p_0)) {
      continue;
    }

    // 4. [Final Test]
    bn_copy(&temp, n);
    // repeatadly divide n by p
    while (true) {
      bn_divmod(&b, &r, &temp, &p_0);
      if (!bn_is_zero(&r)) {
        break;
      }
      bn_copy(&temp, &b);
    }

    if (bn_is_eq_i64(&temp, 1)) {
      verdict = true;
      if (p) {
        bn_copy(p, &p_0);
      }
    } else {
      verdict = false;
    }

    break;
  }

  bn_free_multi(&a, &b, &p_0, &temp, &r, NULL);
  return verdict;
}