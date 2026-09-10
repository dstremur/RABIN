/*
 * bigexp.c
 *
 * bignum exponentiation routines.
 *
 * This file implements integer exponentiation (a^b) and modular
 * exponentiation (a^b mod m), both in plain form and in the Montgomery
 * domain. The modular paths use left-to-right binary exponentiation.
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

#include <assert.h>
#include <ctype.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>

#include "../../include/bignum.h"

void bn_pow(bignum* r, const bignum* a, const bignum* b)
{
  if (bn_is_zero(b)) {
    bn_set_u64(r, 1);
    return;
  }

  bn_set_u64(r, 1);

  for (i64 i = bn_bit_length(b) - 1; i >= 0; i--) {
    bn_mul(r, r, r);
    if (bn_get_bit(b, i)) {
      bn_mul(r, r, a);
    }
  }
}

void bn_mont_exp(bignum* r_bar, const bignum* a_bar, const bignum* d,
                 const bn_mont_ctx* ctx)
{
  i64 bits = bn_bit_length(d);

  // left-to-right binary for short exponents
  // Algorithm 14.79 Handbook of Applied Crypto
  if (bits < 64) {
    bn_copy(r_bar, &ctx->one_mont);

    for (i64 i = bits - 1; i >= 0; i--) {
      bn_mont_mul(r_bar, r_bar, r_bar, ctx);

      if (bn_get_bit(d, i)) {
        bn_mont_mul(r_bar, r_bar, a_bar, ctx);
      }
    }
    return;
  }

  // Fixed window w = 4: precompute the odd powers a^1, a^3, ..., a^15
  // in Montgomery form (tab[1] = a_bar, tab[v] = tab[v-2] * a^2)
  // Algorithm 14.85 Handbook of Applied Crypto
  bignum tab[16];
  for (u64 i = 0; i < 16; i++) bn_init(&tab[i]);

  bn_copy(&tab[1], a_bar);
  bn_mont_mul(&tab[2], a_bar, a_bar, ctx);  // a^2
  for (u64 v = 3; v <= 15; v += 2) {
    bn_mont_mul(&tab[v], &tab[v - 2], &tab[2], ctx);
  }

  // Start with the most significant bit: r = a
  bn_copy(r_bar, a_bar);
  i64 i = bits - 2;

  while (i >= 0) {
    if (!bn_get_bit(d, i)) {
      bn_mont_mul(r_bar, r_bar, r_bar, ctx);
      i--;
      continue;
    }

    // d[i] == 1: window of length L (1..4) ending in a 1-bit
    i64 L = 1;
    while (L < 4 && (i - L) >= 0 && bn_get_bit(d, i - L)) {
      L++;
    }

    // wval = the L-bit value of bits i .. i-L+1 (odd, 1..15)
    u64 wval = 0;
    for (i64 k = i; k >= i - L + 1; k--) {
      wval = (wval << 1) | (u64)bn_get_bit(d, k);
    }

    // r = r^(2^L) * a^wval
    for (i64 k = 0; k < L; k++) {
      bn_mont_mul(r_bar, r_bar, r_bar, ctx);
    }
    bn_mont_mul(r_bar, r_bar, &tab[wval], ctx);

    i = i - L;
  }

  for (u64 i = 0; i < 16; i++) bn_free(&tab[i]);
}

void bn_mod_exp_slow(bignum* r, const bignum* a, const bignum* b,
                     const bignum* m)
{
  bignum base, exp, res, tmp;
  bn_init_multi(&base, &exp, &res, &tmp, NULL);

  bn_copy(&base, a);
  bn_copy(&exp, b);
  bn_set_u64(&res, 1);

  while (!bn_is_zero(&exp)) {
    if (!bn_is_even(&exp)) {
      bn_mul(&tmp, &res, &base);
      bn_mod(&res, &tmp, m);
    }
    bn_mul(&tmp, &base, &base);
    bn_mod(&base, &tmp, m);
    bn_rshift1(&exp);
  }

  bn_copy(r, &res);
  bn_free_multi(&base, &exp, &res, &tmp, NULL);
}

void bn_mod_exp(bignum* r, const bignum* a, const bignum* b, const bignum* m)
{
  if (bn_is_even(m)) {
    bn_mod_exp_slow(r, a, b, m);
    return;
  }

  // fast path
  bn_mont_ctx ctx;
  bn_mont_ctx_init(&ctx, m);
  bn_mod_exp_mont(r, a, b, m, &ctx);
  bn_mont_ctx_free(&ctx);
}

void bn_mod_exp_mont(bignum* r, const bignum* a, const bignum* b,
                     const bignum* m, const bn_mont_ctx* ctx)
{
  assert(!bn_is_zero(m) && !bn_is_even(m));
  assert(bn_cmp(&ctx->n, m) == 0);

  bignum base, result;
  bn_init_multi(&base, &result, NULL);

  bn_mont_in(&base, a, ctx);

  bn_mont_exp(&result, &base, b, ctx);

  bn_mont_out(r, &result, ctx);

  bn_free(&base);
  bn_free(&result);
}
