/*
 * rz_exp.c
 *
 * rz_t exponentiation routines.
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

#include "../../include/rabin.h"

rabin_err_t rz_pow(rz_t* r, const rz_t* a, const rz_t* b)
{
  if (r == NULL || a == NULL || b == NULL) return RABIN_ERR_NULL_PTR;

  if (rz_is_zero(b)) {
    return rz_set_u64(r, 1);
  }

  rabin_err_t err = rz_set_u64(r, 1);
  if (err != RABIN_SUCCESS) return err;

  for (i64 i = rz_bit_length(b) - 1; i >= 0; i--) {
    if ((err = rz_mul(r, r, r)) != RABIN_SUCCESS) return err;
    if (rz_get_bit(b, i)) {
      if ((err = rz_mul(r, r, a)) != RABIN_SUCCESS) return err;
    }
  }
  return RABIN_SUCCESS;
}

rabin_err_t rz_mont_exp(rz_t* r_bar, const rz_t* a_bar, const rz_t* d,
                        rz_mont_ctx* ctx)
{
  if (r_bar == NULL || a_bar == NULL || d == NULL || ctx == NULL)
    return RABIN_ERR_NULL_PTR;

  i64 bits = rz_bit_length(d);
  rabin_err_t err;

  // left-to-right binary for short exponents
  // Algorithm 14.79 Handbook of Applied Crypto
  if (bits < 64) {
    if ((err = rz_copy(r_bar, &ctx->one_mont)) != RABIN_SUCCESS) return err;

    for (i64 i = bits - 1; i >= 0; i--) {
      if ((err = rz_mont_mul(r_bar, r_bar, r_bar, ctx)) != RABIN_SUCCESS)
        return err;

      if (rz_get_bit(d, i)) {
        if ((err = rz_mont_mul(r_bar, r_bar, a_bar, ctx)) != RABIN_SUCCESS)
          return err;
      }
    }
    return RABIN_SUCCESS;
  }

  // Fixed window w = 4: precompute the odd powers a^1, a^3, ..., a^15
  // in Montgomery form (tab[1] = a_bar, tab[v] = tab[v-2] * a^2)
  // Algorithm 14.85 Handbook of Applied Crypto
  rz_t tab[16];
  for (u64 i = 0; i < 16; i++) rz_init(&tab[i]);

  if ((err = rz_copy(&tab[1], a_bar)) != RABIN_SUCCESS) goto cleanup;
  if ((err = rz_mont_mul(&tab[2], a_bar, a_bar, ctx)) != RABIN_SUCCESS)
    goto cleanup;
  for (u64 v = 3; v <= 15; v += 2) {
    if ((err = rz_mont_mul(&tab[v], &tab[v - 2], &tab[2], ctx)) !=
        RABIN_SUCCESS)
      goto cleanup;
  }

  // Start with the most significant bit: r = a
  if ((err = rz_copy(r_bar, a_bar)) != RABIN_SUCCESS) goto cleanup;
  i64 i = bits - 2;

  while (i >= 0) {
    if (!rz_get_bit(d, i)) {
      if ((err = rz_mont_mul(r_bar, r_bar, r_bar, ctx)) != RABIN_SUCCESS)
        goto cleanup;
      i--;
      continue;
    }

    // d[i] == 1: window of length L (1..4) ending in a 1-bit
    i64 L = 1;
    while (L < 4 && (i - L) >= 0 && rz_get_bit(d, i - L)) {
      L++;
    }

    // wval = the L-bit value of bits i .. i-L+1 (odd, 1..15)
    u64 wval = 0;
    for (i64 k = i; k >= i - L + 1; k--) {
      wval = (wval << 1) | (u64)rz_get_bit(d, k);
    }

    // r = r^(2^L) * a^wval
    for (i64 k = 0; k < L; k++) {
      if ((err = rz_mont_mul(r_bar, r_bar, r_bar, ctx)) != RABIN_SUCCESS)
        goto cleanup;
    }
    if ((err = rz_mont_mul(r_bar, r_bar, &tab[wval], ctx)) != RABIN_SUCCESS)
      goto cleanup;

    i = i - L;
  }

  err = RABIN_SUCCESS;
cleanup:
  for (u64 i = 0; i < 16; i++) rz_clear(&tab[i]);
  return err;
}

rabin_err_t rz_mod_exp_slow(rz_t* r, const rz_t* a, const rz_t* b,
                            const rz_t* m)
{
  if (r == NULL || a == NULL || b == NULL || m == NULL)
    return RABIN_ERR_NULL_PTR;

  rz_t base, exp, res, tmp;
  rz_init_multi(&base, &exp, &res, &tmp, NULL);

  rabin_err_t err = rz_copy(&base, a);
  if (err == RABIN_SUCCESS) err = rz_copy(&exp, b);
  if (err == RABIN_SUCCESS) err = rz_set_u64(&res, 1);

  while (err == RABIN_SUCCESS && !rz_is_zero(&exp)) {
    if (!rz_is_even(&exp)) {
      if ((err = rz_mul(&tmp, &res, &base)) != RABIN_SUCCESS) break;
      if ((err = rz_mod(&res, &tmp, m)) != RABIN_SUCCESS) break;
    }
    if ((err = rz_mul(&tmp, &base, &base)) != RABIN_SUCCESS) break;
    if ((err = rz_mod(&base, &tmp, m)) != RABIN_SUCCESS) break;
    if ((err = rz_rshift1(&exp)) != RABIN_SUCCESS) break;
  }

  if (err == RABIN_SUCCESS) err = rz_copy(r, &res);

  rz_clear_multi(&base, &exp, &res, &tmp, NULL);
  return err;
}

rabin_err_t rz_mod_exp(rz_t* r, const rz_t* a, const rz_t* b, const rz_t* m)
{
  if (r == NULL || a == NULL || b == NULL || m == NULL)
    return RABIN_ERR_NULL_PTR;

  if (rz_is_even(m)) {
    return rz_mod_exp_slow(r, a, b, m);
  }

  // fast path
  rz_mont_ctx ctx;
  rabin_err_t err = rz_mont_ctx_init(&ctx, m);
  if (err == RABIN_SUCCESS) {
    err = rz_mod_exp_mont(r, a, b, m, &ctx);
  }
  rz_mont_ctx_clear(&ctx);
  return err;
}

rabin_err_t rz_mod_exp_mont(rz_t* r, const rz_t* a, const rz_t* b,
                            const rz_t* m, rz_mont_ctx* ctx)
{
  if (r == NULL || a == NULL || b == NULL || m == NULL || ctx == NULL)
    return RABIN_ERR_NULL_PTR;

  assert(!rz_is_zero(m) && !rz_is_even(m));
  assert(rz_cmp(&ctx->n, m) == 0);

  rz_t base, result;
  rz_init_multi(&base, &result, NULL);

  rabin_err_t err = rz_mont_in(&base, a, ctx);
  if (err == RABIN_SUCCESS) {
    err = rz_mont_exp(&result, &base, b, ctx);
  }
  if (err == RABIN_SUCCESS) {
    err = rz_mont_out(r, &result, ctx);
  }

  rz_clear(&base);
  rz_clear(&result);
  return err;
}
