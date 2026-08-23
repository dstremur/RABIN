/*
 * bigmont.c
 *
 * Montgomery multiplication for bignums.
 *
 * This file implements the Montgomery context (modulus, R mod n, R^2
 * mod n, and -n^{-1} mod 2^64), the REDC reduction, and the
 * conversion/multiplication operations of the Montgomery domain.
 *
 * The modulus n must be odd and greater than 1. With R = 2^(64 * n->size),
 * a value A_bar in the Montgomery domain represents A = A_bar * R^{-1}
 * mod n.
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

#include "../../include/bighelper.h"
#include "../../include/bignum.h"

/*
 * Initialize a Montgomery context for the modulus n.
 *
 * Let n_l = n->size, measured in 64-bit limbs.
 *
 * Computes the context fields:
 *
 *   ctx->n        = n
 *   ctx->n_inv    = -n^{-1} mod 2^64   (from the low limb of n)
 *   ctx->one_mont = R mod n            (R = 2^(64 * n_l))
 *   ctx->r_square = R^2 mod n
 *
 * and preallocates ctx->tmp with 2 * n_l + 1 limbs of scratch space.
 *
 * Precondition: n is odd and greater than 1.
 *
 * Complexity:
 *   Time: O(n_l^3) - one_mont is computed with bn_mod_exp_slow(), which
 *         performs a full multiply+divide per exponent bit (64 * n_l
 *         bits); r_square costs one extra multiply+divide
 *   Auxiliary memory: O(n_l) limbs
 *   Output memory: O(n_l) limbs per context field
 */
void bn_mont_ctx_init(bn_mont_ctx* ctx, const bignum* n)
{
  bn_init(&ctx->n);
  bn_init(&ctx->one_mont);
  bn_init(&ctx->r_square);
  bn_init(&ctx->tmp);

  bn_copy(&ctx->n, n);

  // mod_inverse_u64 returns -n^{-1} mod 2^64, exactly what REDC needs
  ctx->n_inv = mod_inverse_u64(n->limbs[0]);

  bignum two, exp_r;
  bn_init_multi(&two, &exp_r, NULL);
  bn_set_u64(&two, 2);

  // Calculate 64 * size as a native u64 first
  u64 exp_val = (u64)n->size * 64;
  bn_set_u64(&exp_r, exp_val);

  // one_mont = 2^(64 * size) mod n
  bn_mod_exp_slow(&ctx->one_mont, &two, &exp_r, n);

  // r_square = (one_mont * one_mont) mod n
  bn_mul(&ctx->r_square, &ctx->one_mont, &ctx->one_mont);
  bn_mod(&ctx->r_square, &ctx->r_square, n);

  bn_free_multi(&two, &exp_r, NULL);

  bn_alloc(&ctx->tmp, 2 * n->size + 1);
}

/*
 * Free all bignum storage held by a Montgomery context.
 *
 * After this call the context must not be used until re-initialized.
 *
 * Complexity:
 *   Time: O(1)
 *   Auxiliary memory: O(1)
 *   Output memory: O(1)
 */
void bn_mont_ctx_free(bn_mont_ctx* ctx)
{
  bn_free(&ctx->one_mont);
  bn_free(&ctx->n);
  bn_free(&ctx->r_square);
  bn_free(&ctx->tmp);
}

/*
 * Montgomery reduction (REDC).
 *
 * Let n_l = ctx->n.size, measured in 64-bit limbs.
 *
 * Given t with 0 <= t < n * 2^(64 * n_l), this computes:
 *
 *   r = t * R^{-1} mod n
 *
 * in O(n_l^2) time. For each limb i it forms the multiple
 * m = t[i] * (-n^{-1} mod 2^64) and adds m * n shifted by i limbs to
 * t, which zeroes out limb i (mod 2^64). After n_l steps the lower
 * half of t is zero, so the result is the upper half, followed by a
 * final conditional subtraction of n to bring r into [0, n).
 *
 * t is destroyed (overwritten) and must have at least 2 * n_l + 1
 * limbs of capacity.
 *
 * Complexity:
 *   Time: O(n_l^2) for the reduction loop, plus O(n_l) for the final
 *         correction
 *   Auxiliary memory: O(1)
 *   Output memory: O(n_l) limbs
 */
void bn_mont_redc(bignum* r, bignum* t, bn_mont_ctx* ctx)
{
  u64 size = ctx->n.size;
  u64* n_limbs = ctx->n.limbs;
  u64* t_limbs = t->limbs;

  // Word for word reduction
  for (u64 i = 0; i < size; i++) {
    u64 m = t_limbs[i] * ctx->n_inv;
    u64 carry = 0;

    // do in assembly
    for (u64 j = 0; j < size; j++) {
      unsigned __int128 prod =
          (unsigned __int128)m * n_limbs[j] + t_limbs[i + j] + carry;
      t_limbs[i + j] = (u64)prod;
      carry = (u64)(prod >> 64);
    }

    // Handle the final carry for this row
    u64 k = i + size;
    unsigned char c = adc64(0, t_limbs[k], carry, &t_limbs[k]);
    k++;
    while (c && k < t->size) {
      c = adc64(c, t_limbs[k], 0, &t_limbs[k]);
      k++;
    }
  }

  // Allocate space for size + 1 limbs to prevent upper-limb truncation
  u64 res_size = size + 1;
  bn_alloc(r, res_size);

  // Copy size + 1 limbs from the upper half of the buffer
  memcpy(r->limbs, &t_limbs[size], res_size * sizeof(u64));
  r->size = res_size;

  // Trim leading zeros so bn_cmp works accurately
  bn_trim(r);

  while (bn_cmp(r, &ctx->n) >= 0) {
    bn_sub_abs(r, r, &ctx->n);
  }
}

/*
 * Convert a value from the normal domain into the Montgomery domain.
 *
 * Let n_l = ctx->n.size, measured in 64-bit limbs.
 *
 * This computes:
 *
 *   A_bar = A * R mod n
 *
 * as one Montgomery multiplication of A by R^2 mod n:
 * REDC(A * R^2) = A * R mod n.
 *
 * Complexity:
 *   Time: O(n_l^2)
 *   Auxiliary memory: O(n_l) limbs (ctx->tmp)
 *   Output memory: O(n_l) limbs
 */
void bn_mont_in(bignum* A_bar, const bignum* A, bn_mont_ctx* ctx)
{
  bn_mont_mul(A_bar, A, &ctx->r_square, ctx);
}

/*
 * Convert a value from the Montgomery domain back to the normal domain.
 *
 * Let n_l = ctx->n.size, measured in 64-bit limbs.
 *
 * This computes:
 *
 *   A = A_bar * R^{-1} mod n
 *
 * as one Montgomery multiplication of A_bar by 1.
 *
 * Complexity:
 *   Time: O(n_l^2)
 *   Auxiliary memory: O(n_l) limbs (ctx->tmp and a temporary for 1)
 *   Output memory: O(n_l) limbs
 */
void bn_mont_out(bignum* A, const bignum* A_bar, bn_mont_ctx* ctx)
{
  bignum one;
  bn_init(&one);
  bn_set_u64(&one, 1);

  bn_mont_mul(A, A_bar, &one, ctx);

  bn_free(&one);
}

/*
 * Montgomery multiplication of two values in the Montgomery domain.
 *
 * Let n_l = ctx->n.size, measured in 64-bit limbs.
 *
 * This computes:
 *
 *   r = a_bar * b_bar * R^{-1} mod n
 *
 * so that if a_bar = A * R mod n and b_bar = B * R mod n, then
 * r = A * B * R mod n (the Montgomery form of A * B).
 *
 * Thin wrapper around bn_mont_mul_raw().
 *
 * Complexity:
 *   Time: O(n_l^2)
 *   Auxiliary memory: O(n_l) limbs (ctx->tmp)
 *   Output memory: O(n_l) limbs
 */
void bn_mont_mul(bignum* r, const bignum* a_bar, const bignum* b_bar,
                 bn_mont_ctx* ctx)
{
  bn_mont_mul_raw(r, a_bar, b_bar, ctx);
}

/*
 * Montgomery multiplication using the context's scratch buffer.
 *
 * Let n_l = ctx->n.size, measured in 64-bit limbs.
 *
 * This computes:
 *
 *   result = A_bar * B_bar * R^{-1} mod n
 *
 * by first forming the full product T = A_bar * B_bar (at most
 * 2 * n_l limbs) in ctx->tmp, zero-padding it to 2 * n_l + 1 limbs,
 * and applying bn_mont_redc().
 *
 * Note: uses ctx->tmp as scratch, so it is not reentrant and must not
 * be called with A_bar or B_bar aliasing ctx->tmp.
 *
 * Complexity:
 *   Time: O(n_l^2) - one bignum multiply plus one REDC
 *   Auxiliary memory: O(n_l) limbs (ctx->tmp)
 *   Output memory: O(n_l) limbs
 */
void bn_mont_mul_raw(bignum* result, const bignum* A_bar, const bignum* B_bar,
                     bn_mont_ctx* ctx)
{
  bignum* T = &ctx->tmp;

  u64 req_size = 2 * ctx->n.size + 1;
  // should already be big enough
  if (T->capacity < req_size) {
    bn_alloc(T, req_size);
  }

  bn_mul(T, A_bar, B_bar);

  if (T->size < req_size) {
    for (u64 i = T->size; i < req_size; i++) {
      T->limbs[i] = 0;
    }
    T->size = req_size;
  }

  bn_mont_redc(result, T, ctx);
}
