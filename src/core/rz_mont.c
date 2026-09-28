/*
 * rz_mont.c
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

#include "../../include/rabin.h"
#include "../../include/rzlimb.h"
#include "rz_internal.h"

rabin_err_t rz_mont_ctx_init(rz_mont_ctx* ctx, const rz_t* n)
{
  if (ctx == NULL || n == NULL) return RABIN_ERR_NULL_PTR;
  if (rz_is_zero(n) || rz_is_even(n)) return RABIN_ERR_INVALID_ARG;

  rz_init(&ctx->n);
  rz_init(&ctx->one_mont);
  rz_init(&ctx->r_square);
  rz_init(&ctx->tmp);

  rabin_err_t err = rz_copy(&ctx->n, n);

  if (err == RABIN_SUCCESS) {
    // mod_inverse_u64 returns -n^{-1} mod 2^64, exactly what REDC needs
    ctx->n_inv = mod_inverse_u64(n->limbs[0]);
  }

  if (err == RABIN_SUCCESS) {
    // one_mont = R mod n where R = 2^(64 * size): build R as a rz_t
    // (a single 1 at limb position size) and reduce it with ONE division
    rz_t r;
    rz_init(&r);

    err = rz_set_u64(&r, 1);
    if (err == RABIN_SUCCESS) err = rz_lshift(&r, &r, (int)(n->size * 64));
    if (err == RABIN_SUCCESS) err = rz_mod(&ctx->one_mont, &r, n);
    if (err == RABIN_SUCCESS) {
      // r_square = (one_mont * one_mont) mod n
      err = rz_mul(&ctx->r_square, &ctx->one_mont, &ctx->one_mont);
    }
    if (err == RABIN_SUCCESS) {
      err = rz_mod(&ctx->r_square, &ctx->r_square, n);
    }

    rz_clear(&r);
  }

  if (err == RABIN_SUCCESS) {
    err = rz_alloc(&ctx->tmp, 2 * n->size + 1);
  }

  if (err != RABIN_SUCCESS) {
    rz_mont_ctx_clear(ctx);
  }
  return err;
}

rabin_err_t rz_mont_ctx_clear(rz_mont_ctx* ctx)
{
  if (ctx == NULL) return RABIN_ERR_NULL_PTR;

  rz_clear(&ctx->one_mont);
  rz_clear(&ctx->n);
  rz_clear(&ctx->r_square);
  rz_clear(&ctx->tmp);
  return RABIN_SUCCESS;
}

rabin_err_t rz_mont_redc(rz_t* r, rz_t* t, const rz_mont_ctx* ctx)
{
  if (r == NULL || t == NULL || ctx == NULL) return RABIN_ERR_NULL_PTR;

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
  rabin_err_t err = rz_alloc(r, res_size);
  if (err != RABIN_SUCCESS) return err;

  // Copy size + 1 limbs from the upper half of the buffer
  memcpy(r->limbs, &t_limbs[size], res_size * sizeof(u64));
  r->size = res_size;

  // Trim leading zeros so rz_cmp works accurately
  rz_trim(r);

  if (rz_cmp(r, &ctx->n) >= 0) {
    return rz_sub_abs(r, r, &ctx->n);
  }
  return RABIN_SUCCESS;
}

rabin_err_t rz_mont_in(rz_t* A_bar, const rz_t* A, rz_mont_ctx* ctx)
{
  if (A_bar == NULL || A == NULL || ctx == NULL) return RABIN_ERR_NULL_PTR;

  // REDC requires both operands < n (T < n*R); reduce A if necessary
  if (rz_cmp(A, &ctx->n) >= 0) {
    rz_t work;
    rz_init(&work);

    rabin_err_t err = rz_mod(&work, A, &ctx->n);
    if (err == RABIN_SUCCESS) {
      err = rz_mont_mul(A_bar, &work, &ctx->r_square, ctx);
    }
    rz_clear(&work);
    return err;
  }

  return rz_mont_mul(A_bar, A, &ctx->r_square, ctx);
}

rabin_err_t rz_mont_out(rz_t* A, const rz_t* A_bar, rz_mont_ctx* ctx)
{
  if (A == NULL || A_bar == NULL || ctx == NULL) return RABIN_ERR_NULL_PTR;

  rz_t one;
  rz_init(&one);
  rabin_err_t err = rz_set_u64(&one, 1);
  if (err == RABIN_SUCCESS) {
    err = rz_mont_mul(A, A_bar, &one, ctx);
  }

  rz_clear(&one);
  return err;
}

rabin_err_t rz_mont_mul(rz_t* r, const rz_t* a_bar, const rz_t* b_bar,
                        rz_mont_ctx* ctx)
{
  return rz_mont_mul_raw(r, a_bar, b_bar, ctx);
}

rabin_err_t rz_mont_mul_raw(rz_t* result, const rz_t* A_bar, const rz_t* B_bar,
                            rz_mont_ctx* ctx)
{
  if (result == NULL || A_bar == NULL || B_bar == NULL || ctx == NULL)
    return RABIN_ERR_NULL_PTR;

  rz_t* T = &ctx->tmp;

  u64 req_size = 2 * ctx->n.size + 1;
  // should already be big enough
  if (T->capacity < req_size) {
    rabin_err_t err = rz_alloc(T, req_size);
    if (err != RABIN_SUCCESS) return err;
  }

  rabin_err_t err = rz_mul(T, A_bar, B_bar);
  if (err != RABIN_SUCCESS) return err;

  if (T->size < req_size) {
    for (u64 i = T->size; i < req_size; i++) {
      T->limbs[i] = 0;
    }
    T->size = req_size;
  }

  return rz_mont_redc(result, T, ctx);
}
