/*
 * rzrabin.c
 *
 * Miller-Rabin primality test.
 *
 * This file implements the Miller-Rabin (Rabin) probabilistic
 * primality test for a single base a, both in plain modular arithmetic
 * and in the Montgomery domain.
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

#include <ctype.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>

#include "../core/rz_internal.h"
#include "../include/rabin.h"

bool rz_rabin(const rz_t* n, const rz_t* a)
{
  if (n == NULL || a == NULL) return false;
  if (rz_is_even(n)) return false;

  bool composite = true;

  rz_t one, two, d, n_minus_1, x, tmp;
  if (rz_init_multi(&one, &two, &d, &n_minus_1, &x, &tmp, NULL) !=
      RABIN_SUCCESS)
    return false;

  if (rz_set_u64(&one, 1) != RABIN_SUCCESS) goto cleanup;
  if (rz_set_u64(&two, 2) != RABIN_SUCCESS) goto cleanup;

  if (rz_cmp(n, &one) <= 0) {
    goto cleanup;
  }

  if (rz_cmp(n, &two) == 0) {
    composite = false;
    goto cleanup;
  }

  if (rz_sub(&n_minus_1, n, &one) != RABIN_SUCCESS) goto cleanup;
  if (rz_copy(&d, &n_minus_1) != RABIN_SUCCESS) goto cleanup;

  u64 s = 0;
  while (rz_is_even(&d) && !rz_is_zero(&d)) {
    if (rz_rshift1(&d) != RABIN_SUCCESS) goto cleanup;
    s++;
  }

  if (rz_mod_exp(&x, a, &d, n) != RABIN_SUCCESS) goto cleanup;

  if (rz_cmp(&x, &one) == 0 || rz_cmp(&x, &n_minus_1) == 0) {
    composite = false;
    goto cleanup;
  }

  for (u64 r = 1; r < s; r++) {
    if (rz_mul(&tmp, &x, &x) != RABIN_SUCCESS) goto cleanup;
    if (rz_mod(&x, &tmp, n) != RABIN_SUCCESS) goto cleanup;

    if (rz_cmp(&x, &one) == 0) {
      composite = true;
      break;
    }

    if (rz_cmp(&x, &n_minus_1) == 0) {
      composite = false;
      break;
    }
  }

cleanup:
  rz_clear_multi(&one, &two, &d, &n_minus_1, &x, &tmp, NULL);
  return !composite;
}

bool rz_rabin_mont(const rz_t* n, const rz_t* a)
{
  if (n == NULL || a == NULL) return false;
  if (rz_is_even(n)) return false;

  // Fast handling for small numbers
  if (n->size == 1 && n->limbs[0] <= 3) {
    return n->limbs[0] == 2 || n->limbs[0] == 3;
  }

  bool composite = true;

  // Step 1: Find d and s such that n-1 = d * 2^s
  rz_t d, n_minus_1, one;
  rz_init(&d);
  rz_init(&n_minus_1);
  rz_init(&one);

  // Step 2: Setup Montgomery context
  rz_mont_ctx ctx = {0};

  // Step 3: Precompute Montgomery representations
  rz_t a_bar, x_bar, n_minus_1_mont, tmp;
  rz_init(&a_bar);
  rz_init(&x_bar);
  rz_init(&n_minus_1_mont);
  rz_init(&tmp);

  if (rz_set_u64(&one, 1) != RABIN_SUCCESS) goto cleanup;

  if (rz_sub(&n_minus_1, n, &one) != RABIN_SUCCESS) goto cleanup;  // n-1
  if (rz_copy(&d, &n_minus_1) != RABIN_SUCCESS) goto cleanup;

  u64 s = 0;
  while (rz_is_even(&d) && !rz_is_zero(&d)) {
    if (rz_rshift1(&d) != RABIN_SUCCESS) goto cleanup;
    s++;
  }

  if (rz_mont_ctx_init(&ctx, n) != RABIN_SUCCESS) goto cleanup;

  if (rz_alloc(&tmp, n->size) != RABIN_SUCCESS) goto cleanup;
  if (rz_alloc(&a_bar, n->size) != RABIN_SUCCESS) goto cleanup;
  if (rz_alloc(&x_bar, n->size) != RABIN_SUCCESS) goto cleanup;
  if (rz_alloc(&n_minus_1_mont, n->size) != RABIN_SUCCESS) goto cleanup;

  if (rz_mont_in(&a_bar, a, &ctx) != RABIN_SUCCESS)
    goto cleanup;  // a in
                   // Montgomery
  // (n-1) in Montgomery
  if (rz_mont_in(&n_minus_1_mont, &n_minus_1, &ctx) != RABIN_SUCCESS)
    goto cleanup;

  // Step 4: x = a^d mod n (Montgomery)
  if (rz_mont_exp(&x_bar, &a_bar, &d, &ctx) != RABIN_SUCCESS) goto cleanup;

  // Step 5: Rabin-Miller checks
  if (rz_cmp(&x_bar, &ctx.one_mont) == 0 ||
      rz_cmp(&x_bar, &n_minus_1_mont) == 0) {
    composite = false;
    goto cleanup;
  }

  // #pragma GCC unroll 4
  for (u64 r = 1; r < s; r++) {
    if (rz_mont_mul(&x_bar, &x_bar, &x_bar, &ctx) != RABIN_SUCCESS)
      goto cleanup;  // x = x^2 mod n

    if (rz_cmp(&x_bar, &ctx.one_mont) == 0) {
      composite = true;  // Non-trivial square root of 1
      break;
    }

    if (rz_cmp(&x_bar, &n_minus_1_mont) == 0) {
      composite = false;  // Probably prime
      break;
    }
  }

cleanup:
  // Free all temporaries
  rz_clear(&d);
  rz_clear(&n_minus_1);
  rz_clear(&one);
  rz_clear(&a_bar);
  rz_clear(&x_bar);
  rz_clear(&n_minus_1_mont);
  rz_clear(&tmp);
  rz_mont_ctx_clear(&ctx);

  return !composite;
}
bool rz_rabin_mont_ctx(const rz_t* n, const rz_t* a, rz_mont_ctx* ctx)
{
  if (n == NULL || a == NULL || ctx == NULL) return false;
  if (rz_is_even(n)) return false;

  // Fast handling for small numbers
  if (n->size == 1 && n->limbs[0] <= 3) {
    return n->limbs[0] == 2 || n->limbs[0] == 3;
  }

  bool composite = true;

  // Step 1: Find d and s such that n-1 = d * 2^s
  rz_t d, n_minus_1, one;
  rz_init(&d);
  rz_init(&n_minus_1);
  rz_init(&one);

  // Step 2: Precompute Montgomery representations
  rz_t a_bar, x_bar, n_minus_1_mont, tmp;
  rz_init(&a_bar);
  rz_init(&x_bar);
  rz_init(&n_minus_1_mont);
  rz_init(&tmp);

  if (rz_set_u64(&one, 1) != RABIN_SUCCESS) goto cleanup;

  if (rz_sub(&n_minus_1, n, &one) != RABIN_SUCCESS) goto cleanup;  // n-1
  if (rz_copy(&d, &n_minus_1) != RABIN_SUCCESS) goto cleanup;

  u64 s = 0;
  while (rz_is_even(&d) && !rz_is_zero(&d)) {
    if (rz_rshift1(&d) != RABIN_SUCCESS) goto cleanup;
    s++;
  }

  if (rz_alloc(&tmp, n->size) != RABIN_SUCCESS) goto cleanup;
  if (rz_alloc(&a_bar, n->size) != RABIN_SUCCESS) goto cleanup;
  if (rz_alloc(&x_bar, n->size) != RABIN_SUCCESS) goto cleanup;
  if (rz_alloc(&n_minus_1_mont, n->size) != RABIN_SUCCESS) goto cleanup;

  if (rz_mont_in(&a_bar, a, ctx) != RABIN_SUCCESS)
    goto cleanup;  // a in
                   // Montgomery
  // (n-1) in Montgomery
  if (rz_mont_in(&n_minus_1_mont, &n_minus_1, ctx) != RABIN_SUCCESS)
    goto cleanup;

  // Step 3: x = a^d mod n (Montgomery)
  if (rz_mont_exp(&x_bar, &a_bar, &d, ctx) != RABIN_SUCCESS) goto cleanup;

  // Step 4: Rabin-Miller checks
  if (rz_cmp(&x_bar, &ctx->one_mont) == 0 ||
      rz_cmp(&x_bar, &n_minus_1_mont) == 0) {
    composite = false;
    goto cleanup;
  }

  // #pragma GCC unroll 4
  for (u64 r = 1; r < s; r++) {
    if (rz_mont_mul(&x_bar, &x_bar, &x_bar, ctx) != RABIN_SUCCESS)
      goto cleanup;  // x = x^2 mod n

    if (rz_cmp(&x_bar, &ctx->one_mont) == 0) {
      composite = true;  // Non-trivial square root of 1
      break;
    }

    if (rz_cmp(&x_bar, &n_minus_1_mont) == 0) {
      composite = false;  // Probably prime
      break;
    }
  }

cleanup:
  // Free all temporaries
  rz_clear(&d);
  rz_clear(&n_minus_1);
  rz_clear(&one);
  rz_clear(&a_bar);
  rz_clear(&x_bar);
  rz_clear(&n_minus_1_mont);
  rz_clear(&tmp);

  return !composite;
}
