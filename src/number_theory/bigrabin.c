/*
 * bigrabin.c
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

#include "../include/bignum.h"

/*
 * Miller-Rabin test of n to base a (plain modular arithmetic).
 *
 * Let n_l = n->size, measured in 64-bit limbs.
 *
 * Writes n - 1 = d * 2^s with d odd, computes x = a^d mod n, and
 * checks:
 *
 *   - x == 1 or x == n - 1: probably prime,
 *   - otherwise square x up to s - 1 times; if x becomes n - 1 it is
 *     probably prime, if it becomes 1 (a nontrivial square root of 1)
 *     or never reaches n - 1 it is composite.
 *
 * Returns true if n passes the test (probably prime), false if n is
 * even, n <= 1, or a witness for compositeness was found. n == 2 and
 * n == 3 are reported prime.
 *
 * Complexity:
 *   Time: O(s * n_l^2) - one modular exponentiation plus O(s)
 *         squarings
 *   Auxiliary memory: O(n_l) limbs for temporaries
 *   Output memory: O(1)
 */
bool bn_rabin(const bignum* n, const bignum* a)
{
  if (bn_is_even(n)) return false;

  bignum one, two;
  //  bn_init(&one);
  // bn_init(&two);
  bn_init_multi(&one, &two);
  bn_set_u64(&one, 1);
  bn_set_u64(&two, 2);

  if (bn_cmp(n, &one) <= 0) {
    bn_free(&one);
    bn_free(&two);
    return false;
  }

  if (bn_cmp(n, &two) == 0) {
    bn_free(&one);
    bn_free(&two);
    return true;
  }

  bignum d, n_minus_1;
  bn_init(&d);
  bn_init(&n_minus_1);

  bn_sub(&n_minus_1, n, &one);
  bn_copy(&d, &n_minus_1);

  u64 s = 0;
  while (bn_is_even(&d) && !bn_is_zero(&d)) {
    bn_rshift1(&d);
    s++;
  }

  bignum x, tmp;
  bn_init(&x);
  bn_init(&tmp);

  bn_mod_exp(&x, a, &d, n);

  bool composite = true;

  if (bn_cmp(&x, &one) == 0 || bn_cmp(&x, &n_minus_1) == 0) {
    composite = false;
    goto cleanup;
  }

  for (u64 r = 1; r < s; r++) {
    bn_mul(&tmp, &x, &x);
    bn_mod(&x, &tmp, n);

    if (bn_cmp(&x, &one) == 0) {
      composite = true;
      break;
    }

    if (bn_cmp(&x, &n_minus_1) == 0) {
      composite = false;
      break;
    }
  }

cleanup:
  bn_free(&one);
  bn_free(&two);
  bn_free(&d);
  bn_free(&n_minus_1);
  bn_free(&x);
  bn_free(&tmp);

  return !composite;
}

/*
 * Miller-Rabin test of n to base a (Montgomery domain).
 *
 * Let n_l = n->size, measured in 64-bit limbs.
 *
 * Same test as bn_rabin(), but all modular squarings are Montgomery
 * multiplications, which avoids the division in each reduction. The
 * comparisons are done against one_mont (the Montgomery form of 1) and
 * the Montgomery form of n - 1.
 *
 * Returns true if n passes the test (probably prime), false if n is
 * even, n <= 1, or a witness for compositeness was found. Small
 * single-limb n are handled directly (2 and 3 are prime).
 *
 * Complexity:
 *   Time: O(n_l^3) for the context initialization (see
 *         bn_mont_ctx_init), then O(s * n_l^2) for the test itself
 *   Auxiliary memory: O(n_l) limbs for the context and temporaries
 *   Output memory: O(1)
 */
bool bn_rabin_mont(const bignum* n, const bignum* a)
{
  if (bn_is_even(n)) return false;

  // Fast handling for small numbers
  if (n->size == 1 && n->limbs[0] <= 3) {
    return n->limbs[0] == 2 || n->limbs[0] == 3;
  }

  // Step 1: Find d and s such that n-1 = d * 2^s
  bignum d, n_minus_1, one;
  bn_init(&d);
  bn_init(&n_minus_1);
  bn_init(&one);
  bn_set_u64(&one, 1);

  bn_sub(&n_minus_1, n, &one);  // n-1
  bn_copy(&d, &n_minus_1);

  u64 s = 0;
  while (bn_is_even(&d) && !bn_is_zero(&d)) {
    bn_rshift1(&d);
    s++;
  }

  // Step 2: Setup Montgomery context
  bn_mont_ctx ctx;
  bn_mont_ctx_init(&ctx, n);

  // Step 3: Precompute Montgomery representations
  bignum a_bar, x_bar, n_minus_1_mont, tmp;
  bn_init(&a_bar);
  bn_init(&x_bar);
  bn_init(&n_minus_1_mont);
  bn_init(&tmp);

  bn_alloc(&tmp, n->size);
  bn_alloc(&a_bar, n->size);
  bn_alloc(&x_bar, n->size);
  bn_alloc(&n_minus_1_mont, n->size);

  bn_mont_in(&a_bar, a, &ctx);                    // a in Montgomery
  bn_mont_in(&n_minus_1_mont, &n_minus_1, &ctx);  // (n-1) in Montgomery

  // Step 4: x = a^d mod n (Montgomery)
  bn_mont_exp(&x_bar, &a_bar, &d, &ctx);

  bool composite = true;

  // Step 5: Rabin-Miller checks
  if (bn_cmp(&x_bar, &ctx.one_mont) == 0 ||
      bn_cmp(&x_bar, &n_minus_1_mont) == 0) {
    composite = false;
    goto cleanup;
  }

  // #pragma GCC unroll 4
  for (u64 r = 1; r < s; r++) {
    bn_mont_mul(&x_bar, &x_bar, &x_bar, &ctx);  // x = x^2 mod n
    // bn_copy(&x_bar, &tmp);

    if (bn_cmp(&x_bar, &ctx.one_mont) == 0) {
      composite = true;  // Non-trivial square root of 1
      break;
    }

    if (bn_cmp(&x_bar, &n_minus_1_mont) == 0) {
      composite = false;  // Probably prime
      break;
    }
  }

cleanup:
  // Free all temporaries
  bn_free(&d);
  bn_free(&n_minus_1);
  bn_free(&one);
  bn_free(&a_bar);
  bn_free(&x_bar);
  bn_free(&n_minus_1_mont);
  bn_free(&tmp);
  bn_mont_ctx_free(&ctx);

  return !composite;
}
