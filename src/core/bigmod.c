/*
 * bigmod.c
 *
 * bignum modulo routines.
 *
 * This file implements reduction of a bignum modulo a 64-bit divisor
 * (with and without quotient output) and the modular multiplicative
 * inverse via the binary extended GCD.
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

#include "../../include/bigcore.h"
#include "../../include/bighelper.h"

uint64_t bn_mod_u64(const bignum* a, uint64_t d)
{
  unsigned __int128 rem = 0;

  // Purely mathematical reduction with zero memory allocation
  for (i64 i = a->size - 1; i >= 0; i--) {
    unsigned __int128 cur = (rem << 64) | a->limbs[i];
    rem = cur % d;
  }

  u64 res = (u64)rem;

  if (a->is_neg && res != 0) {
    return d - res;
  }

  return res;
}

uint64_t bn_divmod_u64(bignum* q, const bignum* a, uint64_t d)
{
  // aliasing
  if (q == a) {
    bignum tmp;
    bn_init(&tmp);
    uint64_t rem = bn_divmod_u64(&tmp, a, d);
    bn_copy(q, &tmp);
    bn_free(&tmp);
    return rem;
  }
  bn_alloc(q, a->size);
  q->size = a->size;

  unsigned __int128 rem = 0;

  for (i64 i = a->size - 1; i >= 0; i--) {
    unsigned __int128 cur = (rem << 64) | a->limbs[i];

    q->limbs[i] = (u64)(cur / d);
    rem = cur % d;
  }

  bn_trim(q);
  return (u64)rem;
}

bool bn_mod_inverse(bignum* res, const bignum* a, const bignum* m)
{
  // If a is 0 or modulus is <= 1, no inverse exists
  if (bn_is_zero(a) || bn_is_zero(m) || bn_cmp(m, &BN_ONE) == 0) {
    return false;
  }

  bignum u, v, x1, x2;
  bn_init(&u);
  bn_init(&v);
  bn_init(&x1);
  bn_init(&x2);

  bn_copy(&u, a);
  bn_copy(&v, m);
  bn_set_u64(&x1, 1);
  bn_set_u64(&x2, 0);

  while (!bn_is_zero(&u) && !bn_is_zero(&v)) {
    // Eliminate powers of 2 in u
    while (bn_is_even(&u)) {
      bn_rshift1(&u);  // u = u / 2
      if (bn_is_even(&x1)) {
        bn_rshift1(&x1);
      } else {
        bn_add(&x1, &x1, m);
        bn_rshift1(&x1);  // x1 = (x1 + m) / 2
      }
    }

    // Eliminate powers of 2 in v
    while (bn_is_even(&v)) {
      bn_rshift1(&v);  // v = v / 2
      if (bn_is_even(&x2)) {
        bn_rshift1(&x2);
      } else {
        bn_add(&x2, &x2, m);
        bn_rshift1(&x2);  // x2 = (x2 + m) / 2
      }
    }

    // Step-down subtraction
    if (bn_cmp(&u, &v) >= 0) {
      bn_sub(&u, &u, &v);
      // Simulating signed subtraction under unsigned bignum bounds
      if (bn_cmp(&x1, &x2) < 0) {
        bn_add(&x1, &x1, m);
      }
      bn_sub(&x1, &x1, &x2);
    } else {
      bn_sub(&v, &v, &u);
      if (bn_cmp(&x2, &x1) < 0) {
        bn_add(&x2, &x2, m);
      }
      bn_sub(&x2, &x2, &x1);
    }
  }

  bool success = false;
  // If gcd is 1, the matching variable holds the modular inverse
  if (bn_cmp(&u, &BN_ONE) == 0) {
    bn_copy(res, &x1);
    success = true;
  } else if (bn_cmp(&v, &BN_ONE) == 0) {
    bn_copy(res, &x2);
    success = true;
  }

  bn_free(&u);
  bn_free(&v);
  bn_free(&x1);
  bn_free(&x2);

  return success;
}
