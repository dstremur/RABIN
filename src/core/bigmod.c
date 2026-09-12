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
#include "../../include/bigmath.h"

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
  q->is_neg = a->is_neg;

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
  bignum u, v, d, one;
  bool success = false;

  bn_init_multi(&u, &v, &d, &one, NULL);
  bn_set_u64(&one, 1);

  bn_gcd_extended_lehmer(&u, &v, &d, a, m);

  if (bn_cmp(&d, &one) != 0) {
    bn_set_u64(res, 0);
    success = false;
  } else {
    bn_mod(res, &u, m);
    if (res->is_neg) {
      bn_add(res, res, m);
    }
    success = true;
  }

  bn_free_multi(&u, &v, &d, &one, NULL);

  return success;
}
