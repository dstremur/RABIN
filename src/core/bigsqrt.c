/*
 * bigsqrt.c
 *
 * Integer square root routines.
 *
 * This file implements the integer (floor) square root of a bignum using
 * Heron's method (Newton's method applied to x^2 - a = 0).
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

#include "../../include/bignum.h"

/*
 * Calculate the integer square root of a using Heron's method.
 *
 * Let n = a->size, measured in 64-bit limbs.
 *
 * This computes:
 *
 *   r = floor(sqrt(a))
 *
 * The iteration starts at x0 = 2^(ceil(k/2) + 1) where k is the bit
 * length of a, and repeatedly applies:
 *
 *   x_{i+1} = (x_i + a / x_i) / 2
 *
 * until the sequence stops decreasing. Each iteration roughly doubles
 * the number of correct bits, so O(log n) iterations suffice.
 *
 * If a is negative, r is left unchanged (no real square root exists).
 *
 * Complexity:
 *   Time: O(n^2 log n) - O(log n) iterations, each dominated by a
 *         division of an n-limb value by an n-limb value
 *   Auxiliary memory: O(n) limbs for temporaries
 *   Output memory: O(n) limbs (at most n/2 + 1 limbs)
 */
void bn_isqrt_heron(bignum* r, bignum* a)
{
  if (a->is_neg) return;

  if (bn_is_zero(a)) {
    bn_set_u64(r, 0);
    return;
  }

  bignum xn, xnext, tmp;
  bn_init_multi(&xn, &xnext, &tmp, NULL);

  // set x0 = 2^{log_2(a) / 2 + 1}
  u64 k = bn_bit_length(a);
  bn_set_u64(&xn, 1);

  bn_lshift(&xn, &xn, (k + 1) / 2 + 1);

  while (true) {
    bn_div(&tmp, a, &xn);

    bn_add(&xnext, &xn, &tmp);

    bn_rshift1(&xnext);
    if (bn_cmp(&xnext, &xn) >= 0) {
      break;
    }

    bn_swap(&xn, &xnext);
  }
  bn_copy(r, &xn);

  bn_free_multi(&xn, &xnext, &tmp, NULL);
}

/*
 * Calculate the integer square root of a.
 *
 * Let n = a->size, measured in 64-bit limbs.
 *
 * This computes:
 *
 *   r = floor(sqrt(a))
 *
 * Currently a thin wrapper around bn_isqrt_heron().
 *
 * Complexity:
 *   Time: O(n^2 log n), see bn_isqrt_heron()
 *   Auxiliary memory: O(n) limbs
 *   Output memory: O(n) limbs
 */
void bn_isqrt(bignum* r, bignum* a) { bn_isqrt_heron(r, a); }
