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

void bn_isqrt_heron(bignum* r, const bignum* a)
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

void bn_isqrt(bignum* r, const bignum* a) { bn_isqrt_heron(r, a); }
