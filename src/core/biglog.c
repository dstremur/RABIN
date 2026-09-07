/*
 * biglog.c
 *
 * Integer logarithm routines.
 *
 * This file implements integer logarithms for bignums: the base-2
 * logarithm (exact, via the bit length) and a fixed-point approximation
 * of the natural logarithm.
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

void bn_ln(bignum* r, bignum* a)
{
  if (bn_is_zero(a)) {
    return;
  }

  // ln(2) = 0.6931471805599453094172321214

  bignum log2, tmp, tmp2;
  bn_init_multi(&log2, &tmp, &tmp2, NULL);

  bn_log_2(&log2, a);

  bn_set_u64(&tmp, 69314718);
  bn_set_u64(&tmp2, 100000000);

  // r = (k * 69314718) / 100000000
  bn_mul(r, &log2, &tmp);
  bn_div(r, r, &tmp2);

  bn_free_multi(&log2, &tmp, &tmp2, NULL);
}

void bn_log_2(bignum* r, bignum* a)
{
  if (bn_is_zero(a)) {
    return;
  }

  u64 k = bn_bit_length(a) - 1;

  bn_set_u64(r, k);
}
