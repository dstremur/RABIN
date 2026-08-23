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

/**
 * @brief Calculate the integer natural logarithm of a, truncated.
 *
 * Let n = a->size, measured in 64-bit limbs.
 *
 * This computes:
 *
 *   r = floor(ln(a))
 *
 * using the identity ln(a) = ln(2) * log_2(a). The base-2 logarithm is
 * exact (see bn_log_2), and ln(2) is approximated by the fixed-point
 * rational 69314718 / 100000000 (8 decimal digits), so the result is
 * accurate to within roughly 1 for large a.
 *
 * If a is zero, r is left unchanged (ln(0) is undefined).
 *
 * Complexity:
 *   Time: O(n) for the bit length, plus O(n^2) for the bn_mul() and
 *         O(n^2) for the bn_div() of two n-limb values
 *   Auxiliary memory: O(n) limbs for temporaries
 *   Output memory: O(1) limbs (the result fits in a few limbs)
 *
 * @param[out] r Result storing floor(ln(a)).
 * @param[in]  a Value to take the natural logarithm of.
 */
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

/**
 * @brief Calculate the integer base-2 logarithm of a.
 *
 * Let n = a->size, measured in 64-bit limbs.
 *
 * This computes:
 *
 *   r = floor(log_2(a))
 *
 * which is exact and equals the bit length of a minus one.
 *
 * If a is zero, r is left unchanged (log_2(0) is undefined).
 *
 * Complexity:
 *   Time: O(1) (only the most significant limb is inspected)
 *   Auxiliary memory: O(1)
 *   Output memory: O(1) limbs
 *
 * @param[out] r Result storing floor(log_2(a)).
 * @param[in]  a Value to take the base-2 logarithm of.
 */
void bn_log_2(bignum* r, bignum* a)
{
  if (bn_is_zero(a)) {
    return;
  }

  u64 k = bn_bit_length(a) - 1;

  bn_set_u64(r, k);
}
