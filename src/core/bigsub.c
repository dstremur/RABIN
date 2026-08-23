/*
 * bigsub.c
 *
 * bignum subtraction routines.
 *
 * This file implements magnitude subtraction and signed subtraction for
 * the bignum library.
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

#include "../../include/bighelper.h"
#include "../../include/bignum.h"

/*
 * Subtract the absolute value of b from the absolute value of a into r.
 *
 * Let n = a->size, measured in 64-bit limbs.
 *
 * This computes:
 *
 *   r = |a| - |b|
 *
 * and requires |a| >= |b|; if |b| > |a| the result wraps around
 * (two's-complement style) and is meaningless. It does not interpret or
 * modify the sign of the operands or the result. The caller is
 * responsible for setting r->is_neg.
 *
 * r may alias a or b.
 *
 * Complexity:
 *   Time: O(n)
 *   Auxiliary memory: O(1) in the normal case,
 *                     O(n) if r aliases a or b and a temporary is used
 *   Output memory: O(n) limbs
 */
void bn_sub_abs(bignum* r, const bignum* a, const bignum* b)
{
  // aliasing
  if (r == a || r == b) {
    bignum tmp;
    bn_init(&tmp);
    bn_sub_abs(&tmp, a, b);
    bn_copy(r, &tmp);
    bn_free(&tmp);
    return;
  }
  bn_alloc(r, a->size);

  u64 borrow = 0;

  for (u64 i = 0; i < a->size; i++) {
    u64 av = a->limbs[i];
    u64 bv = (i < b->size) ? b->limbs[i] : 0;

    unsigned __int128 diff = (unsigned __int128)av - bv - borrow;

    r->limbs[i] = (u64)diff;
    borrow = (diff >> 127) & 1;  // detect underflow
  }

  r->size = a->size;
  bn_trim(r);
}

/*
 * Subtract two signed bignums: r = a - b.
 *
 * Let n = max(a->size, b->size), measured in 64-bit limbs.
 *
 * Subtraction is reduced to magnitude addition or subtraction:
 *
 *   a - (-b) = a + b          (opposite signs, b negative)
 *   (-a) - b = -(a + b)       (opposite signs, a negative)
 *   same signs: subtract the smaller magnitude from the larger one and
 *               take the sign of the operand with the larger magnitude
 *
 * If the magnitudes are equal, the result is zero.
 *
 * Complexity:
 *   Time: O(n)
 *   Auxiliary memory: O(1), except for any temporary storage used by
 *                     bn_add_abs(), bn_sub_abs(), or bn_alloc()
 *   Output memory: O(n) limbs
 */
void bn_sub(bignum* r, const bignum* a, const bignum* b)
{
  // a - (-b) = a + b
  if (!a->is_neg && b->is_neg) {
    bn_add_abs(r, a, b);
    r->is_neg = false;
    return;
  }

  // (-a) - b = -(a + b)
  if (a->is_neg && !b->is_neg) {
    bn_add_abs(r, a, b);
    r->is_neg = true;
    return;
  }

  i64 cmp = bn_cmp_abs(a, b);

  if (cmp == 0) {
    bn_set_u64(r, 0);
    return;
  }

  if (cmp > 0) {
    // |a| > |b|
    bn_sub_abs(r, a, b);
    r->is_neg = a->is_neg;
  } else {
    // |b| > |a|
    bn_sub_abs(r, b, a);
    r->is_neg = !a->is_neg;
  }
}
