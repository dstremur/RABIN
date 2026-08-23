/*
 * bigadd.c
 *
 * bignum addition routines.
 *
 * This file implements magnitude addition, signed addition, and helper
 * addition operations for the bignum library.
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
 * Add the absolute values of a and b into r.
 *
 * Let n = max(a->size, b->size), measured in 64-bit limbs.
 *
 * This function only adds magnitudes. It does not interpret or modify
 * the sign of the operands or the result. The caller is responsible for
 * setting r->is_neg.
 *
 * r may alias a or b.
 *
 * Complexity:
 *   Time: O(n)
 *   Auxiliary memory: O(1) in the normal case,
 *                     O(n) if r aliases a or b and a temporary is used
 *   Output memory: O(n) limbs, at most n + 1 limbs
 */
void bn_add_abs(bignum* r, const bignum* a, const bignum* b)
{
  // aliasing
  if (r == a || r == b) {
    bignum tmp;
    bn_init(&tmp);
    bn_add_abs(&tmp, a, b);
    bn_copy(r, &tmp);
    bn_free(&tmp);
    return;
  }

  // ensure a is always larger
  if (a->size < b->size) {
    const bignum* tmp = a;
    a = b;
    b = tmp;
  }

  if (r->capacity < a->size + 1) {
    bn_alloc(r, a->size + 1);
  }

  u64 carry = bn_add_inner(r->limbs, a->limbs, a->size, b->limbs, b->size);

  r->limbs[a->size] = carry;
  r->size = a->size + carry;
}

/*
 * Add two signed bignums.
 *
 * Let n = max(a->size, b->size), measured in 64-bit limbs.
 *
 * If a and b have the same sign, their magnitudes are added.
 *
 * If a and b have different signs, the smaller magnitude is subtracted
 * from the larger magnitude, and the result takes the sign of the operand
 * with the larger magnitude.
 *
 * If the magnitudes are equal, the result is normalized to positive zero.
 *
 * Complexity:
 *   Time: O(n)
 *   Auxiliary memory: O(1), except for any temporary storage used by
 *                     bn_add_abs(), bn_sub_abs(), bn_copy(), or bn_alloc()
 *   Output memory: O(n) limbs
 */
void bn_add(bignum* r, const bignum* a, const bignum* b)
{
  if (a->is_neg == b->is_neg) {
    bn_add_abs(r, a, b);
    r->is_neg = a->is_neg;
    return;
  }

  i64 cmp = bn_cmp_abs(a, b);

  if (cmp == 0) {
    bn_set_u64(r, 0);
    r->is_neg = false;
    return;
  }

  if (cmp > 0) {
    bn_sub_abs(r, a, b);
    r->is_neg = a->is_neg;
  } else {
    bn_sub_abs(r, b, a);
    r->is_neg = b->is_neg;
  }
}

/*
 * Add an unsigned 64-bit value to a signed bignum.
 *
 * Let n = a->size, measured in 64-bit limbs.
 *
 * This computes:
 *
 *   r = a + b
 *
 * where b is a nonnegative u64.
 *
 * If a is positive, this is ordinary magnitude addition.
 *
 * If a is negative, this is magnitude subtraction:
 *
 *   r = -|a| + b
 *
 * which is equivalent to subtracting b from |a| and preserving the
 * correct sign.
 *
 * r may alias a, assuming bn_copy() supports aliasing.
 *
 * Complexity:
 *   Time: O(n) worst case.
 *         The carry/borrow loop may stop early if propagation ends
 *         before reaching the most significant limb.
 *   Auxiliary memory: O(1), except for any temporary storage used by
 *                     bn_copy() or bn_alloc()
 *   Output memory: O(n) limbs, or O(n + 1) if a new carry limb is created
 */
void bn_add_u64(bignum* r, const bignum* a, u64 b)
{
  if (b == 0) {
    bn_copy(r, a);
    return;
  }

  if (!a->is_neg) {
    // POSITIVE + POSITIVE
    bn_copy(r, a);
    u64 carry = b;
    for (u64 i = 0; i < r->size; i++) {
      u64 old = r->limbs[i];
      r->limbs[i] += carry;
      if (r->limbs[i] >= old) {
        carry = 0;
        break;
      }  // Optimization: exit early
      carry = 1;
    }
    if (carry) {
      if (r->capacity < r->size + 1) bn_alloc(r, r->size + 1);
      r->limbs[r->size++] = 1;
    }
  } else {
    // NEGATIVE + POSITIVE (e.g., -80 + 16)
    // This is Magnitude Subtraction: |a| - b
    if (a->size == 1 && a->limbs[0] <= b) {
      bn_set_u64(r, b - a->limbs[0]);
      r->is_neg = false;
    } else {
      bn_copy(r, a);
      u64 borrow = b;
      for (u64 i = 0; i < r->size && borrow; i++) {
        u64 old = r->limbs[i];
        r->limbs[i] -= borrow;
        borrow = (old < borrow) ? 1 : 0;
      }
      r->is_neg = true;
      bn_trim(r);
      if (r->size == 0) r->is_neg = false;
    }
  }
}

/*
 * Add a to r at the given limb offset.
 *
 * Let n = a->size, measured in 64-bit limbs.
 * Let o = offset.
 * Let m = o + n + 1 be the maximum resulting size in limbs.
 *
 * This effectively performs:
 *
 *   r = r + (a << (offset * 64))
 *
 * for unsigned magnitudes.
 *
 * The addition loop itself only touches about n + 1 limbs, but allocation
 * and trimming may need to consider the full output range up to m limbs.
 *
 * Complexity:
 *   Time: O(m) worst case, where m = offset + a->size + 1.
 *         The inner addition loop is O(n), but bn_alloc() and bn_trim()
 *         may make the total worst case proportional to the output size.
 *   Auxiliary memory: O(1)
 *   Output memory: O(m) limbs
 *
 * Note:
 *   The current implementation does not explicitly handle r == a safely
 *   for all nonzero offsets. If aliasing is not supported, document that
 *   here and enforce it in the API.
 */
void bn_add_at_offset(bignum* r, const bignum* a, u64 offset)
{
  if (a->size == 0) return;

  u64 carry = 0;
  // Ensure r has enough limbs to hold a + carry at this offset
  u64 required_size = offset + a->size + 1;
  if (r->size < required_size) {
    bn_alloc(r, required_size);
    r->size = required_size;
  }

  for (u64 i = 0; (i < a->size || carry > 0); i++) {
    u64 r_idx = offset + i;
    u64 av = (i < a->size) ? a->limbs[i] : 0;
    u64 rv = r->limbs[r_idx];

    unsigned __int128 sum = (unsigned __int128)av + rv + carry;
    r->limbs[r_idx] = (u64)sum;
    carry = (u64)(sum >> 64);
  }
  bn_trim(r);
}
