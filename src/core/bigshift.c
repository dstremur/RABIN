/*
 * bigshift.c
 *
 * bignum shift routines.
 *
 * This file implements bit-level and limb-level left and right shifts
 * for bignums, both in place (single-bit shifts) and out of place
 * (arbitrary shift amounts).
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

#include <stddef.h>
#include <string.h>

#include "../include/bignum.h"

/**
 * @brief Shift r left by one bit, in place: r = r * 2.
 *
 * Let n = r->size, measured in 64-bit limbs.
 *
 * The bits are propagated from the least significant limb to the most
 * significant one; if the top bit overflows, a new limb is appended.
 * The sign is preserved.
 *
 * Complexity:
 *   Time: O(n)
 *   Auxiliary memory: O(1)
 *   Output memory: O(n) limbs, or O(n + 1) if a carry limb is created
 *
 * @param[in,out] r Value to shift left by one bit (modified in place).
 */
void bn_lshift1(bignum* r)
{
  if (r->size == 0) return;

  u64 carry = 0;

  for (u64 i = 0; i < r->size; i++) {
    u64 new_carry = r->limbs[i] >> 63;
    r->limbs[i] = (r->limbs[i] << 1) | carry;
    carry = new_carry;
  }

  if (carry) {
    bn_alloc(r, r->size + 1);
    r->limbs[r->size++] = carry;
  }
}

/**
 * @brief Shift r right by one bit, in place: r = r / 2 (truncated).
 *
 * Let n = r->size, measured in 64-bit limbs.
 *
 * The bits are propagated from the most significant limb to the least
 * significant one; leading zero limbs are trimmed afterwards. The sign
 * is preserved.
 *
 * Complexity:
 *   Time: O(n)
 *   Auxiliary memory: O(1)
 *   Output memory: O(n) limbs
 *
 * @param[in,out] r Value to shift right by one bit (modified in place).
 */
void bn_rshift1(bignum* r)
{
  if (r->size == 0) return;
  u64 carry = 0;
  for (i64 i = r->size - 1; i >= 0; i--) {
    u64 next = (r->limbs[i] & 1) << 63;
    r->limbs[i] = (r->limbs[i] >> 1) | carry;
    carry = next;
  }

  bn_trim(r);
}

/**
 * @brief Shift a left by shift bits into r: r = a << shift.
 *
 * Let n = a->size, measured in 64-bit limbs.
 * Let w = shift / 64 (whole limbs) and b = shift % 64 (remaining bits).
 *
 * Each limb of a is shifted left by b bits and placed w limbs higher
 * in r, with the overflow bits carried into the next limb. The sign is
 * preserved.
 *
 * r may alias a.
 *
 * Complexity:
 *   Time: O(n)
 *   Auxiliary memory: O(1) in the normal case,
 *                     O(n) if r aliases a and a temporary is used
 *   Output memory: O(n) limbs, at most n + w + 1 limbs
 *
 * @param[out] r     Result of a << shift.
 * @param[in]  a     Value to shift.
 * @param[in]  shift Number of bits to shift left.
 */
void bn_lshift(bignum* r, const bignum* a, int shift)
{
  // aliasing
  if (r == a) {
    bignum tmp;
    bn_init(&tmp);

    bn_lshift(&tmp, a, shift);
    bn_copy(r, &tmp);

    bn_free(&tmp);
    return;
  }

  // early exit
  if (shift == 0) {
    bn_copy(r, a);
    return;
  }

  u64 words = shift / 64;
  u64 bits = shift % 64;

  u64 max_size = a->size + words + 1;
  // if alloc fails just return 0
  if (!bn_alloc(r, max_size)) {
    bn_set_u64(r, 0);
    return;
  }

  // zero the limbs
  memset(r->limbs, 0, max_size * sizeof(u64));

  r->size = max_size;
  r->is_neg = a->is_neg;

  u64 carry = 0;

  // shift left limb by limb
  for (u64 i = 0; i < a->size; i++) {
    u64 src = a->limbs[i];
    u64 dst = i + words;

    // use 128 bit sum, maybe optimize with assembly
    unsigned __int128 sum = (unsigned __int128)src << bits;
    sum += carry;

    // top 64 bits
    r->limbs[dst] = (u64)sum;
    // lower 64 bits go to the carry
    carry = (u64)(sum >> 64);
  }

  // handle remaining carry
  if (carry) {
    r->limbs[a->size + words] = carry;
    r->size = a->size + words + 1;
  } else {
    r->size = a->size + words;
  }

  bn_trim(r);
}

/**
 * @brief Shift a right by shift bits into r: r = a >> shift (truncated).
 *
 * Let n = a->size, measured in 64-bit limbs.
 * Let w = shift / 64 (whole limbs) and b = shift % 64 (remaining bits).
 *
 * The top w limbs are dropped and the remaining limbs are shifted
 * right by b bits, pulling in the low bits of the next limb. If the
 * shift is at least the size of a, the result is zero. The sign is
 * preserved.
 *
 * r may alias a.
 *
 * Complexity:
 *   Time: O(n)
 *   Auxiliary memory: O(1) in the normal case,
 *                     O(n) if r aliases a and a temporary is used
 *   Output memory: O(n) limbs, at most n - w limbs
 *
 * @param[out] r     Result of a >> shift.
 * @param[in]  a     Value to shift.
 * @param[in]  shift Number of bits to shift right.
 */
void bn_rshift(bignum* r, const bignum* a, int shift)
{
  // aliasing
  if (r == a) {
    bignum tmp;
    bn_init(&tmp);

    bn_rshift(&tmp, a, shift);
    bn_copy(r, &tmp);

    bn_free(&tmp);
    return;
  }

  // early exit
  if (shift == 0) {
    bn_copy(r, a);
    return;
  }

  u64 words = shift / 64;
  u64 bits = shift % 64;

  // we shift more than the number of limbs
  if (words >= a->size) {
    bn_set_u64(r, 0);
    return;
  }

  u64 new_size = a->size - words;

  if (!bn_alloc(r, new_size)) {
    bn_set_u64(r, 0);
    return;
  }

  u64* new_limbs = r->limbs;

  if (bits == 0) {
    for (u64 i = 0; i < new_size; i++) {
      new_limbs[i] = a->limbs[i + words];
    }
  } else {
    for (u64 i = 0; i < new_size; i++) {
      u64 curr = a->limbs[i + words];
      u64 next = (i + words + 1 < a->size) ? a->limbs[i + words + 1] : 0;

      r->limbs[i] = (curr >> bits) | (next << (64 - bits));
    }
  }

  r->size = new_size;

  if (r->capacity > r->size) {
    memset(r->limbs + r->size, 0, (r->capacity - r->size) * sizeof(u64));
  }
  bn_trim(r);
}

/**
 * @brief Shift r left by one bit and add bit (0 or 1), in place:
 * r = (r << 1) + bit.
 *
 * Let n = r->size, measured in 64-bit limbs.
 *
 * This is the primitive used when building numbers bit by bit (e.g.
 * during parsing or in exponentiation loops). The added bit enters at
 * the least significant position and the carry propagates upward; if
 * the top bit overflows, a new limb is appended. The sign is
 * preserved.
 *
 * Complexity:
 *   Time: O(n)
 *   Auxiliary memory: O(1)
 *   Output memory: O(n) limbs, or O(n + 1) if a carry limb is created
 *
 * @param[in,out] r   Value to shift and add into (modified in place).
 * @param[in]     bit Bit (0 or 1) to add at the least significant position.
 */
void bn_lshift1_add(bignum* r, int bit)
{
  u64 carry = (bit != 0);
  for (u64 i = 0; i < r->size; i++) {
    u64 tmp = r->limbs[i] >> 63;
    r->limbs[i] = (r->limbs[i] << 1) | carry;
    carry = tmp;
  }

  if (carry) {
    bn_alloc(r, r->size + 1);
    r->limbs[r->size++] = carry;
  }
}
