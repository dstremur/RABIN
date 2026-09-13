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

  // truncation toward zero keeps the sign of the dividend
  r->is_neg = a->is_neg && !bn_is_zero(r);
}

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

// void bn_lshift_avx512(bignum* r, const bignum* a, int shift)
// {
//   // aliasing
//   if (r == a) {
//     bignum tmp;
//     bn_init(&tmp);
//     bn_lshift_avx512(&tmp, a, shift);
//     bn_copy(r, &tmp);
//     bn_free(&tmp);
//     return;
//   }

//   // early exit for 0
//   if (shift == 0) {
//     bn_copy(r, a);
//     return;
//   }

//   // handle zero-valued bignum
//   if (a->size == 0) {
//     bn_set_u64(r, 0);
//     return;
//   }

//   u64 words = shift / 64;
//   u64 bits = shift % 64;

//   u64 max_size = a->size + words + 1;
//   // if alloc fails just return 0
//   if (!bn_alloc(r, max_size)) {
//     bn_set_u64(r, 0);
//     return;
//   }

//   // zero the limbs
//   memset(r->limbs, 0, max_size * sizeof(u64));

//   r->is_neg = a->is_neg;

//   // Fallback to memcpy for exact word multiples
//   if (bits == 0) {
//     memcpy(r->limbs + words, a->limbs, a->size * sizeof(u64));
//     r->size = a->size + words;
//     bn_trim(r);
//     return;
//   }

//   u64 i = 0;
//   __m512i prev_carry_vec = _mm512_setzero_si512();

//   // Process 8 limbs (512 bits) per iteration
//   for (; i + 8 <= a->size; i += 8) {
//     __m512i vec = _mm512_loadu_si512((const void*)&a->limbs[i]);

//     // Shift limbs left and calculate bits crossing over
//     __m512i shifted = _mm512_slli_epi64(vec, bits);
//     __m512i carries = _mm512_srli_epi64(vec, 64 - bits);

//     // VALIGNQ logic: Shift the 'carries' vector right across elements by 1
//     // and inject the highest limb from the previous vector into lane 0
//     __m512i aligned_carries = _mm512_alignr_epi64(carries, prev_carry_vec,
//     7);

//     // Combine left-shifted bits with incoming carry bits
//     __m512i result = _mm512_or_si512(shifted, aligned_carries);

//     _mm512_storeu_si512((void*)&r->limbs[i + words], result);
//     prev_carry_vec = carries;
//   }

//   // Handle remaining partial vector using AVX-512 masked operations
//   if (i < a->size) {
//     __mmask8 mask = (1U << (a->size - i)) - 1;

//     __m512i vec = _mm512_maskz_loadu_epi64(mask, &a->limbs[i]);

//     __m512i shifted = _mm512_slli_epi64(vec, bits);
//     __m512i carries = _mm512_srli_epi64(vec, 64 - bits);
//     __m512i aligned_carries = _mm512_alignr_epi64(carries, prev_carry_vec,
//     7);

//     __m512i result = _mm512_or_si512(shifted, aligned_carries);

//     _mm512_mask_storeu_epi64(&r->limbs[i + words], mask, result);
//   }

//   // The final scalar carry out of the most significant limb
//   u64 final_carry = a->limbs[a->size - 1] >> (64 - bits);

//   if (final_carry) {
//     r->limbs[a->size + words] = final_carry;
//     r->size = a->size + words + 1;
//   } else {
//     r->size = a->size + words;
//   }

//   bn_trim(r);
// }
