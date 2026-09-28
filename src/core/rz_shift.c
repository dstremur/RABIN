/*
 * rz_shift.c
 *
 * Bit shift routines for rz_t.
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

#include "../include/rabin.h"
#include "rz_internal.h"

rabin_err_t rz_lshift1(rz_t* r)
{
  if (r == NULL) return RABIN_ERR_NULL_PTR;
  if (r->size == 0) return RABIN_SUCCESS;

  u64 carry = 0;

  for (u64 i = 0; i < r->size; i++) {
    u64 new_carry = r->limbs[i] >> 63;
    r->limbs[i] = (r->limbs[i] << 1) | carry;
    carry = new_carry;
  }

  if (carry) {
    rabin_err_t err = rz_alloc(r, r->size + 1);
    if (err != RABIN_SUCCESS) return err;
    r->limbs[r->size++] = carry;
  }
  return RABIN_SUCCESS;
}

rabin_err_t rz_rshift1(rz_t* r)
{
  if (r == NULL) return RABIN_ERR_NULL_PTR;
  if (r->size == 0) return RABIN_SUCCESS;

  u64 carry = 0;
  for (i64 i = r->size - 1; i >= 0; i--) {
    u64 next = (r->limbs[i] & 1) << 63;
    r->limbs[i] = (r->limbs[i] >> 1) | carry;
    carry = next;
  }

  rz_trim(r);
  return RABIN_SUCCESS;
}

rabin_err_t rz_lshift(rz_t* r, const rz_t* a, int shift)
{
  if (r == NULL || a == NULL) return RABIN_ERR_NULL_PTR;

  // aliasing
  if (r == a) {
    rz_t tmp;
    rz_init(&tmp);

    rabin_err_t err = rz_lshift(&tmp, a, shift);
    if (err != RABIN_SUCCESS) {
      rz_clear(&tmp);
      return err;
    }
    err = rz_copy(r, &tmp);
    rz_clear(&tmp);
    return err;
  }

  // early exit
  if (shift == 0) {
    return rz_copy(r, a);
  }

  u64 words = (u64)shift / 64;
  u64 bits = (u64)shift % 64;

  u64 max_size = a->size + words + 1;
  rabin_err_t err = rz_alloc(r, max_size);
  if (err != RABIN_SUCCESS) return err;

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

  rz_trim(r);
  return RABIN_SUCCESS;
}

rabin_err_t rz_rshift(rz_t* r, const rz_t* a, int shift)
{
  if (r == NULL || a == NULL) return RABIN_ERR_NULL_PTR;

  // aliasing
  if (r == a) {
    rz_t tmp;
    rz_init(&tmp);

    rabin_err_t err = rz_rshift(&tmp, a, shift);
    if (err != RABIN_SUCCESS) {
      rz_clear(&tmp);
      return err;
    }
    err = rz_copy(r, &tmp);
    rz_clear(&tmp);
    return err;
  }

  // early exit
  if (shift == 0) {
    return rz_copy(r, a);
  }

  u64 words = (u64)shift / 64;
  u64 bits = (u64)shift % 64;

  // we shift more than the number of limbs
  if (words >= a->size) {
    return rz_set_u64(r, 0);
  }

  u64 new_size = a->size - words;

  rabin_err_t err = rz_alloc(r, new_size);
  if (err != RABIN_SUCCESS) return err;

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
  rz_trim(r);

  // truncation toward zero keeps the sign of the dividend
  r->is_neg = a->is_neg && !rz_is_zero(r);
  return RABIN_SUCCESS;
}

rabin_err_t rz_lshift1_add(rz_t* r, int bit)
{
  if (r == NULL) return RABIN_ERR_NULL_PTR;

  u64 carry = (bit != 0);
  for (u64 i = 0; i < r->size; i++) {
    u64 tmp = r->limbs[i] >> 63;
    r->limbs[i] = (r->limbs[i] << 1) | carry;
    carry = tmp;
  }

  if (carry) {
    rabin_err_t err = rz_alloc(r, r->size + 1);
    if (err != RABIN_SUCCESS) return err;
    r->limbs[r->size++] = carry;
  }
  return RABIN_SUCCESS;
}
