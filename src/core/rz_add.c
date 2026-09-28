/*
 * rz_add.c
 *
 * Integer addition routines.
 *
 * This file implements magnitude addition, signed addition, and helper
 * addition operations for the RABIN integer library.
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

#include "../../include/rabin.h"
#include "../../include/rzlimb.h"
#include "rz_internal.h"

rabin_err_t rz_add_abs(rz_t* r, const rz_t* a, const rz_t* b)
{
  if (r == NULL || a == NULL || b == NULL) return RABIN_ERR_NULL_PTR;

  // aliasing
  if (r == a || r == b) {
    rz_t tmp;
    rz_init(&tmp);
    rabin_err_t err = rz_add_abs(&tmp, a, b);
    if (err != RABIN_SUCCESS) {
      rz_clear(&tmp);
      return err;
    }
    err = rz_copy(r, &tmp);
    rz_clear(&tmp);
    return err;
  }

  // ensure a is always larger
  if (a->size < b->size) {
    const rz_t* tmp = a;
    a = b;
    b = tmp;
  }

  if (r->capacity < a->size + 1) {
    rabin_err_t err = rz_alloc(r, a->size + 1);
    if (err != RABIN_SUCCESS) return err;
  }

  u64 carry = rz_add_inner(r->limbs, a->limbs, a->size, b->limbs, b->size);

  r->limbs[a->size] = carry;
  r->size = a->size + carry;
  return RABIN_SUCCESS;
}

rabin_err_t rz_add(rz_t* r, const rz_t* a, const rz_t* b)
{
  if (r == NULL || a == NULL || b == NULL) return RABIN_ERR_NULL_PTR;

  if (r == a || r == b) {
    rz_t tmp;
    rz_init(&tmp);
    rabin_err_t err = rz_add(&tmp, a, b);
    if (err != RABIN_SUCCESS) {
      rz_clear(&tmp);
      return err;
    }
    err = rz_copy(r, &tmp);
    rz_clear(&tmp);
    return err;
  }

  if (a->is_neg == b->is_neg) {
    rabin_err_t err = rz_add_abs(r, a, b);
    if (err != RABIN_SUCCESS) return err;
    r->is_neg = a->is_neg;
    return RABIN_SUCCESS;
  }

  i64 cmp = rz_cmp_abs(a, b);

  if (cmp == 0) {
    rabin_err_t err = rz_set_u64(r, 0);
    if (err != RABIN_SUCCESS) return err;
    r->is_neg = false;
    return RABIN_SUCCESS;
  }

  rabin_err_t err;
  if (cmp > 0) {
    err = rz_sub_abs(r, a, b);
    if (err != RABIN_SUCCESS) return err;
    r->is_neg = a->is_neg;
  } else {
    err = rz_sub_abs(r, b, a);
    if (err != RABIN_SUCCESS) return err;
    r->is_neg = b->is_neg;
  }
  return RABIN_SUCCESS;
}

rabin_err_t rz_add_u64(rz_t* r, const rz_t* a, u64 b)
{
  if (r == NULL || a == NULL) return RABIN_ERR_NULL_PTR;

  if (b == 0) {
    return rz_copy(r, a);
  }

  if (!a->is_neg) {
    // POSITIVE + POSITIVE
    rabin_err_t err = rz_copy(r, a);
    if (err != RABIN_SUCCESS) return err;
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
      if (r->capacity < r->size + 1) {
        err = rz_alloc(r, r->size + 1);
        if (err != RABIN_SUCCESS) return err;
      }
      r->limbs[r->size++] = 1;
    }
    return RABIN_SUCCESS;
  } else {
    // NEGATIVE + POSITIVE (e.g., -80 + 16)
    // This is Magnitude Subtraction: |a| - b
    if (a->size == 1 && a->limbs[0] <= b) {
      return rz_set_u64(r, b - a->limbs[0]);
    } else {
      rabin_err_t err = rz_copy(r, a);
      if (err != RABIN_SUCCESS) return err;
      u64 borrow = b;
      for (u64 i = 0; i < r->size && borrow; i++) {
        u64 old = r->limbs[i];
        r->limbs[i] -= borrow;
        borrow = (old < borrow) ? 1 : 0;
      }
      r->is_neg = true;
      rz_trim(r);
      if (r->size == 0) r->is_neg = false;
      return RABIN_SUCCESS;
    }
  }
}

rabin_err_t rz_add_at_offset(rz_t* r, const rz_t* a, u64 offset)
{
  if (r == NULL || a == NULL) return RABIN_ERR_NULL_PTR;
  if (a->size == 0) return RABIN_SUCCESS;

  u64 carry = 0;
  // Ensure r has enough limbs to hold a + carry at this offset
  if (offset > UINT64_MAX - a->size - 1) return RABIN_ERR_OVERFLOW;
  u64 required_size = offset + a->size + 1;
  if (r->size < required_size) {
    rabin_err_t err = rz_alloc(r, required_size);
    if (err != RABIN_SUCCESS) return err;
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
  rz_trim(r);
  return RABIN_SUCCESS;
}

rabin_err_t rz_neg(rz_t* r, const rz_t* a)
{
  if (r == NULL || a == NULL) {
    return RABIN_ERR_NULL_PTR;
  }

  // Copy if destination is a different object
  if (r != a) {
    rabin_err_t err = rz_copy(r, a);
    if (err != RABIN_SUCCESS) return err;
  }

  // Zero must always remain positive (is_neg = false)
  if (rz_is_zero(r)) {
    r->is_neg = false;
  } else {
    r->is_neg = !r->is_neg;
  }
  return RABIN_SUCCESS;
}
