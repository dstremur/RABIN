/*
 * rz_sub.c
 *
 * rz_t subtraction routines.
 *
 * This file implements magnitude subtraction and signed subtraction for
 * the rz_t library.
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

rabin_err_t rz_sub_abs(rz_t* r, const rz_t* a, const rz_t* b)
{
  if (r == NULL || a == NULL || b == NULL) return RABIN_ERR_NULL_PTR;

  // aliasing
  if (r == a || r == b) {
    rz_t tmp;
    rz_init(&tmp);
    rabin_err_t err = rz_sub_abs(&tmp, a, b);
    if (err != RABIN_SUCCESS) {
      rz_clear(&tmp);
      return err;
    }
    err = rz_copy(r, &tmp);
    rz_clear(&tmp);
    return err;
  }

  rabin_err_t err = rz_alloc(r, a->size);
  if (err != RABIN_SUCCESS) return err;

  u64 borrow = 0;

  for (u64 i = 0; i < a->size; i++) {
    u64 av = a->limbs[i];
    u64 bv = (i < b->size) ? b->limbs[i] : 0;

    unsigned __int128 diff = (unsigned __int128)av - bv - borrow;

    r->limbs[i] = (u64)diff;
    borrow = (diff >> 127) & 1;  // detect underflow
  }

  r->size = a->size;
  rz_trim(r);
  return RABIN_SUCCESS;
}

rabin_err_t rz_sub(rz_t* r, const rz_t* a, const rz_t* b)
{
  if (r == NULL || a == NULL || b == NULL) return RABIN_ERR_NULL_PTR;

  if (r == a || r == b) {
    rz_t tmp;
    rz_init(&tmp);
    rabin_err_t err = rz_sub(&tmp, a, b);
    if (err != RABIN_SUCCESS) {
      rz_clear(&tmp);
      return err;
    }
    err = rz_copy(r, &tmp);
    rz_clear(&tmp);
    return err;
  }

  // a - (-b) = a + b
  if (!a->is_neg && b->is_neg) {
    rabin_err_t err = rz_add_abs(r, a, b);
    if (err != RABIN_SUCCESS) return err;
    r->is_neg = false;
    return RABIN_SUCCESS;
  }

  // (-a) - b = -(a + b)
  if (a->is_neg && !b->is_neg) {
    rabin_err_t err = rz_add_abs(r, a, b);
    if (err != RABIN_SUCCESS) return err;
    r->is_neg = true;
    return RABIN_SUCCESS;
  }

  i64 cmp = rz_cmp_abs(a, b);

  if (cmp == 0) {
    return rz_set_u64(r, 0);
  }

  rabin_err_t err;
  if (cmp > 0) {
    // |a| > |b|
    err = rz_sub_abs(r, a, b);
    if (err != RABIN_SUCCESS) return err;
    r->is_neg = a->is_neg;
  } else {
    // |b| > |a|
    err = rz_sub_abs(r, b, a);
    if (err != RABIN_SUCCESS) return err;
    r->is_neg = !a->is_neg;
  }
  return RABIN_SUCCESS;
}
