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

void bn_sub(bignum* r, const bignum* a, const bignum* b)
{
  if (r == a || r == b) {
    bignum tmp;
    bn_init(&tmp);
    bn_sub(&tmp, a, b);
    bn_copy(r, &tmp);
    bn_free(&tmp);
    return;
  }

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
