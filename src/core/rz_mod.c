/*
 * rz_mod.c
 *
 * rz_t modulo routines.
 *
 * This file implements reduction of a rz_t modulo a 64-bit divisor
 * (with and without quotient output) and the modular multiplicative
 * inverse via the binary extended GCD.
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

#include "../../include/rz.h"
#include "../../include/rzlimb.h"
#include "../../include/rzmath.h"
#include "rz_internal.h"

uint64_t rz_mod_u64(const rz_t* a, uint64_t d)
{
  unsigned __int128 rem = 0;

  // Purely mathematical reduction with zero memory allocation
  for (i64 i = a->size - 1; i >= 0; i--) {
    unsigned __int128 cur = (rem << 64) | a->limbs[i];
    rem = cur % d;
  }

  u64 res = (u64)rem;

  if (a->is_neg && res != 0) {
    return d - res;
  }

  return res;
}

uint64_t rz_divmod_u64(rz_t* q, const rz_t* a, uint64_t d)
{
  // aliasing
  if (q == a) {
    rz_t tmp;
    rz_init(&tmp);
    uint64_t rem = rz_divmod_u64(&tmp, a, d);
    rz_copy(q, &tmp);
    rz_clear(&tmp);
    return rem;
  }
  rz_alloc(q, a->size);
  q->size = a->size;
  q->is_neg = a->is_neg;

  unsigned __int128 rem = 0;

  for (i64 i = a->size - 1; i >= 0; i--) {
    unsigned __int128 cur = (rem << 64) | a->limbs[i];

    q->limbs[i] = (u64)(cur / d);
    rem = cur % d;
  }

  rz_trim(q);
  return (u64)rem;
}

rabin_err_t rz_mod_inverse(rz_t* res, const rz_t* a, const rz_t* m)
{
  if (res == NULL || a == NULL || m == NULL) return RABIN_ERR_NULL_PTR;

  rz_t u, v, d, one;
  if (rz_init_multi(&u, &v, &d, &one, NULL) != RABIN_SUCCESS)
    return RABIN_ERR_OUT_OF_MEMORY;

  rabin_err_t err = rz_set_u64(&one, 1);
  if (err == RABIN_SUCCESS) {
    err = rz_gcd_extended_lehmer(&u, &v, &d, a, m);
  }

  if (err == RABIN_SUCCESS && rz_cmp(&d, &one) != 0) {
    // gcd(a, m) != 1: no modular inverse exists
    rabin_err_t e = rz_set_u64(res, 0);
    err = (e == RABIN_SUCCESS) ? RABIN_ERR_INVALID_ARG : e;
  }

  if (err == RABIN_SUCCESS) {
    err = rz_mod(res, &u, m);
    if (err == RABIN_SUCCESS && res->is_neg) {
      err = rz_add(res, res, m);
    }
  }

  rz_clear_multi(&u, &v, &d, &one, NULL);

  return err;
}
