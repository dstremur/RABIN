/*
 * rz_sqrt.c
 *
 * Integer square root routines.
 *
 * This file implements the integer (floor) square root of a rz_t using
 * Heron's method (Newton's method applied to x^2 - a = 0).
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

#include "../../include/rabin.h"

rabin_err_t rz_isqrt_heron(rz_t* r, const rz_t* a)
{
  if (r == NULL || a == NULL) return RABIN_ERR_NULL_PTR;
  if (a->is_neg) return RABIN_SUCCESS; /* no-op: r left unchanged */

  if (rz_is_zero(a)) {
    return rz_set_u64(r, 0);
  }

  rz_t xn, xnext, tmp;
  rz_init_multi(&xn, &xnext, &tmp, NULL);

  rabin_err_t err;

  // set x0 = 2^{log_2(a) / 2 + 1}
  u64 k = rz_bit_length(a);
  if ((err = rz_set_u64(&xn, 1)) != RABIN_SUCCESS) goto out;

  if ((err = rz_lshift(&xn, &xn, (k + 1) / 2 + 1)) != RABIN_SUCCESS) goto out;

  while (true) {
    if ((err = rz_div(&tmp, a, &xn)) != RABIN_SUCCESS) goto out;

    if ((err = rz_add(&xnext, &xn, &tmp)) != RABIN_SUCCESS) goto out;

    if ((err = rz_rshift1(&xnext)) != RABIN_SUCCESS) goto out;
    if (rz_cmp(&xnext, &xn) >= 0) {
      break;
    }

    if ((err = rz_swap(&xn, &xnext)) != RABIN_SUCCESS) goto out;
  }
  if ((err = rz_copy(r, &xn)) != RABIN_SUCCESS) goto out;

  err = RABIN_SUCCESS;
out:
  rz_clear_multi(&xn, &xnext, &tmp, NULL);
  return err;
}

rabin_err_t rz_isqrt(rz_t* r, const rz_t* a) { return rz_isqrt_heron(r, a); }
