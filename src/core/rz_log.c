/*
 * rz_log.c
 *
 * Integer logarithm routines.
 *
 * This file implements integer logarithms for bignums: the base-2
 * logarithm (exact, via the bit length) and a fixed-point approximation
 * of the natural logarithm.
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

rabin_err_t rz_ln(rz_t* r, rz_t* a)
{
  if (r == NULL || a == NULL) return RABIN_ERR_NULL_PTR;
  if (rz_is_zero(a)) return RABIN_SUCCESS; /* no-op: r left unchanged */

  // ln(2) = 0.6931471805599453094172321214

  rz_t log2, tmp, tmp2;
  rz_init_multi(&log2, &tmp, &tmp2, NULL);

  rabin_err_t err = rz_log_2(&log2, a);
  if (err == RABIN_SUCCESS) err = rz_set_u64(&tmp, 69314718);
  if (err == RABIN_SUCCESS) err = rz_set_u64(&tmp2, 100000000);

  // r = (k * 69314718) / 100000000
  if (err == RABIN_SUCCESS) err = rz_mul(r, &log2, &tmp);
  if (err == RABIN_SUCCESS) err = rz_div(r, r, &tmp2);

  rz_clear_multi(&log2, &tmp, &tmp2, NULL);
  return err;
}

rabin_err_t rz_log_2(rz_t* r, rz_t* a)
{
  if (r == NULL || a == NULL) return RABIN_ERR_NULL_PTR;
  if (rz_is_zero(a)) return RABIN_SUCCESS; /* no-op: r left unchanged */

  u64 k = rz_bit_length(a) - 1;

  return rz_set_u64(r, k);
}
