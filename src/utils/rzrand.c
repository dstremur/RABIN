/*
 * rzrand.c
 *
 * Random rz_t generation.
 *
 * This file implements random rz_t generation from /dev/urandom:
 * random odd numbers of an exact bit length, and random numbers in a
 * closed interval [low, high].
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
#include <fcntl.h>
#include <math.h>
#include <stdarg.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "../../include/rabin.h"
#include "../core/rz_internal.h"

rabin_err_t rz_gen_random(rz_t* r, u64 bits)
{
  if (r == NULL) return RABIN_ERR_NULL_PTR;

  int fd = open("/dev/urandom", O_RDONLY);
  if (fd < 0) return RABIN_ERR_INVALID_ARG;

  rabin_err_t err = rz_gen_random_with_fd(r, bits, fd);
  close(fd);

  return err;
}

rabin_err_t rz_gen_random_odd_with_fd(rz_t* r, u64 bits, int fd)
{
  if (r == NULL) return RABIN_ERR_NULL_PTR;
  if (bits == 0) return RABIN_ERR_INVALID_ARG;

  u64 limbs_needed = (bits + 63) / 64;

  rabin_err_t err = rz_alloc(r, limbs_needed);
  if (err != RABIN_SUCCESS) return err;
  r->size = limbs_needed;

  // Read random bytes directly using the open file descriptor
  if (read(fd, r->limbs, limbs_needed * sizeof(u64)) !=
      (ssize_t)(limbs_needed * sizeof(u64))) {
    return RABIN_ERR_INVALID_ARG;
  }

  // Mask the top limb to fit the exact bit length
  u64 top_bits = bits % 64;
  if (top_bits != 0) {
    u64 mask = ((u64)1 << top_bits) - 1;
    r->limbs[r->size - 1] &= mask;
  }

  // Ensure it's exactly 'bits' long by setting the MSB
  r->limbs[r->size - 1] |= ((u64)1 << ((bits - 1) % 64));

  // Ensure it's odd by setting the LSB
  r->limbs[0] |= 1;

  rz_trim(r);
  return RABIN_SUCCESS;
}

rabin_err_t rz_gen_random_with_fd(rz_t* r, u64 bits, int fd)
{
  if (r == NULL) return RABIN_ERR_NULL_PTR;
  if (bits == 0) {
    return rz_set_u64(r, 0);
  }

  u64 limbs_needed = (bits + 63) / 64;

  rabin_err_t err = rz_alloc(r, limbs_needed);
  if (err != RABIN_SUCCESS) return err;
  r->size = limbs_needed;

  // Read random bytes directly using the open file descriptor
  if (read(fd, r->limbs, limbs_needed * sizeof(u64)) !=
      (ssize_t)(limbs_needed * sizeof(u64))) {
    return RABIN_ERR_INVALID_ARG;
  }

  // Mask the top limb to fit the exact bit length
  u64 top_bits = bits % 64;
  if (top_bits != 0) {
    u64 mask = ((u64)1 << top_bits) - 1;
    r->limbs[r->size - 1] &= mask;
  }

  rz_trim(r);
  return RABIN_SUCCESS;
}

rabin_err_t rz_gen_random_range(rz_t* r, const rz_t* low, const rz_t* high)
{
  if (r == NULL || low == NULL || high == NULL) return RABIN_ERR_NULL_PTR;
  if (rz_cmp(low, high) > 0) return RABIN_ERR_INVALID_ARG;

  rz_t range;
  rz_init(&range);
  rabin_err_t err = rz_sub(&range, high, low);
  if (err != RABIN_SUCCESS) {
    rz_clear(&range);
    return err;
  }

  if (rz_is_zero(&range)) {
    err = rz_copy(r, low);
    rz_clear(&range);
    return err;
  }

  u64 bits = (u64)rz_bit_length(&range);

  do {
    if ((err = rz_gen_random(r, bits)) != RABIN_SUCCESS) break;
  } while (rz_cmp(r, &range) > 0);

  if (err == RABIN_SUCCESS) {
    err = rz_add(r, r, low);
  }

  rz_clear(&range);
  return err;
}
