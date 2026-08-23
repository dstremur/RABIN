/*
 * bigrand.c
 *
 * Random bignum generation.
 *
 * This file implements random bignum generation from /dev/urandom:
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

#include "../../include/bignum.h"

/**
 * @brief Generate a random bits-long odd number.
 *
 * Opens /dev/urandom, delegates to bn_gen_random_with_fd(), and
 * closes the descriptor.
 *
 * Returns true on success, false if /dev/urandom cannot be opened or
 * the read fails.
 *
 * Complexity:
 *   Time: O(bits / 64)
 *   Auxiliary memory: O(1)
 *   Output memory: O(bits / 64) limbs
 *
 * @param[out] r    Result storing the random odd number.
 * @param[in]  bits Desired bit length of the number.
 *
 * @return true  On success.
 * @return false If /dev/urandom cannot be opened or the read fails.
 */
bool bn_gen_random(bignum* r, u64 bits)
{
  int fd = open("/dev/urandom", O_RDONLY);
  if (fd < 0) return false;

  bn_gen_random_with_fd(r, bits, fd);
  close(fd);

  return true;
}

/**
 * @brief Generate a random bits-long odd number using a given urandom fd.
 *
 * Reads ceil(bits / 64) random limbs from fd, masks the top limb to
 * the exact bit length, sets the MSB so the number is exactly
 * bits long, and sets the LSB so it is odd.
 *
 * Returns true on success, false on allocation or read failure.
 *
 * Complexity:
 *   Time: O(bits / 64)
 *   Auxiliary memory: O(1)
 *   Output memory: O(bits / 64) limbs
 *
 * @param[out] r    Result storing the random odd number.
 * @param[in]  bits Desired bit length of the number.
 * @param[in]  fd   Open file descriptor for /dev/urandom.
 *
 * @return true  On success.
 * @return false On allocation or read failure.
 */
bool bn_gen_random_with_fd(bignum* r, u64 bits, int fd)
{
  int limbs_needed = (bits + 63) / 64;

  if (!bn_alloc(r, limbs_needed)) return false;
  r->size = limbs_needed;

  // Read random bytes directly using the open file descriptor
  if (read(fd, r->limbs, limbs_needed * sizeof(u64)) !=
      (ssize_t)(limbs_needed * sizeof(u64))) {
    return false;
  }

  // Mask the top limb to fit the exact bit length
  int top_bits = bits % 64;
  if (top_bits != 0) {
    u64 mask = ((u64)1 << top_bits) - 1;
    r->limbs[r->size - 1] &= mask;
  }

  // Ensure it's exactly 'bits' long by setting the MSB
  r->limbs[r->size - 1] |= ((u64)1 << ((bits - 1) % 64));

  // Ensure it's odd by setting the LSB
  r->limbs[0] |= 1;

  return true;
}

/**
 * @brief Generate a random number in the closed range [low, high].
 *
 * Let b = bit length of (high - low).
 *
 * Rejection sampling: draws random b-bit numbers until one is <=
 * (high - low), then adds low. If low == high, copies low directly.
 *
 * Complexity:
 *   Time: O(b / 64) per draw, expected O(b / 64) overall (geometric
 *         number of draws), plus O(b) for the range computation
 *   Auxiliary memory: O(b) limbs for the range
 *   Output memory: O(b) limbs
 *
 * @param[out] r    Result storing the random number in [low, high].
 * @param[in]  low  Lower bound (inclusive).
 * @param[in]  high Upper bound (inclusive).
 */
void bn_gen_random_range(bignum* r, const bignum* low, const bignum* high)
{
  bignum range;
  bn_init(&range);
  bn_sub(&range, high, low);

  if (bn_is_zero(&range)) {
    bn_copy(r, low);
    bn_free(&range);
    return;
  }

  int bits = bn_bit_length(&range);

  do {
    bn_gen_random(r, bits);
  } while (bn_cmp(r, &range) > 0);

  bn_add(r, r, low);
  bn_free(&range);
}
