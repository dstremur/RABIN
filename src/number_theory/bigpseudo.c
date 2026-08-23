/*
 * bigpseudo.c
 *
 * Pseudoprime generation.
 *
 * This file generates strong pseudoprimes to base 2: composite numbers
 * that pass the Miller-Rabin test for base 2. Useful for testing
 * primality-test implementations.
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

#include <fcntl.h>
#include <stdarg.h>
#include <stddef.h>
#include <unistd.h>

#include "../include/bignum.h"

/*
 * Generate a k-bit strong pseudoprime to base 2 into p.
 *
 * Repeatedly draws random k-bit candidates from /dev/urandom until one
 * is found that is composite (fails BPSW) yet passes the base-2
 * Miller-Rabin test (bn_rabin_mont). Such numbers are strong
 * pseudoprimes to base 2.
 *
 * Returns true on success, false if /dev/urandom cannot be opened.
 *
 * Complexity:
 *   Time: unbounded expected - each attempt costs one BPSW test plus
 *         one base-2 Miller-Rabin test; the hit rate depends on k
 *   Auxiliary memory: O(k/64) limbs
 *   Output memory: O(k/64) limbs
 */
bool bn_gen_strps(bignum* p, u64 k)
{
  bignum two;
  bn_init(&two);
  bn_set_u64(&two, 2);
  int fd = open("/dev/urandom", O_RDONLY);
  if (fd < 0) return false;

  do {
    bn_gen_random_with_fd(p, k, fd);

  } while (bn_bpsw(p) || !bn_rabin_mont(p, &two));

  close(fd);
  bn_free(&two);
  return true;
}
