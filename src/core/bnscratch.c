/*
 * bnscratch.c
 *
 * Thread-local scratch arena for hot-path bignum operations.
 *
 * Provides a per-thread, grow-only bump allocator so that
 * multiplication, division, and NTT routines can obtain large scratch
 * buffers without per-call malloc/free traffic.
 *
 * Usage discipline: a caller obtains at most one block per call
 * (bn_scratch_get) and releases it before returning
 * (bn_scratch_release). Scratch must not be held across calls into the
 * library, because a later get() may regrow the arena and invalidate
 * outstanding pointers.
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
#include <stdlib.h>

#include "../../include/bignum.h"

static __thread u64* tls_buf = NULL;
static __thread u64 tls_cap = 0; /* capacity in u64 units */
static __thread u64 tls_off = 0; /* current bump offset */

/**
 * @brief Reserve n u64s of thread-local scratch space.
 *
 * Returns a pointer to at least n u64s. The memory is NOT zeroed.
 * The block stays valid until bn_scratch_release() is called.
 *
 * Complexity:
 *   Time: O(1) amortized, O(1) realloc on arena growth
 *   Auxiliary memory: O(n) u64s (grow-only per thread)
 *   Output memory: O(n) u64s
 *
 * @param[in] n Number of u64s to reserve.
 *
 * @return A pointer to at least n u64s of scratch (not zeroed).
 */
u64* bn_scratch_get(u64 n)
{
  if (n == 0) n = 1;

  if (tls_off + n <= tls_cap) {
    u64* p = tls_buf + tls_off;
    tls_off += n;
    return p;
  }

  /* Grow the arena. Per the usage discipline a grow only happens with
   * no outstanding block (tls_off == 0), so rewinding is safe. */
  u64 new_cap = tls_cap ? tls_cap : 1024;
  while (new_cap < n) new_cap <<= 1;

  u64* new_buf = realloc(tls_buf, new_cap * sizeof(u64));
  if (new_buf) {
    tls_buf = new_buf;
    tls_cap = new_cap;
    tls_off = n;
    return new_buf;
  }

  /* Fallback for a failed realloc: one-off heap allocation */
  return malloc(n * sizeof(u64));
}

/**
 * @brief Release the current thread's scratch block(s) by rewinding the bump
 * pointer. The arena memory itself is kept for reuse.
 *
 * Complexity:
 *   Time: O(1)
 *   Auxiliary memory: O(1)
 *   Output memory: O(1)
 */
void bn_scratch_release(void) { tls_off = 0; }
