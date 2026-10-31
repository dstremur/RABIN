/*
 * u64.c
 *
 * 64-bit modular arithmetic primitives.
 *
 * This file implements the single-limb (u64) building blocks used by
 * the NTT and RNS code paths: modular add/subtract/multiply, modular
 * exponentiation, modular inverses, Barrett reduction, and the
 * Goldilocks field (p = 2^64 - 2^32 + 1) reduction.
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

#include "../../include/u64.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../../include/rns.h"

typedef unsigned __int128 u128;

inline u64 u64_mod_add(u64 a, u64 b, u64 p)
{
  u64 res = a + b;

  u64 mask = -(u64)((res >= p) | (res < a));
  return res - (p & mask);
}

u64 u64_mod_sub(u64 a, u64 b, u64 p)
{
  u64 res = a - b;
  // If a < b, a borrow occurred.
  // (a < b) evaluates to 1 or 0. Negating it creates a mask of all 1s or all
  // 0s.
  u64 mask = -(u64)(a < b);
  return res + (p & mask);
}

inline u64 u64_mod_mul(u64 a, u64 b, u64 p)
{
  unsigned __int128 res = (unsigned __int128)a * b;
  return (u64)(res % p);
}

u64 u64_mod_pow(u64 base, u64 exp, u64 p)
{
  u64 res = 1;
  base %= p;
  while (exp > 0) {
    if (exp % 2 == 1) res = u64_mod_mul(res, base, p);
    base = u64_mod_mul(base, base, p);
    exp /= 2;
  }
  return res;
}

u64 u64_mod_inverse_euclid(u64 a, u64 p)
{
  if (a == 0) return 0;  // Should not happen with primes

  __int128 t = 0, newt = 1;
  __int128 r = p, newr = a;

  while (newr != 0) {
    __int128 q = r / newr;

    __int128 tmp_t = newt;
    newt = t - q * newt;
    t = tmp_t;

    __int128 tmp_r = newr;
    newr = r - q * newr;
    r = tmp_r;
  }

  if (t < 0) t += p;

  return (u64)t;
}

u64 u64_mod_inverse(u64 n, u64 p) { return u64_mod_pow(n, p - 2, p); }

u128 u64_compute_mu(u64 q)
{
  if (q == 0) return 0;
  return ~((u128)0) / q;
}

u64 u64_barrett(u128 c, u64 q, u128 mu)
{
  if (q == 0) return 0;

  u64 c_hi = (u64)(c >> 64), c_lo = (u64)c;
  u64 m_hi = (u64)(mu >> 64), m_lo = (u64)mu;

  // top 128 bits of the 256-bit product c * mu:
  //   c*mu = (c_hi*m_hi)*2^128 + (c_hi*m_lo + c_lo*m_hi)*2^64 + c_lo*m_lo
  u128 a = (u128)c_hi * m_lo;
  u128 mid = a + (u128)c_lo * m_hi;
  u64 mid_carry = mid < a;
  u128 lo_carry = ((u128)c_lo * m_lo) >> 64;
  mid += lo_carry;
  if (mid < lo_carry) mid_carry++;

  u128 q_est = (u128)c_hi * m_hi + (mid >> 64) + ((u128)mid_carry << 64);

  // q_est <= floor(c/q) <= q_est + 1, so one correction subtraction suffices
  u128 r = c - q_est * q;
  while (r >= (u128)q) {
    r -= q;
  }

  return (u64)r;
}

u64 u64_goldilock_red(u128 c)
{
  const u64 p = 0xFFFFFFFF00000001ULL;

  /* Fold once in the 128-bit domain: with c = y*2^64 + x and
     2^64 = 2^32 - 1 (mod p),  c = x + y*(2^32 - 1)  (mod p).
     The folded value is < 2^96 + 2^64, so at most 32+1 bits remain
     above the low word. */
  u64 x = (u64)c, y = (u64)(c >> 64);
  u128 t = (u128)x + (u128)y * 0xFFFFFFFFULL;
  u64 lo = (u64)t, hi = (u64)(t >> 64); /* hi <= 2^32 - 1 */

  /* Second fold in 64-bit arithmetic: hi*2^64 = hi*(2^32 - 1). */
  u64 s = lo + ((hi << 32) - hi);
  if (s < lo) s += 0xFFFFFFFFULL; /* fold the single 2^64 wrap */

  /* Now 0 <= s < 2^64 < 2p, so one conditional subtraction canonicalizes. */
  if (s >= p) s -= p;
  return s;
}

u64 u64_mul_goldilocks(u64 a, u64 b)
{
  const u64 p = 0xFFFFFFFF00000001ULL;

  /* c = lo + hi*2^64. Splitting hi = hl + hh*2^32 and using
     2^64 = 2^32 - 1, 2^96 = -1 (mod p) gives the exact identity
     c = lo - hh + hl*(2^32 - 1) (mod p), where every term fits in 64 bits
     and hl*(2^32 - 1) needs no further folding. q tracks the net 2^64
     overflow of the running sum (q in -1..1, kept two's complement).
     The whole tail is branchless: masked adds/subtracts only, so the
     transform butterflies never hit a mispredicted branch. */
  u128 x = (u128)a * b;
  u64 lo = (u64)x, hi = (u64)(x >> 64);
  u64 hl = hi & 0xFFFFFFFFULL, hh = hi >> 32;
  u64 hlf = (hl << 32) - hl; /* < 2^64 */

  u64 s = lo - hh;
  u64 q = 0 - (lo < hh);
  u64 s0 = s + hlf;
  q += (s0 < s);

  /* fold q*(2^32-1) = (q<<32)-q. A wrap with q = -1 is the intended
     mod-2^64 form of the debt (adding p); only q = +1 wraps carry a
     surplus 2^64 that needs one more (2^32-1). */
  s = s0 + ((q << 32) - q);
  u64 w = -(u64)(s < s0);
  u64 pos = ((q - 1) <= 0); /* q == +1 only */
  s += 0xFFFFFFFFULL & w & -(u64)pos;

  /* final representative is < 2p: one masked subtract canonicalizes */
  s -= p & -(u64)(s >= p);
  return s;
}
