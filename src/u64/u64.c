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

#include "../include/u64.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../../include/bigrns.h"

typedef unsigned __int128 u128;

/**
 * @brief Modular addition: (a + b) mod p, with 0 <= a, b < p.
 *
 * Adds in u64 and conditionally subtracts p, using a branch-free
 * overflow/carry mask.
 *
 * Complexity:
 *   Time: O(1)
 *   Auxiliary memory: O(1)
 *   Output memory: O(1)
 *
 * @param[in] a First operand.
 * @param[in] b Second operand.
 * @param[in] p Modulus.
 *
 * @return (a + b) mod p.
 */
inline u64 mod_add(u64 a, u64 b, u64 p)
{
  u64 res = a + b;

  u64 mask = -(u64)((res >= p) | (res < a));
  return res - (p & mask);
}

/**
 * @brief Modular subtraction: (a - b) mod p, with 0 <= a, b < p.
 *
 * Subtracts in u64 and conditionally adds p when a borrow occurred,
 * using a branch-free mask.
 *
 * Complexity:
 *   Time: O(1)
 *   Auxiliary memory: O(1)
 *   Output memory: O(1)
 *
 * @param[in] a First operand.
 * @param[in] b Second operand.
 * @param[in] p Modulus.
 *
 * @return (a - b) mod p.
 */
u64 mod_sub(u64 a, u64 b, u64 p)
{
  u64 res = a - b;
  // If a < b, a borrow occurred.
  // (a < b) evaluates to 1 or 0. Negating it creates a mask of all 1s or all
  // 0s.
  u64 mask = -(u64)(a < b);
  return res + (p & mask);
}

/**
 * @brief Modular multiplication: (a * b) mod p, with 0 <= a, b < p.
 *
 * Multiplies in 128 bits and reduces with a hardware 128/64 division.
 *
 * Complexity:
 *   Time: O(1)
 *   Auxiliary memory: O(1)
 *   Output memory: O(1)
 *
 * @param[in] a First operand.
 * @param[in] b Second operand.
 * @param[in] p Modulus.
 *
 * @return (a * b) mod p.
 */
inline u64 mod_mul(u64 a, u64 b, u64 p)
{
  unsigned __int128 res = (unsigned __int128)a * b;
  return (u64)(res % p);
}

/**
 * @brief Modular exponentiation: base^exp mod p.
 *
 * Right-to-left binary exponentiation (square-and-multiply).
 *
 * Complexity:
 *   Time: O(log exp) modular multiplications
 *   Auxiliary memory: O(1)
 *   Output memory: O(1)
 *
 * @param[in] base Base.
 * @param[in]  exp Exponent.
 * @param[in]    p Modulus.
 *
 * @return base^exp mod p.
 */
u64 mod_pow(u64 base, u64 exp, u64 p)
{
  u64 res = 1;
  base %= p;
  while (exp > 0) {
    if (exp % 2 == 1) res = mod_mul(res, base, p);
    base = mod_mul(base, base, p);
    exp /= 2;
  }
  return res;
}

/**
 * @brief Modular inverse via the extended Euclidean algorithm:
 * a^(-1) mod p.
 *
 * p needs to be prime (more generally, gcd(a, p) must be 1). Returns
 * 0 if a is 0 (should not happen with primes).
 *
 * Complexity:
 *   Time: O(log p)
 *   Auxiliary memory: O(1)
 *   Output memory: O(1)
 *
 * @param[in] a Value to invert.
 * @param[in] p Modulus (prime).
 *
 * @return a^(-1) mod p, or 0 if a is 0.
 */
u64 mod_inverse_euclid(u64 a, u64 p)
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

/**
 * @brief Modular inverse via Fermat's little theorem: a^(-1) = a^(p-2) mod p.
 *
 * p needs to be prime.
 *
 * Complexity:
 *   Time: O(log p) modular multiplications
 *   Auxiliary memory: O(1)
 *   Output memory: O(1)
 *
 * @param[in] n Value to invert.
 * @param[in] p Modulus (prime).
 *
 * @return n^(-1) mod p.
 */
u64 mod_inverse(u64 n, u64 p) { return mod_pow(n, p - 2, p); }

/**
 * @brief Compute the Barrett reduction constant mu = floor((2^128 - 1) / q).
 *
 * Complexity:
 *   Time: O(1)
 *   Auxiliary memory: O(1)
 *   Output memory: O(1)
 *
 * @param[in] q Modulus.
 *
 * @return The Barrett reduction constant mu.
 */
u64 compute_mu(u64 q)
{
  u128 dividend = ~((u128)0);
  return (u64)(dividend / q);
}

/**
 * @brief Barrett reduction: c mod q, where mu = floor((2^128 - 1) / q).
 *
 * Estimates the quotient as q_est = floor(c * mu / 2^128) using the
 * high 128 bits of the 192-bit product, computes r = c - q_est * q,
 * and subtracts q at most twice to correct the (rare) overestimate.
 *
 * Complexity:
 *   Time: O(1)
 *   Auxiliary memory: O(1)
 *   Output memory: O(1)
 *
 * @param[in]  c  Value to reduce (128-bit).
 * @param[in]  q  Modulus.
 * @param[in] mu Barrett reduction constant.
 *
 * @return c mod q.
 */
u64 barrett_reduction(u128 c, u64 q, u64 mu)
{
  unsigned __int128 q_est = (unsigned __int128)((c * mu) >> 128);

  u64 r = (u64)(c - q_est * q);

  while (r >= q) {
    r -= q;
  }

  return r;
}

/**
 * @brief Goldilocks field reduction: c mod p with p = 2^64 - 2^32 + 1.
 *
 * Splits the 128-bit value c into 32-bit chunks and uses the field
 * relation 2^32 = 1 (mod p) to fold the upper chunks down, followed by
 * a single conditional subtraction.
 *
 * Complexity:
 *   Time: O(1)
 *   Auxiliary memory: O(1)
 *   Output memory: O(1)
 *
 * @param[in] c Value to reduce (128-bit).
 *
 * @return c mod p, where p = 2^64 - 2^32 + 1.
 */
u64 goldilock_red(u128 c)
{
  u128 X_3 = c >> 96;
  u128 X_2 = (c >> 64) & ((1ULL << 32) - 1);
  u128 X_1 = c & ((1ULL >> 64) - 1);
  u128 C_out = X_1 + (X_2 * (1ULL << 32) - 1) - X_3;
  if (C_out >= 18446744069414584321ULL) {
    C_out -= 18446744069414584321ULL;
  }

  return (u64)C_out;
}
