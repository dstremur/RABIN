/*
 * u64_mont.c
 *
 * Single-limb (u64) Montgomery multiplication.
 *
 * This file implements the Montgomery context and operations for a
 * 64-bit prime modulus p with R = 2^64: context initialization, the
 * REDC reduction, multiplication, domain conversion, and inversion.
 *
 * A value a_hat in the Montgomery domain represents a = a_hat * R^{-1}
 * mod p.
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

#include "../../include/bighelper.h"
#include "../../include/u64.h"

/*
 * Initialize a single-limb Montgomery context for the prime p.
 *
 * Computes:
 *
 *   ctx->p        = p
 *   ctx->p_inv    = -p^{-1} mod 2^64
 *   ctx->r2_mod_p = R^2 mod p, with R = 2^64
 *
 * Precondition: p is odd (in practice, prime).
 *
 * Complexity:
 *   Time: O(1)
 *   Auxiliary memory: O(1)
 *   Output memory: O(1)
 */
void mont_init(mont_ctx* ctx, u64 p)
{
  ctx->p = p;

  ctx->p_inv = mod_inverse_u64(p);

  // Calculate 2^128 mod p
  unsigned __int128 r2 = ((unsigned __int128)1 << 64) % p;
  r2 = (r2 * r2) % p;
  ctx->r2_mod_p = (u64)r2;
}

/*
 * Single-limb Montgomery reduction (REDC).
 *
 * Given T with 0 <= T < p * 2^64, this computes:
 *
 *   t = T * R^{-1} mod p,   R = 2^64
 *
 * by choosing m = (T mod 2^64) * (-p^{-1} mod 2^64) so that
 * T + m*p is divisible by 2^64, shifting right by 64 bits, and
 * applying one final conditional subtraction of p.
 *
 * Complexity:
 *   Time: O(1)
 *   Auxiliary memory: O(1)
 *   Output memory: O(1)
 */
u64 mont_redc(unsigned __int128 T, mont_ctx* ctx)
{
  u64 p = ctx->p;

  // T mod R
  u64 T_modR = (u64)T;

  u64 m = T_modR * ctx->p_inv;

  unsigned __int128 t_wide = (T + (unsigned __int128)m * p);
  u64 t = (u64)(t_wide >> 64);

  return (t >= p) ? (t - p) : t;
}

/*
 * Montgomery multiplication: (a * b * R^{-1}) mod p.
 *
 * If a and b are in the Montgomery domain, the result is the
 * Montgomery form of their represented values' product.
 *
 * Complexity:
 *   Time: O(1)
 *   Auxiliary memory: O(1)
 *   Output memory: O(1)
 */
u64 mont_mul(u64 a, u64 b, mont_ctx* ctx)
{
  unsigned __int128 R = (unsigned __int128)a * b;

  return mont_redc(R, ctx);
}

/*
 * Convert a value from the normal domain into the Montgomery domain:
 * a_hat = a * R mod p.
 *
 * Implemented as one Montgomery multiplication by R^2 mod p.
 *
 * Complexity:
 *   Time: O(1)
 *   Auxiliary memory: O(1)
 *   Output memory: O(1)
 */
u64 mont_in(u64 a, mont_ctx* ctx) { return mont_mul(a, ctx->r2_mod_p, ctx); }

/*
 * Convert a value from the Montgomery domain back to the normal
 * domain: a = a_hat * R^{-1} mod p.
 *
 * Implemented as one Montgomery multiplication by 1.
 *
 * Complexity:
 *   Time: O(1)
 *   Auxiliary memory: O(1)
 *   Output memory: O(1)
 */
u64 mont_out(u64 a_hat, mont_ctx* ctx) { return mont_mul(a_hat, 1, ctx); }

/*
 * Modular inverse of a value given in the Montgomery domain.
 *
 * Converts out of the Montgomery domain, inverts with the extended
 * Euclidean algorithm, and converts the result back in.
 *
 * Complexity:
 *   Time: O(log p)
 *   Auxiliary memory: O(1)
 *   Output memory: O(1)
 */
u64 mont_inverse(u64 a_mont, const mont_ctx* ctx)
{
  u64 a = mont_out(a_mont, ctx);

  u64 inv = mod_inverse_euclid(a, ctx->p);

  return mont_in(inv, ctx);
}
