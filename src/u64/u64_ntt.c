/*
 * u64_ntt.c
 *
 * Single-limb (u64) Number Theoretic Transform.
 *
 * This file implements a cyclic (Cooley-Tukey) NTT and its inverse for
 * u64 arrays over a prime field, with all butterflies performed in the
 * Montgomery domain for speed. It also provides context initialization
 * for a generic prime (with caller-supplied roots of unity) and for
 * the Goldilocks field p = 5 * 2^55 + 1.
 *
 * The transform length is n = 2^k for k <= 54.
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

#include <pthread.h>

#include "../../include/u64.h"

/*
 * Forward cyclic NTT: a_hat = NTT(a).
 *
 * Let n = ctx->n = 2^k.
 *
 * Iterative Cooley-Tukey (decimation in time):
 *
 *   1. bit-reversal permutation of the input, converting each entry
 *      into the Montgomery domain on the fly,
 *   2. log2(n) stages of butterflies, each using precomputed twiddle
 *      factors (omega powers, already in Montgomery form) and
 *      Montgomery multiplication.
 *
 * The output a_hat is in the Montgomery domain; pair it with
 * ntt_u64_cyclic_inverse() (or the montgomery_in variant) to recover
 * standard integers.
 *
 * a_hat must not alias a.
 *
 * Complexity:
 *   Time: O(n log n)
 *   Auxiliary memory: O(1)
 *   Output memory: O(n)
 */
void ntt_u64_cyclic_forward(u64* a_hat, const u64* a, ntt_ctx_u64* ctx)
{
  u64 n = ctx->n;
  u64 q = ctx->q;

  // 1. Bit-reversal permutation AND entering Montgomery space

  for (u64 i = 0; i < n; i++) {
    u64 rev = ctx->bit_rev_indices[i];
    // Convert to Montgomery form right as we load the data
    a_hat[rev] = mont_in(a[i], &ctx->mctx);
  }

  // 2. Cooley-Tukey Butterfly
  for (u64 len = 2; len <= n; len <<= 1) {
    u64 half = len >> 1;
    u64 step = n / len;

    for (u64 i = 0; i < n; i += len) {
      for (u64 j = 0; j < half; j++) {
        u64 twiddle =
            ctx->omega_powers[j * step];  // Already in Montgomery form

        u64 even = i + j;
        u64 odd = i + j + half;

        // t = a_hat[odd] * twiddle (Montgomery multiplication)
        u64 t = mont_mul(a_hat[odd], twiddle, &ctx->mctx);
        u64 u = a_hat[even];

        // Montgomery form preserves addition and subtraction natively
        a_hat[even] = mod_add(u, t, q);
        a_hat[odd] = mod_sub(u, t, q);
      }
    }
  }
}

/*
 * Initialize an NTT context for the Goldilocks field p = 5 * 2^55 + 1.
 *
 * Let k = log2 of the transform length (n = 2^k, k <= 54).
 *
 * Derives the roots of unity from the primitive root g = 3:
 *
 *   psi   = g^(5 * 2^(54-k)) mod p   (primitive (k+1)-th root)
 *   omega = psi^2 mod p              (primitive k-th root)
 *
 * and delegates to ntt_ctx_u64_init().
 *
 * Returns false (and leaves ctx unchanged) if k > 54.
 *
 * Complexity:
 *   Time: O(n) for the precomputed tables, plus O(log p) for the root
 *         derivation
 *   Auxiliary memory: O(1)
 *   Output memory: O(n) per precomputed table (3 tables)
 */
bool ntt_ctx_u64_init_golden(ntt_ctx_u64* ctx, u64 k)
{
  if (k > 54) return false;

  // p = 5 * 2^55 + 1
  u64 p = 180143985094819841ULL;
  u64 g = 3;

  // Calculate the exponent for psi: c * 2^(55 - (k+1))
  // We use 5ULL to ensure it shifts as a 64-bit integer
  u64 exp = 5ULL << (55 - (k + 1));

  // psi = g^exp mod p
  u64 psi = mod_pow(g, exp, p);

  // omega = psi^2 mod p
  u64 omega = mod_mul(psi, psi, p);

  // Initialize the Montgomery NTT context
  return ntt_ctx_u64_init(ctx, p, k, omega, psi);
}

/*
 * Cached Goldilocks NTT contexts, one per transform size k (0..54).
 *
 * The omega power tables and bit-reversal indices are pure functions
 * of k, so each context is built lazily on first use and reused
 * forever after. Contexts are read-only after initialization, so they
 * can be shared across threads.
 */
static ntt_ctx_u64 golden_ctxs[55];
static bool golden_ctx_ready[55] = {false};
static pthread_mutex_t golden_ctx_lock = PTHREAD_MUTEX_INITIALIZER;

/*
 * Return the shared Goldilocks NTT context for transform length 2^k
 * (k <= 54), initializing it on first use.
 *
 * Returns NULL if k > 54 or initialization fails. The returned context
 * is owned by the library and must not be freed or modified by the
 * caller.
 *
 * Complexity:
 *   Time: O(1) after the first call for a given k, O(2^k) the first
 *         time
 *   Auxiliary memory: O(1)
 *   Output memory: O(2^k) u64s per table, once per k
 */
ntt_ctx_u64* ntt_ctx_u64_golden_cached(u64 k)
{
  if (k > 54) return NULL;

  pthread_mutex_lock(&golden_ctx_lock);
  if (!golden_ctx_ready[k]) {
    if (!ntt_ctx_u64_init_golden(&golden_ctxs[k], k)) {
      pthread_mutex_unlock(&golden_ctx_lock);
      return NULL;
    }
    golden_ctx_ready[k] = true;
  }
  pthread_mutex_unlock(&golden_ctx_lock);

  return &golden_ctxs[k];
}

/*
 * Initialize a u64 NTT context for a prime p and transform length
 * n = 2^k.
 *
 * Precomputes:
 *
 *   - the Montgomery context for p,
 *   - n_inv = n^{-1} mod p in Montgomery form (for the inverse
 *     transform scaling),
 *   - the omega and omega^{-1} power tables (in Montgomery form),
 *   - the bit-reversal index table.
 *
 * The caller must supply omega, a primitive k-th root of unity mod p
 * (psi is accepted for interface compatibility but only omega is used
 * for the tables).
 *
 * Complexity:
 *   Time: O(n)
 *   Auxiliary memory: O(1)
 *   Output memory: O(n) per precomputed table (3 tables)
 */
bool ntt_ctx_u64_init(ntt_ctx_u64* ctx, u64 p, u64 k, u64 omega, u64 psi)
{
  u64 n = (1ULL << k);
  ctx->n = n;
  ctx->k = k;
  ctx->q = p;

  mont_init(&ctx->mctx, p);

  ctx->omega_powers = malloc(sizeof(u64) * n);
  ctx->omega_inv_powers = malloc(sizeof(u64) * n);
  ctx->bit_rev_indices = malloc(sizeof(u64) * n);

  // Calculate n^-1 mod q and put it in Montgomery form
  u64 n_inv_standard = mod_inverse_euclid(n, p);
  ctx->n_inv = mont_in(n_inv_standard, &ctx->mctx);

  // Precompute powers and move to Montgomery form immediately
  u64 current_omega = 1;
  u64 omega_inv = mod_inverse_euclid(omega, p);
  u64 current_omega_inv = 1;

  for (u64 i = 0; i < n; i++) {
    ctx->omega_powers[i] = mont_in(current_omega, &ctx->mctx);
    ctx->omega_inv_powers[i] = mont_in(current_omega_inv, &ctx->mctx);

    current_omega =
        mod_mul(current_omega, omega, p);  // Standard mod_mul for precalc
    current_omega_inv = mod_mul(current_omega_inv, omega_inv, p);
  }

  // Precompute Bit-Reversal
  u64 bits = 0;
  while (((u64)1 << bits) < n) bits++;
  for (u64 i = 0; i < n; i++) {
    u64 rev = 0, temp = i;
    for (u64 j = 0; j < bits; j++) {
      rev = (rev << 1) | (temp & 1);
      temp >>= 1;
    }
    ctx->bit_rev_indices[i] = rev;
  }

  return true;
}

/*
 * Free the precomputed tables of a u64 NTT context.
 *
 * Complexity:
 *   Time: O(1)
 *   Auxiliary memory: O(1)
 *   Output memory: O(1)
 */
void ntt_ctx_u64_free(ntt_ctx_u64* ctx)
{
  if (ctx->omega_powers) free(ctx->omega_powers);
  if (ctx->omega_inv_powers) free(ctx->omega_inv_powers);
  if (ctx->bit_rev_indices) free(ctx->bit_rev_indices);
}

/*
 * Inverse cyclic NTT: a_hat = NTT^{-1}(a).
 *
 * Let n = ctx->n = 2^k.
 *
 * Same structure as the forward transform, but with the inverse
 * twiddle factors, followed by scaling with n^{-1} mod p and
 * conversion out of the Montgomery domain. The output a_hat contains
 * standard (non-Montgomery) integers.
 *
 * a_hat must not alias a.
 *
 * Complexity:
 *   Time: O(n log n)
 *   Auxiliary memory: O(1)
 *   Output memory: O(n)
 */
void ntt_u64_cyclic_inverse(u64* a_hat, const u64* a, ntt_ctx_u64* ctx)
{
  u64 n = ctx->n;
  u64 q = ctx->q;

  // 1. Bit-reversal permutation AND entering Montgomery space
  for (u64 i = 0; i < n; i++) {
    u64 rev = ctx->bit_rev_indices[i];
    a_hat[rev] = mont_in(a[i], &ctx->mctx);
  }

  // 2. Cooley-Tukey Butterfly
  for (u64 len = 2; len <= n; len <<= 1) {
    u64 half = len >> 1;
    u64 step = n / len;

    for (u64 i = 0; i < n; i += len) {
      for (u64 j = 0; j < half; j++) {
        u64 twiddle =
            ctx->omega_inv_powers[j * step];  // Already in Montgomery form

        u64 even = i + j;
        u64 odd = i + j + half;

        u64 t = mont_mul(a_hat[odd], twiddle, &ctx->mctx);
        u64 u = a_hat[even];

        a_hat[even] = mod_add(u, t, q);
        a_hat[odd] = mod_sub(u, t, q);
      }
    }
  }

  // 3. Scale by n^-1 mod q AND exit Montgomery space
  for (u64 i = 0; i < n; i++) {
    // Multiply by n_inv (which is in Montgomery form)
    a_hat[i] = mont_mul(a_hat[i], ctx->n_inv, &ctx->mctx);

    // Exit Montgomery space to get the final standard integer
    a_hat[i] = mont_out(a_hat[i], &ctx->mctx);
  }
}

/*
 * Inverse cyclic NTT with Montgomery-domain input:
 * a_hat = NTT^{-1}(a), where a is already in the Montgomery domain.
 *
 * Let n = ctx->n = 2^k.
 *
 * Identical to ntt_u64_cyclic_inverse() except that step 1 performs
 * only the bit-reversal permutation (no mont_in), so the input must
 * already be in the Montgomery domain. This saves one Montgomery
 * multiplication per entry when chaining forward and inverse
 * transforms. The output contains standard (non-Montgomery) integers.
 *
 * a_hat must not alias a.
 *
 * Complexity:
 *   Time: O(n log n)
 *   Auxiliary memory: O(1)
 *   Output memory: O(n)
 */
void ntt_u64_cyclic_inverse_montgomery_in(u64* a_hat, const u64* a,
                                          ntt_ctx_u64* ctx)
{
  u64 n = ctx->n;
  u64 q = ctx->q;

  // 1. Bit-reversal permutation AND entering Montgomery space
  for (u64 i = 0; i < n; i++) {
    u64 rev = ctx->bit_rev_indices[i];
    a_hat[rev] = a[i];
  }

  // 2. Cooley-Tukey Butterfly
  for (u64 len = 2; len <= n; len <<= 1) {
    u64 half = len >> 1;
    u64 step = n / len;

    for (u64 i = 0; i < n; i += len) {
      for (u64 j = 0; j < half; j++) {
        u64 twiddle =
            ctx->omega_inv_powers[j * step];  // Already in Montgomery form

        u64 even = i + j;
        u64 odd = i + j + half;

        u64 t = mont_mul(a_hat[odd], twiddle, &ctx->mctx);
        u64 u = a_hat[even];

        a_hat[even] = mod_add(u, t, q);
        a_hat[odd] = mod_sub(u, t, q);
      }
    }
  }

  // 3. Scale by n^-1 mod q AND exit Montgomery space
  for (u64 i = 0; i < n; i++) {
    // Multiply by n_inv (which is in Montgomery form)
    a_hat[i] = mont_mul(a_hat[i], ctx->n_inv, &ctx->mctx);

    // Exit Montgomery space to get the final standard integer
    a_hat[i] = mont_out(a_hat[i], &ctx->mctx);
  }
}
