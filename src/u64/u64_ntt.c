/*
 * u64_ntt.c
 *
 * Single-limb (u64) Number Theoretic Transform.
 *
 * This file implements a cyclic (Cooley-Tukey) NTT and its inverse for
 * u64 arrays over a prime field, with all butterflies performed in the
 * Montgomery domain for speed. It also provides context initialization
 * for a generic prime (with caller-supplied roots of unity) and for
 * the golden prime p = 5 * 2^55 + 1.
 *
 * A second transform family operates on the Goldilocks prime
 * p = 2^64 - 2^32 + 1 in the canonical (non-Montgomery) domain: p > 2^63
 * breaks the 2q butterfly bounds, so all reductions there use 2^64 =
 * 2^32 - 1 fold arithmetic (u64_mul_goldilocks()).
 *
 * The transform length is n = 2^k for k <= 54 (golden prime) or
 * k <= 31 (Goldilocks).
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

/*
 * TODO: switch to Gentleman-Sande NTT for the forward pass to eliminate bit
 * reversal use optimized mod with bit ops only reduce every 4th step preorder
 * twiddle factors
 *
 */

#include <pthread.h>

#include "../../include/u64.h"

static inline u64 reduce_2q(u64 val, u64 two_q)
{
  u64 mask = -(u64)(val >= two_q);
  return val - (two_q & mask);
}

// Harvey Gentleman-Sande (DIF) Butterfly
static inline void harvey_dif_butterfly(u64* u_ptr, u64* v_ptr, u64 twiddle,
                                        u64 q, u64 two_q,
                                        const u64_mont_ctx_t* mctx)
{
  u64 u = *u_ptr;
  u64 v = *v_ptr;

  // Sum term: bounded in [0, 2q)
  *u_ptr = reduce_2q(u + v, two_q);

  // Diff term: (u + 2q - v) is strictly positive and < 4q
  u64 diff = u + two_q - v;

  // Montgomery multiplication reduces output back to < q
  *v_ptr = u64_mont_mul(diff, twiddle, mctx);
}

// Harvey Cooley-Tukey (DIT) Butterfly
static inline void harvey_dit_butterfly(u64* u_ptr, u64* v_ptr, u64 twiddle,
                                        u64 q, u64 two_q,
                                        const u64_mont_ctx_t* mctx)
{
  u64 u = *u_ptr;
  u64 v = *v_ptr;

  // t = v * twiddle < q
  u64 t = u64_mont_mul(v, twiddle, mctx);

  *u_ptr = reduce_2q(u + t, two_q);
  /* reduce_2q on the difference term too: an unreduced (u + 2q - t)
     grows the lazy domain by 2q every stage (to (2s+1)q after s
     stages), which overflows u64 for large q at deep transforms. One
     masked subtract here keeps the invariant [0, 2q) for all stages,
     leaving only the 4q < 2^64 butterfly bound (q < 2^62). */
  *v_ptr = reduce_2q(u + two_q - t, two_q);
}

void u64_ntt_cyclic_forward(u64* a_hat, const u64* a, const u64_ntt_ctx_t* ctx)
{
  u64 n = ctx->n;
  u64 q = ctx->q;
  u64 two_q = q << 1;
  u64 tw_offset = 0;
  u64 half = n >> 1;
  const u64* twiddles = &ctx->twiddle_forward[tw_offset];

  // 1. Entering Montgomery space and computing first pass
  for (u64 i = 0; i < half; i++) {
    u64 u = u64_mont_in(a[i], &ctx->mctx);
    u64 v = u64_mont_in(a[i + half], &ctx->mctx);

    a_hat[i] = reduce_2q(u + v, two_q);
    u64 diff = u + two_q - v;
    a_hat[i + half] = u64_mont_mul(diff, twiddles[i], &ctx->mctx);
  }

  tw_offset += half;

  // 2. Gentleman-Sande (DIF) Butterfly
  for (u64 len = n >> 1; len >= 2; len >>= 1) {
    u64 half = len >> 1;
    u64 step = n / len;
    twiddles = &ctx->twiddle_forward[tw_offset];
    for (u64 i = 0; i < n; i += len) {
      u64 j = 0;
      // loop unrolling for simd
      for (; j + 3 < half; j += 4) {
        harvey_dif_butterfly(&a_hat[i + j + 0], &a_hat[i + j + 0 + half],
                             twiddles[j + 0], q, two_q, &ctx->mctx);
        harvey_dif_butterfly(&a_hat[i + j + 1], &a_hat[i + j + 1 + half],
                             twiddles[j + 1], q, two_q, &ctx->mctx);
        harvey_dif_butterfly(&a_hat[i + j + 2], &a_hat[i + j + 2 + half],
                             twiddles[j + 2], q, two_q, &ctx->mctx);
        harvey_dif_butterfly(&a_hat[i + j + 3], &a_hat[i + j + 3 + half],
                             twiddles[j + 3], q, two_q, &ctx->mctx);
      }
      for (; j < half; j++) {
        harvey_dif_butterfly(&a_hat[i + j], &a_hat[i + j + half], twiddles[j],
                             q, two_q, &ctx->mctx);
      }
    }
    tw_offset += half;
  }
}

rabin_err_t u64_ntt_ctx_init_golden(u64_ntt_ctx_t* ctx, u64 k)
{
  if (ctx == NULL) return RABIN_ERR_NULL_PTR;
  if (k > 54) return RABIN_ERR_INVALID_ARG;

  // p = 5 * 2^55 + 1
  u64 p = 180143985094819841ULL;
  u64 g = 3;

  // Calculate the exponent for psi: c * 2^(55 - (k+1))
  // We use 5ULL to ensure it shifts as a 64-bit integer
  u64 exp = 5ULL << (55 - (k + 1));

  // psi = g^exp mod p
  u64 psi = u64_mod_pow(g, exp, p);

  // omega = psi^2 mod p
  u64 omega = u64_mod_mul(psi, psi, p);

  // Initialize the Montgomery NTT context
  return u64_ntt_ctx_init(ctx, p, k, omega, psi);
}

rabin_err_t u64_ntt_ctx_init_goldilocks(u64_ntt_ctx_t* ctx, u64 k)
{
  if (ctx == NULL) return RABIN_ERR_NULL_PTR;

  /* p - 1 = (2^32 - 1) * 2^32 has 2-adicity 32; table size caps k at 31. */
  if (k > 31) return RABIN_ERR_INVALID_ARG;

  const u64 p = 0xFFFFFFFF00000001ULL;  // 2^64 - 2^32 + 1
  const u64 g = 7ULL;                   // minimal primitive root mod p
  u64 n = 1ULL << k;

  /* omega = g^((p-1)/n) is a primitive n-th root of unity:
     (p - 1) / 2^k = (2^32 - 1) << (32 - k), a valid shift for k <= 31. */
  u64 omega = u64_mod_pow(g, 0xFFFFFFFFULL << (32 - k), p);

  *ctx = (u64_ntt_ctx_t){0};
  ctx->n = n;
  ctx->k = k;
  ctx->q = p;

  ctx->omega_powers = malloc(sizeof(u64) * n);
  ctx->omega_inv_powers = malloc(sizeof(u64) * n);
  ctx->twiddle_forward = malloc(sizeof(u64) * n);
  ctx->twiddle_inverse = malloc(sizeof(u64) * n);

  if (ctx->omega_powers == NULL || ctx->omega_inv_powers == NULL ||
      ctx->twiddle_forward == NULL || ctx->twiddle_inverse == NULL) {
    free(ctx->omega_powers);
    free(ctx->omega_inv_powers);
    free(ctx->bit_rev_indices);
    free(ctx->twiddle_forward);
    free(ctx->twiddle_inverse);
    ctx->omega_powers = NULL;
    ctx->omega_inv_powers = NULL;
    ctx->bit_rev_indices = NULL;
    ctx->twiddle_forward = NULL;
    ctx->twiddle_inverse = NULL;
    return RABIN_ERR_OUT_OF_MEMORY;
  }

  /* All tables canonical (Montgomery is unusable for p > 2^63). */
  ctx->n_inv = u64_mod_inverse_euclid(n, p);

  u64 current_omega = 1;
  u64 omega_inv = u64_mod_inverse_euclid(omega, p);
  u64 current_omega_inv = 1;
  for (u64 i = 0; i < n; i++) {
    ctx->omega_powers[i] = current_omega;
    ctx->omega_inv_powers[i] = current_omega_inv;

    current_omega = u64_mul_goldilocks(current_omega, omega);
    current_omega_inv = u64_mul_goldilocks(current_omega_inv, omega_inv);
  }

  /* Pre-ordered twiddle tables, same layout as u64_ntt_ctx_init(). */
  u64 fwd_idx = 0;
  for (u64 len = n; len >= 2; len >>= 1) {
    u64 half = len >> 1;
    u64 step = n / len;
    for (u64 j = 0; j < half; j++)
      ctx->twiddle_forward[fwd_idx++] = ctx->omega_powers[j * step];
  }

  u64 inv_idx = 0;
  for (u64 len = 2; len <= n; len <<= 1) {
    u64 half = len >> 1;
    u64 step = n / len;
    for (u64 j = 0; j < half; j++)
      ctx->twiddle_inverse[inv_idx++] = ctx->omega_inv_powers[j * step];
  }

  return RABIN_SUCCESS;
}

rabin_err_t u64_ntt_ctx_init_proth2(u64_ntt_ctx_t* ctx, u64 k)
{
  if (ctx == NULL) return RABIN_ERR_NULL_PTR;
  if (k > 54) return RABIN_ERR_INVALID_ARG;

  /* Second RNS prime q2 = 69 * 2^55 + 1 (same 2-adicity and analogous
     generator as the golden prime). q2 < 2^62 is required: the lazy
     butterfly domain keeps values in [0, 2q) and the DIF diff term
     reaches 4q, so REDC needs 4q^2 < q * 2^64, i.e. q < 2^62. With
     q1 = 5*2^55+1 the product q1*q2 ~ 2^118 keeps W=32 chunk
     convolutions exact for every admissible transform length. */
  u64 p = 2485986994308513793ULL;
  u64 g = 5;

  // Exponent for psi: c * 2^(55 - (k+1)); 69 * 2^55 < 2^62 fits in u64
  u64 exp = 69ULL << (55 - (k + 1));

  u64 psi = u64_mod_pow(g, exp, p);
  u64 omega = u64_mod_mul(psi, psi, p);

  return u64_ntt_ctx_init(ctx, p, k, omega, psi);
}

/*
 * Cached golden-prime NTT contexts, one per transform size k (0..54).
 *
 * The omega power tables and bit-reversal indices are pure functions
 * of k, so each context is built lazily on first use and reused
 * forever after. Contexts are read-only after initialization, so they
 * can be shared across threads.
 */
static u64_ntt_ctx_t golden_ctxs[55];
static bool golden_ctx_ready[55] = {false};
static pthread_mutex_t golden_ctx_lock = PTHREAD_MUTEX_INITIALIZER;

u64_ntt_ctx_t* u64_ntt_ctx_golden_cached(u64 k)
{
  if (k > 54) return NULL;

  pthread_mutex_lock(&golden_ctx_lock);
  if (!golden_ctx_ready[k]) {
    if (u64_ntt_ctx_init_golden(&golden_ctxs[k], k) != RABIN_SUCCESS) {
      pthread_mutex_unlock(&golden_ctx_lock);
      return NULL;
    }
    golden_ctx_ready[k] = true;
  }
  pthread_mutex_unlock(&golden_ctx_lock);

  return &golden_ctxs[k];
}

/*
 * Cached q2 = 131*2^55+1 NTT contexts, one per transform size k (0..54).
 * Same lazy, read-only-after-init scheme as the golden-prime cache.
 */
static u64_ntt_ctx_t proth2_ctxs[55];
static bool proth2_ctx_ready[55] = {false};
static pthread_mutex_t proth2_ctx_lock = PTHREAD_MUTEX_INITIALIZER;

u64_ntt_ctx_t* u64_ntt_ctx_proth2_cached(u64 k)
{
  if (k > 54) return NULL;

  pthread_mutex_lock(&proth2_ctx_lock);
  if (!proth2_ctx_ready[k]) {
    if (u64_ntt_ctx_init_proth2(&proth2_ctxs[k], k) != RABIN_SUCCESS) {
      pthread_mutex_unlock(&proth2_ctx_lock);
      return NULL;
    }
    proth2_ctx_ready[k] = true;
  }
  pthread_mutex_unlock(&proth2_ctx_lock);

  return &proth2_ctxs[k];
}

/*
 * Cached Goldilocks NTT contexts, one per transform size k (0..31).
 * Same lazy, read-only-after-init scheme as the golden-prime cache.
 */
static u64_ntt_ctx_t goldilocks_ctxs[32];
static bool goldilocks_ctx_ready[32] = {false};
static pthread_mutex_t goldilocks_ctx_lock = PTHREAD_MUTEX_INITIALIZER;

u64_ntt_ctx_t* u64_ntt_ctx_goldilocks_cached(u64 k)
{
  if (k > 31) return NULL;

  pthread_mutex_lock(&goldilocks_ctx_lock);
  if (!goldilocks_ctx_ready[k]) {
    if (u64_ntt_ctx_init_goldilocks(&goldilocks_ctxs[k], k) != RABIN_SUCCESS) {
      pthread_mutex_unlock(&goldilocks_ctx_lock);
      return NULL;
    }
    goldilocks_ctx_ready[k] = true;
  }
  pthread_mutex_unlock(&goldilocks_ctx_lock);

  return &goldilocks_ctxs[k];
}

/* ---- Goldilocks fold arithmetic (canonical domain, p = 2^64-2^32+1) ----
 * p > 2^63 makes the 2q Montgomery butterfly bounds unusable, so every
 * value stays reduced in [0, p) using 2^64 = 2^32 - 1 fold identities:
 * a wrapped add carries a surplus 2^64 (fix by adding 2^32 - 1), a
 * borrowed sub carries a debt of 2^64 (fix by subtracting 2^32 - 1;
 * equivalent to adding p mod 2^64). Property-tested against u128
 * modulo in tests/test_ntt.c. */
#define GL_P 0xFFFFFFFF00000001ULL
#define GL_FOLD 0xFFFFFFFFULL /* 2^32 - 1 */

static inline u64 gl_add(u64 u, u64 v)
{
  u64 s = u + v;
  u64 w = -(u64)(s < u); /* wrapped: surplus 2^64 -> add 2^32 - 1 */
  s += GL_FOLD & w;
  s -= GL_P & ~w & -(u64)(s >= GL_P);
  return s;
}

static inline u64 gl_sub(u64 u, u64 v)
{
  u64 d = u - v;
  return d - (GL_FOLD & -(u64)(u < v)); /* borrowed: debt 2^64 -> sub 2^32-1 */
}

/* Gentleman-Sande (DIF) butterfly, canonical Goldilocks domain. */
static inline void gl_dif_butterfly(u64* u_ptr, u64* v_ptr, u64 twiddle)
{
  u64 u = *u_ptr;
  u64 v = *v_ptr;

  *u_ptr = gl_add(u, v);
  *v_ptr = u64_mul_goldilocks(gl_sub(u, v), twiddle);
}

/* Cooley-Tukey (DIT) butterfly, canonical Goldilocks domain. */
static inline void gl_dit_butterfly(u64* u_ptr, u64* v_ptr, u64 twiddle)
{
  u64 u = *u_ptr;
  u64 t = u64_mul_goldilocks(*v_ptr, twiddle);

  *u_ptr = gl_add(u, t);
  *v_ptr = gl_sub(u, t);
}

void u64_ntt_goldilocks_forward(u64* a_hat, const u64* a,
                                const u64_ntt_ctx_t* ctx)
{
  u64 n = ctx->n;
  u64 tw_offset = 0;
  u64 half = n >> 1;
  const u64* twiddles = &ctx->twiddle_forward[tw_offset];

  /* 1. First pass (len = n), reading straight from the input array. */
  for (u64 i = 0; i < half; i++) {
    u64 u = a[i], v = a[i + half];

    a_hat[i] = gl_add(u, v);
    a_hat[i + half] = u64_mul_goldilocks(gl_sub(u, v), twiddles[i]);
  }

  tw_offset += half;

  /* 2. Gentleman-Sande (DIF) butterfly. */
  for (u64 len = n >> 1; len >= 2; len >>= 1) {
    u64 half = len >> 1;
    u64 j;
    twiddles = &ctx->twiddle_forward[tw_offset];
    for (u64 i = 0; i < n; i += len) {
      /* loop unrolling for simd */
      for (j = 0; j + 3 < half; j += 4) {
        gl_dif_butterfly(&a_hat[i + j + 0], &a_hat[i + j + 0 + half],
                         twiddles[j + 0]);
        gl_dif_butterfly(&a_hat[i + j + 1], &a_hat[i + j + 1 + half],
                         twiddles[j + 1]);
        gl_dif_butterfly(&a_hat[i + j + 2], &a_hat[i + j + 2 + half],
                         twiddles[j + 2]);
        gl_dif_butterfly(&a_hat[i + j + 3], &a_hat[i + j + 3 + half],
                         twiddles[j + 3]);
      }
      for (; j < half; j++)
        gl_dif_butterfly(&a_hat[i + j], &a_hat[i + j + half], twiddles[j]);
    }
    tw_offset += half;
  }
}

void u64_ntt_goldilocks_inverse(u64* a_hat, const u64* a,
                                const u64_ntt_ctx_t* ctx)
{
  u64 n = ctx->n;
  u64 tw_offset = 0;

  for (u64 i = 0; i < n; i++) a_hat[i] = a[i];

  /* 1. Cooley-Tukey (DIT) butterfly. */
  for (u64 len = 2; len < n; len <<= 1) {
    u64 half = len >> 1;
    u64 j;
    const u64* twiddles = &ctx->twiddle_inverse[tw_offset];

    for (u64 i = 0; i < n; i += len) {
      for (j = 0; j + 3 < half; j += 4) {
        gl_dit_butterfly(&a_hat[i + j + 0], &a_hat[i + j + 0 + half],
                         twiddles[j + 0]);
        gl_dit_butterfly(&a_hat[i + j + 1], &a_hat[i + j + 1 + half],
                         twiddles[j + 1]);
        gl_dit_butterfly(&a_hat[i + j + 2], &a_hat[i + j + 2 + half],
                         twiddles[j + 2]);
        gl_dit_butterfly(&a_hat[i + j + 3], &a_hat[i + j + 3 + half],
                         twiddles[j + 3]);
      }
      for (; j < half; j++)
        gl_dit_butterfly(&a_hat[i + j], &a_hat[i + j + half], twiddles[j]);
    }
    tw_offset += half;
  }

  /* 2. Final pass fused with the n^{-1} scaling. */
  u64 half = n >> 1;
  const u64* twiddles = &ctx->twiddle_inverse[tw_offset];
  u64 j = 0;

#define FUSED_INV_STEP(idx)                                             \
  do {                                                                  \
    u64 u = a_hat[(idx)];                                               \
    u64 v = a_hat[(idx) + half];                                        \
    u64 t = u64_mul_goldilocks(v, twiddles[(idx)]);                     \
    a_hat[(idx)] = u64_mul_goldilocks(gl_add(u, t), ctx->n_inv);        \
    a_hat[(idx) + half] = u64_mul_goldilocks(gl_sub(u, t), ctx->n_inv); \
  } while (0)

  for (; j + 3 < half; j += 4) {
    FUSED_INV_STEP(j + 0);
    FUSED_INV_STEP(j + 1);
    FUSED_INV_STEP(j + 2);
    FUSED_INV_STEP(j + 3);
  }
  for (; j < half; j++) FUSED_INV_STEP(j);

#undef FUSED_INV_STEP
}

rabin_err_t u64_ntt_ctx_init(u64_ntt_ctx_t* ctx, u64 p, u64 k, u64 omega,
                             u64 psi)
{
  (void)psi;  // Accepted for interface compatibility; only omega is used.
  if (ctx == NULL) return RABIN_ERR_NULL_PTR;
  if (k >= 64) return RABIN_ERR_INVALID_ARG;

  u64 n = (1ULL << k);
  *ctx = (u64_ntt_ctx_t){0};
  ctx->n = n;
  ctx->k = k;
  ctx->q = p;

  u64_mont_init(&ctx->mctx, p);

  ctx->omega_powers = malloc(sizeof(u64) * n);
  ctx->omega_inv_powers = malloc(sizeof(u64) * n);
  ctx->bit_rev_indices = malloc(sizeof(u64) * n);

  ctx->twiddle_forward = malloc(sizeof(u64) * n);
  ctx->twiddle_inverse = malloc(sizeof(u64) * n);

  if (ctx->omega_powers == NULL || ctx->omega_inv_powers == NULL ||
      ctx->bit_rev_indices == NULL || ctx->twiddle_forward == NULL ||
      ctx->twiddle_inverse == NULL) {
    free(ctx->omega_powers);
    free(ctx->omega_inv_powers);
    free(ctx->bit_rev_indices);
    free(ctx->twiddle_forward);
    free(ctx->twiddle_inverse);
    ctx->omega_powers = NULL;
    ctx->omega_inv_powers = NULL;
    ctx->bit_rev_indices = NULL;
    ctx->twiddle_forward = NULL;
    ctx->twiddle_inverse = NULL;
    // u64 ntt ct clear
    return RABIN_ERR_OUT_OF_MEMORY;
  }

  // Calculate n^-1 mod q and put it in Montgomery form
  u64 n_inv_standard = u64_mod_inverse_euclid(n, p);
  ctx->n_inv = u64_mont_in(n_inv_standard, &ctx->mctx);

  // Precompute powers and move to Montgomery form immediately
  u64 current_omega = 1;
  u64 omega_inv = u64_mod_inverse_euclid(omega, p);
  u64 current_omega_inv = 1;

  for (u64 i = 0; i < n; i++) {
    ctx->omega_powers[i] = u64_mont_in(current_omega, &ctx->mctx);
    ctx->omega_inv_powers[i] = u64_mont_in(current_omega_inv, &ctx->mctx);

    current_omega = u64_mod_mul(current_omega, omega,
                                p);  // Standard u64_mod_mul for precalc
    current_omega_inv = u64_mod_mul(current_omega_inv, omega_inv, p);
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

  // --- PRE-ORDER TWIDDLE TABLES (Unit Stride Access) ---

  // 1. Forward Pass Twiddles (DIF)
  u64 fwd_idx = 0;
  for (u64 len = n; len >= 2; len >>= 1) {
    u64 half = len >> 1;
    u64 step = n / len;
    for (u64 j = 0; j < half; j++) {
      ctx->twiddle_forward[fwd_idx++] = ctx->omega_powers[j * step];
    }
  }

  // 2. Inverse Pass Twiddles (DIT)
  u64 inv_idx = 0;
  for (u64 len = 2; len <= n; len <<= 1) {
    u64 half = len >> 1;
    u64 step = n / len;
    for (u64 j = 0; j < half; j++) {
      ctx->twiddle_inverse[inv_idx++] = ctx->omega_inv_powers[j * step];
    }
  }

  return RABIN_SUCCESS;
}

void u64_ntt_ctx_clear(u64_ntt_ctx_t* ctx)
{
  if (ctx == NULL) return;

  free(ctx->omega_powers);
  free(ctx->omega_inv_powers);
  free(ctx->bit_rev_indices);
  free(ctx->twiddle_forward);
  free(ctx->twiddle_inverse);
  ctx->omega_powers = NULL;
  ctx->omega_inv_powers = NULL;
  ctx->bit_rev_indices = NULL;
  ctx->twiddle_forward = NULL;
  ctx->twiddle_inverse = NULL;
}

void u64_ntt_cyclic_inverse(u64* a_hat, const u64* a, const u64_ntt_ctx_t* ctx)
{
  u64 n = ctx->n;
  u64 q = ctx->q;
  u64 two_q = q << 1;
  u64 tw_offset = 0;

  // 1. Entering Montgomery space
  for (u64 i = 0; i < n; i++) {
    a_hat[i] = u64_mont_in(a[i], &ctx->mctx);
  }

  // 2. Cooley-Tukey Butterfly
  for (u64 len = 2; len < n; len <<= 1) {
    u64 half = len >> 1;
    const u64* twiddles = &ctx->twiddle_inverse[tw_offset];

    for (u64 i = 0; i < n; i += len) {
      u64 j = 0;
      for (; j + 3 < half; j += 4) {
        harvey_dit_butterfly(&a_hat[i + j + 0], &a_hat[i + j + 0 + half],
                             twiddles[j + 0], q, two_q, &ctx->mctx);
        harvey_dit_butterfly(&a_hat[i + j + 1], &a_hat[i + j + 1 + half],
                             twiddles[j + 1], q, two_q, &ctx->mctx);
        harvey_dit_butterfly(&a_hat[i + j + 2], &a_hat[i + j + 2 + half],
                             twiddles[j + 2], q, two_q, &ctx->mctx);
        harvey_dit_butterfly(&a_hat[i + j + 3], &a_hat[i + j + 3 + half],
                             twiddles[j + 3], q, two_q, &ctx->mctx);
      }

      for (; j < half; j++) {
        harvey_dit_butterfly(&a_hat[i + j], &a_hat[i + j + half], twiddles[j],
                             q, two_q, &ctx->mctx);
      }
    }
    tw_offset += half;
  }

  u64 half = n >> 1;
  const u64* twiddles = &ctx->twiddle_inverse[tw_offset];
  u64 j = 0;

#define FUSED_INV_STEP(idx)                                                    \
  do {                                                                         \
    u64 u = a_hat[(idx)];                                                      \
    u64 v = a_hat[(idx) + half];                                               \
    u64 t = u64_mont_mul(v, twiddles[(idx)], &ctx->mctx);                      \
    u64 res_u = reduce_2q(u + t, two_q);                                       \
    u64 res_v = u + two_q - t;                                                 \
    a_hat[(idx)] =                                                             \
        u64_mont_out(u64_mont_mul(res_u, ctx->n_inv, &ctx->mctx), &ctx->mctx); \
    a_hat[(idx) + half] =                                                      \
        u64_mont_out(u64_mont_mul(res_v, ctx->n_inv, &ctx->mctx), &ctx->mctx); \
  } while (0)

  // Unrolled final fused pass
  for (; j + 3 < half; j += 4) {
    FUSED_INV_STEP(j + 0);
    FUSED_INV_STEP(j + 1);
    FUSED_INV_STEP(j + 2);
    FUSED_INV_STEP(j + 3);
  }
  for (; j < half; j++) {
    FUSED_INV_STEP(j);
  }

#undef FUSED_INV_STEP
}

void u64_ntt_cyclic_inverse_montgomery_in(u64* a_hat, const u64* a,
                                          const u64_ntt_ctx_t* ctx)
{
  u64 n = ctx->n;
  u64 q = ctx->q;
  u64 two_q = q << 1;
  u64 tw_offset = 0;

  for (u64 i = 0; i < n; i++) {
    a_hat[i] = a[i];
  }

  // 2. Cooley-Tukey Butterfly
  for (u64 len = 2; len < n; len <<= 1) {
    u64 half = len >> 1;
    const u64* twiddles = &ctx->twiddle_inverse[tw_offset];

    for (u64 i = 0; i < n; i += len) {
      u64 j = 0;
      for (; j + 3 < half; j += 4) {
        harvey_dit_butterfly(&a_hat[i + j + 0], &a_hat[i + j + 0 + half],
                             twiddles[j + 0], q, two_q, &ctx->mctx);
        harvey_dit_butterfly(&a_hat[i + j + 1], &a_hat[i + j + 1 + half],
                             twiddles[j + 1], q, two_q, &ctx->mctx);
        harvey_dit_butterfly(&a_hat[i + j + 2], &a_hat[i + j + 2 + half],
                             twiddles[j + 2], q, two_q, &ctx->mctx);
        harvey_dit_butterfly(&a_hat[i + j + 3], &a_hat[i + j + 3 + half],
                             twiddles[j + 3], q, two_q, &ctx->mctx);
      }

      for (; j < half; j++) {
        harvey_dit_butterfly(&a_hat[i + j], &a_hat[i + j + half], twiddles[j],
                             q, two_q, &ctx->mctx);
      }
    }
    tw_offset += half;
  }

  u64 half = n >> 1;
  const u64* twiddles = &ctx->twiddle_inverse[tw_offset];
  u64 j = 0;

#define FUSED_INV_STEP(idx)                                                    \
  do {                                                                         \
    u64 u = a_hat[(idx)];                                                      \
    u64 v = a_hat[(idx) + half];                                               \
    u64 t = u64_mont_mul(v, twiddles[(idx)], &ctx->mctx);                      \
    u64 res_u = reduce_2q(u + t, two_q);                                       \
    u64 res_v = u + two_q - t;                                                 \
    a_hat[(idx)] =                                                             \
        u64_mont_out(u64_mont_mul(res_u, ctx->n_inv, &ctx->mctx), &ctx->mctx); \
    a_hat[(idx) + half] =                                                      \
        u64_mont_out(u64_mont_mul(res_v, ctx->n_inv, &ctx->mctx), &ctx->mctx); \
  } while (0)

  // Unrolled final fused pass
  for (; j + 3 < half; j += 4) {
    FUSED_INV_STEP(j + 0);
    FUSED_INV_STEP(j + 1);
    FUSED_INV_STEP(j + 2);
    FUSED_INV_STEP(j + 3);
  }
  for (; j < half; j++) {
    FUSED_INV_STEP(j);
  }

#undef FUSED_INV_STEP
}
