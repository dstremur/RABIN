/*
 * bigntt.c
 *
 * Bignum Number Theoretic Transform.
 *
 * This file implements a cyclic (Cooley-Tukey) NTT over bignums, with
 * all butterflies performed in the Montgomery domain. It also provides
 * the context management (precomputed root-of-unity power tables and
 * bit-reversal indices), generator search for finite fields, and
 * Proth-prime / Goldilocks-field context constructors used by the
 * NTT-based multiplication path.
 *
 * The transform length is n = 2^k.
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

#include "../../include/bigntt.h"

#include <stdio.h>

#include "../../include/bigvector.h"
#include "../../include/u64.h"

/**
 * @brief Find a generator g of the multiplicative group F_p^* given the
 * prime factorization of p - 1.
 *
 * Let n_l = p->size, measured in 64-bit limbs.
 *
 * Repeatedly picks a random a in [2, p-2] and checks that
 * a^((p-1)/r) != 1 mod p for every prime factor r of p - 1; the first
 * candidate that passes all checks is a primitive root, stored in g.
 *
 * Complexity:
 *   Time: O(pi(p-1) * n_l^2) expected - one modular exponentiation per
 *         factor per candidate, O(1) candidates expected
 *   Auxiliary memory: O(n_l) limbs for temporaries
 *   Output memory: O(n_l) limbs
 *
 * @param[out]      g       Result storing the generator (primitive root).
 * @param[in]       p       Prime modulus.
 * @param[in]       factors Vector of prime factors of p - 1.
 */
void bn_find_gen(bignum* g, bignum* p, bigvector* factors)
{
  bignum a, b, exp, n_min_1;
  bn_init_multi(&a, &b, &exp, &n_min_1, NULL);
  bn_sub(&n_min_1, p, &BN_ONE);

start:
  bn_gen_random_range(&a, &BN_TWO, &n_min_1);

  for (u64 i = 0; i < factors->size; i++) {
    bn_div(&exp, &n_min_1, &factors->data[i]);
    bn_mod_exp(&b, &a, &exp, p);
    if (bn_is_eq_i64(&b, 1)) {
      goto start;
    }
  }

  bn_copy(g, &a);

  bn_free_multi(&a, &b, &exp, &n_min_1, NULL);
}

/**
 * @brief Find a generator g of the multiplicative group F_p^*.
 *
 * Let n_l = p->size, measured in 64-bit limbs.
 *
 * Factorizes p - 1 with bn_factorize() and delegates to
 * bn_find_gen().
 *
 * Complexity:
 *   Time: factorization of p - 1 (subexponential in practice) plus
 *         O(n_l^2) for the generator search
 *   Auxiliary memory: O(n_l) limbs
 *   Output memory: O(n_l) limbs
 *
 * @param[out] g Result storing the generator (primitive root).
 * @param[in]  p Prime modulus.
 */
void bn_find_gen_fp(bignum* g, bignum* p)
{
  bigvector factors;
  bigvector_init_dynamic(&factors);

  bignum n_min_1;
  bn_init_multi(&n_min_1, NULL);
  bn_sub(&n_min_1, p, &BN_ONE);

  bn_factorize(&factors, &n_min_1);

  bn_find_gen(g, p, &factors);

  bigvector_free(&factors);
  bn_free_multi(&n_min_1, NULL);
}

/**
 * @brief Find a generator g of F_p^* for a Proth prime p = c * 2^k + 1.
 *
 * Let n_l = p->size, measured in 64-bit limbs.
 *
 * The factorization of p - 1 = c * 2^k is known up to the factorization
 * of c: it is {2} union factorization(c). The generator search then
 * proceeds as in bn_find_gen().
 *
 * Complexity:
 *   Time: factorization of c plus O(n_l^2) for the generator search
 *   Auxiliary memory: O(n_l) limbs
 *   Output memory: O(n_l) limbs
 *
 * @param[out] g Result storing the generator (primitive root).
 * @param[in]  p Proth prime modulus.
 * @param[in]  c Odd multiplier of the Proth prime (p = c * 2^k + 1).
 */
void bn_find_gen_proth(bignum* g, bignum* p, bignum* c)
{
  bigvector factors;
  bigvector_init_dynamic(&factors);

  // 2 is always a prime factor
  bigvector_append(&factors, &BN_TWO);

  // factorize c
  bn_factorize(&factors, c);

  printf("p - 1 factors ");
  // bigvector_println(&factors);

  // find a generator
  bn_find_gen(g, p, &factors);

  bigvector_free(&factors);
}

/**
 * @brief Generate a Proth prime p = c * 2^k + 1 (with c starting at the given
 * odd value) together with a generator g and the roots psi, omega.
 *
 * Let k = exponent of the power of two.
 *
 * Scans odd c upward until c * 2^k + 1 passes BPSW, then factorizes
 * p - 1 = c * 2^k (as {2} union factorization(c)), finds a generator
 * g, and sets:
 *
 *   psi   = g^c mod p      (primitive (k+1)-th root of unity)
 *   omega = psi^2 mod p    (primitive k-th root of unity)
 *
 * Complexity:
 *   Time: O(k^3) expected per candidate for the BPSW test, plus
 *         factorization of c and O(k^2) for the generator search
 *   Auxiliary memory: O(k/64) limbs
 *   Output memory: O(k/64) limbs per output
 *
 * @param[out]     g     Result storing the generator.
 * @param[out]     p     Result storing the Proth prime.
 * @param[out]     omega Result storing the primitive k-th root of unity.
 * @param[out]     psi   Result storing the primitive (k+1)-th root of unity.
 * @param[in]      k     Exponent of the power of two.
 * @param[in]      c     Starting value for the odd multiplier c.
 */
void bn_gen_proth_ntt(bignum* g, bignum* p, bignum* omega, bignum* psi, u64 k,
                      u64 c)
{
  bignum p_bn, c_bn, two_k;
  bn_init_multi(&p_bn, &c_bn, &two_k, NULL);

  // calc 2^k
  bn_lshift(&two_k, &BN_ONE, k);

  // make c odd
  u64 curr_c = (c & 1) ? c : c + 1;

  while (1) {
    // p = curr_c * 2^k
    bn_set_u64(&c_bn, curr_c);
    bn_mul(&p_bn, &two_k, &c_bn);

    // p = p + 1
    bn_add(&p_bn, &p_bn, &BN_ONE);

    if (bn_bpsw(&p_bn)) {
      break;
    }

    curr_c += 2;
  }

  // copy prime
  bn_copy(p, &p_bn);

  // factorize p - 1 = c * 2^k
  bigvector factors;
  bigvector_init_dynamic(&factors);

  // 2 always divides p - 1
  bigvector_append(&factors, &BN_TWO);

  // factorize c
  bn_factorize(&factors, &c_bn);

  // bigvector_println(&factors);

  // now find a generator
  bn_find_gen(g, &p_bn, &factors);

  // psi = g^c mod p
  bn_mod_exp(psi, g, &c_bn, p);

  // omega = psi^2 mod p
  bn_mul(&two_k, psi, psi);
  bn_mod(omega, &two_k, p);

  bn_free_multi(&p_bn, &c_bn, &two_k, NULL);
  bigvector_free(&factors);
}

/**
 * @brief Initialize a bignum NTT context by generating a fresh Proth prime.
 *
 * Let k = log2 of the transform length (n = 2^k).
 *
 * Generates a Proth prime p = c' * 2^(k+1) + 1 (starting from the
 * given c) with its generator and roots via bn_gen_proth_ntt(), then
 * delegates to bigntt_ctx_init().
 *
 * Returns true on success, false on allocation failure.
 *
 * Complexity:
 *   Time: O(k^3) expected for the prime generation, plus O(n * k^2)
 *         for the context tables
 *   Auxiliary memory: O(k/64) limbs
 *   Output memory: O(n) bignums per table (4 tables)
 *
 * @param[out] ctx NTT context to initialize.
 * @param[in]    k log2 of the transform length (n = 2^k).
 * @param[in]    c Starting value for the odd multiplier of the Proth prime.
 *
 * @return true  On success.
 * @return false On allocation failure.
 */
bool bigntt_ctx_init_simple(ntt_ctx* ctx, u64 k, u64 c)
{
  bignum p, g, omega, psi, tmp;
  bn_init_multi(&p, &g, &omega, &psi, &tmp, NULL);

  bn_gen_proth_ntt(&g, &p, &omega, &psi, k + 1, c);

  if (!bigntt_ctx_init(ctx, &p, &g, &omega, &psi, k, c)) {
    bn_free_multi(&p, &g, &omega, &psi, &tmp, NULL);
    return false;
  }

  bn_free_multi(&p, &g, &omega, &psi, &tmp, NULL);
  return true;
}

/**
 * @brief Initialize a bignum NTT context for the Goldilocks field
 * p = 5 * 2^55 + 1.
 *
 * Let k = log2 of the transform length (n = 2^k, k <= 54).
 *
 * Uses the known primitive root g = 3 and derives:
 *
 *   psi   = g^(5 * 2^(54-k)) mod p   (primitive (k+1)-th root)
 *   omega = psi^2 mod p              (primitive k-th root)
 *
 * then delegates to bigntt_ctx_init().
 *
 * Returns false (and leaves ctx unchanged) if k > 54 or on allocation
 * failure.
 *
 * Complexity:
 *   Time: O(n * k^2) for the context tables, plus O(k^2) for the root
 *         derivation
 *   Auxiliary memory: O(k/64) limbs
 *   Output memory: O(n) bignums per table (4 tables)
 *
 * @param[out] ctx NTT context to initialize.
 * @param[in]    k log2 of the transform length (n = 2^k, k <= 54).
 *
 * @return true  On success.
 * @return false If k > 54 or on allocation failure.
 */
bool bigntt_ctx_init_golden(ntt_ctx* ctx, u64 k)
{
  if (k > 54) return false;

  bignum p, g, omega, psi, tmp, c_bn;
  bn_init_multi(&p, &g, &omega, &psi, &tmp, &c_bn, NULL);

  // p = 5* 2^55 + 1
  bn_set_u64(&p, 180143985094819841);
  bn_set_u64(&c_bn, 5);

  // factorize p - 1 = 5 * 2^55
  bigvector factors;
  bigvector_init_dynamic(&factors);

  // 2,5 always divide p - 1
  bigvector_append(&factors, &BN_TWO);
  bigvector_append(&factors, &c_bn);

  bn_set_u64(&g, 3);
  // psi = g^(p - 1 / 2^k+1) = c * 2^(55 - (k+1)) mod p
  bn_set_u64(&tmp, 5);
  bn_lshift(&tmp, &tmp, 55 - (k + 1));
  bn_mod_exp(&psi, &g, &tmp, &p);

  // omega = psi^2 mod p
  bn_mul(&tmp, &psi, &psi);
  bn_mod(&omega, &tmp, &p);

  bigvector_free(&factors);

  if (!bigntt_ctx_init(ctx, &p, &g, &omega, &psi, k, 5)) {
    bn_free_multi(&p, &g, &omega, &psi, &tmp, &c_bn, NULL);
    return false;
  }

  bn_free_multi(&p, &g, &omega, &psi, &tmp, &c_bn, NULL);
  return true;
}

/**
 * @brief Check that an NTT of the given size cannot wrap around modulo p.
 *
 * Let ntt_size = N (number of coefficients) and bit_width = W.
 *
 * The largest coefficient of the true (unreduced) convolution of two
 * degree < N polynomials with coefficients < 2^W is bounded by
 * N * (2^W - 1)^2. This returns true iff that bound is strictly less
 * than p, i.e. the cyclic convolution modulo p is exact (no
 * wrap-around).
 *
 * Complexity:
 *   Time: O(k^2) where k is the size of p in limbs
 *   Auxiliary memory: O(k) limbs for temporaries
 *   Output memory: O(1)
 *
 * @param[in]  ntt_size Number of coefficients (transform length N).
 * @param[in] bit_width Bit width W of the input coefficients.
 * @param[in]         p Modulus prime.
 *
 * @return true  If the convolution is exact (no wrap-around modulo p).
 * @return false If the convolution may wrap around modulo p.
 */
bool bn_check_ntt_safety(u64 ntt_size, u64 bit_width, const bignum* p)
{
  // Max value of a resulting coefficient is roughly ntt_size * (2^bit_width -
  // 1)^2 For bit_width = 16 and ntt_size = 2^20, max_coeff is ~2^52. Our prime
  // is ~2^{57.3}, so this is safe.
  bignum max_val, base_minus_1;
  bn_init_multi(&max_val, &base_minus_1, NULL);

  bn_lshift(&base_minus_1, &BN_ONE, bit_width);
  bn_sub(&base_minus_1, &base_minus_1, &BN_ONE);  // (2^W - 1)

  bn_mul(&max_val, &base_minus_1, &base_minus_1);  // (2^W - 1)^2

  bignum n_bn;
  bn_init(&n_bn);
  bn_set_u64(&n_bn, ntt_size);
  bn_mul(&max_val, &max_val, &n_bn);  // N * (2^W - 1)^2

  bool safe = (bn_cmp(&max_val, p) < 0);

  bn_free_multi(&max_val, &base_minus_1, &n_bn, NULL);
  return safe;
}

/**
 * @brief Initialize a bignum NTT context for the prime p and transform length
 * n = 2^k.
 *
 * Precomputes:
 *
 *   - the Montgomery context for p,
 *   - the omega^i and omega^{-i} tables (i = 0 .. n-1),
 *   - the psi^i and psi^{-i} tables (i = 0 .. n-1),
 *   - n_inv = n^{-1} mod p (for the inverse transform scaling),
 *   - the bit-reversal index table.
 *
 * All tables are converted into the Montgomery domain at the end.
 * The caller must supply omega, a primitive k-th root of unity mod p
 * (and psi, a primitive (k+1)-th root; g and c are accepted for
 * interface compatibility).
 *
 * Returns true on success. On failure (allocation failure or a
 * non-invertible root) the context is freed and false is returned.
 *
 * Complexity:
 *   Time: O(n * k^2) - n modular multiplications per table (4 tables)
 *   Auxiliary memory: O(k) limbs for temporaries
 *   Output memory: O(n) bignums per table (4 tables) plus O(n) u64
 *                  indices
 *
 * @param[out]  ctx   NTT context to initialize.
 * @param[in]     p   Prime modulus.
 * @param[in]     g   Generator (accepted for interface compatibility).
 * @param[in] omega Primitive k-th root of unity mod p.
 * @param[in]   psi   Primitive (k+1)-th root of unity mod p.
 * @param[in]     k   log2 of the transform length (n = 2^k).
 * @param[in]     c   Odd multiplier (accepted for interface compatibility).
 *
 * @return true  On success.
 * @return false On allocation failure or a non-invertible root.
 */
bool bigntt_ctx_init(ntt_ctx* ctx, bignum* p, bignum* g, bignum* omega,
                     bignum* psi, u64 k, u64 c)
{
  u64 n = (1ULL << k);
  ctx->k = k;
  ctx->n = n;

  bn_init(&ctx->q);
  bn_copy(&ctx->q, p);

  // init montgomery context
  bn_mont_ctx_init(&ctx->mctx, p);

  ctx->omega_powers = NULL;
  ctx->omega_inv_powers = NULL;
  ctx->psi_powers = NULL;
  ctx->psi_inv_powers = NULL;
  ctx->bit_rev_indices = NULL;
  bn_init(&ctx->n_inv);

  ctx->omega_powers = malloc(sizeof(bignum) * n);
  ctx->omega_inv_powers = malloc(sizeof(bignum) * n);
  ctx->bit_rev_indices = malloc(sizeof(u64) * n);

  if (!ctx->omega_powers || !ctx->omega_inv_powers || !ctx->bit_rev_indices) {
    goto fail;
  }

  for (u64 i = 0; i < n; i++) {
    bn_init(&ctx->omega_powers[i]);
    bn_init(&ctx->omega_inv_powers[i]);
  }

  // precompute w^i mod q
  bn_set_u64(&ctx->omega_powers[0], 1);

  for (u64 i = 1; i < n; i++) {
    bn_mul(&ctx->omega_powers[i], &ctx->omega_powers[i - 1], omega);
    bn_mod(&ctx->omega_powers[i], &ctx->omega_powers[i], &ctx->q);
  }

  bignum omega_inv;
  bn_init(&omega_inv);
  if (!bn_mod_inverse(&omega_inv, omega, &ctx->q)) {
    bn_free(&omega_inv);
    goto fail;
  }

  bn_set_u64(&ctx->omega_inv_powers[0], 1);
  for (u64 i = 1; i < n; i++) {
    bn_mul(&ctx->omega_inv_powers[i], &ctx->omega_inv_powers[i - 1],
           &omega_inv);
    bn_mod(&ctx->omega_inv_powers[i], &ctx->omega_inv_powers[i], &ctx->q);
  }

  bn_free(&omega_inv);

  ctx->psi_powers = malloc(sizeof(bignum) * n);
  ctx->psi_inv_powers = malloc(sizeof(bignum) * n);
  if (!ctx->psi_powers || !ctx->psi_inv_powers) {
    goto fail;
  }

  for (u64 i = 0; i < n; i++) {
    bn_init(&ctx->psi_powers[i]);
    bn_init(&ctx->psi_inv_powers[i]);
  }

  // compute psi^i mod q
  bn_set_u64(&ctx->psi_powers[0], 1);
  for (u64 i = 1; i < n; i++) {
    bn_mul(&ctx->psi_powers[i], &ctx->psi_powers[i - 1], psi);
    bn_mod(&ctx->psi_powers[i], &ctx->psi_powers[i], &ctx->q);
  }

  bignum psi_inv;
  bn_init(&psi_inv);
  if (!bn_mod_inverse(&psi_inv, psi, &ctx->q)) {
    bn_free(&psi_inv);
    goto fail;
  }

  bn_set_u64(&ctx->psi_inv_powers[0], 1);
  for (u64 i = 1; i < n; i++) {
    bn_mul(&ctx->psi_inv_powers[i], &ctx->psi_inv_powers[i - 1], &psi_inv);
    bn_mod(&ctx->psi_inv_powers[i], &ctx->psi_inv_powers[i], &ctx->q);
  }

  bn_free(&psi_inv);

  // 6. Compute n_inv = n^-1 mod q
  bignum bn_n;
  bn_init(&bn_n);
  bn_set_u64(&bn_n, n);
  if (!bn_mod_inverse(&ctx->n_inv, &bn_n, &ctx->q)) {
    bn_free(&bn_n);
    goto fail;
  }
  bn_free(&bn_n);

  for (u64 i = 0; i < n; i++) {
    bn_mont_in(&ctx->omega_powers[i], &ctx->omega_powers[i], &ctx->mctx);
    bn_mont_in(&ctx->omega_inv_powers[i], &ctx->omega_inv_powers[i],
               &ctx->mctx);
    bn_mont_in(&ctx->psi_powers[i], &ctx->psi_powers[i], &ctx->mctx);
    bn_mont_in(&ctx->psi_inv_powers[i], &ctx->psi_inv_powers[i], &ctx->mctx);
  }
  bn_mont_in(&ctx->n_inv, &ctx->n_inv, &ctx->mctx);

  // 7. Precompute Bit-Reversal Permutation Indices
  u64 bits = 0;
  while (((u64)1 << bits) < n) {
    bits++;
  }

  for (u64 i = 0; i < n; i++) {
    u64 rev = 0;
    u64 temp = i;
    for (u64 j = 0; j < bits; j++) {
      rev = (rev << 1) | (temp & 1);
      temp >>= 1;
    }
    ctx->bit_rev_indices[i] = rev;
  }

  return true;

fail:
  bigntt_ctx_free(ctx);
  return false;
}

/**
 * @brief Free all storage held by a bignum NTT context.
 *
 * Frees the four power tables (each entry is a bignum), the
 * bit-reversal table, the modulus, n_inv, and the Montgomery context.
 * Safe to call on a partially initialized context.
 *
 * Complexity:
 *   Time: O(n)
 *   Auxiliary memory: O(1)
 *   Output memory: O(1)
 *
 * @param[in,out] ctx NTT context to free.
 */
void bigntt_ctx_free(ntt_ctx* ctx)
{
  if (!ctx) return;

  if (ctx->omega_powers) {
    for (u64 i = 0; i < ctx->n; i++) bn_free(&ctx->omega_powers[i]);
    free(ctx->omega_powers);
    ctx->omega_powers = NULL;
  }
  if (ctx->omega_inv_powers) {
    for (u64 i = 0; i < ctx->n; i++) bn_free(&ctx->omega_inv_powers[i]);
    free(ctx->omega_inv_powers);
    ctx->omega_inv_powers = NULL;
  }
  if (ctx->psi_powers) {
    for (u64 i = 0; i < ctx->n; i++) bn_free(&ctx->psi_powers[i]);
    free(ctx->psi_powers);
    ctx->psi_powers = NULL;
  }
  if (ctx->psi_inv_powers) {
    for (u64 i = 0; i < ctx->n; i++) bn_free(&ctx->psi_inv_powers[i]);
    free(ctx->psi_inv_powers);
    ctx->psi_inv_powers = NULL;
  }

  if (ctx->bit_rev_indices) {
    free(ctx->bit_rev_indices);
    ctx->bit_rev_indices = NULL;
  }

  bn_free(&ctx->q);
  bn_free(&ctx->n_inv);
  bn_mont_ctx_free(&ctx->mctx);
}

/**
 * @brief Forward cyclic NTT of a bignum polynomial: a_hat = NTT(a).
 *
 * Let n = ctx->n = 2^k and k_l = the size of ctx->q in limbs.
 *
 * Iterative Cooley-Tukey (decimation in time):
 *
 *   1. bit-reversal permutation of the coefficients, converting each
 *      into the Montgomery domain on the fly (coefficients beyond
 *      a->deg are zero),
 *   2. log2(n) stages of butterflies with precomputed twiddle factors
 *      (omega powers, in Montgomery form); additions and subtractions
 *      are reduced with a single conditional add/subtract of q.
 *
 * The output a_hat is in the Montgomery domain with a_hat->deg = n - 1.
 *
 * Complexity:
 *   Time: O(n log n * k_l^2) - each butterfly costs a few k_l-limb
 *         operations
 *   Auxiliary memory: O(n) bignums for the working copy
 *   Output memory: O(n) bignums
 *
 * @param[out] a_hat Result storing the forward NTT (Montgomery domain).
 * @param[in]      a Input polynomial.
 * @param[in]    ctx NTT context.
 */
void bigntt_cyclic_forward(bigpoly* a_hat, bigpoly* a, ntt_ctx* ctx)
{
  u64 n = ctx->n;
  bigpoly res;
  bigpoly_init(&res);
  bigpoly_alloc(&res, n);

  // copy and reverse bits and into montgomery form
  for (u64 i = 0; i < n; i++) {
    u64 rev = ctx->bit_rev_indices[i];
    if (i <= a->deg) {
      bn_mont_in(&res.coeff[rev], &a->coeff[i], &ctx->mctx);
    } else {
      bn_set_u64(&res.coeff[rev], 0);
    }
  }

  // cooley tukey butterfly

  bignum t, u;
  bn_init_multi(&t, &u, NULL);

  for (u64 len = 2; len <= n; len <<= 1) {
    u64 half = len >> 1;
    u64 step = n / len;

    for (u64 i = 0; i < n; i += len) {
      for (u64 j = 0; j < half; j++) {
        u64 twiddle = j * step;
        bignum* w = &ctx->omega_powers[twiddle];

        u64 even = i + j;
        u64 odd = i + j + half;

        // t = res[odd] * w mod q
        bn_mont_mul(&t, &res.coeff[odd], w, &ctx->mctx);
        // bn_mod(&t, &t, &ctx->q);

        // u = res.data[even]
        bn_copy(&u, &res.coeff[even]);

        // res[even] = (u + t) mod q
        bn_add(&res.coeff[even], &u, &t);
        // bn_mod(&res.coeff[even], &res.coeff[even], &ctx->q);

        if (bn_cmp(&res.coeff[even], &ctx->q) >= 0) {
          bn_sub(&res.coeff[even], &res.coeff[even], &ctx->q);
        }

        // res[odd] = (u - t) mod q
        bn_sub(&res.coeff[odd], &u, &t);
        bn_add(&res.coeff[odd], &res.coeff[odd], &ctx->q);
        // bn_mod(&res.coeff[odd], &res.coeff[odd], &ctx->q);
        if (bn_cmp(&res.coeff[odd], &ctx->q) >= 0) {
          bn_sub(&res.coeff[odd], &res.coeff[odd], &ctx->q);
        }
      }
    }
  }

  bn_free_multi(&t, &u, NULL);

  bigpoly_free(a_hat);
  bigpoly_init(a_hat);
  bigpoly_alloc(a_hat, n);

  // bigpoly_copy(a_hat, &res);

  for (u64 i = 0; i < n; i++) {
    bn_copy(&a_hat->coeff[i], &res.coeff[i]);
  }

  a_hat->deg = n - 1;

  bigpoly_free(&res);
}

/**
 * @brief Inverse cyclic NTT of a bignum polynomial: a_hat = NTT^{-1}(a).
 *
 * Let n = ctx->n = 2^k and k_l = the size of ctx->q in limbs.
 *
 * Thin wrapper around bigntt_cyclic_inverse_mont_in() for inputs in
 * the normal domain. The output contains standard (non-Montgomery)
 * integers with a_hat->deg = n - 1.
 *
 * Complexity:
 *   Time: O(n log n * k_l^2)
 *   Auxiliary memory: O(n) bignums for the working copy
 *   Output memory: O(n) bignums
 *
 * @param[out] a_hat Result storing the inverse NTT (normal domain).
 * @param[in]      a Input polynomial (normal domain).
 * @param[in]    ctx NTT context.
 */
void bigntt_cyclic_inverse(bigpoly* a_hat, bigpoly* a, ntt_ctx* ctx)
{
  bigntt_cyclic_inverse_mont_in(a_hat, a, ctx);
}

/**
 * @brief Inverse cyclic NTT with Montgomery-domain input:
 * a_hat = NTT^{-1}(a), where a is already in the Montgomery domain.
 *
 * Let n = ctx->n = 2^k and k_l = the size of ctx->q in limbs.
 *
 * Same structure as the forward transform, but with the inverse
 * twiddle factors (omega^{-i}), followed by scaling with n^{-1} mod q
 * and conversion out of the Montgomery domain. This variant skips the
 * mont_in of step 1, so the input must already be in the Montgomery
 * domain (e.g. the output of bigntt_cyclic_forward()).
 *
 * The output contains standard (non-Montgomery) integers with
 * a_hat->deg = n - 1.
 *
 * Complexity:
 *   Time: O(n log n * k_l^2)
 *   Auxiliary memory: O(n) bignums for the working copy
 *   Output memory: O(n) bignums
 *
 * @param[out] a_hat Result storing the inverse NTT (normal domain).
 * @param[in]      a Input polynomial (Montgomery domain).
 * @param[in]    ctx NTT context.
 */
void bigntt_cyclic_inverse_mont_in(bigpoly* a_hat, bigpoly* a, ntt_ctx* ctx)
{
  u64 n = ctx->n;
  bigpoly res;
  bigpoly_init(&res);
  bigpoly_alloc(&res, n);

  // copy and reverse bits
  for (u64 i = 0; i < n; i++) {
    u64 rev = ctx->bit_rev_indices[i];
    if (i <= a->deg) {
      bn_copy(&res.coeff[rev], &a->coeff[i]);
    } else {
      bn_set_u64(&res.coeff[rev], 0);
    }
  }
  // cooley tukey butterfly

  bignum t, u;
  bn_init_multi(&t, &u, NULL);

  for (u64 len = 2; len <= n; len <<= 1) {
    u64 half = len >> 1;
    u64 step = n / len;

    for (u64 i = 0; i < n; i += len) {
      for (u64 j = 0; j < half; j++) {
        u64 twiddle = j * step;

        bignum* w = &ctx->omega_inv_powers[twiddle];

        u64 even = i + j;
        u64 odd = i + j + half;

        // t = res[odd] * w mod q
        bn_mont_mul(&t, &res.coeff[odd], w, &ctx->mctx);
        // bn_mod(&t, &t, &ctx->q);

        // u = res.data[even]
        bn_copy(&u, &res.coeff[even]);

        // res[even] = (u + t) mod q
        bn_add(&res.coeff[even], &u, &t);
        // bn_mod(&res.coeff[even], &res.coeff[even], &ctx->q);

        if (bn_cmp(&res.coeff[even], &ctx->q) >= 0) {
          bn_sub(&res.coeff[even], &res.coeff[even], &ctx->q);
        }

        // res[odd] = (u - t) mod q
        bn_sub(&res.coeff[odd], &u, &t);
        bn_add(&res.coeff[odd], &res.coeff[odd], &ctx->q);
        // bn_mod(&res.coeff[odd], &res.coeff[odd], &ctx->q);
        if (bn_cmp(&res.coeff[odd], &ctx->q) >= 0) {
          bn_sub(&res.coeff[odd], &res.coeff[odd], &ctx->q);
        }
      }
    }
  }

  // scale by n^-1 mod q
  for (u64 i = 0; i < n; i++) {
    bn_mont_mul(&res.coeff[i], &res.coeff[i], &ctx->n_inv, &ctx->mctx);
    // bn_mod(&res.coeff[i], &res.coeff[i], &ctx->q);

    bn_mont_out(&res.coeff[i], &res.coeff[i], &ctx->mctx);
  }

  bn_free_multi(&t, &u, NULL);

  bigpoly_free(a_hat);
  bigpoly_init(a_hat);
  bigpoly_alloc(a_hat, n);

  for (u64 i = 0; i < n; i++) {
    bn_copy(&a_hat->coeff[i], &res.coeff[i]);
  }
  // bigpoly_copy(a_hat, &res);

  a_hat->deg = n - 1;

  bigpoly_free(&res);
}
