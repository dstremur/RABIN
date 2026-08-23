/*
 * bigrns.c
 *
 * Residue Number System (RNS) arithmetic for bignums.
 *
 * This file implements a mixed-radix RNS over a set of 64-bit primes:
 * context management (CRT weights, Garner weights, per-prime
 * Montgomery contexts), conversion between bignums and residue
 * vectors, component-wise addition, and a parallel RNS-based matrix
 * determinant (determinants mod each prime computed with the u64
 * matrix routines, then recombined with the CRT).
 *
 * A ctx_rns holds the primes, their product M, the CRT weights
 * w_i = (M / p_i) * (M / p_i)^{-1} mod p_i, the Garner weights, and
 * one Montgomery context per prime.
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

#include "../../include/bigrns.h"

#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../../include/u64.h"

/*
 * Estimate how many 62-bit primes are needed to represent a bignum.
 *
 * Let k = bit length of a.
 *
 * Returns ceil(k / 62), i.e. the number of 62-bit primes whose product
 * has at least k bits. (The header declares this as
 * bigrns_estimate_primes().)
 *
 * Complexity:
 *   Time: O(n) where n is the size of a in limbs
 *   Auxiliary memory: O(1)
 *   Output memory: O(1)
 */
u64 rns_estimate_primes(const bignum* a)
{
  u64 k = bn_bit_length(a);
  u64 res = (k + 62) / 62;

  return res;
}

/*
 * Initialize an RNS context from a list of primes.
 *
 * Let c = count.
 *
 * Stores a copy of the primes, computes their product M, the CRT
 * weights w_i = (M / p_i) * (M / p_i)^{-1} mod p_i (as bignums), one
 * Montgomery context per prime, and the Garner weights
 * g_i = (p_0 * ... * p_{i-1})^{-1} mod p_i.
 *
 * Complexity:
 *   Time: O(c^2 * n) where n is the size of M in limbs (bignum
 *         multiplications and divisions by small primes)
 *   Auxiliary memory: O(n) limbs for temporaries
 *   Output memory: O(c) bignums for the CRT weights, O(c) u64s for
 *         the primes and Garner weights, O(c) mont_ctxs
 */
void rns_context_init(ctx_rns* ctx, const u64* primes, u64 count)
{
  ctx->primes = malloc(sizeof(u64) * count);
  ctx->crt_weights = malloc(sizeof(bignum) * count);
  memcpy(ctx->primes, primes, sizeof(u64) * count);
  ctx->count = count;

  bignum tmp;
  bn_init_multi(&tmp, NULL);
  bn_init(&ctx->prod);
  bn_set_u64(&ctx->prod, 1);

  // calculate product M = p_0 * ... * p_{c-1}
  for (u64 i = 0; i < count; i++) {
    bn_set_u64(&tmp, primes[i]);
    bn_mul(&ctx->prod, &ctx->prod, &tmp);
  }

  // calculate CRT weights w_i = (M / p_i) * (M / p_i)^{-1} mod p_i
  for (u64 i = 0; i < count; i++) {
    bn_init(&ctx->crt_weights[i]);

    bignum M_div_tmp, inv_bn;
    bn_init_multi(&M_div_tmp, &inv_bn, NULL);

    bn_set_u64(&tmp, primes[i]);
    bn_div(&M_div_tmp, &ctx->prod, &tmp);

    // find inverse of (M / p_i) mod p_i
    u64 m_mod_p = bn_mod_u64(&M_div_tmp, primes[i]);
    u64 inv = mod_inverse_euclid(m_mod_p, primes[i]);

    bn_set_u64(&inv_bn, inv);

    bn_mul(&ctx->crt_weights[i], &M_div_tmp, &inv_bn);

    bn_free_multi(&M_div_tmp, &inv_bn, NULL);
  }

  // one Montgomery context per prime
  ctx->m_ctxs = malloc(sizeof(mont_ctx) * count);
  for (u64 i = 0; i < count; i++) {
    mont_init(&ctx->m_ctxs[i], primes[i]);
  }

  // calculate Garner weights g_i = (p_0 * ... * p_{i-1})^{-1} mod p_i
  ctx->garner_weights = malloc(sizeof(u64) * count);
  ctx->garner_weights[0] = 1;

  bignum prev;
  bn_init_multi(&prev, NULL);
  bn_set_u64(&prev, 1);

  for (u64 i = 1; i < count; i++) {
    u64 m_prev_mod = bn_mod_u64(&prev, primes[i]);
    ctx->garner_weights[i] = mod_inverse_euclid(m_prev_mod, primes[i]);

    bn_set_u64(&tmp, primes[i - 1]);
    bn_mul(&prev, &prev, &tmp);
  }

  bn_free_multi(&tmp, &prev, NULL);
}

/*
 * Free all storage owned by an RNS context.
 *
 * Complexity:
 *   Time: O(count)
 *   Auxiliary memory: O(1)
 *   Output memory: O(1)
 */
void rns_context_free(ctx_rns* ctx)
{
  if (!ctx) return;

  // 1. Free the individual bignum weights
  if (ctx->crt_weights) {
    for (u64 i = 0; i < ctx->count; i++) {
      bn_free(&ctx->crt_weights[i]);
    }
    free(ctx->crt_weights);
    ctx->crt_weights = NULL;
  }

  // 2. Free the big product
  bn_free(&ctx->prod);

  // 3. Free the primes array
  if (ctx->primes) {
    free(ctx->primes);
    ctx->primes = NULL;
  }

  // free the montgomery contexts
  if (ctx->m_ctxs) {
    free(ctx->m_ctxs);
    ctx->m_ctxs = NULL;
  }

  // free garner weights
  if (ctx->garner_weights) {
    free(ctx->garner_weights);
    ctx->garner_weights = NULL;
  }

  ctx->count = 0;
}

/*
 * Convert a bignum to its RNS residue vector: r_i = a mod p_i.
 *
 * Let c = ctx->count.
 *
 * Allocates the residue array on first use, then reduces a modulo
 * each prime.
 *
 * Complexity:
 *   Time: O(c * n) where n is the size of a in limbs
 *   Auxiliary memory: O(1)
 *   Output memory: O(c) u64s
 */
void bignum_to_rns(rns_num* r, const bignum* a, ctx_rns* ctx)
{
  if (!r->residues) {
    r->residues = malloc(sizeof(u64) * ctx->count);
  }
  r->size = ctx->count;

  for (u64 i = 0; i < r->size; i++) {
    r->residues[i] = bn_mod_u64(a, ctx->primes[i]);
  }
}

/*
 * Component-wise addition of two RNS residue vectors: r = a + b.
 *
 * Let c = r->size.
 *
 * Adds the residues modulo the corresponding prime for each
 * component.
 *
 * Complexity:
 *   Time: O(c)
 *   Auxiliary memory: O(1)
 *   Output memory: O(c) u64s
 */
void rns_add(rns_num* r, const rns_num* a, const rns_num* b, const ctx_rns* ctx)
{
  for (u64 i = 0; i < r->size; i++) {
    r->residues[i] = mod_add(a->residues[i], b->residues[i], ctx->primes[i]);
  }
}

/*
 * Reconstruct a bignum from its RNS residue vector via the CRT.
 *
 * Let c = r->size and n = size of the product M in limbs.
 *
 * Computes a = sum_i r_i * w_i mod M using the precomputed CRT
 * weights. The result is the unique value in [0, M) congruent to the
 * residues.
 *
 * Note: the start/end prints are leftover debug output.
 *
 * Complexity:
 *   Time: O(c * n^2) for the bignum multiplications, plus O(n^2) for
 *         the final reduction
 *   Auxiliary memory: O(n) limbs for temporaries
 *   Output memory: O(n) limbs
 */
void rns_to_bignum(bignum* a, const rns_num* r, ctx_rns* ctx)
{
  printf("start \n");
  bignum sum, tmp;
  bn_init_multi(&sum, &tmp, NULL);

  bn_set_u64(a, 0);

  for (u64 i = 0; i < r->size; i++) {
    bn_set_u64(&tmp, r->residues[i]);
    bn_mul(&sum, &ctx->crt_weights[i], &tmp);

    bn_add(a, a, &sum);
  }

  bn_mod(a, a, &ctx->prod);

  bn_free_multi(&sum, &tmp, NULL);

  printf("end \n");
}

/*
 * Estimate the number of primes needed for an RNS matrix determinant.
 *
 * Let n = A->r_size.
 *
 * Uses the Hadamard bound on |det(A)|, doubles it to cover the range
 * [-M, M] (so the sign can be recovered from the CRT result), and
 * returns the number of 62-bit primes whose product exceeds that
 * bound.
 *
 * Complexity:
 *   Time: O(n^3 * k^2) where k is the size of the entries in limbs
 *         (Hadamard bound), plus O(n) for the estimate
 *   Auxiliary memory: O(k) limbs
 *   Output memory: O(1)
 */
u64 rns_estimate_determinant(bigmatrix* A)
{
  bignum res;
  bn_init(&res);

  // hadamard bound
  bigmatrix_hadamard(&res, A);

  // double to get to range [-M, M]
  bn_lshift1(&res);

  u64 k = rns_estimate_primes(&res);

  bn_free(&res);

  return k;
}

/*
 * Compute the determinant of a bignum matrix via RNS.
 *
 * Let n = A->r_size and c = ctx->count.
 *
 * For each prime p_i, reduces the whole matrix mod p_i and computes
 * the determinant in the u64 Montgomery domain (parallelized over the
 * primes with OpenMP), then recombines the residue vector with the
 * CRT. If the result exceeds M / 2 it is interpreted as negative
 * (two's-complement style) and M is subtracted.
 *
 * The primes must be chosen so that 2 * |det(A)| < M (see
 * rns_estimate_determinant()).
 *
 * Complexity:
 *   Time: O(c * n^3) for the u64 determinants (parallel over c), plus
 *         O(c * n^2) for the reductions and O(c * n_b^2) for the CRT
 *         where n_b is the size of M in limbs
 *   Auxiliary memory: O(n^2) u64s per thread for the reduced matrix
 *   Output memory: O(n_b) limbs
 */
void bigmatrix_det_rns(bignum* det, const bigmatrix* A, const ctx_rns* ctx)
{
  // create array of n matrices mod p_i
  u64 n = A->r_size;
  u64* residues = malloc(sizeof(u64) * ctx->count);

#pragma omp parallel
  {
    u64 raw_size = sizeof(u64) * n * n;
    u64 aligned_size = (raw_size + 63) & ~63;
    u64* reduced_data = aligned_alloc(64, aligned_size);

#pragma omp for schedule(dynamic)
    for (u64 k = 0; k < ctx->count; k++) {
      u64 p = ctx->primes[k];
      for (u64 i = 0; i < n; i++) {
        for (u64 j = 0; j < n; j++) {
          reduced_data[i * n + j] = bn_mod_u64(GET(A, i, j), p);
        }
      }
      residues[k] = matrix_u64_det_optimized(reduced_data, n, &ctx->m_ctxs[k]);
    }
    free(reduced_data);
  }

  rns_num r;
  r.residues = residues;
  r.size = ctx->count;

  // reconstruct x in [0, M - 1]
  rns_to_bignum(det, &r, ctx);

  // if d > M / 2 det is negative
  bignum M_half;
  bn_init_multi(&M_half, NULL);
  bn_copy(&M_half, &ctx->prod);
  bn_rshift1(&M_half);

  if (bn_cmp(det, &M_half) > 0) {
    bn_sub(det, det, &ctx->prod);
  }

  bn_free(&M_half);
  free(residues);
}
