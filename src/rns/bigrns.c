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

u64 rns_estimate_primes(const bignum* a)
{
  u64 k = bn_bit_length(a);

  // every prime in RNS_PRIMES and RNS_PRIMES2 is >= 2^59, so the product
  // of k of them has > 59k bits; budgeting fewer than 62 bits per prime
  // would let the product fall short of a (e.g. 2 primes cover 120 bits,
  // not 124) and the CRT reconstruction would silently wrap
  u64 res = (k + 58) / 59;

  return res;
}

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
    bn_set_u64(&tmp, primes[i - 1]);
    bn_mul(&prev, &prev, &tmp);

    u64 m_prev_mod = bn_mod_u64(&prev, primes[i]);
    ctx->garner_weights[i] = mod_inverse_euclid(m_prev_mod, primes[i]);
  }

  bn_free_multi(&tmp, &prev, NULL);
}

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

void rns_add(rns_num* r, const rns_num* a, const rns_num* b, const ctx_rns* ctx)
{
  for (u64 i = 0; i < r->size; i++) {
    r->residues[i] = mod_add(a->residues[i], b->residues[i], ctx->primes[i]);
  }
}

void rns_to_bignum(bignum* a, rns_num* v, const ctx_rns* ctx)
{
  u64 t = ctx->count;
  if (t == 0) {
    bn_init(a);
    bn_set_u64(a, 0);
    return;
  }

  // Step 2: Initialize x with v_1 (using 0-based indexing: v_0)
  bn_init(a);
  bn_set_u64(a, v->residues[0]);

  bignum prod_prev, term, bn_u, bn_p;
  bn_init_multi(&prod_prev, &term, &bn_u, &bn_p, NULL);

  // Initialize the running product of primes (m_1 in the algorithm)
  bn_set_u64(&prod_prev, ctx->primes[0]);

  // Step 3: For i from 2 to t (1 to t-1 in 0-based indexing)
  for (u64 i = 1; i < t; i++) {
    // Calculate (v_i - x) mod m_i
    u64 x_mod = bn_mod_u64(a, ctx->primes[i]);
    u64 diff;
    if (v->residues[i] >= x_mod) {
      diff = v->residues[i] - x_mod;
    } else {
      diff = (v->residues[i] + ctx->primes[i]) - x_mod;
    }

    // u = (v_i - x) * C_i mod m_i
    // (Cast to __int128 to prevent overflow if primes are up to 64-bit)
    u64 u = (u64)(((unsigned __int128)diff * ctx->garner_weights[i]) %
                  ctx->primes[i]);

    // x = x + u * prod_{j=1}^{i-1} m_j
    bn_set_u64(&bn_u, u);
    bn_mul(&term, &bn_u, &prod_prev);
    bn_add(a, a, &term);

    // Update the running product for the next iteration: prod_prev *= m_i
    bn_set_u64(&bn_p, ctx->primes[i]);
    bn_mul(&prod_prev, &prod_prev, &bn_p);
  }

  bn_free_multi(&prod_prev, &term, &bn_u, &bn_p, NULL);
}

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
