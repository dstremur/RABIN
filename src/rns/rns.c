/*
 * rns.c
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
 * A rns_ctx_t holds the primes, their product M, the CRT weights
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

#include "../../include/rns.h"

#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../../include/u64.h"

u64 rns_estimate_primes(const rz_t* a)
{
  if (a == NULL) return 0;
  u64 k = (u64)rz_bit_length(a);

  // every prime in RNS_PRIMES and RNS_PRIMES2 is >= 2^59, so the product
  // of k of them has > 59k bits; budgeting fewer than 62 bits per prime
  // would let the product fall short of a (e.g. 2 primes cover 120 bits,
  // not 124) and the CRT reconstruction would silently wrap
  u64 res = (k + 58) / 59;

  return res;
}

rabin_err_t rns_ctx_init(rns_ctx_t* ctx, const u64* primes, u64 count)
{
  if (ctx == NULL || primes == NULL) return RABIN_ERR_NULL_PTR;
  if (count == 0) return RABIN_ERR_INVALID_ARG;

  // zero the context so rns_ctx_clear() on the failure path never touches
  // uninitialized storage
  *ctx = (rns_ctx_t){0};

  if (count > SIZE_MAX / sizeof(u64) || count > SIZE_MAX / sizeof(rz_t) ||
      count > SIZE_MAX / sizeof(u64_mont_ctx_t)) {
    return RABIN_ERR_OVERFLOW;
  }

  ctx->primes = malloc(sizeof(u64) * count);
  ctx->crt_weights = malloc(sizeof(rz_t) * count);
  if (ctx->primes == NULL || ctx->crt_weights == NULL)
    return RABIN_ERR_OUT_OF_MEMORY;
  memcpy(ctx->primes, primes, sizeof(u64) * count);
  ctx->count = count;

  rz_t tmp, prev;
  rz_init_multi(&tmp, &prev, NULL);
  rz_init(&ctx->prod);
  rabin_err_t err = rz_set_u64(&ctx->prod, 1);
  if (err != RABIN_SUCCESS) goto out;

  // calculate product M = p_0 * ... * p_{c-1}
  for (u64 i = 0; i < count; i++) {
    if ((err = rz_set_u64(&tmp, primes[i])) != RABIN_SUCCESS) goto out;
    if ((err = rz_mul(&ctx->prod, &ctx->prod, &tmp)) != RABIN_SUCCESS) goto out;
  }

  // calculate CRT weights w_i = (M / p_i) * (M / p_i)^{-1} mod p_i
  rz_t M_div_tmp, inv_bn;
  rz_init_multi(&M_div_tmp, &inv_bn, NULL);

  for (u64 i = 0; i < count; i++) {
    rz_init(&ctx->crt_weights[i]);

    if ((err = rz_set_u64(&tmp, primes[i])) != RABIN_SUCCESS) goto crt_fail;
    if ((err = rz_div(&M_div_tmp, &ctx->prod, &tmp)) != RABIN_SUCCESS)
      goto crt_fail;

    // find inverse of (M / p_i) mod p_i
    u64 m_mod_p = rz_mod_u64(&M_div_tmp, primes[i]);
    u64 inv = u64_mod_inverse_euclid(m_mod_p, primes[i]);
    if (inv == 0) {
      err = RABIN_ERR_INVALID_ARG;
      goto crt_fail;
    }

    if ((err = rz_set_u64(&inv_bn, inv)) != RABIN_SUCCESS) goto crt_fail;
    if ((err = rz_mul(&ctx->crt_weights[i], &M_div_tmp, &inv_bn)) !=
        RABIN_SUCCESS)
      goto crt_fail;
  }

  rz_clear_multi(&M_div_tmp, &inv_bn, NULL);

  // one Montgomery context per prime
  ctx->m_ctxs = malloc(sizeof(u64_mont_ctx_t) * count);
  if (ctx->m_ctxs == NULL) {
    err = RABIN_ERR_OUT_OF_MEMORY;
    goto out;
  }
  for (u64 i = 0; i < count; i++) {
    u64_mont_init(&ctx->m_ctxs[i], primes[i]);
  }

  // calculate Garner weights g_i = (p_0 * ... * p_{i-1})^{-1} mod p_i
  ctx->garner_weights = malloc(sizeof(u64) * count);
  if (ctx->garner_weights == NULL) {
    err = RABIN_ERR_OUT_OF_MEMORY;
    goto out;
  }
  ctx->garner_weights[0] = 1;

  if ((err = rz_set_u64(&prev, 1)) != RABIN_SUCCESS) goto out;

  for (u64 i = 1; i < count; i++) {
    if ((err = rz_set_u64(&tmp, primes[i - 1])) != RABIN_SUCCESS) goto out;
    if ((err = rz_mul(&prev, &prev, &tmp)) != RABIN_SUCCESS) goto out;

    u64 m_prev_mod = rz_mod_u64(&prev, primes[i]);
    ctx->garner_weights[i] = u64_mod_inverse_euclid(m_prev_mod, primes[i]);
  }

  rz_clear_multi(&M_div_tmp, &inv_bn, NULL);
  rz_clear_multi(&tmp, &prev, NULL);
  return RABIN_SUCCESS;

crt_fail:
  rz_clear_multi(&M_div_tmp, &inv_bn, NULL);
out:
  rz_clear_multi(&tmp, &prev, NULL);
  rns_ctx_clear(ctx);
  return err;
}

rabin_err_t rns_ctx_clear(rns_ctx_t* ctx)
{
  if (!ctx) return RABIN_ERR_NULL_PTR;

  // 1. Free the individual rz_t weights
  if (ctx->crt_weights) {
    for (u64 i = 0; i < ctx->count; i++) {
      rz_clear(&ctx->crt_weights[i]);
    }
    free(ctx->crt_weights);
    ctx->crt_weights = NULL;
  }

  // 2. Free the big product
  rz_clear(&ctx->prod);

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
  return RABIN_SUCCESS;
}

rabin_err_t rns_import(rns_num_t* r, const rz_t* a, rns_ctx_t* ctx)
{
  if (r == NULL || a == NULL || ctx == NULL) return RABIN_ERR_NULL_PTR;

  if (r->residues == NULL) {
    if (ctx->count > SIZE_MAX / sizeof(u64)) return RABIN_ERR_OVERFLOW;
    r->residues = malloc(sizeof(u64) * ctx->count);
    if (r->residues == NULL) return RABIN_ERR_OUT_OF_MEMORY;
  }
  r->size = ctx->count;

  for (u64 i = 0; i < r->size; i++) {
    r->residues[i] = rz_mod_u64(a, ctx->primes[i]);
  }
  return RABIN_SUCCESS;
}

rabin_err_t rns_add(rns_num_t* r, const rns_num_t* a, const rns_num_t* b,
                    const rns_ctx_t* ctx)
{
  if (r == NULL || a == NULL || b == NULL || ctx == NULL)
    return RABIN_ERR_NULL_PTR;
  if (r->size != ctx->count || a->size != ctx->count || b->size != ctx->count)
    return RABIN_ERR_MATRIX_DIM;

  for (u64 i = 0; i < r->size; i++) {
    r->residues[i] =
        u64_mod_add(a->residues[i], b->residues[i], ctx->primes[i]);
  }
  return RABIN_SUCCESS;
}

rabin_err_t rns_export(rz_t* a, rns_num_t* v, const rns_ctx_t* ctx)
{
  if (a == NULL || v == NULL || ctx == NULL) return RABIN_ERR_NULL_PTR;

  u64 t = ctx->count;
  if (t == 0) {
    rabin_err_t err = rz_init(a);
    if (err == RABIN_SUCCESS) err = rz_set_u64(a, 0);
    return err;
  }

  // Step 2: Initialize x with v_1 (using 0-based indexing: v_0)
  rabin_err_t err = rz_init(a);
  if (err == RABIN_SUCCESS) err = rz_set_u64(a, v->residues[0]);
  if (err != RABIN_SUCCESS) return err;

  rz_t prod_prev, term, rz_u, rz_p;
  rz_init_multi(&prod_prev, &term, &rz_u, &rz_p, NULL);

  // Initialize the running product of primes (m_1 in the algorithm)
  if ((err = rz_set_u64(&prod_prev, ctx->primes[0])) != RABIN_SUCCESS) goto out;

  // Step 3: For i from 2 to t (1 to t-1 in 0-based indexing)
  for (u64 i = 1; i < t; i++) {
    // Calculate (v_i - x) mod m_i
    u64 x_mod = rz_mod_u64(a, ctx->primes[i]);
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
    if ((err = rz_set_u64(&rz_u, u)) != RABIN_SUCCESS) goto out;
    if ((err = rz_mul(&term, &rz_u, &prod_prev)) != RABIN_SUCCESS) goto out;
    if ((err = rz_add(a, a, &term)) != RABIN_SUCCESS) goto out;

    // Update the running product for the next iteration: prod_prev *= m_i
    if ((err = rz_set_u64(&rz_p, ctx->primes[i])) != RABIN_SUCCESS) goto out;
    if ((err = rz_mul(&prod_prev, &prod_prev, &rz_p)) != RABIN_SUCCESS)
      goto out;
  }

  err = RABIN_SUCCESS;
out:
  rz_clear_multi(&prod_prev, &term, &rz_u, &rz_p, NULL);
  return err;
}

u64 rns_estimate_det(const rmat_t* A)
{
  if (A == NULL) return 0;

  rz_t res;
  rz_init(&res);

  // hadamard bound
  if (rmat_hadamard(&res, A) != RABIN_SUCCESS) {
    rz_clear(&res);
    return 0;
  }

  // double to get to range [-M, M]
  if (rz_lshift1(&res) != RABIN_SUCCESS) {
    rz_clear(&res);
    return 0;
  }

  u64 k = rns_estimate_primes(&res);

  rz_clear(&res);

  return k;
}

rabin_err_t rmat_det_rns(rz_t* det, const rmat_t* A, const rns_ctx_t* ctx)
{
  if (det == NULL || A == NULL || ctx == NULL) return RABIN_ERR_NULL_PTR;
  if (A->rows != A->cols) return RABIN_ERR_MATRIX_DIM;
  if (ctx->count == 0) return RABIN_ERR_INVALID_ARG;

  // create array of n matrices mod p_i
  u64 n = A->rows;
  if (ctx->count > SIZE_MAX / sizeof(u64)) return RABIN_ERR_OVERFLOW;
  u64* residues = malloc(sizeof(u64) * ctx->count);
  if (residues == NULL) return RABIN_ERR_OUT_OF_MEMORY;

  rabin_err_t err = RABIN_SUCCESS;

  // per-thread scratch for the reduced matrix; allocation failures inside
  // the parallel region cannot be reported, so sizes are pre-checked
  if (n > 0 && n > SIZE_MAX / sizeof(u64) / n) {
    free(residues);
    return RABIN_ERR_OVERFLOW;
  }

#pragma omp parallel
  {
    u64 raw_size = sizeof(u64) * n * n;
    u64 aligned_size = (raw_size + 63) & ~63;
    // aligned_alloc() requires a positive multiple of the alignment
    u64* reduced_data = aligned_size ? aligned_alloc(64, aligned_size) : NULL;

#pragma omp for schedule(dynamic)
    for (u64 k = 0; k < ctx->count; k++) {
      u64 p = ctx->primes[k];
      for (u64 i = 0; i < n; i++) {
        for (u64 j = 0; j < n; j++) {
          reduced_data[i * n + j] = rz_mod_u64(RMAT_GET(A, i, j), p);
        }
      }
      residues[k] = u64_mat_det_optimized(reduced_data, n, &ctx->m_ctxs[k]);
    }
    free(reduced_data);
  }

  rns_num_t r;
  r.residues = residues;
  r.size = ctx->count;

  // reconstruct x in [0, M - 1]
  if ((err = rns_export(det, &r, ctx)) != RABIN_SUCCESS) goto out;

  // if d > M / 2 det is negative
  rz_t M_half;
  rz_init_multi(&M_half, NULL);
  if ((err = rz_copy(&M_half, &ctx->prod)) != RABIN_SUCCESS) {
    rz_clear(&M_half);
    goto out;
  }
  if ((err = rz_rshift1(&M_half)) != RABIN_SUCCESS) {
    rz_clear(&M_half);
    goto out;
  }

  if (rz_cmp(det, &M_half) > 0) {
    if ((err = rz_sub(det, det, &ctx->prod)) != RABIN_SUCCESS) {
      rz_clear(&M_half);
      goto out;
    }
  }

  rz_clear(&M_half);
out:
  free(residues);
  return err;
}
