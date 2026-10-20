/*
 * rntt.c
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

#include "../../include/rntt.h"

#include <stdio.h>
#include <stdlib.h>

#include "../../include/rvec.h"
#include "../../include/rzfactor.h"
#include "../../include/rzprime.h"
#include "../../include/rzrand.h"
#include "../../include/u64.h"

rabin_err_t rz_find_gen(rz_t* g, const rz_t* p, const rvec_t* factors)
{
  if (g == NULL || p == NULL || factors == NULL) return RABIN_ERR_NULL_PTR;

  rz_t a, b, exp, n_min_1;
  rz_init_multi(&a, &b, &exp, &n_min_1, NULL);
  rabin_err_t err = rz_sub(&n_min_1, p, &RZ_ONE);
  if (err != RABIN_SUCCESS) goto cleanup;

start:
  TRY(rz_gen_random_range(&a, &RZ_TWO, &n_min_1));

  for (u64 i = 0; i < factors->size; i++) {
    TRY(rz_div(&exp, &n_min_1, &factors->data[i]));
    TRY(rz_mod_exp(&b, &a, &exp, p));
    if (rz_is_eq_i64(&b, 1)) {
      goto start;
    }
  }

  TRY(rz_copy(g, &a));

  err = RABIN_SUCCESS;
cleanup:
  rz_clear_multi(&a, &b, &exp, &n_min_1, NULL);
  return err;
}

rabin_err_t rz_find_gen_fp(rz_t* g, const rz_t* p)
{
  if (g == NULL || p == NULL) return RABIN_ERR_NULL_PTR;

  rvec_t factors = {0};
  rabin_err_t err = rvec_init_dynamic(&factors);
  if (err != RABIN_SUCCESS) return err;

  rz_t n_min_1;
  rz_init_multi(&n_min_1, NULL);
  TRY(rz_sub(&n_min_1, p, &RZ_ONE));

  TRY(rz_factorize(&factors, &n_min_1));

  TRY(rz_find_gen(g, p, &factors));

  err = RABIN_SUCCESS;
cleanup:
  rvec_clear(&factors);
  rz_clear_multi(&n_min_1, NULL);
  return err;
}

rabin_err_t rz_find_gen_proth(rz_t* g, const rz_t* p, rz_t* c)
{
  if (g == NULL || p == NULL || c == NULL) return RABIN_ERR_NULL_PTR;

  rvec_t factors = {0};
  rabin_err_t err = rvec_init_dynamic(&factors);
  if (err != RABIN_SUCCESS) return err;

  // 2 is always a prime factor
  TRY(rvec_append(&factors, &RZ_TWO));

  // factorize c
  TRY(rz_factorize(&factors, c));

  // find a generator
  TRY(rz_find_gen(g, p, &factors));

  err = RABIN_SUCCESS;
cleanup:
  rvec_clear(&factors);
  return err;
}

rabin_err_t rz_gen_proth_ntt(rz_t* g, rz_t* p, rz_t* omega, rz_t* psi, u64 k,
                             u64 c)
{
  if (g == NULL || p == NULL || omega == NULL || psi == NULL)
    return RABIN_ERR_NULL_PTR;
  if (k >= 64) return RABIN_ERR_INVALID_ARG;

  rz_t p_bn, c_bn, two_k;
  rz_init_multi(&p_bn, &c_bn, &two_k, NULL);
  rvec_t factors = {0};

  rabin_err_t err = rz_lshift(&two_k, &RZ_ONE, k);
  if (err != RABIN_SUCCESS) goto cleanup;

  // make c odd
  u64 curr_c = (c & 1) ? c : c + 1;

  while (1) {
    // p = curr_c * 2^k
    TRY(rz_set_u64(&c_bn, curr_c));
    TRY(rz_mul(&p_bn, &two_k, &c_bn));

    // p = p + 1
    TRY(rz_add(&p_bn, &p_bn, &RZ_ONE));

    if (rz_bpsw(&p_bn)) {
      break;
    }

    curr_c += 2;
  }

  // copy prime
  TRY(rz_copy(p, &p_bn));

  // factorize p - 1 = c * 2^k
  TRY(rvec_init_dynamic(&factors));

  // 2 always divides p - 1
  TRY(rvec_append(&factors, &RZ_TWO));

  // factorize c
  TRY(rz_factorize(&factors, &c_bn));

  // now find a generator
  TRY(rz_find_gen(g, &p_bn, &factors));

  // psi = g^c mod p
  TRY(rz_mod_exp(psi, g, &c_bn, p));

  // omega = psi^2 mod p
  TRY(rz_mul(&two_k, psi, psi));
  TRY(rz_mod(omega, &two_k, p));

cleanup:
  rvec_clear(&factors);
  rz_clear_multi(&p_bn, &c_bn, &two_k, NULL);
  return err;
}

rabin_err_t rntt_ctx_init_simple(rntt_ctx_t* ctx, u64 k, u64 c)
{
  if (ctx == NULL) return RABIN_ERR_NULL_PTR;

  rz_t p, g, omega, psi, tmp;
  rz_init_multi(&p, &g, &omega, &psi, &tmp, NULL);

  rabin_err_t err = rz_gen_proth_ntt(&g, &p, &omega, &psi, k + 1, c);
  if (err != RABIN_SUCCESS) goto cleanup;

  TRY(rntt_ctx_init(ctx, &p, &g, &omega, &psi, k, c));

  err = RABIN_SUCCESS;
cleanup:
  rz_clear_multi(&p, &g, &omega, &psi, &tmp, NULL);
  return err;
}

rabin_err_t rntt_ctx_init_golden(rntt_ctx_t* ctx, u64 k)
{
  if (ctx == NULL) return RABIN_ERR_NULL_PTR;
  if (k > 54) return RABIN_ERR_INVALID_ARG;

  rz_t p, g, omega, psi, tmp, c_bn;
  rz_init_multi(&p, &g, &omega, &psi, &tmp, &c_bn, NULL);

  // p = 5* 2^55 + 1
  rabin_err_t err = rz_set_u64(&p, 180143985094819841);
  if (err != RABIN_SUCCESS) goto cleanup;
  TRY(rz_set_u64(&c_bn, 5));

  // factorize p - 1 = 5 * 2^55
  rvec_t factors = {0};
  TRY(rvec_init_dynamic(&factors));

  // 2,5 always divide p - 1
  TRY(rvec_append(&factors, &RZ_TWO));
  TRY(rvec_append(&factors, &c_bn));

  TRY(rz_set_u64(&g, 3));
  // psi = g^(p - 1 / 2^k+1) = c * 2^(55 - (k+1)) mod p
  TRY(rz_set_u64(&tmp, 5));
  TRY(rz_lshift(&tmp, &tmp, 55 - (k + 1)));
  TRY(rz_mod_exp(&psi, &g, &tmp, &p));

  // omega = psi^2 mod p
  TRY(rz_mul(&tmp, &psi, &psi));
  TRY(rz_mod(&omega, &tmp, &p));

  TRY(rntt_ctx_init(ctx, &p, &g, &omega, &psi, k, 5));

cleanup:
  rvec_clear(&factors);
  rz_clear_multi(&p, &g, &omega, &psi, &tmp, &c_bn, NULL);
  return err;
}

bool rz_check_ntt_safety(u64 ntt_size, u64 bit_width, const rz_t* p)
{
  // Max value of a resulting coefficient is roughly ntt_size * (2^bit_width -
  // 1)^2 For bit_width = 16 and ntt_size = 2^20, max_coeff is ~2^52. Our prime
  // is ~2^{57.3}, so this is safe.
  rz_t max_val, base_minus_1;
  rz_init_multi(&max_val, &base_minus_1, NULL);

  rz_lshift(&base_minus_1, &RZ_ONE, bit_width);
  rz_sub(&base_minus_1, &base_minus_1, &RZ_ONE);  // (2^W - 1)

  rz_mul(&max_val, &base_minus_1, &base_minus_1);  // (2^W - 1)^2

  rz_t n_bn;
  rz_init(&n_bn);
  rz_set_u64(&n_bn, ntt_size);
  rz_mul(&max_val, &max_val, &n_bn);  // N * (2^W - 1)^2

  bool safe = (rz_cmp(&max_val, p) < 0);

  rz_clear_multi(&max_val, &base_minus_1, &n_bn, NULL);
  return safe;
}

rabin_err_t rntt_ctx_init(rntt_ctx_t* ctx, const rz_t* p, const rz_t* g,
                          const rz_t* omega, const rz_t* psi, u64 k, u64 c)
{
  (void)g;  // accepted for interface compatibility
  (void)c;  // accepted for interface compatibility
  if (ctx == NULL || p == NULL || omega == NULL || psi == NULL)
    return RABIN_ERR_NULL_PTR;
  if (k == 0 || k > 54) return RABIN_ERR_INVALID_ARG;

  // zero the context so rntt_ctx_clear() on the failure path never touches
  // uninitialized storage
  *ctx = (rntt_ctx_t){0};

  u64 n = (1ULL << k);
  ctx->k = k;
  ctx->n = n;

  rz_t omega_inv, psi_inv, rz_n;
  rz_init_multi(&omega_inv, &psi_inv, &rz_n, NULL);
  rabin_err_t err = RABIN_SUCCESS;

  rz_init(&ctx->q);
  TRY(rz_copy(&ctx->q, p));

  // init montgomery context
  TRY(rz_mont_ctx_init(&ctx->mctx, p));

  rz_init(&ctx->n_inv);

  if (n > SIZE_MAX / sizeof(rz_t) || n > SIZE_MAX / sizeof(u64)) {
    err = RABIN_ERR_OVERFLOW;
    goto cleanup;
  }
  ctx->omega_powers = malloc(sizeof(rz_t) * n);
  ctx->omega_inv_powers = malloc(sizeof(rz_t) * n);
  ctx->bit_rev_indices = malloc(sizeof(u64) * n);

  if (!ctx->omega_powers || !ctx->omega_inv_powers || !ctx->bit_rev_indices) {
    err = RABIN_ERR_OUT_OF_MEMORY;
    goto cleanup;
  }

  for (u64 i = 0; i < n; i++) {
    rz_init(&ctx->omega_powers[i]);
    rz_init(&ctx->omega_inv_powers[i]);
  }

  // precompute w^i mod q
  TRY(rz_set_u64(&ctx->omega_powers[0], 1));

  for (u64 i = 1; i < n; i++) {
    TRY(rz_mul(&ctx->omega_powers[i], &ctx->omega_powers[i - 1], omega));
    TRY(rz_mod(&ctx->omega_powers[i], &ctx->omega_powers[i], &ctx->q));
  }

  TRY(rz_mod_inverse(&omega_inv, omega, &ctx->q));

  TRY(rz_set_u64(&ctx->omega_inv_powers[0], 1));
  for (u64 i = 1; i < n; i++) {
    TRY(rz_mul(&ctx->omega_inv_powers[i], &ctx->omega_inv_powers[i - 1],
               &omega_inv));
    TRY(rz_mod(&ctx->omega_inv_powers[i], &ctx->omega_inv_powers[i], &ctx->q));
  }

  rz_clear(&omega_inv);

  ctx->psi_powers = malloc(sizeof(rz_t) * n);
  ctx->psi_inv_powers = malloc(sizeof(rz_t) * n);
  if (!ctx->psi_powers || !ctx->psi_inv_powers) {
    err = RABIN_ERR_OUT_OF_MEMORY;
    goto cleanup;
  }

  for (u64 i = 0; i < n; i++) {
    rz_init(&ctx->psi_powers[i]);
    rz_init(&ctx->psi_inv_powers[i]);
  }

  // compute psi^i mod q
  TRY(rz_set_u64(&ctx->psi_powers[0], 1));
  for (u64 i = 1; i < n; i++) {
    TRY(rz_mul(&ctx->psi_powers[i], &ctx->psi_powers[i - 1], psi));
    TRY(rz_mod(&ctx->psi_powers[i], &ctx->psi_powers[i], &ctx->q));
  }

  TRY(rz_mod_inverse(&psi_inv, psi, &ctx->q));

  TRY(rz_set_u64(&ctx->psi_inv_powers[0], 1));
  for (u64 i = 1; i < n; i++) {
    TRY(rz_mul(&ctx->psi_inv_powers[i], &ctx->psi_inv_powers[i - 1], &psi_inv));
    TRY(rz_mod(&ctx->psi_inv_powers[i], &ctx->psi_inv_powers[i], &ctx->q));
  }

  rz_clear(&psi_inv);

  // 6. Compute n_inv = n^-1 mod q
  TRY(rz_set_u64(&rz_n, n));
  TRY(rz_mod_inverse(&ctx->n_inv, &rz_n, &ctx->q));
  rz_clear(&rz_n);

  for (u64 i = 0; i < n; i++) {
    TRY(rz_mont_in(&ctx->omega_powers[i], &ctx->omega_powers[i], &ctx->mctx));
    TRY(rz_mont_in(&ctx->omega_inv_powers[i], &ctx->omega_inv_powers[i],
                   &ctx->mctx));
    TRY(rz_mont_in(&ctx->psi_powers[i], &ctx->psi_powers[i], &ctx->mctx));
    TRY(rz_mont_in(&ctx->psi_inv_powers[i], &ctx->psi_inv_powers[i],
                   &ctx->mctx));
  }
  TRY(rz_mont_in(&ctx->n_inv, &ctx->n_inv, &ctx->mctx));

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

cleanup:
  rz_clear_multi(&omega_inv, &psi_inv, &rz_n, NULL);
  if (err != RABIN_SUCCESS) rntt_ctx_clear(ctx);
  return err;
}

rabin_err_t rntt_ctx_clear(rntt_ctx_t* ctx)
{
  if (ctx == NULL) return RABIN_ERR_NULL_PTR;

  if (ctx->omega_powers) {
    for (u64 i = 0; i < ctx->n; i++) rz_clear(&ctx->omega_powers[i]);
    free(ctx->omega_powers);
    ctx->omega_powers = NULL;
  }
  if (ctx->omega_inv_powers) {
    for (u64 i = 0; i < ctx->n; i++) rz_clear(&ctx->omega_inv_powers[i]);
    free(ctx->omega_inv_powers);
    ctx->omega_inv_powers = NULL;
  }
  if (ctx->psi_powers) {
    for (u64 i = 0; i < ctx->n; i++) rz_clear(&ctx->psi_powers[i]);
    free(ctx->psi_powers);
    ctx->psi_powers = NULL;
  }
  if (ctx->psi_inv_powers) {
    for (u64 i = 0; i < ctx->n; i++) rz_clear(&ctx->psi_inv_powers[i]);
    free(ctx->psi_inv_powers);
    ctx->psi_inv_powers = NULL;
  }

  if (ctx->bit_rev_indices) {
    free(ctx->bit_rev_indices);
    ctx->bit_rev_indices = NULL;
  }

  rz_clear(&ctx->q);
  rz_clear(&ctx->n_inv);
  rz_mont_ctx_clear(&ctx->mctx);
  return RABIN_SUCCESS;
}

rabin_err_t rntt_cyclic_forward(rpol_t* a_hat, const rpol_t* a, rntt_ctx_t* ctx)
{
  if (a_hat == NULL || a == NULL || ctx == NULL) return RABIN_ERR_NULL_PTR;

  u64 n = ctx->n;
  rpol_t res = {0};
  rz_t t, u;
  rz_init_multi(&t, &u, NULL);
  rabin_err_t err = rpol_init(&res);
  if (err != RABIN_SUCCESS) {
    rpol_clear(&res);
    return err;
  }
  TRY(rpol_alloc(&res, n));

  // copy and reverse bits and into montgomery form
  for (u64 i = 0; i < n; i++) {
    u64 rev = ctx->bit_rev_indices[i];
    if (i <= a->deg) {
      if ((err = rz_mont_in(&res.coeff[rev], &a->coeff[i], &ctx->mctx)) !=
          RABIN_SUCCESS)
        goto cleanup;
    } else {
      if ((err = rz_set_u64(&res.coeff[rev], 0)) != RABIN_SUCCESS) goto cleanup;
    }
  }

  // cooley tukey butterfly

  for (u64 len = 2; len <= n; len <<= 1) {
    u64 half = len >> 1;
    u64 step = n / len;

    for (u64 i = 0; i < n; i += len) {
      for (u64 j = 0; j < half; j++) {
        u64 twiddle = j * step;
        rz_t* w = &ctx->omega_powers[twiddle];

        u64 even = i + j;
        u64 odd = i + j + half;

        // t = res[odd] * w mod q
        TRY(rz_mont_mul(&t, &res.coeff[odd], w, &ctx->mctx));

        // u = res.data[even]
        TRY(rz_copy(&u, &res.coeff[even]));

        // res[even] = (u + t) mod q
        TRY(rz_add(&res.coeff[even], &u, &t));
        // rz_mod(&res.coeff[even], &res.coeff[even], &ctx->q);

        if (rz_cmp(&res.coeff[even], &ctx->q) >= 0) {
          TRY(rz_sub(&res.coeff[even], &res.coeff[even], &ctx->q));
        }

        // res[odd] = (u - t) mod q
        TRY(rz_sub(&res.coeff[odd], &u, &t));
        TRY(rz_add(&res.coeff[odd], &res.coeff[odd], &ctx->q));
        // rz_mod(&res.coeff[odd], &res.coeff[odd], &ctx->q);
        if (rz_cmp(&res.coeff[odd], &ctx->q) >= 0) {
          TRY(rz_sub(&res.coeff[odd], &res.coeff[odd], &ctx->q));
        }
      }
    }
  }

cleanup:
  rz_clear_multi(&t, &u, NULL);
  if (err == RABIN_SUCCESS) {
    rpol_clear(a_hat);
    rabin_err_t terr = rpol_init(a_hat);
    if (terr == RABIN_SUCCESS) terr = rpol_alloc(a_hat, n);
    for (u64 i = 0; terr == RABIN_SUCCESS && i < n; i++)
      terr = rz_copy(&a_hat->coeff[i], &res.coeff[i]);
    if (terr == RABIN_SUCCESS) a_hat->deg = n - 1;
    err = terr;
  }
  rpol_clear(&res);
  return err;
}

rabin_err_t rntt_cyclic_inverse(rpol_t* a_hat, const rpol_t* a, rntt_ctx_t* ctx)
{
  return rntt_cyclic_inverse_mont_in(a_hat, a, ctx);
}

rabin_err_t rntt_cyclic_inverse_mont_in(rpol_t* a_hat, const rpol_t* a,
                                        rntt_ctx_t* ctx)
{
  if (a_hat == NULL || a == NULL || ctx == NULL) return RABIN_ERR_NULL_PTR;

  u64 n = ctx->n;
  rpol_t res = {0};
  rz_t t, u;
  rz_init_multi(&t, &u, NULL);
  rabin_err_t err = rpol_init(&res);
  if (err != RABIN_SUCCESS) {
    rpol_clear(&res);
    return err;
  }
  TRY(rpol_alloc(&res, n));

  // copy and reverse bits
  for (u64 i = 0; i < n; i++) {
    u64 rev = ctx->bit_rev_indices[i];
    if (i <= a->deg) {
      if ((err = rz_copy(&res.coeff[rev], &a->coeff[i])) != RABIN_SUCCESS)
        goto cleanup;
    } else {
      if ((err = rz_set_u64(&res.coeff[rev], 0)) != RABIN_SUCCESS) goto cleanup;
    }
  }
  // cooley tukey butterfly

  for (u64 len = 2; len <= n; len <<= 1) {
    u64 half = len >> 1;
    u64 step = n / len;

    for (u64 i = 0; i < n; i += len) {
      for (u64 j = 0; j < half; j++) {
        u64 twiddle = j * step;

        rz_t* w = &ctx->omega_inv_powers[twiddle];

        u64 even = i + j;
        u64 odd = i + j + half;

        // t = res[odd] * w mod q
        TRY(rz_mont_mul(&t, &res.coeff[odd], w, &ctx->mctx));

        // u = res.data[even]
        TRY(rz_copy(&u, &res.coeff[even]));

        // res[even] = (u + t) mod q
        TRY(rz_add(&res.coeff[even], &u, &t));
        // rz_mod(&res.coeff[even], &res.coeff[even], &ctx->q);

        if (rz_cmp(&res.coeff[even], &ctx->q) >= 0) {
          TRY(rz_sub(&res.coeff[even], &res.coeff[even], &ctx->q));
        }

        // res[odd] = (u - t) mod q
        TRY(rz_sub(&res.coeff[odd], &u, &t));
        TRY(rz_add(&res.coeff[odd], &res.coeff[odd], &ctx->q));
        // rz_mod(&res.coeff[odd], &res.coeff[odd], &ctx->q);
        if (rz_cmp(&res.coeff[odd], &ctx->q) >= 0) {
          TRY(rz_sub(&res.coeff[odd], &res.coeff[odd], &ctx->q));
        }
      }
    }
  }

  // scale by n^-1 mod q
  for (u64 i = 0; i < n; i++) {
    TRY(rz_mont_mul(&res.coeff[i], &res.coeff[i], &ctx->n_inv, &ctx->mctx));
    // rz_mod(&res.coeff[i], &res.coeff[i], &ctx->q);

    TRY(rz_mont_out(&res.coeff[i], &res.coeff[i], &ctx->mctx));
  }

cleanup:
  rz_clear_multi(&t, &u, NULL);
  if (err == RABIN_SUCCESS) {
    rpol_clear(a_hat);
    rabin_err_t terr = rpol_init(a_hat);
    if (terr == RABIN_SUCCESS) terr = rpol_alloc(a_hat, n);
    for (u64 i = 0; terr == RABIN_SUCCESS && i < n; i++)
      terr = rz_copy(&a_hat->coeff[i], &res.coeff[i]);
    if (terr == RABIN_SUCCESS) a_hat->deg = n - 1;
    err = terr;
  }
  rpol_clear(&res);
  return err;
}
