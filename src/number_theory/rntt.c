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
  if (err != RABIN_SUCCESS) goto out;

start:
  if ((err = rz_gen_random_range(&a, &RZ_TWO, &n_min_1)) != RABIN_SUCCESS)
    goto out;

  for (u64 i = 0; i < factors->size; i++) {
    if ((err = rz_div(&exp, &n_min_1, &factors->data[i])) != RABIN_SUCCESS)
      goto out;
    if ((err = rz_mod_exp(&b, &a, &exp, p)) != RABIN_SUCCESS) goto out;
    if (rz_is_eq_i64(&b, 1)) {
      goto start;
    }
  }

  if ((err = rz_copy(g, &a)) != RABIN_SUCCESS) goto out;

  err = RABIN_SUCCESS;
out:
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
  if ((err = rz_sub(&n_min_1, p, &RZ_ONE)) != RABIN_SUCCESS) goto out;

  if ((err = rz_factorize(&factors, &n_min_1)) != RABIN_SUCCESS) goto out;

  if ((err = rz_find_gen(g, p, &factors)) != RABIN_SUCCESS) goto out;

  err = RABIN_SUCCESS;
out:
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
  if ((err = rvec_append(&factors, &RZ_TWO)) != RABIN_SUCCESS) goto out;

  // factorize c
  if ((err = rz_factorize(&factors, c)) != RABIN_SUCCESS) goto out;

  // find a generator
  if ((err = rz_find_gen(g, p, &factors)) != RABIN_SUCCESS) goto out;

  err = RABIN_SUCCESS;
out:
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

  rabin_err_t err = rz_lshift(&two_k, &RZ_ONE, k);
  if (err != RABIN_SUCCESS) goto out;

  // make c odd
  u64 curr_c = (c & 1) ? c : c + 1;

  while (1) {
    // p = curr_c * 2^k
    if ((err = rz_set_u64(&c_bn, curr_c)) != RABIN_SUCCESS) goto out;
    if ((err = rz_mul(&p_bn, &two_k, &c_bn)) != RABIN_SUCCESS) goto out;

    // p = p + 1
    if ((err = rz_add(&p_bn, &p_bn, &RZ_ONE)) != RABIN_SUCCESS) goto out;

    if (rz_bpsw(&p_bn)) {
      break;
    }

    curr_c += 2;
  }

  // copy prime
  if ((err = rz_copy(p, &p_bn)) != RABIN_SUCCESS) goto out;

  // factorize p - 1 = c * 2^k
  rvec_t factors = {0};
  if ((err = rvec_init_dynamic(&factors)) != RABIN_SUCCESS) goto out;

  // 2 always divides p - 1
  if ((err = rvec_append(&factors, &RZ_TWO)) != RABIN_SUCCESS) goto out_f;

  // factorize c
  if ((err = rz_factorize(&factors, &c_bn)) != RABIN_SUCCESS) goto out_f;

  // now find a generator
  if ((err = rz_find_gen(g, &p_bn, &factors)) != RABIN_SUCCESS) goto out_f;

  // psi = g^c mod p
  if ((err = rz_mod_exp(psi, g, &c_bn, p)) != RABIN_SUCCESS) goto out_f;

  // omega = psi^2 mod p
  if ((err = rz_mul(&two_k, psi, psi)) != RABIN_SUCCESS) goto out_f;
  if ((err = rz_mod(omega, &two_k, p)) != RABIN_SUCCESS) goto out_f;

  err = RABIN_SUCCESS;
out_f:
  rvec_clear(&factors);
out:
  rz_clear_multi(&p_bn, &c_bn, &two_k, NULL);
  return err;
}

rabin_err_t rntt_ctx_init_simple(rntt_ctx_t* ctx, u64 k, u64 c)
{
  if (ctx == NULL) return RABIN_ERR_NULL_PTR;

  rz_t p, g, omega, psi, tmp;
  rz_init_multi(&p, &g, &omega, &psi, &tmp, NULL);

  rabin_err_t err = rz_gen_proth_ntt(&g, &p, &omega, &psi, k + 1, c);
  if (err != RABIN_SUCCESS) goto out;

  if ((err = rntt_ctx_init(ctx, &p, &g, &omega, &psi, k, c)) != RABIN_SUCCESS)
    goto out;

  err = RABIN_SUCCESS;
out:
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
  if (err != RABIN_SUCCESS) goto out;
  if ((err = rz_set_u64(&c_bn, 5)) != RABIN_SUCCESS) goto out;

  // factorize p - 1 = 5 * 2^55
  rvec_t factors = {0};
  if ((err = rvec_init_dynamic(&factors)) != RABIN_SUCCESS) goto out;

  // 2,5 always divide p - 1
  if ((err = rvec_append(&factors, &RZ_TWO)) != RABIN_SUCCESS) goto out_f;
  if ((err = rvec_append(&factors, &c_bn)) != RABIN_SUCCESS) goto out_f;

  if ((err = rz_set_u64(&g, 3)) != RABIN_SUCCESS) goto out_f;
  // psi = g^(p - 1 / 2^k+1) = c * 2^(55 - (k+1)) mod p
  if ((err = rz_set_u64(&tmp, 5)) != RABIN_SUCCESS) goto out_f;
  if ((err = rz_lshift(&tmp, &tmp, 55 - (k + 1))) != RABIN_SUCCESS) goto out_f;
  if ((err = rz_mod_exp(&psi, &g, &tmp, &p)) != RABIN_SUCCESS) goto out_f;

  // omega = psi^2 mod p
  if ((err = rz_mul(&tmp, &psi, &psi)) != RABIN_SUCCESS) goto out_f;
  if ((err = rz_mod(&omega, &tmp, &p)) != RABIN_SUCCESS) goto out_f;

  if ((err = rntt_ctx_init(ctx, &p, &g, &omega, &psi, k, 5)) != RABIN_SUCCESS)
    goto out_f;

  err = RABIN_SUCCESS;
out_f:
  rvec_clear(&factors);
out:
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

  rz_init(&ctx->q);
  rabin_err_t err = rz_copy(&ctx->q, p);
  if (err != RABIN_SUCCESS) goto fail;

  // init montgomery context
  if ((err = rz_mont_ctx_init(&ctx->mctx, p)) != RABIN_SUCCESS) goto fail;

  rz_init(&ctx->n_inv);

  if (n > SIZE_MAX / sizeof(rz_t) || n > SIZE_MAX / sizeof(u64)) {
    err = RABIN_ERR_OVERFLOW;
    goto fail;
  }
  ctx->omega_powers = malloc(sizeof(rz_t) * n);
  ctx->omega_inv_powers = malloc(sizeof(rz_t) * n);
  ctx->bit_rev_indices = malloc(sizeof(u64) * n);

  if (!ctx->omega_powers || !ctx->omega_inv_powers || !ctx->bit_rev_indices) {
    err = RABIN_ERR_OUT_OF_MEMORY;
    goto fail;
  }

  for (u64 i = 0; i < n; i++) {
    rz_init(&ctx->omega_powers[i]);
    rz_init(&ctx->omega_inv_powers[i]);
  }

  // precompute w^i mod q
  if ((err = rz_set_u64(&ctx->omega_powers[0], 1)) != RABIN_SUCCESS) goto fail;

  for (u64 i = 1; i < n; i++) {
    if ((err = rz_mul(&ctx->omega_powers[i], &ctx->omega_powers[i - 1],
                      omega)) != RABIN_SUCCESS)
      goto fail;
    if ((err = rz_mod(&ctx->omega_powers[i], &ctx->omega_powers[i], &ctx->q)) !=
        RABIN_SUCCESS)
      goto fail;
  }

  rz_t omega_inv;
  rz_init(&omega_inv);
  if ((err = rz_mod_inverse(&omega_inv, omega, &ctx->q)) != RABIN_SUCCESS) {
    rz_clear(&omega_inv);
    goto fail;
  }

  if ((err = rz_set_u64(&ctx->omega_inv_powers[0], 1)) != RABIN_SUCCESS)
    goto fail_inv;
  for (u64 i = 1; i < n; i++) {
    if ((err = rz_mul(&ctx->omega_inv_powers[i], &ctx->omega_inv_powers[i - 1],
                      &omega_inv)) != RABIN_SUCCESS)
      goto fail_inv;
    if ((err = rz_mod(&ctx->omega_inv_powers[i], &ctx->omega_inv_powers[i],
                      &ctx->q)) != RABIN_SUCCESS)
      goto fail_inv;
  }

  rz_clear(&omega_inv);

  ctx->psi_powers = malloc(sizeof(rz_t) * n);
  ctx->psi_inv_powers = malloc(sizeof(rz_t) * n);
  if (!ctx->psi_powers || !ctx->psi_inv_powers) {
    err = RABIN_ERR_OUT_OF_MEMORY;
    goto fail;
  }

  for (u64 i = 0; i < n; i++) {
    rz_init(&ctx->psi_powers[i]);
    rz_init(&ctx->psi_inv_powers[i]);
  }

  // compute psi^i mod q
  if ((err = rz_set_u64(&ctx->psi_powers[0], 1)) != RABIN_SUCCESS) goto fail;
  for (u64 i = 1; i < n; i++) {
    if ((err = rz_mul(&ctx->psi_powers[i], &ctx->psi_powers[i - 1], psi)) !=
        RABIN_SUCCESS)
      goto fail;
    if ((err = rz_mod(&ctx->psi_powers[i], &ctx->psi_powers[i], &ctx->q)) !=
        RABIN_SUCCESS)
      goto fail;
  }

  rz_t psi_inv;
  rz_init(&psi_inv);
  if ((err = rz_mod_inverse(&psi_inv, psi, &ctx->q)) != RABIN_SUCCESS) {
    rz_clear(&psi_inv);
    goto fail;
  }

  if ((err = rz_set_u64(&ctx->psi_inv_powers[0], 1)) != RABIN_SUCCESS)
    goto fail_psi_inv;
  for (u64 i = 1; i < n; i++) {
    if ((err = rz_mul(&ctx->psi_inv_powers[i], &ctx->psi_inv_powers[i - 1],
                      &psi_inv)) != RABIN_SUCCESS)
      goto fail_psi_inv;
    if ((err = rz_mod(&ctx->psi_inv_powers[i], &ctx->psi_inv_powers[i],
                      &ctx->q)) != RABIN_SUCCESS)
      goto fail_psi_inv;
  }

  rz_clear(&psi_inv);

  // 6. Compute n_inv = n^-1 mod q
  rz_t rz_n;
  rz_init(&rz_n);
  if ((err = rz_set_u64(&rz_n, n)) != RABIN_SUCCESS) {
    rz_clear(&rz_n);
    goto fail;
  }
  if ((err = rz_mod_inverse(&ctx->n_inv, &rz_n, &ctx->q)) != RABIN_SUCCESS) {
    rz_clear(&rz_n);
    goto fail;
  }
  rz_clear(&rz_n);

  for (u64 i = 0; i < n; i++) {
    if ((err = rz_mont_in(&ctx->omega_powers[i], &ctx->omega_powers[i],
                          &ctx->mctx)) != RABIN_SUCCESS)
      goto fail;
    if ((err = rz_mont_in(&ctx->omega_inv_powers[i], &ctx->omega_inv_powers[i],
                          &ctx->mctx)) != RABIN_SUCCESS)
      goto fail;
    if ((err = rz_mont_in(&ctx->psi_powers[i], &ctx->psi_powers[i],
                          &ctx->mctx)) != RABIN_SUCCESS)
      goto fail;
    if ((err = rz_mont_in(&ctx->psi_inv_powers[i], &ctx->psi_inv_powers[i],
                          &ctx->mctx)) != RABIN_SUCCESS)
      goto fail;
  }
  if ((err = rz_mont_in(&ctx->n_inv, &ctx->n_inv, &ctx->mctx)) != RABIN_SUCCESS)
    goto fail;

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

  err = RABIN_SUCCESS;
  return err;

fail_inv:
  rz_clear(&omega_inv);
  goto fail;
fail_psi_inv:
  rz_clear(&psi_inv);
fail:
  rntt_ctx_clear(ctx);
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
  rabin_err_t err = rpol_init(&res);
  if (err != RABIN_SUCCESS) return err;
  if ((err = rpol_alloc(&res, n)) != RABIN_SUCCESS) goto out_res;

  // copy and reverse bits and into montgomery form
  for (u64 i = 0; i < n; i++) {
    u64 rev = ctx->bit_rev_indices[i];
    if (i <= a->deg) {
      if ((err = rz_mont_in(&res.coeff[rev], &a->coeff[i], &ctx->mctx)) !=
          RABIN_SUCCESS)
        goto out_res;
    } else {
      if ((err = rz_set_u64(&res.coeff[rev], 0)) != RABIN_SUCCESS) goto out_res;
    }
  }

  // cooley tukey butterfly

  rz_t t, u;
  rz_init_multi(&t, &u, NULL);

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
        if ((err = rz_mont_mul(&t, &res.coeff[odd], w, &ctx->mctx)) !=
            RABIN_SUCCESS)
          goto out;

        // u = res.data[even]
        if ((err = rz_copy(&u, &res.coeff[even])) != RABIN_SUCCESS) goto out;

        // res[even] = (u + t) mod q
        if ((err = rz_add(&res.coeff[even], &u, &t)) != RABIN_SUCCESS) goto out;
        // rz_mod(&res.coeff[even], &res.coeff[even], &ctx->q);

        if (rz_cmp(&res.coeff[even], &ctx->q) >= 0) {
          if ((err = rz_sub(&res.coeff[even], &res.coeff[even], &ctx->q)) !=
              RABIN_SUCCESS)
            goto out;
        }

        // res[odd] = (u - t) mod q
        if ((err = rz_sub(&res.coeff[odd], &u, &t)) != RABIN_SUCCESS) goto out;
        if ((err = rz_add(&res.coeff[odd], &res.coeff[odd], &ctx->q)) !=
            RABIN_SUCCESS)
          goto out;
        // rz_mod(&res.coeff[odd], &res.coeff[odd], &ctx->q);
        if (rz_cmp(&res.coeff[odd], &ctx->q) >= 0) {
          if ((err = rz_sub(&res.coeff[odd], &res.coeff[odd], &ctx->q)) !=
              RABIN_SUCCESS)
            goto out;
        }
      }
    }
  }

out:
  rz_clear_multi(&t, &u, NULL);
  rpol_clear(a_hat);
  if (err != RABIN_SUCCESS) goto out_res;
  if ((err = rpol_init(a_hat)) != RABIN_SUCCESS) goto out_res;
  if ((err = rpol_alloc(a_hat, n)) != RABIN_SUCCESS) goto out_res;

  // rpol_copy(a_hat, &res);

  for (u64 i = 0; i < n; i++) {
    if ((err = rz_copy(&a_hat->coeff[i], &res.coeff[i])) != RABIN_SUCCESS)
      goto out_res;
  }

  a_hat->deg = n - 1;

  err = RABIN_SUCCESS;
out_res:
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
  rabin_err_t err = rpol_init(&res);
  if (err != RABIN_SUCCESS) return err;
  if ((err = rpol_alloc(&res, n)) != RABIN_SUCCESS) goto out_res;

  // copy and reverse bits
  for (u64 i = 0; i < n; i++) {
    u64 rev = ctx->bit_rev_indices[i];
    if (i <= a->deg) {
      if ((err = rz_copy(&res.coeff[rev], &a->coeff[i])) != RABIN_SUCCESS)
        goto out_res;
    } else {
      if ((err = rz_set_u64(&res.coeff[rev], 0)) != RABIN_SUCCESS) goto out_res;
    }
  }
  // cooley tukey butterfly

  rz_t t, u;
  rz_init_multi(&t, &u, NULL);

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
        if ((err = rz_mont_mul(&t, &res.coeff[odd], w, &ctx->mctx)) !=
            RABIN_SUCCESS)
          goto out;

        // u = res.data[even]
        if ((err = rz_copy(&u, &res.coeff[even])) != RABIN_SUCCESS) goto out;

        // res[even] = (u + t) mod q
        if ((err = rz_add(&res.coeff[even], &u, &t)) != RABIN_SUCCESS) goto out;
        // rz_mod(&res.coeff[even], &res.coeff[even], &ctx->q);

        if (rz_cmp(&res.coeff[even], &ctx->q) >= 0) {
          if ((err = rz_sub(&res.coeff[even], &res.coeff[even], &ctx->q)) !=
              RABIN_SUCCESS)
            goto out;
        }

        // res[odd] = (u - t) mod q
        if ((err = rz_sub(&res.coeff[odd], &u, &t)) != RABIN_SUCCESS) goto out;
        if ((err = rz_add(&res.coeff[odd], &res.coeff[odd], &ctx->q)) !=
            RABIN_SUCCESS)
          goto out;
        // rz_mod(&res.coeff[odd], &res.coeff[odd], &ctx->q);
        if (rz_cmp(&res.coeff[odd], &ctx->q) >= 0) {
          if ((err = rz_sub(&res.coeff[odd], &res.coeff[odd], &ctx->q)) !=
              RABIN_SUCCESS)
            goto out;
        }
      }
    }
  }

  // scale by n^-1 mod q
  for (u64 i = 0; i < n; i++) {
    if ((err = rz_mont_mul(&res.coeff[i], &res.coeff[i], &ctx->n_inv,
                           &ctx->mctx)) != RABIN_SUCCESS)
      goto out;
    // rz_mod(&res.coeff[i], &res.coeff[i], &ctx->q);

    if ((err = rz_mont_out(&res.coeff[i], &res.coeff[i], &ctx->mctx)) !=
        RABIN_SUCCESS)
      goto out;
  }

out:
  rz_clear_multi(&t, &u, NULL);
  rpol_clear(a_hat);
  if (err != RABIN_SUCCESS) goto out_res;
  if ((err = rpol_init(a_hat)) != RABIN_SUCCESS) goto out_res;
  if ((err = rpol_alloc(a_hat, n)) != RABIN_SUCCESS) goto out_res;

  for (u64 i = 0; i < n; i++) {
    if ((err = rz_copy(&a_hat->coeff[i], &res.coeff[i])) != RABIN_SUCCESS)
      goto out_res;
  }
  // rpol_copy(a_hat, &res);

  a_hat->deg = n - 1;

  err = RABIN_SUCCESS;
out_res:
  rpol_clear(&res);
  return err;
}
