#include "../../include/bigntt.h"

#include <stdio.h>

#include "../../include/bigvector.h"
#include "../../include/u64.h"

// finds a generator g of Fp given the factorization of p - 1
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

// find a generator g of Fp
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

// Find a generator g for a Proth prime p = c * 2^k + 1
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

// use 1 * 2^61 + 1
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

// output in montgomery space
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

void bigntt_cyclic_inverse(bigpoly* a_hat, bigpoly* a, ntt_ctx* ctx)
{
  bigntt_cyclic_inverse_mont_in(a_hat, a, ctx);
}

// input in montgomery form
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
