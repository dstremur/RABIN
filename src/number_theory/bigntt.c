#include "../../include/bigntt.h"

#include <stdio.h>

#include "../../include/bigvector.h"
#include "../../include/u64.h"

static inline bool is_pow_2(u64 n) { return n && !(n & (n - 1)); }

static inline u64 ilog2_u64(u64 x)
{
  u64 r = 0;
  while (x >>= 1) r++;
  return r;
}

static u64 reverse_bits(u64 x, u64 bits)
{
  u64 r = 0;

  for (u64 i = 0; i < bits; i++) {
    r <<= 1;
    r |= (x & 1);
    x >>= 1;
  }

  return r;
}

bool bn_mod_inverse(bignum* res, const bignum* a, const bignum* m)
{
  // If a is 0 or modulus is <= 1, no inverse exists
  if (bn_is_zero(a) || bn_is_zero(m) || bn_cmp(m, &BN_ONE) == 0) {
    return false;
  }

  bignum u, v, x1, x2;
  bn_init(&u);
  bn_init(&v);
  bn_init(&x1);
  bn_init(&x2);

  bn_copy(&u, a);
  bn_copy(&v, m);
  bn_set_u64(&x1, 1);
  bn_set_u64(&x2, 0);

  while (!bn_is_zero(&u) && !bn_is_zero(&v)) {
    // Eliminate powers of 2 in u
    while (bn_is_even(&u)) {
      bn_rshift1(&u);  // u = u / 2
      if (bn_is_even(&x1)) {
        bn_rshift1(&x1);
      } else {
        bn_add(&x1, &x1, m);
        bn_rshift1(&x1);  // x1 = (x1 + m) / 2
      }
    }

    // Eliminate powers of 2 in v
    while (bn_is_even(&v)) {
      bn_rshift1(&v);  // v = v / 2
      if (bn_is_even(&x2)) {
        bn_rshift1(&x2);
      } else {
        bn_add(&x2, &x2, m);
        bn_rshift1(&x2);  // x2 = (x2 + m) / 2
      }
    }

    // Step-down subtraction
    if (bn_cmp(&u, &v) >= 0) {
      bn_sub(&u, &u, &v);
      // Simulating signed subtraction under unsigned bignum bounds
      if (bn_cmp(&x1, &x2) < 0) {
        bn_add(&x1, &x1, m);
      }
      bn_sub(&x1, &x1, &x2);
    } else {
      bn_sub(&v, &v, &u);
      if (bn_cmp(&x2, &x1) < 0) {
        bn_add(&x2, &x2, m);
      }
      bn_sub(&x2, &x2, &x1);
    }
  }

  bool success = false;
  // If gcd is 1, the matching variable holds the modular inverse
  if (bn_cmp(&u, &BN_ONE) == 0) {
    bn_copy(res, &x1);
    success = true;
  } else if (bn_cmp(&v, &BN_ONE) == 0) {
    bn_copy(res, &x2);
    success = true;
  }

  bn_free(&u);
  bn_free(&v);
  bn_free(&x1);
  bn_free(&x2);

  return success;
}

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
  bigvector_println(&factors);

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

  bigvector_println(&factors);

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

  u64 n = (1ULL << k);
  ctx->k = k;
  ctx->n = n;
  bn_copy(&ctx->q, &p);

  ctx->omega_powers = NULL;
  ctx->omega_inv_powers = NULL;
  ctx->psi_powers = NULL;
  ctx->psi_inv_powers = NULL;
  ctx->bit_rev_indices = NULL;
  bn_init(&ctx->n_inv);

  ctx->omega_powers = malloc(sizeof(bignum) * n);
  ctx->omega_inv_powers = malloc(sizeof(bignum) * n);
  ctx->bit_rev_indices = malloc(sizeof(u64) * n);

  for (u64 i = 0; i < n; i++) {
    bn_init(&ctx->omega_powers[i]);
    bn_init(&ctx->omega_inv_powers[i]);
  }

  // precompute w^i mod q
  bn_set_u64(&ctx->omega_powers[0], 1);

  for (u64 i = 1; i < n; i++) {
    bn_mul(&ctx->omega_powers[i], &ctx->omega_powers[i - 1], &omega);
    bn_mod(&ctx->omega_powers[i], &ctx->omega_powers[i], &ctx->q);
  }

  bignum omega_inv;
  bn_init(&omega_inv);
  if (!bn_mod_inverse(&omega_inv, &omega, &ctx->q)) {
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
    bn_mul(&ctx->psi_powers[i], &ctx->psi_powers[i - 1], &psi);
    bn_mod(&ctx->psi_powers[i], &ctx->psi_powers[i], &ctx->q);
  }

  bignum psi_inv;
  bn_init(&psi_inv);
  if (!bn_mod_inverse(&psi_inv, &psi, &ctx->q)) {
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

  bn_free_multi(&p, &g, &omega, &psi, &tmp, NULL);

  return true;

fail:
  bigntt_ctx_free(ctx);
  bn_free_multi(&p, &g, &omega, &psi, &tmp, NULL);
  return false;
}

bool bigntt_ctx_init(ntt_ctx* ctx, u64 n, const bignum* q, const bignum* omega,
                     const bignum* psi)
{
  if (!ctx || !q || !omega) return false;

  if (!is_pow_2(n)) {
    printf("transform size must be power of 2\n");
    return false;
  }

  ctx->n = n;
  bn_init(&ctx->q);
  bn_copy(&ctx->q, q);

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

  if (psi != NULL) {
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
      bn_mod(&ctx->psi_inv_powers[i], &ctx->psi_inv_powers[i], &psi_inv);
    }

    bn_free(&psi_inv);
  }

  // 6. Compute n_inv = n^-1 mod q
  bignum bn_n;
  bn_init(&bn_n);
  bn_set_u64(&bn_n, n);
  if (!bn_mod_inverse(&ctx->n_inv, &bn_n, &ctx->q)) {
    bn_free(&bn_n);
    goto fail;
  }
  bn_free(&bn_n);

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

  free(ctx->bit_rev_indices);
  ctx->bit_rev_indices = NULL;

  bn_free(&ctx->q);
  bn_free(&ctx->n_inv);
}

bool bigntt_find_prime(bignum* q, u64 n, u64 bits)
{
  // n must be a power of 2
  if (!is_pow_2(n)) {
    return false;
  }

  // q needs to be of the form
  // q = k * (2n) + 1
  u64 stride = 2 * n;

  bignum cand;
  bn_init_multi(&cand, NULL);

  // cand = k * 2^(bits - 1) + 1
  // k = ceil(2^(bits-1) / 2n)
  bn_lshift(&cand, &BN_ONE, bits - 1);

  u64 rem = bn_mod_u64(&cand, stride);

  if (rem != 1) {
    u64 delta = (stride + 1 - rem) % stride;
    bn_add_u64(&cand, &cand, delta);
  }
  while (true) {
    if (bn_bpsw(&cand)) {
      bn_copy(q, &cand);
      bn_free(&cand);
      return true;
    }

    bn_add_u64(&cand, &cand, stride);
  }

  bn_free(&cand);
  return false;
}
