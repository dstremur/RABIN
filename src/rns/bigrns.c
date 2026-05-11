#include "../../include/bigrns.h"

#include <stdlib.h>
#include <string.h>

u64 mod_add(u64 a, u64 b, u64 p)
{
  u64 res = a + b;
  if (res >= p || res < a) {
    res -= p;
  }
  return res;
}

u64 mod_sub(u64 a, u64 b, u64 p)
{
  if (a < b) {
    return a + p - b;
  }
  return a - b;
}

u64 mod_mul(u64 a, u64 b, u64 p)
{
  unsigned __int128 res = (unsigned __int128)a * b;
  return (u64)(res % p);
}

u64 mod_pow(u64 base, u64 exp, u64 p)
{
  u64 res = 1;
  base %= p;
  while (exp > 0) {
    if (exp % 2 == 1) res = mod_mul(res, base, p);
    base = mod_mul(base, base, p);
    exp /= 2;
  }
  return res;
}

// p needs to be prime
u64 mod_inverse(u64 n, u64 p) { return mod_pow(n, p - 2, p); }

// estimate how my primes are needed
u64 bigrns_estimate_primes(const bignum* a)
{
  u64 k = bn_bit_length(a);
  u64 res = (k + 1) / 62;

  return res;
}

void rns_context_init(ctx_rns* ctx, const u64* primes, u64 count)
{
  ctx->primes = malloc(sizeof(u64) * count);
  memcpy(ctx->primes, primes, sizeof(u64) * count);
  ctx->count = count;

  bignum tmp;
  bn_init_multi(&tmp, NULL);
  bn_init(&ctx->prod);
  bn_set_u64(&ctx->prod, 1);

  // calculate produt
  for (u64 i = 0; i < count; i++) {
    bn_set_u64(&tmp, primes[i]);
    bn_mul(&ctx->prod, &ctx->prod, &tmp);
  }

  // calculate crt weights
  for (u64 i = 0; i < count; i++) {
    bn_init(&ctx->crt_weights[i]);

    bignum M_div_tmp, inv_bn;
    bn_init_multi(&M_div_tmp, &inv_bn, NULL);

    bn_set_u64(&tmp, primes[i]);
    bn_div(&M_div_tmp, &ctx->prod, &tmp);

    // find inverse
    u64 m_mod_p = bn_mod_u64(&M_div_tmp, primes[i]);
    u64 inv = mod_inverse(m_mod_p, primes[i]);

    bn_set_u64(&inv_bn, inv);

    bn_mul(&ctx->crt_weights[i], &M_div_tmp, &inv_bn);

    bn_free_multi(&M_div_tmp, &inv_bn, NULL);
  }

  bn_free_multi(&tmp, NULL);
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
    r->residues[i] = mod_mul(a->residues[i], b->residues[i], ctx->primes[i]);
  }
}

void rns_to_bignum(bignum* a, const rns_num* r, ctx_rns* ctx)
{
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
}

// Estimates primes for a simple product or single number
u64 estimate_for_value(const bignum* max_val)
{
  u64 bits = bn_bit_length(max_val);
  // Add 1 bit for safety/sign
  return (bits + 1 + 61) / 62;
}

// Estimates primes for a matrix determinant (The most robust way)
u64 estimate_for_determinant(u64 n, const bignum* max_element)
{
  double log_v = (double)bn_bit_length(max_element);
  double log_n = bn_log_2((double)n);

  // Hadamard's: (n/2 * log2(n)) + (n * log2(V))
  double total_bits = ((double)n / 2.0 * log_n) + ((double)n * log_v);

  // If result can be negative, we need 1 extra bit for the range
  total_bits += 1.0;

  return (u64)(total_bits / 62.0) + 1;
}

void rns_context_init_for_matrix(ctx_rns* ctx, const bigmatrix* A)
{
  bignum max_val;
  bn_init(&max_val);
  bigmatrix_get_max_element(&max_val, A);  // Helper to find largest entry

  u64 num_primes = estimate_for_determinant(A->rows, &max_val);

  // Generate or fetch 'num_primes' distinct 62-bit primes
  u64* prime_list = generate_primes(num_primes);

  rns_context_init(ctx, prime_list, num_primes);

  free(prime_list);
  bn_free(&max_val);
}
