#include "../../include/bigrns.h"
#include "../../include/u64.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>



// estimate how my primes are needed
u64 rns_estimate_primes(const bignum* a)
{
  u64 k = bn_bit_length(a);
  u64 res = (k + 62) / 62;

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
    u64 inv = mod_inverse_euclid(m_mod_p, primes[i]);

    bn_set_u64(&inv_bn, inv);

    bn_mul(&ctx->crt_weights[i], &M_div_tmp, &inv_bn);

    bn_free_multi(&M_div_tmp, &inv_bn, NULL);
  }

  bn_free_multi(&tmp, NULL);
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

// Estimates primes for a matrix determinant (The most robust way)
u64 rns_estimate_determinant(bigmatrix* A)
{
  bignum res;
  bn_init(&res);

  // hadamard bound
  bigmatrix_hadamard(&res, A);

  // double to get to range [-M, M]
  bn_lshift1(&res);

  return rns_estimate_primes(&res);
}

void bigmatrix_reduce(matrix_u64* R, const bigmatrix* A, u64 m)
{
  R->data = malloc(sizeof(u64) * A->c_size * A->r_size);
  R->c_size = A->c_size;
  R->r_size = A->r_size;
  R->modulus = m;

  for (u64 i = 0; i < A->r_size; i++) {
    for (u64 j = 0; j < A->c_size; j++) {
      R->data[i * R->c_size + j] = bn_mod_u64(GET(A, i, j), m);
    }
  }
}



void bigmatrix_det_rns(bignum* det, const bigmatrix* A, const ctx_rns* ctx)
{
  // create array of n matrices mod p_i

  u64 n = A->r_size;
  u64* residues = malloc(sizeof(u64) * ctx->count);

#pragma omp parallel for schedule(static)
  for (u64 i = 0; i < ctx->count; i++) {
    u64 p = ctx->primes[i];

    barrett_ctx b_ctx;
    barrett_init(&b_ctx, p);


    u64* reduced_data = malloc(sizeof(u64) * n * n);

#pragma omp simd
    for (u64 r = 0; r < n; r++) {
#pragma omp simd
      for (u64 c = 0; c < n; c++) {
        reduced_data[r * n + c] = bn_mod_u64(GET(A, r, c), p);
      }
    }
    residues[i] = matrix_u64_det_optimized(reduced_data, n, &b_ctx);
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
