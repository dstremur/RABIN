#include "../../include/bigrns.h"

#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../../include/u64.h"

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

  ctx->m_ctxs = malloc(sizeof(mont_ctx) * count);
  for (u64 i = 0; i < count; i++) {
    mont_init(&ctx->m_ctxs[i], primes[i]);
  }

  // calculate garner weights
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

// Estimates primes for a matrix determinant
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
