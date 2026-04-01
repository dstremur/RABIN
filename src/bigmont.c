#include <assert.h>
#include <ctype.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>

#include "../include/bignum.h"

void bn_mont_ctx_init(bn_mont_ctx* ctx, bignum* n)
{
  bn_init(&ctx->n);
  bn_init(&ctx->one_mont);
  bn_init(&ctx->r_square);

  bn_alloc(&ctx->n, n->size);

  bn_copy(&ctx->n, n);

  ctx->n_inv = mod_inverse_u64(n->limbs[0]);

  u64 k = n->size;

  bn_alloc(&ctx->one_mont, k);
  ctx->one_mont.limbs[0] = 1;
  ctx->one_mont.size = 1;

  u64 bits = k * 64;

  for (u64 i = 0; i < bits; i++) {
    bn_lshift1(&ctx->one_mont);
    if (bn_cmp(&ctx->one_mont, n) >= 0) {
      bn_sub(&ctx->one_mont, &ctx->one_mont, n);
    }
  }

  bn_init(&ctx->r_square);
  bn_alloc(&ctx->r_square, k);
  bn_copy(&ctx->r_square, &ctx->one_mont);

  for (u64 i = 0; i < bits; i++) {
    bn_lshift1(&ctx->r_square);
    if (bn_cmp(&ctx->r_square, n) >= 0) {
      bn_sub(&ctx->r_square, &ctx->r_square, n);
    }
  }
}

void bn_mont_redc(bignum* r, bignum* t, bn_mont_ctx* ctx)
{
  u64 size = ctx->n.size;
  u64 inv = ctx->n_inv;
  u64* n_limbs = ctx->n.limbs;
  u64* t_limbs = t->limbs;

  // Word for word reduction
  for (u64 i = 0; i < size; i++) {
    u64 m = t->limbs[i] * ctx->n_inv;

    u64 carry = 0;
    for (u64 j = 0; j < size; j++) {
      unsigned __int128 product = (unsigned __int128)m * n_limbs[j];
      unsigned __int128 sum =
          (unsigned __int128)t_limbs[i + j] + carry + (uint64_t)product;

      t_limbs[i + j] = (uint64_t)sum;
      carry = (uint64_t)(sum >> 64) + (uint64_t)(product >> 64);
    }

    size_t k = i + size;
    while (carry > 0 && k < t->size) {
      unsigned __int128 sum = (unsigned __int128)t_limbs[k] + carry;
      t_limbs[k] = (uint64_t)sum;
      carry = (uint64_t)(sum >> 64);
      k++;
    }
  }

  for (u64 i = 0; i < size; i++) {
    r->limbs[i] = t_limbs[i + size];
  }
  r->size = size;

  if (bn_cmp(r, &ctx->n) >= 0) {
    bn_sub(r, r, &ctx->n);
  }
}

void bn_mont_in(bignum* A_bar, bignum* A, bn_mont_ctx* ctx)
{
  bignum T;

  bn_init(&T);
  bn_alloc(&T, 2 * ctx->n.size);

  // T = A * R^2
  bn_mul(&T, A, &ctx->r_square);

  // REDC(A * R^2) = (A * R^2) / R mod N = A * R mod N
  bn_mont_redc(A_bar, &T, ctx);

  bn_free(&T);
}

void bn_mont_out(bignum* A, bignum* A_bar, bn_mont_ctx* ctx)
{
  bignum T;
  bn_init(&T);
  bn_alloc(&T, ctx->n.size * 2);

  memset(T.limbs, 0, sizeof(uint64_t) * ctx->n.size * 2);

  for (size_t i = 0; i < A_bar->size; i++) {
    T.limbs[i] = A_bar->limbs[i];
  }
  T.size = ctx->n.size * 2;

  bn_mont_redc(A, &T, ctx);

  bn_free(&T);
}

void bn_mont_mul(bignum* result, bignum* A_bar, bignum* B_bar, bn_mont_ctx* ctx)
{
  // Use a static-ish approach or a workspace to debug
  bignum T;
  bn_init(&T);

  bn_alloc(&T, A_bar->size + B_bar->size);

  bn_mul(&T, A_bar, B_bar);

  // Ensure result is big enough to hold the modulus size
  if (result->capacity < ctx->n.size) {
    bn_alloc(result, ctx->n.size);
  }

  bn_mont_redc(result, &T, ctx);

  bn_free(&T);
}
