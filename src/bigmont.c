#include <assert.h>
#include <ctype.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>

#include "../include/bignum.h"

void bn_mont_ctx_init(bn_mont_ctx* ctx, bignum* n) {
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

void bn_mont_redc(bignum* r, bignum* t, bn_mont_ctx* ctx) {
  u64 size = ctx->n.size;
  u64* n_limbs = ctx->n.limbs;
  u64* t_limbs = t->limbs;

  // Word for word reduction
  for (u64 i = 0; i < size; i++) {
    u64 m = t_limbs[i] * ctx->n_inv;

    u64 carry = 0;
    for (u64 j = 0; j < size; j++) {
      unsigned __int128 product = (unsigned __int128)m * n_limbs[j];
      unsigned __int128 sum =
          (unsigned __int128)t_limbs[i + j] + carry + product;

      t_limbs[i + j] = (uint64_t)sum;
      carry = (uint64_t)(sum >> 64);
    }

    size_t k = i + size;
    while (carry > 0 && k < t->size) {
      unsigned __int128 sum = (unsigned __int128)t_limbs[k] + carry;
      t_limbs[k] = (uint64_t)sum;
      carry = (uint64_t)(sum >> 64);
      k++;
    }
  }

  bn_alloc(r, size);
  for (u64 i = 0; i < size; i++) {
    r->limbs[i] = t_limbs[i + size];
  }
  r->size = size;

  // Trim r before comparing so bn_cmp works accurately
  bn_trim(r);

  if (bn_cmp(r, &ctx->n) >= 0) {
    bn_sub(r, r, &ctx->n);
  }
}
void bn_mont_in(bignum* A_bar, bignum* A, bn_mont_ctx* ctx) {
  bignum T;
  bn_init(&T);
  bn_mul(&T, A, &ctx->r_square);

  // CRITICAL FIX: Pad T to exactly 2*N + 1 limbs
  u64 req_size = 2 * ctx->n.size + 1;
  if (T.capacity < req_size) {
    bn_alloc(&T, req_size);
  }
  for (u64 i = T.size; i < req_size; i++) {
    T.limbs[i] = 0;
  }
  T.size = req_size;

  bn_mont_redc(A_bar, &T, ctx);
  bn_free(&T);
}
void bn_mont_out(bignum* A, bignum* A_bar, bn_mont_ctx* ctx) {
  bignum T;
  bn_init(&T);

  // Ensure starting size is exactly 2*N + 1
  u64 req_size = 2 * ctx->n.size + 1;
  bn_alloc(&T, req_size);
  memset(T.limbs, 0, req_size * sizeof(u64));

  for (size_t i = 0; i < A_bar->size; i++) {
    T.limbs[i] = A_bar->limbs[i];
  }
  T.size = req_size;

  bn_mont_redc(A, &T, ctx);
  bn_free(&T);
}

void bn_mont_mul(bignum* r, bignum* a_bar, bignum* b_bar, bn_mont_ctx* ctx) {
  if (r == a_bar || r == b_bar) {
    bignum tmp;
    bn_init(&tmp);
    bn_copy(&tmp, r);
    bn_mont_mul(&tmp, a_bar, b_bar, ctx);
    bn_copy(r, &tmp);
    bn_free(&tmp);
  } else {
    bn_mont_mul_raw(r, a_bar, b_bar, ctx);
  }
}

void bn_mont_mul_raw(bignum* result, bignum* A_bar, bignum* B_bar,
                     bn_mont_ctx* ctx) {
  bignum T;
  bn_init(&T);

  bn_mul(&T, A_bar, B_bar);

  // CRITICAL FIX: Pad T to exactly 2*N + 1 limbs
  u64 req_size = 2 * ctx->n.size + 1;
  if (T.capacity < req_size) {
    bn_alloc(&T, req_size);
  }
  for (u64 i = T.size; i < req_size; i++) {
    T.limbs[i] = 0;
  }
  T.size = req_size;

  bn_mont_redc(result, &T, ctx);
  bn_free(&T);
}
