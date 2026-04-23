#include <assert.h>
#include <ctype.h>
#include <immintrin.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>

#include "../../include/bignum.h"

void bn_mont_ctx_init(bn_mont_ctx* ctx, const bignum* n)
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

  bn_init(&ctx->tmp);
  bn_alloc(&ctx->tmp, 2 * n->size + 1);
}

void bn_mont_ctx_free(bn_mont_ctx* ctx)
{
  bn_free(&ctx->one_mont);
  bn_free(&ctx->n);
  bn_free(&ctx->r_square);
  bn_free(&ctx->tmp);
}

void bn_mont_redc(bignum* r, bignum* t, bn_mont_ctx* ctx)
{
  u64 size = ctx->n.size;
  u64* n_limbs = ctx->n.limbs;
  u64* t_limbs = t->limbs;

  // Word for word reduction
  for (u64 i = 0; i < size; i++) {
    u64 m = t_limbs[i] * ctx->n_inv;
    u64 carry = 0;

    // do in assembly
    for (u64 j = 0; j < size; j++) {
      unsigned __int128 prod =
          (unsigned __int128)m * n_limbs[j] + t_limbs[i + j] + carry;
      t_limbs[i + j] = (u64)prod;
      carry = (u64)(prod >> 64);
    }

    // Handle the final carry for this row
    u64 k = i + size;
    unsigned char c =
        _addcarry_u64(0, t_limbs[k], carry, (unsigned long long*)&t_limbs[k]);
    k++;
    while (c && k < t->size) {
      c = _addcarry_u64(c, t_limbs[k], 0, (unsigned long long*)&t_limbs[k]);
      k++;
    }
  }

  if (r->capacity < size) bn_alloc(r, size);
  memcpy(r->limbs, &t_limbs[size], size * sizeof(u64));

  r->size = size;

  // Trim r before comparing so bn_cmp works accurately
  bn_trim(r);

  if (bn_cmp(r, &ctx->n) >= 0) {
    bn_sub(r, r, &ctx->n);
  }
}
void bn_mont_in(bignum* A_bar, const bignum* A, bn_mont_ctx* ctx)
{
  bn_mont_mul(A_bar, A, &ctx->r_square, ctx);
}
void bn_mont_out(bignum* A, const bignum* A_bar, bn_mont_ctx* ctx)
{
  bn_mont_mul(A, A_bar, &ctx->one_mont, ctx);
}

void bn_mont_mul(bignum* r, const bignum* a_bar, const bignum* b_bar, bn_mont_ctx* ctx)
{
  bn_mont_mul_raw(r, a_bar, b_bar, ctx);
}

void bn_mont_mul_raw(bignum* result, const bignum* A_bar, const bignum* B_bar,
                     bn_mont_ctx* ctx)
{
  bignum* T = &ctx->tmp;

  u64 req_size = 2 * ctx->n.size + 1;
  // should already be big enough
  if (T->capacity < req_size) {
    bn_alloc(T, req_size);
  }

  bn_mul(T, A_bar, B_bar);

  if (T->size < req_size) {
    for (u64 i = T->size; i < req_size; i++) {
      T->limbs[i] = 0;
    }
    T->size = req_size;
  }

  bn_mont_redc(result, T, ctx);
}
