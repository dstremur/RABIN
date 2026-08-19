#include <assert.h>
#include <ctype.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>

#include "../../include/bighelper.h"
#include "../../include/bignum.h"

void bn_mont_ctx_init(bn_mont_ctx* ctx, const bignum* n)
{
  bn_init(&ctx->n);
  bn_init(&ctx->one_mont);
  bn_init(&ctx->r_square);
  bn_init(&ctx->tmp);

  bn_copy(&ctx->n, n);

  // Your mod_inverse_u64 already returns -inv, so use it directly
  ctx->n_inv = mod_inverse_u64(n->limbs[0]);

  bignum two, exp_r;
  bn_init_multi(&two, &exp_r, NULL);
  bn_set_u64(&two, 2);

  // Calculate 64 * size as a native u64 first
  u64 exp_val = (u64)n->size * 64;
  bn_set_u64(&exp_r, exp_val);

  // one_mont = 2^(64 * size) mod n
  bn_mod_exp_slow(&ctx->one_mont, &two, &exp_r, n);

  // r_square = (one_mont * one_mont) mod n
  bn_mul(&ctx->r_square, &ctx->one_mont, &ctx->one_mont);
  bn_mod(&ctx->r_square, &ctx->r_square, n);

  bn_free_multi(&two, &exp_r, NULL);

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
    unsigned char c = adc64(0, t_limbs[k], carry, &t_limbs[k]);
    k++;
    while (c && k < t->size) {
      c = adc64(c, t_limbs[k], 0, &t_limbs[k]);
      k++;
    }
  }

  // Allocate space for size + 1 limbs to prevent upper-limb truncation
  u64 res_size = size + 1;
  bn_alloc(r, res_size);

  // Copy size + 1 limbs from the upper half of the buffer
  memcpy(r->limbs, &t_limbs[size], res_size * sizeof(u64));
  r->size = res_size;

  // Trim leading zeros so bn_cmp works accurately
  bn_trim(r);

  while (bn_cmp(r, &ctx->n) >= 0) {
    bn_sub_abs(r, r, &ctx->n);
  }
}

void bn_mont_in(bignum* A_bar, const bignum* A, bn_mont_ctx* ctx)
{
  bn_mont_mul(A_bar, A, &ctx->r_square, ctx);
}
void bn_mont_out(bignum* A, const bignum* A_bar, bn_mont_ctx* ctx)
{
  bignum one;
  bn_init(&one);
  bn_set_u64(&one, 1);

  bn_mont_mul(A, A_bar, &one, ctx);

  bn_free(&one);
}

void bn_mont_mul(bignum* r, const bignum* a_bar, const bignum* b_bar,
                 bn_mont_ctx* ctx)
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
