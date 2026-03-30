#include <assert.h>
#include <ctype.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>

#include "../include/bignum.h"

// calculates a^b into r
// binary exponentiation
void bn_pow(bignum* r, bignum* a, bignum* b) {
  if (bn_is_zero(b)) {
    bn_set_u64(r, 1);
    return;
  }

  bignum base, exp, two;
  bn_init(&base);
  bn_init(&exp);
  bn_init(&two);
  bn_set_u64(&two, 2);
  bn_copy(&base, a);
  bn_copy(&exp, b);

  bn_set_u64(r, 1);

  while (!bn_is_zero(&exp)) {
    if (!bn_is_even(&exp)) {
      bn_mul(r, r, &base);
    }
    bn_mul(&base, &base, &base);

    bn_rshift1(&exp);
  }

  bn_free(&base);
  bn_free(&exp);
  bn_free(&two);
}

void bn_mod_exp(bignum* r, bignum* a, bignum* b, bignum* m) {
  bignum base, exp, res, tmp;
  bn_init(&base);
  bn_init(&exp);
  bn_init(&res);
  bn_init(&tmp);

  bn_copy(&base, a);
  bn_copy(&exp, b);
  bn_set_u64(&res, 1);

  while (!bn_is_zero(&exp)) {
    if (!bn_is_even(&exp)) {
      bn_mul(&tmp, &res, &base);
      bn_mod(&res, &tmp, m);
    }

    bn_mul(&tmp, &base, &base);
    bn_mod(&base, &tmp, m);

    bn_rshift1(&exp);
  }

  bn_copy(r, &res);

  bn_free(&base);
  bn_free(&exp);
  bn_free(&res);
  bn_free(&tmp);
}

void bn_mod_exp_mont(bignum* r, const bignum* a, const bignum* b,
                     const bignum* m, bn_mont_ctx* ctx) {
  assert(!bn_is_zero(m) && !bn_is_even(m));

  assert(bn_cmp(&ctx->n, m) == 0);

  bignum base, result, tmp;
  bn_init_multi(&base, &result, &tmp);

  bn_mont_in(ctx, a, &base);

  bn_copy(&result, &ctx->one_mont);

  bignum exp;
  bn_init(&exp);
  bn_copy(&exp, b);

  while (!bn_is_zero(&exp)) {
    if (!bn_is_even(&exp)) {
      bn_mont_mul(ctx, &result, &base, &tmp);
      bn_copy(&result, &tmp);
    }

    bn_mont_mul(ctx, &base, &base, &tmp);
    bn_copy(&base, &tmp);

    bn_rshift1(&exp);
  }
  bn_free(&exp);

  bn_mont_out(ctx, &result, r);

  bn_free(&base);
  bn_free(&result);
  bn_free(&tmp);
}
