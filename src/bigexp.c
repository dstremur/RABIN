#include <assert.h>
#include <ctype.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>

#include "../include/bignum.h"

// calculates a^b into r
// binary exponentiation
void bn_pow(bignum* r, bignum* a, bignum* b)
{
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

// calculates a^d int r in the montgomery context
void bn_mont_exp(bignum* r_bar, bignum* a_bar, bignum* d, bn_mont_ctx* ctx)
{
  bignum base, exponent;

  bn_init_multi(&base, &exponent, NULL);
  bn_alloc(&base, ctx->n.size);
  bn_alloc(&exponent, d->size);

  bn_copy(&base, a_bar);
  bn_copy(&exponent, d);
  bn_copy(r_bar, &ctx->one_mont);

  while (!bn_is_zero(&exponent)) {
    if (exponent.limbs[0] & 1) {
      bn_mont_mul(r_bar, r_bar, &base, ctx);
    }

    bn_mont_mul(&base, &base, &base, ctx);

    bn_rshift1(&exponent);
  }

  bn_free(&base);
  bn_free(&exponent);
}

// calculates a^b mod m into r
void bn_mod_exp(bignum* r, bignum* a, bignum* b, bignum* m)
{
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

// calculates a^b mod m into r in the montgomery domain
void bn_mod_exp_mont(bignum* r, bignum* a, bignum* b, bignum* m,
                     bn_mont_ctx* ctx)
{
  assert(!bn_is_zero(m) && !bn_is_even(m));

  assert(bn_cmp(&ctx->n, m) == 0);

  bignum base, result;
  bn_init_multi(&base, &result, NULL);

  bn_mont_in(&base, a, ctx);

  bn_copy(&result, &ctx->one_mont);

  bignum exp;
  bn_init(&exp);
  bn_copy(&exp, b);

  while (!bn_is_zero(&exp)) {
    if (!bn_is_even(&exp)) {
      bn_mont_mul(&result, &result, &base, ctx);
    }

    bn_mont_mul(&base, &base, &base, ctx);

    bn_rshift1(&exp);
  }
  bn_free(&exp);

  bn_mont_out(r, &result, ctx);

  bn_free(&base);
  bn_free(&result);
}
