#include <assert.h>
#include <ctype.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>

#include "../../include/bignum.h"

// calculates a^b into r
// binary exponentiation
void bn_pow(bignum* r, const bignum* a, const bignum* b)
{
  if (bn_is_zero(b)) {
    bn_set_u64(r, 1);
    return;
  }

  bn_set_u64(r, 1);

  for (i64 i = bn_bit_length(b) - 1; i >= 0; i--) {
    bn_mul(r, r, r);
    if (bn_get_bit(b, i)) {
      bn_mul(r, r, a);
    }
  }
}

// calculates a^d int r in the montgomery context
void bn_mont_exp(bignum* r_bar, const bignum* a_bar, const bignum* d,
                 bn_mont_ctx* ctx)
{
  bn_copy(r_bar, &ctx->one_mont);

  for (i64 i = bn_bit_length(d) - 1; i >= 0; i--) {
    bn_mont_mul(r_bar, r_bar, r_bar, ctx);

    if (bn_get_bit(d, i)) {
      bn_mont_mul(r_bar, r_bar, a_bar, ctx);
    }
  }
}

void bn_mod_exp_slow(bignum* r, const bignum* a, const bignum* b,
                     const bignum* m)
{
  bignum base, exp, res, tmp;
  bn_init_multi(&base, &exp, &res, &tmp, NULL);

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
  bn_free_multi(&base, &exp, &res, &tmp, NULL);
}

// calculates a^b mod m into r
void bn_mod_exp(bignum* r, const bignum* a, const bignum* b, const bignum* m)
{
  if (bn_is_even(m)) {
    bn_mod_exp_slow(r, a, b, m);
    return;
  }

  // fast path
  bn_mont_ctx ctx;
  bn_mont_ctx_init(&ctx, m);
  bn_mod_exp_mont(r, a, b, m, &ctx);
  bn_mont_ctx_free(&ctx);
}

// calculates a^b mod m into r in the montgomery domain
void bn_mod_exp_mont(bignum* r, const bignum* a, const bignum* b,
                     const bignum* m, bn_mont_ctx* ctx)
{
  assert(!bn_is_zero(m) && !bn_is_even(m));
  assert(bn_cmp(&ctx->n, m) == 0);

  bignum base, result;
  bn_init_multi(&base, &result, NULL);

  bn_mont_in(&base, a, ctx);

  bn_mont_exp(&result, &base, b, ctx);

  bn_mont_out(r, &result, ctx);

  bn_free(&base);
  bn_free(&result);
}
