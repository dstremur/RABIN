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
