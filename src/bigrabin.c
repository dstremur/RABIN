#include <ctype.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>

#include "../include/bignum.h"

bool bn_rabin(bignum* n, bignum* a) {
  if (bn_is_even(n)) return false;

  bignum one, two;
  //  bn_init(&one);
  // bn_init(&two);
  bn_init_multi(&one, &two);
  bn_set_u64(&one, 1);
  bn_set_u64(&two, 2);

  if (bn_cmp(n, &one) <= 0) {
    bn_free(&one);
    bn_free(&two);
    return false;
  }

  if (bn_cmp(n, &two) == 0) {
    bn_free(&one);
    bn_free(&two);
    return true;
  }

  bignum d, n_minus_1;
  bn_init(&d);
  bn_init(&n_minus_1);

  bn_sub(&n_minus_1, n, &one);
  bn_copy(&d, &n_minus_1);

  u64 s = 0;
  while (bn_is_even(&d) && !bn_is_zero(&d)) {
    bn_rshift1(&d);
    s++;
  }

  bignum x, tmp;
  bn_init(&x);
  bn_init(&tmp);

  bn_mod_exp(&x, a, &d, n);

  bool composite = true;

  if (bn_cmp(&x, &one) == 0 || bn_cmp(&x, &n_minus_1) == 0) {
    composite = false;
    goto cleanup;
  }

  for (u64 r = 1; r < s; r++) {
    bn_mul(&tmp, &x, &x);
    bn_mod(&x, &tmp, n);

    if (bn_cmp(&x, &one) == 0) {
      composite = true;
      break;
    }

    if (bn_cmp(&x, &n_minus_1) == 0) {
      composite = false;
      break;
    }
  }

cleanup:
  bn_free(&one);
  bn_free(&two);
  bn_free(&d);
  bn_free(&n_minus_1);
  bn_free(&x);
  bn_free(&tmp);

  return !composite;
}
