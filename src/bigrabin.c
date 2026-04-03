#include <ctype.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>

#include "../include/bignum.h"

bool bn_rabin(bignum* n, bignum* a)
{
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

bool bn_rabin_mont(bignum* n, bignum* a)
{
  if (bn_is_even(n)) return false;

  // Fast handling for small numbers
  if (n->size == 1 && n->limbs[0] <= 3) {
    return n->limbs[0] == 2 || n->limbs[0] == 3;
  }

  // Step 1: Find d and s such that n-1 = d * 2^s
  bignum d, n_minus_1, one;
  bn_init(&d);
  bn_init(&n_minus_1);
  bn_init(&one);
  bn_set_u64(&one, 1);

  bn_sub(&n_minus_1, n, &one);  // n-1
  bn_copy(&d, &n_minus_1);

  u64 s = 0;
  while (bn_is_even(&d) && !bn_is_zero(&d)) {
    bn_rshift1(&d);
    s++;
  }

  // Step 2: Setup Montgomery context
  bn_mont_ctx ctx;
  bn_mont_ctx_init(&ctx, n);

  // Step 3: Precompute Montgomery representations
  bignum a_bar, x_bar, n_minus_1_mont, tmp;
  bn_init(&a_bar);
  bn_init(&x_bar);
  bn_init(&n_minus_1_mont);
  bn_init(&tmp);

  bn_alloc(&tmp, n->size);
  bn_alloc(&a_bar, n->size);
  bn_alloc(&x_bar, n->size);
  bn_alloc(&n_minus_1_mont, n->size);

  bn_mont_in(&a_bar, a, &ctx);                    // a in Montgomery
  bn_mont_in(&n_minus_1_mont, &n_minus_1, &ctx);  // (n-1) in Montgomery

  // Step 4: x = a^d mod n (Montgomery)
  bn_mont_exp(&x_bar, &a_bar, &d, &ctx);

  bool composite = true;

  // Step 5: Rabin-Miller checks
  if (bn_cmp(&x_bar, &ctx.one_mont) == 0 ||
      bn_cmp(&x_bar, &n_minus_1_mont) == 0) {
    composite = false;
    goto cleanup;
  }

  for (u64 r = 1; r < s; r++) {
    bn_mont_mul(&x_bar, &x_bar, &x_bar, &ctx);  // x = x^2 mod n
    // bn_copy(&x_bar, &tmp);

    if (bn_cmp(&x_bar, &ctx.one_mont) == 0) {
      composite = true;  // Non-trivial square root of 1
      break;
    }

    if (bn_cmp(&x_bar, &n_minus_1_mont) == 0) {
      composite = false;  // Probably prime
      break;
    }
  }

cleanup:
  // Free all temporaries
  bn_free(&d);
  bn_free(&n_minus_1);
  bn_free(&one);
  bn_free(&a_bar);
  bn_free(&x_bar);
  bn_free(&n_minus_1_mont);
  bn_free(&ctx.n);
  bn_free(&ctx.one_mont);
  bn_free(&ctx.r_square);

  return !composite;
}
