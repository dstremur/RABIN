#include <stddef.h>
#include <string.h>

#include "../include/bignum.h"

void bn_newton_div(bignum* q, const bignum* a, const bignum* d)
{
  if (bn_is_zero(d)) return;
  if (bn_cmp(a, d) < 0) {
    bn_set_u64(q, 0);
    return;
  }

  // estimate the reciprocal of d using newton iterations
  // find roots of f(x) = 1/x - d => y = x + x * (1 - x * d)

  bignum x, two_p, tmp, error, correction;
  bn_init_multi(&x, &two_p, &tmp, &error, &correction, NULL);

  // choose initial estimate
  u64 a_bits = bn_bit_length(a);
  u64 d_bits = bn_bit_length(d);
  // tune if precision is lacking
  u64 P = a_bits + 32;

  bn_set_u64(&x, 1);
  bn_lshift(&x, &x, P - d_bits);

  for (u64 i = 0; i < 32; i++) {
    // x = x + (x * (2^(P) - d * x)) >> P

    // tmp = b * x
    bn_mul(&tmp, &x, d);

    // two_p = 2^{P}
    bn_set_u64(&two_p, 1);
    bn_lshift(&two_p, &two_p, P);

    // Calculate (2^P - d*x) and update x
    if (bn_cmp(&two_p, &tmp) >= 0) {
      bn_sub(&error, &two_p, &tmp);

      bn_mul(&correction, &x, &error);
      bn_rshift(&correction, &correction, P);

      // x = x + correction
      bn_add(&x, &x, &correction);
    } else {
      // x was slightly too large, subtract correction instead
      bn_sub(&error, &tmp, &two_p);
      bn_mul(&correction, &x, &error);
      bn_rshift(&correction, &correction, P);
      bn_sub(&x, &x, &correction);
    }

    if (bn_is_zero(&correction)) break;
  }

  bn_mul(q, a, &x);
  bn_rshift(q, q, P);

  bignum r;
  bn_init(&r);
  bn_mul(&tmp, q, d);
  bn_sub(&r, a, &tmp);

  while (bn_cmp(&r, d) >= 0) {
    bn_add_u64(q, q, 1);
    bn_sub(&r, &r, d);
  }

  bn_free_multi(&x, &tmp, &two_p, &r, &error, &correction, NULL);
}

void bn_div(bignum* q, const bignum* a, const bignum* b)
{
  if (b->size == 0 || (b->size == 1 && b->limbs[0] == 0)) {
    return;
  }
  // aliasing
  if (q == a || q == b) {
    bignum tmp;
    bn_init(&tmp);
    bn_div(&tmp, a, b);
    bn_copy(q, &tmp);
    bn_free(&tmp);
    return;
  }

  bignum r;
  bn_init(&r);

  bn_alloc(q, a->size);
  memset(q->limbs, 0, q->size * sizeof(u64));
  q->size = a->size;
  q->is_neg = a->is_neg ^ b->is_neg;

  i64 nbits = bn_bit_length(a);

  // binary long division
  for (i64 i = nbits - 1; i >= 0; i--) {
    bn_lshift1(&r);
    if (bn_get_bit(a, i)) {
      bn_set_bit(&r, 0);
    }

    if (bn_cmp(&r, b) >= 0) {
      bn_sub(&r, &r, b);
      bn_set_bit(q, i);
    }
  }
  bn_trim(q);
  bn_free(&r);
}
