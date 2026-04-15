#include <stddef.h>
#include <string.h>

#include "../include/bignum.h"

#define KARATSUBA_LIMIT 128

void bn_mul(bignum* r, const bignum* a, const bignum* b)
{
  // aliasing
  if (r == a || r == b) {
    bignum tmp;
    bn_init(&tmp);
    bn_copy(&tmp, r);
    bn_mul_raw(&tmp, a, b);
    bn_copy(r, &tmp);
    bn_free(&tmp);
  } else {
    bn_mul_raw(r, a, b);
  }
}

// now use karatsuba
void bn_mul_raw(bignum* r, const bignum* a, const bignum* b)
{
  // 1. Basic checks
  if (a->size == 0 || b->size == 0) {
    r->size = 0;
    return;
  }

  if (a->size >= KARATSUBA_LIMIT && b->size >= KARATSUBA_LIMIT) {
    bn_mul_karatsuba(r, a, b);
  } else {
    bn_mul_school(r, a, b);
  }

  r->is_neg = a->is_neg ^ b->is_neg;
}

extern void bn_mul_add_inner(uint64_t* r, const uint64_t* b, uint64_t a_limb, uint64_t len);

void bn_mul_school(bignum* r, const bignum* a, const bignum* b)
{
  u64 max = a->size + b->size;
  bn_alloc(r, max);

  memset(r->limbs, 0, max * sizeof(u64));
  r->size = max;
  r->is_neg = a->is_neg ^ b->is_neg;

  for (u64 i = 0; i < a->size; i++) {
    if (a->limbs[i] == 0) continue;
/*
	u64 carry = 0;

    for (u64 j = 0; j < b->size; j++) {
      u64 idx = i + j;

      unsigned __int128 prod =
          (unsigned __int128)a->limbs[i] * b->limbs[j] + r->limbs[idx] + carry;

      r->limbs[idx] = (u64)prod;
      carry = (u64)(prod >> 64);
    }

    // ripple carry
    u64 k = i + b->size;
    unsigned __int128 ripple = (unsigned __int128)r->limbs[k] + carry;
    r->limbs[k] = (u64)ripple;
    u64 extra = (u64)(ripple >> 64);

    while (extra && ++k < r->size) {
      ripple = (unsigned __int128)r->limbs[k] + extra;
      r->limbs[k] = (u64)ripple;
      extra = (u64)(ripple >> 64);
    }
*/
	bn_mul_add_inner(&r->limbs[i], b->limbs, a->limbs[i], b->size);
  }

  bn_trim(r);
}

void bn_mul_karatsuba(bignum* r, const bignum* a, const bignum* b)
{
  if (r == a || r == b) {
    bignum tmp;
    bn_init(&tmp);
    bn_mul_karatsuba(&tmp, a, b);
    bn_copy(r, &tmp);
    bn_free(&tmp);
    return;
  }

  if (a->size < KARATSUBA_LIMIT || b->size < KARATSUBA_LIMIT) {
    bn_mul(r, a, b);
    return;
  }

  u64 max = MAX(a->size, b->size);

  u64 m = max / 2;

  bignum a0, a1, b0, b1;

  bn_init(&a0);
  bn_init(&a1);
  bn_init(&b0);
  bn_init(&b1);

  if (a->size > m) {
    bn_alloc(&a0, m);
    memcpy(a0.limbs, a->limbs, m * sizeof(u64));
    a0.size = m;
    bn_trim(&a0);

    bn_alloc(&a1, a->size - m);
    memcpy(a1.limbs, a->limbs + m, (a->size - m) * sizeof(u64));
    a1.size = a->size - m;
    bn_trim(&a1);
  } else {
    bn_copy(&a0, a);
  }

  a0.is_neg = 0;
  a1.is_neg = 0;

  if (b->size > m) {
    bn_alloc(&b0, m);
    memcpy(b0.limbs, b->limbs, m * sizeof(u64));
    b0.size = m;
    bn_trim(&b0);

    bn_alloc(&b1, b->size - m);
    memcpy(b1.limbs, b->limbs + m, (b->size - m) * sizeof(u64));
    b1.size = b->size - m;
    bn_trim(&b1);
  } else {
    bn_copy(&b0, b);
  }

  b0.is_neg = 0;
  b1.is_neg = 0;

  bignum z0, z1, z2, s_a, s_b, tmp;

  bn_init(&z0);
  bn_init(&z1);
  bn_init(&z2);
  bn_init(&s_a);
  bn_init(&s_b);
  bn_init(&tmp);

  bn_mul(&z0, &a0, &b0);

  bn_mul(&z2, &a1, &b1);

  bn_add(&s_a, &a0, &a1);
  bn_add(&s_b, &b0, &b1);
  bn_mul(&z1, &s_a, &s_b);

  bn_sub(&tmp, &z1, &z0);
  bn_sub(&z1, &tmp, &z2);

  u64 max_len = a->size + b->size;
  bn_alloc(r, max_len);
  memset(r->limbs, 0, max_len * sizeof(u64));
  r->size = max_len;
  r->is_neg = a->is_neg ^ b->is_neg;

  bn_add_at_offset(r, &z0, 0);
  bn_add_at_offset(r, &z1, m);
  bn_add_at_offset(r, &z2, 2 * m);

  bn_free(&a0);
  bn_free(&a1);
  bn_free(&b0);
  bn_free(&b1);
  bn_free(&z0);
  bn_free(&z1);
  bn_free(&z2);
  bn_free(&s_a);
  bn_free(&s_b);
  bn_free(&tmp);

  bn_trim(r);
}
