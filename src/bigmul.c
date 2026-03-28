#include <ctype.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>

#include "../include/bignum.h"

void bn_mul(bignum* r, const bignum* a, const bignum* b) {
  // Handle aliasing
  if (r == a || r == b) {
    bignum tmp;
    bn_init(&tmp);
    bn_mul(&tmp, a, b);
    bn_copy(r, &tmp);
    bn_free(&tmp);
    return;
  }

  u64 max = a->size + b->size;
  bn_alloc(r, max);

  memset(r->limbs, 0, max * sizeof(u64));
  r->size = max;
  r->is_neg = a->is_neg ^ b->is_neg;

  for (u64 i = 0; i < a->size; i++) {
    if (a->limbs[i] == 0) continue;
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
    __int128 ripple = (__int128)r->limbs[k] + carry;
    r->limbs[k] = (u64)ripple;
    u64 extra = (u64)(ripple >> 64);

    while (extra && ++k < r->size) {
      ripple = (__int128)r->limbs[k] + extra;
      r->limbs[k] = (u64)ripple;
      extra = (u64)(ripple >> 64);
    }
  }

  bn_trim(r);
}
