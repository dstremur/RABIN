#include "../include/bignum.h"

#include <ctype.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>

void bn_sub_abs(bignum* r, const bignum* a, const bignum* b) {
  r->limbs = calloc(a->size, sizeof(u64));

  u64 borrow = 0;

  for (u64 i = 0; i < a->size; i++) {
    u64 av = a->limbs[i];
    u64 bv = (i < b->size) ? b->limbs[i] : 0;

    unsigned __int128 diff = (unsigned __int128)av - bv - borrow;

    r->limbs[i] = (u64)diff;
    borrow = (diff >> 127) & 1;  // detect underflow
  }

  r->size = a->size;
  bn_trim(r);
}
// Subtract: r = a - b
void bn_sub(bignum* r, const bignum* a, const bignum* b) {
  // a - (-b) = a + b
  if (!a->is_neg && b->is_neg) {
    bn_add(r, a, b);
    r->is_neg = false;
    return;
  }

  // (-a) - b = -(a + b)
  if (a->is_neg && !b->is_neg) {
    bn_add(r, a, b);
    r->is_neg = true;
    return;
  }

  i64 cmp = bn_cmp(a, b);

  if (cmp == 0) {
    bn_init(r);
    return;
  }

  if (cmp > 0) {
    // |a| > |b|
    bn_sub_abs(r, a, b);
    r->is_neg = a->is_neg;
  } else {
    // |b| > |a|
    bn_sub_abs(r, b, a);
    r->is_neg = !a->is_neg;
  }
}