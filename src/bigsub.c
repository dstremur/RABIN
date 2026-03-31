#include <ctype.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>

#include "../include/bignum.h"

void bn_sub_abs(bignum* r, const bignum* a, const bignum* b)
{
  // aliasing

  if (r == a || r == b) {
    bignum tmp;
    bn_init(&tmp);
    bn_sub_abs(&tmp, a, b);
    bn_copy(r, &tmp);
    bn_free(&tmp);
    return;
  }
  bn_alloc(r, a->size);

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
void bn_sub(bignum* r, const bignum* a, const bignum* b)
{
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
    bn_set_u64(r, 0);
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
