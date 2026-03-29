#include <ctype.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>

#include "../include/bignum.h"

void bn_add(bignum* r, const bignum* a, const bignum* b) {
  // aliasing
  //
  if (r == a || r == b) {
    bignum tmp;
    bn_init(&tmp);
    bn_add(&tmp, a, b);
    bn_copy(r, &tmp);
    bn_free(&tmp);
    return;
  }

  u64 max = MAX(a->size, b->size);

  bn_alloc(r, max + 1);

  u64 carry = 0;

  for (u64 i = 0; i < max; i++) {
    u64 av = (i < a->size) ? a->limbs[i] : 0;
    u64 bv = (i < b->size) ? b->limbs[i] : 0;

    unsigned __int128 sum = (unsigned __int128)av + bv + carry;

    r->limbs[i] = (u64)sum;    // low 64 bits
    carry = (u64)(sum >> 64);  // high 64 bits
  }

  r->limbs[max] = carry;
  r->size = max + (carry ? 1 : 0);
  bn_trim(r);
}

void bn_add_u64(bignum* r, const bignum* a, u64 b) {
  if (b == 0) {
    bn_copy(r, a);
    return;
  }

  bn_copy(r, a);

  u64 carry = b;

  for (u64 i = 0; i < r->size && carry > 0; i++) {
    u64 old = r->limbs[i];
    r->limbs[i] += carry;

    if (r->limbs[i] < old) {
      carry = 1;
    } else {
      carry = 0;
    }
  }

  if (carry) {
    if (bn_alloc(r, r->size + 1)) {
      r->limbs[r->size] = carry;
      r->size++;
    }
  }
}
