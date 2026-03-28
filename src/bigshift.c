#include <ctype.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>

#include "../include/bignum.h"

void bn_lshift1(bignum* r) {
  u64 carry = 0;

  for (u64 i = 0; i < r->size; i++) {
    u64 new_carry = r->limbs[i] >> 63;
    r->limbs[i] = (r->limbs[i] << 1) | carry;
    carry = new_carry;
  }

  if (carry) {
    r->limbs = realloc(r->limbs, (r->size + 1) * sizeof(u64));
    r->limbs[r->size++] = carry;
  }
}

void bn_rshift(bignum* r, const bignum* a, int shift) {
  if (shift == 0) {
    bn_copy(r, a);
    return;
  }

  u64 words = shift / 64;
  u64 bits = shift % 64;

  if (words >= a->size) {
    bn_set_u64(r, 0);
    return;
  }

  u64 new_size = a->size - words;

  bn_alloc(r, new_size);

  u64* new_limbs = r->limbs;

  if (bits == 0) {
    for (u64 i = 0; i < new_size; i++) {
      new_limbs[i] = a->limbs[i + words];
    }
  } else {
    for (u64 i = 0; i < new_size; i++) {
      u64 curr = a->limbs[i + words];
      u64 next = (i + words + 1 < a->size) ? a->limbs[i + words + 1] : 0;

      r->limbs[i] = (curr >> bits) | (next << (64 - bits));
    }
  }

  // r->limbs = new_limbs;
  r->size = new_size;

  if (r->capacity > r->size) {
    memset(r->limbs + r->size, 0, (r->capacity - r->size) * sizeof(u64));
  }
  bn_trim(r);
}

// Shift r left by 1 and add 0 or 1
void bn_lshift1_add(bignum* r, int bit) {
  u64 carry = (bit != 0);
  for (u64 i = 0; i < r->size; i++) {
    u64 tmp = r->limbs[i] >> 63;
    r->limbs[i] = (r->limbs[i] << 1) | carry;
    carry = tmp;
  }

  if (carry) {
    bn_alloc(r, r->size + 1);
    r->limbs[r->size++] = carry;
  }
}
