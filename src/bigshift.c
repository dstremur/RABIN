#include "../include/bignum.h"

#include <ctype.h>

#include <stddef.h>
#include <stdio.h>
#include <string.h>

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

// Shift r left by 1 and add 0 or 1
void bn_lshift1_add(bignum* r, int bit) {
  u64 carry = bit;
  for (u64 i = 0; i < r->size; i++) {
    __uint128_t tmp = ((__uint128_t)r->limbs[i] << 1) | carry;
    r->limbs[i] = (u64)tmp;
    carry = tmp >> 64;
  }

  if (carry) {
    r->limbs = realloc(r->limbs, (r->size + 1) * sizeof(u64));
    r->limbs[r->size++] = carry;
  }
}