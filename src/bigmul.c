#include "../include/bignum.h"

#include <ctype.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>


void bn_mul(bignum* r, const bignum* a, const bignum* b) {
  u64 max = a->size + b->size;
  u64* temp = calloc(max, sizeof(u64));

  for (u64 i = 0; i < a->size; i++) {
    u64 carry = 0;

    for (u64 j = 0; j < b->size; j++) {
      u64 idx = i + j;

      unsigned __int128 prod =
          (unsigned __int128)a->limbs[i] * b->limbs[j] + temp[idx] + carry;

      temp[idx] = (u64)prod;
      carry = (u64)(prod >> 64);
    }

    temp[i + b->size] += carry;
  }

  free(r->limbs);
  r->limbs = temp;
  r->size = max;
  bn_trim(r);
}