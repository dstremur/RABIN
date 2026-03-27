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

void bn_rshift(bignum *r, const bignum *a, int shift) {
	if (shift == 0) {
		bn_copy(r, a);
		return;
	}

	u64 words = shift / 64; 
	u64 bits = shift % 64; 

	if (words >= a->size) {
		bn_set_u64(a, 0);
		return;
	}

	u64 new_size = a->size - words;
	u64* new_limbs = calloc(new_size, sizeof(u64)); 

	if (bits == 0) {
		for (u64 i = 0; i < new_size; i++) {
			new_limbs[i] = a->limbs[i + words];
		}
	} else {
		for (u64 i = 0; i < new_size; i++) {
			u64 curr = a->limbs[i + words];
			u64 next = (i + words + 1 < a->size) ? a->limbs[i + words + 1] : 0;

			new_limbs[i] = (curr >> bits) | (next << (64 - bits)); 
		}
	}

	free(r->limbs);

	r->limbs = new_limbs;
	r->size = new_size;
	bn_trim(r); 
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
