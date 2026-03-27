#include "../include/bignum.h"

#include <ctype.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>

void bn_mod(bignum* r, const bignum* a, const bignum* b) {
  bn_init(r);
  r->limbs = calloc(1, sizeof(u64));
  r->size = 1;

  int nbits = bn_bit_length(a);

  for (int i = nbits - 1; i >= 0; i--) {
    bn_lshift1_add(r, bn_get_bit(a, i));

    if (bn_cmp(r, b) >= 0) {
      bignum tmp;
      bn_init(&tmp);
      bn_sub_abs(&tmp, r, b);  // safe subtraction
      bn_free(r);
      *r = tmp;
    }
  }

  bn_trim(r);
}



uint64_t bn_divmod_u64(bignum* q, const bignum* a, uint64_t d) {
  q->limbs = calloc(a->size, sizeof(u64));
  q->size = a->size;

  __int128 rem = 0;

  for (i64 i = a->size - 1; i >= 0; i--) {
    __int128 cur = (rem << 64) | a->limbs[i];

    q->limbs[i] = (u64)(cur / d);
    rem = cur % d;
  }

  bn_trim(q);
  return (uint64_t)rem;
}
