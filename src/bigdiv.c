#include <ctype.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>

#include "../include/bignum.h"

void bn_div(bignum* q, const bignum* a, const bignum* b)
{
  if (b->size == 0 || (b->size == 1 && b->limbs[0] == 0)) {
    return;
  }

  if (q == a || q == b) {
    bignum tmp;
    bn_init(&tmp);
    bn_div(&tmp, a, b);
    bn_copy(q, &tmp);
    bn_free(&tmp);
    return;
  }

  bignum r;
  bn_init(&r);

  bn_alloc(q, a->size);
  memset(q->limbs, 0, q->size * sizeof(u64));
  q->size = a->size;
  q->is_neg = a->is_neg ^ b->is_neg;

  i64 nbits = bn_bit_length(a);

  for (i64 i = nbits - 1; i >= 0; i--) {
    bn_lshift1_add(&r, bn_get_bit(a, i));

    if (bn_cmp(&r, b) >= 0) {
      bn_sub_abs(&r, &r, b);
      bn_set_bit(q, i);
    }
  }
  bn_trim(q);
  bn_free(&r);
}
