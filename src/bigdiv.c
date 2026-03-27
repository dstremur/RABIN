#include "../include/bignum.h"

#include <ctype.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>





void bn_div(bignum* q, const bignum* a, const bignum* b) {
	if (q == a || q == b) {
    bignum tmp;
    bn_init(&tmp);
    bn_div(&tmp, a, b);
    bn_copy(q, &tmp);
    bn_free(&tmp);
    return;
}

	bignum r;
  bn_init(q);
  bn_init(&r);

  q->limbs = calloc(a->size, sizeof(u64));
  q->size = a->size; 

  i64 nbits = bn_bit_length(a);

  for (i64 i = nbits - 1; i >= 0; i--) {
	  bn_lshift1_add(&r, bn_get_bit(a, i));

	  if (bn_cmp(&r, b) >= 0) {
		  // TODO overwrite directly no temp
		  bignum tmp;
		  bn_init(&tmp);
		  bn_sub_abs(&tmp, &r, b);
		  bn_copy(&r, &tmp);
		  bn_free(&tmp);

		  bn_set_bit(q, i);
	  }
  }
  bn_trim(q);
  bn_free(&r);

}
