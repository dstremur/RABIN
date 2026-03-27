#include "../include/bignum.h"

#include <ctype.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>





void bn_div(bignum* q, const bignum* a, const bignum* b) {
  bignum* r = malloc(sizeof(bignum));
  bn_init(q);
  bn_init(r);
  r->limbs = calloc(1, sizeof(u64)); 
  r->size = 1;

  q->limbs = calloc(a->size, sizeof(u64));
  q->size = a->size; 

  i64 nbits = bn_bit_length(a);

  for (i64 i = nbits - 1; i >= 0; i--) {
	  bn_lshift1_add(r, bn_get_bit(a, i));

	  if (bn_cmp(r, b) >= 0) {
		  // TODO overwrite directly no temp
		  bignum tmp;
		  bn_init(&tmp);
		  bn_sub_abs(&tmp, r, b);
		  bn_free(r);
		  *r = tmp;

		  bn_set_bit(q, i);
	  }
  }
  bn_trim(q);
  bn_free(r);

}
