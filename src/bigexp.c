#include "../include/bignum.h"

#include <ctype.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>

void bn_rshift1(bignum* r) {
    if (r->size == 0) return;

    for (i64 i = 0; i < r->size; i++) {
        r->limbs[i] >>= 1;
    	if (i + 1 < r->size) {
			if (r->limbs[i  + 1] & 1ULL) {
				r->limbs[i] |= (1ULL << 63);
			}
		}
	}

    bn_trim(r);
}
// calculates a^b into r
// binary exponentiation
void bn_pow(bignum* r, bignum* a, bignum* b){

	if (bn_is_zero(b)) {
		bn_set_u64(r, 1);
		return;
	}

	bignum base, exp, two;
	bn_init(&base);
	bn_init(&exp);
	bn_init(&two);
	bn_set_u64(&two, 2);
	bn_copy(&base, a);
	bn_copy(&exp, b); 

	bn_set_u64(r, 1);


    while(!bn_is_zero(&exp)){
        if (!bn_is_even(&exp)){
            bn_mul(r, r, &base);
        }       
        bn_mul(&base, &base, &base);

        bn_rshift1(&exp); 
    }

	bn_free(&base);
	bn_free(&exp);
	bn_free(&two);
}
