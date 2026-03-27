#include "../include/bignum.h"

#include <ctype.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>


// calculates a^b into r
// binary exponentiation
void bn_pow(bignum* r, bignum* a, bignum* b){
    bn_init(r);

    while(!bn_is_zero(b)){
        if (!bn_is_even(b)){
            bn_mul(r, r, a);
        }       
        bn_mul(a, a, a);
        bn_rshift(a,a,1);
    }

}

/* 
void bn_pow(bignum* r, const bignum* a, const bignum* b) {
    // 1. Handle the base case: a^0 = 1
    if (bn_is_zero(b)) {
        bn_set_u64(r, 1); 
        return;
    }

    // 2. Create temps so we don't destroy the user's input variables
    bignum base, exp;
    bn_copy(&base, a);
    bn_copy(&exp, b);

    // 3. Initialize result to 1
    bn_set_u64(r, 1);

    while (!bn_is_zero(&exp)) {
        // If exponent is odd, multiply result by current base
        if (!bn_is_even(&exp)) {
            bn_mul(r, r, &base);
        }
        
        // Square the base
        bn_mul(&base, &base, &base);
        
        // Shift the EXPONENT, not the base
        bn_rshift(&exp, &exp, 1);
    }

    // Clean up temps
    bn_free(&base);
    bn_free(&exp);
}
*/