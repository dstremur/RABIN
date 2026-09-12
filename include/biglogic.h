#ifndef BIGLOGIC_H
#define BIGLOGIC_H

#include "bigcore.h"

void bn_and(bignum* r, const bignum* a, const bignum* m);

void bn_or(bignum* r, const bignum* a, const bignum* m);

void bn_xor(bignum* r, const bignum* a, const bignum* m);

void bn_not(bignum* r, const bignum* a);

#endif