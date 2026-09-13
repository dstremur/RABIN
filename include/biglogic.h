#ifndef BIGLOGIC_H
#define BIGLOGIC_H

#include <stdbool.h>

#include "bigcore.h"

bool bn_supports_avx512(void);

void bn_and(bignum* r, const bignum* a, const bignum* m);

void bn_or(bignum* r, const bignum* a, const bignum* m);

void bn_xor(bignum* r, const bignum* a, const bignum* m);

void bn_not(bignum* r, const bignum* a);

#endif