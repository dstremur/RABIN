#ifndef BIGRAT_H
#define BIGRAT_H

#include "bignum.h"

/*
Invariants:
den > 0
gcd(num, den) = 1
*/

typedef struct bigrat {
  bignum num;
  bignum den;
} bigrat;

void br_add(bigrat* r, const bigrat* a, const bigrat* b);

void br_sub(bigrat* r, bigrat* a, bigrat* b);

void br_mul(bigrat* r, const bigrat* a, const bigrat* b);

void br_div(bigrat* r, bigrat* a, bigrat* b);

bool br_normalize(bigrat* r);

void br_print(const bigrat* r);

void br_neg(bigrat* r);

void br_init(bigrat* r);

void br_free(bigrat* r);

#endif