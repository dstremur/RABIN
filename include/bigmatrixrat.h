#ifndef BIGMATRIXRAT_H
#define BIGMATRIXRAT_H

#include "bignum.h"
#include "bigrat.h"

typedef struct {
  bigmatrix M;
  bignum den;
} bigmatrixrat;

#endif
