#ifndef BIGMATRIX_H
#define BIGMATRIX_H

#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>

#include "bignum.h"
#include "bigvector.h"
// use a flat array structure
// index = (row * nr_cols) + col
typedef struct bigmatrix {
  bignum* data;
  u64 r_size;
  u64 c_size;
} bigmatrix;

#define GET(M, r, c) (&(M)->data[(r) * (M)->c_size + (c)])

void bigmatrix_init(bigmatrix* A, u64 r, u64 c);
void bigmatrix_free(bigmatrix* A);
void bigmatrix_add(bigmatrix* R, bigmatrix* A, bigmatrix* B);
void bigmatrix_get(bignum* R, bigmatrix* A, u64 r, u64 c);
void bigmatrix_set(bigmatrix* A, bignum* a, u64 r, u64 c);
void bigmatrix_mul(bigmatrix* R, bigmatrix* A, bigmatrix* B);
void bigmatrix_copy(bigmatrix* R, bigmatrix* A);
void bigmatrix_det(bignum* d, bigmatrix* A); 
void bigmatrix_print(bigmatrix* A); 
void bigmatrix_get_col(bigvector* c, bigmatrix* A, u64 col);
void bigmatrix_get_row(bigvector* r, bigmatrix* A, u64 row);
void bigmatrix_hadamard(bignum* r, bigmatrix* A);
void bigmatrix_mv(bigvector* r, bigmatrix* A, bigvector* a);
void bigmatrix_vm(bigvector* r, bigmatrix* A, bigvector* a);
#endif
