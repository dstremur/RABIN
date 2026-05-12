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
void bigmatrix_add(bigmatrix* R, const bigmatrix* A, const bigmatrix* B);
void bigmatrix_get(bignum* R, const bigmatrix* A, u64 r, u64 c);
void bigmatrix_set(bigmatrix* A, const bignum* a, u64 r, u64 c);
void bigmatrix_mul(bigmatrix* R, const bigmatrix* A, const bigmatrix* B);
void bigmatrix_copy(bigmatrix* R, bigmatrix* A);
void bigmatrix_det(bignum* d, const bigmatrix* A); 
void bigmatrix_print(const bigmatrix* A); 
void bigmatrix_get_col(bigvector* c, const bigmatrix* A, u64 col);
void bigmatrix_get_row(bigvector* r, const bigmatrix* A, u64 row);
void bigmatrix_hadamard(bignum* r, const bigmatrix* A);
void bigmatrix_mv(bigvector* r, const bigmatrix* A, bigvector* a);
void bigmatrix_vm(bigvector* r, const bigmatrix* A, bigvector* a);
void bigmatrix_print_python(const bigmatrix* A);


#endif
