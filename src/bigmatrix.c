#include "../include/bignum.h"

void bigmatrix_init(bigmatrix* M) {
  M->cols = NULL;
  M->rows = NULL;
  M->c_size = 0;
  M->r_size = 0;
}
