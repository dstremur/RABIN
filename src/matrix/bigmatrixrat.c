
#include "../../include/bigmatrixrat.h"

void bigmatrixrat_init(bigmatrixrat* A, u64 m, u64 n)
{
  bigmatrix_init(&A->M, m, n);
  bn_init(&A->den);
  bn_set_u64(&A->den, 1);
}

void bigmatrixrat_free(bigmatrixrat* A)
{
  bigmatrix_free(&A->M);
  bn_free(&A->den);
}
void bigmatrixrat_copy(bigmatrixrat* R, bigmatrixrat* A)
{
  bigmatrix_copy(&R->M, &A->M);
  bn_copy(&R->den, &A->den);
}

void bigmatrixrat_copyZ(bigmatrixrat* R, bigmatrix* A)
{
  bigmatrix_copy(&R->M, A);
  bn_set_u64(&R->den, 1);
}

void bigmatrixrat_add(bigmatrixrat* R, const bigmatrixrat* A,
                      const bigmatrixrat* B)
{
  bigmatrix TA, TB;
  bigmatrix_init(&TA, A->M.r_size, A->M.c_size);
  bigmatrix_init(&TB, A->M.r_size, A->M.c_size);

  bignum l, sA, sB;
  bn_init_multi(&l, &sA, &sB);

  bn_lcm(&l, &A->den, &B->den);

  // find scaling factor
  bn_div_exact(&sA, &l, &A->den);
  bn_div_exact(&sB, &l, &B->den);

  // scale the matrices
  bigmatrix_scalar(&TA, &A->M, &sA);
  bigmatrix_scalar(&TB, &B->M, &sB);

  bigmatrix_add(&R->M, &TA, &TB);
  bn_copy(&R->den, &l);

  bigmatrix_free(&TA);
  bigmatrix_free(&TB);
  bn_free_multi(&l, &sA, &sB);
}

void bigmatrixrat_get_element(bigrat* out, const bigmatrixrat* A, u64 i, u64 j)
{
  bigmatrix_get(&out->num, &A->M, i, j);
  bn_copy(&out->den, &A->den);
  br_normalize(out);
}

void bigmatrixrat_print(const bigmatrixrat* A)
{
  bignum temp;
  bn_init(&temp);

  for (u64 i = 0; i < A->M.r_size; i++) {
    for (u64 j = 0; j < A->M.c_size; j++) {
      bigmatrix_get(&temp, &A->M, i, j);
      bn_print(&temp);
      printf("/");
      bn_print(&A->den);
      printf(" ");
    }
    printf("\n");
  }

  bn_free(&temp);
}

void bigmatrixrat_sub(bigmatrixrat* R, const bigmatrixrat* A,
                      const bigmatrixrat* B)
{
  bigmatrix TA, TB;
  bigmatrix_init(&TA, A->M.r_size, A->M.c_size);
  bigmatrix_init(&TB, A->M.r_size, A->M.c_size);

  bignum l, sA, sB;
  bn_init_multi(&l, &sA, &sB);

  bn_lcm(&l, &A->den, &B->den);

  // find scaling factor
  bn_div_exact(&sA, &l, &A->den);
  bn_div_exact(&sB, &l, &B->den);

  // scale the matrices
  bigmatrix_scalar(&TA, &A->M, &sA);
  bigmatrix_scalar(&TB, &B->M, &sB);

  bigmatrix_sub(&R->M, &TA, &TB);
  bn_copy(&R->den, &l);

  bigmatrix_free(&TA);
  bigmatrix_free(&TB);
  bn_free_multi(&l, &sA, &sB);
}

void bigmatrixrat_normalize(bigmatrixrat* M)
{
  bignum d;
  bn_init(&d);

  // Compute GCD of m->den and all entries in m->num
  bigmatrix_gcd_all(&d, &M->M);
  bn_gcd(&d, &d, &M->den);

  // If GCD > 1, divide both numerator matrix entries and denominator by g
  if (!bn_is_eq_i64(&d, 1)) {
    bigmatrix_div_exact_scalar(&M->M, &M->M, &d);
    bn_div_exact(&M->den, &M->den, &d);
  }

  bn_free(&d);
}
