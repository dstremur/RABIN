#include "../../include/rmatq.h"

rabin_err_t rmatq_init(rmatq_t* A, u64 m, u64 n)
{
  if (A == NULL) return RABIN_ERR_NULL_PTR;

  rabin_err_t err = rmat_init(&A->M, m, n);
  if (err != RABIN_SUCCESS) return err;

  rz_init(&A->den);
  return rz_set_u64(&A->den, 1);
}

rabin_err_t rmatq_clear(rmatq_t* A)
{
  if (A == NULL) return RABIN_ERR_NULL_PTR;

  rmat_clear(&A->M);
  rz_clear(&A->den);
  return RABIN_SUCCESS;
}

rabin_err_t rmatq_copy(rmatq_t* R, rmatq_t* A)
{
  if (R == NULL || A == NULL) return RABIN_ERR_NULL_PTR;

  rabin_err_t err = rmat_copy(&R->M, &A->M);
  if (err != RABIN_SUCCESS) return err;

  return rz_copy(&R->den, &A->den);
}

rabin_err_t rmatq_copyZ(rmatq_t* R, rmat_t* A)
{
  if (R == NULL || A == NULL) return RABIN_ERR_NULL_PTR;

  rabin_err_t err = rmat_copy(&R->M, A);
  if (err != RABIN_SUCCESS) return err;

  return rz_set_u64(&R->den, 1);
}

rabin_err_t rmatq_add(rmatq_t* R, const rmatq_t* A, const rmatq_t* B)
{
  if (R == NULL || A == NULL || B == NULL) return RABIN_ERR_NULL_PTR;
  if (A->M.rows != B->M.rows || A->M.cols != B->M.cols)
    return RABIN_ERR_MATRIX_DIM;
  if (R->M.rows != A->M.rows || R->M.cols != A->M.cols)
    return RABIN_ERR_MATRIX_DIM;

  rmat_t TA = {0}, TB = {0};
  rz_t l, sA, sB;
  rz_init_multi(&l, &sA, &sB, NULL);
  rabin_err_t err = RABIN_SUCCESS;

  TRY(rmat_init(&TA, A->M.rows, A->M.cols));
  TRY(rmat_init(&TB, A->M.rows, A->M.cols));

  TRY(rz_lcm(&l, &A->den, &B->den));

  // find scaling factor
  TRY(rz_div_exact(&sA, &l, &A->den));
  TRY(rz_div_exact(&sB, &l, &B->den));

  // scale the matrices
  TRY(rmat_scalar(&TA, &A->M, &sA));
  TRY(rmat_scalar(&TB, &B->M, &sB));

  TRY(rmat_add(&R->M, &TA, &TB));
  TRY(rz_copy(&R->den, &l));

cleanup:
  rmat_clear(&TB);
  rmat_clear(&TA);
  rz_clear_multi(&l, &sA, &sB, NULL);
  return err;
}

#if 0
/* WIP: unfinished element accessors, not declared in rmatq.h and unused.
 * rmatq_set() also has the copy direction inverted and neither returns a
 * status. Kept for reference, excluded from the build. */
rabin_err_t rmatq_get(rq_t* a, rmatq_t* A, u64 r, u64 c)
{
  if (!A || !a) return RABIN_ERR_INVALID_ARG;

  rz_copy(&a->den, &A->den);

  rz_copy(&a->num, RMAT_GET(&A->M, r, c));

  rq_normalize(a);
}

rabin_err_t rmatq_set(rmatq_t* A, rq_t* a, u64 r, u64 c)
{
  if (!A || !a) return RABIN_ERR_INVALID_ARG;

  rz_copy(&a->num, RMAT_GET(&A->M, r, c));

  rq_normalize(a);
}
#endif

rabin_err_t rmatq_print(const rmatq_t* A)
{
  if (A == NULL) return RABIN_ERR_NULL_PTR;

  rz_t temp;
  rz_init(&temp);

  for (u64 i = 0; i < A->M.rows; i++) {
    for (u64 j = 0; j < A->M.cols; j++) {
      rmat_get(&temp, &A->M, i, j);
      rz_print(&temp);
      printf("/");
      rz_print(&A->den);
      printf(" ");
    }
    printf("\n");
  }

  rz_clear(&temp);
  return RABIN_SUCCESS;
}

rabin_err_t rmatq_sub(rmatq_t* R, const rmatq_t* A, const rmatq_t* B)
{
  if (R == NULL || A == NULL || B == NULL) return RABIN_ERR_NULL_PTR;
  if (A->M.rows != B->M.rows || A->M.cols != B->M.cols)
    return RABIN_ERR_MATRIX_DIM;
  if (R->M.rows != A->M.rows || R->M.cols != A->M.cols)
    return RABIN_ERR_MATRIX_DIM;

  rmat_t TA = {0}, TB = {0};
  rz_t l, sA, sB;
  rz_init_multi(&l, &sA, &sB, NULL);
  rabin_err_t err = RABIN_SUCCESS;

  TRY(rmat_init(&TA, A->M.rows, A->M.cols));
  TRY(rmat_init(&TB, A->M.rows, A->M.cols));

  TRY(rz_lcm(&l, &A->den, &B->den));

  // find scaling factor
  TRY(rz_div_exact(&sA, &l, &A->den));
  TRY(rz_div_exact(&sB, &l, &B->den));

  // scale the matrices
  TRY(rmat_scalar(&TA, &A->M, &sA));
  TRY(rmat_scalar(&TB, &B->M, &sB));

  TRY(rmat_sub(&R->M, &TA, &TB));
  TRY(rz_copy(&R->den, &l));

cleanup:
  rmat_clear(&TB);
  rmat_clear(&TA);
  rz_clear_multi(&l, &sA, &sB, NULL);
  return err;
}

rabin_err_t rmatq_normalize(rmatq_t* M)
{
  if (M == NULL) return RABIN_ERR_NULL_PTR;

  rz_t d;
  rz_init(&d);

  // Compute GCD of m->den and all entries in m->num
  rabin_err_t err = RABIN_SUCCESS;
  TRY(rmat_gcd_all(&d, &M->M));
  TRY(rz_gcd(&d, &d, &M->den));

  // If GCD > 1, divide both numerator matrix entries and denominator by g
  if (!rz_is_eq_i64(&d, 1)) {
    TRY(rmat_div_exact_scalar(&M->M, &M->M, &d));
    TRY(rz_div_exact(&M->den, &M->den, &d));
  }

  err = RABIN_SUCCESS;
cleanup:
  rz_clear(&d);
  return err;
}
