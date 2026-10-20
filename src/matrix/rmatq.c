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

rabin_err_t rmatq_get(rq_t* a, const rmatq_t* A, u64 r, u64 c)
{
  if (a == NULL || A == NULL) return RABIN_ERR_NULL_PTR;
  if (r >= A->M.rows || c >= A->M.cols) return RABIN_ERR_INVALID_ARG;

  rabin_err_t err = rz_copy(&a->num, RMAT_GET(&A->M, r, c));
  if (err != RABIN_SUCCESS) return err;

  err = rz_copy(&a->den, &A->den);
  if (err != RABIN_SUCCESS) return err;

  return rq_normalize(a);
}

rabin_err_t rmatq_set(rmatq_t* A, const rq_t* a, u64 r, u64 c)
{
  if (A == NULL || a == NULL) return RABIN_ERR_NULL_PTR;
  if (r >= A->M.rows || c >= A->M.cols) return RABIN_ERR_INVALID_ARG;
  if (rz_is_zero(&a->den)) return RABIN_ERR_INVALID_ARG;

  rz_t l, s, t;
  rz_init_multi(&l, &s, &t, NULL);
  rabin_err_t err = RABIN_SUCCESS;

  // Grow the common denominator so the new entry fits without rounding.
  TRY(rz_lcm(&l, &A->den, &a->den));

  if (rz_cmp(&l, &A->den) != 0) {
    TRY(rz_div_exact(&s, &l, &A->den));

    for (u64 k = 0; k < A->M.rows * A->M.cols; k++) {
      TRY(rz_mul(&t, &A->M.data[k], &s));
      TRY(rz_copy(&A->M.data[k], &t));
    }

    TRY(rz_copy(&A->den, &l));
  }

  // Store a = a.num / a.den as (a.num * den / a.den) / den.
  TRY(rz_div_exact(&s, &A->den, &a->den));
  TRY(rz_mul(&t, &a->num, &s));
  TRY(rz_copy(RMAT_GET(&A->M, r, c), &t));

  TRY(rmatq_normalize(A));

cleanup:
  rz_clear_multi(&l, &s, &t, NULL);
  return err;
}

rabin_err_t rmatq_set_zero(rmatq_t* A, u64 r, u64 c)
{
  if (A == NULL) return RABIN_ERR_NULL_PTR;
  if (r >= A->M.rows || c >= A->M.cols) return RABIN_ERR_INVALID_ARG;

  return rz_set_u64(RMAT_GET(&A->M, r, c), 0);
}

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

/**
 * Exact rational square root: sets \f$r = \sqrt a\f$ when both \f$a.num\f$
 * and \f$a.den\f$ are perfect squares, and fails otherwise. Internal to
 * rmatq_cholesky(), which never calls it with a negative or zero value.
 */
static rabin_err_t rq_sqrt_exact(rq_t* r, const rq_t* a)
{
  if (r == NULL || a == NULL) return RABIN_ERR_NULL_PTR;
  if (rz_is_zero(&a->den) || a->num.is_neg) return RABIN_ERR_INVALID_ARG;

  rz_t sn, sd;
  rz_init_multi(&sn, &sd, NULL);

  rabin_err_t err = RABIN_SUCCESS;

  if (!rz_is_square(&sn, &a->num) || !rz_is_square(&sd, &a->den)) {
    err = RABIN_ERR_INVALID_ARG;
    goto cleanup;
  }

  TRY(rz_copy(&r->num, &sn));
  TRY(rz_copy(&r->den, &sd));

cleanup:
  rz_clear_multi(&sn, &sd, NULL);
  return err;
}

rabin_err_t rmatq_cholesky(rmatq_t* C, const rmatq_t* A)
{
  if (C == NULL || A == NULL) return RABIN_ERR_NULL_PTR;
  if (A->M.rows != A->M.cols) return RABIN_ERR_MATRIX_DIM;
  if (C->M.rows != A->M.rows || C->M.cols != A->M.cols)
    return RABIN_ERR_MATRIX_DIM;

  const u64 n = A->M.rows;
  if (n > SIZE_MAX / n) return RABIN_ERR_OVERFLOW;
  if (n == 0) return rz_set_u64(&C->den, 1);

  // Scratch upper triangle Q of rationals, row-major n x n.
  rq_t* Q = malloc(n * n * sizeof(rq_t));
  if (Q == NULL) return RABIN_ERR_OUT_OF_MEMORY;

  rq_t tmp;
  rz_t D, s, t;
  rq_init(&tmp);
  rz_init_multi(&D, &s, &t, NULL);

  rabin_err_t err = RABIN_SUCCESS;
  u64 inited = 0;

  for (u64 k = 0; k < n * n; k++) {
    TRY(rq_init(&Q[k]));
    inited = k + 1;
  }

  // 1. [Initialize]
  for (u64 i = 0; i < n; i++) {
    for (u64 j = i; j < n; j++) {
      TRY(rmatq_get(&Q[i * n + j], A, i, j));
    }
  }

  // 2. [Loop on i]
  for (u64 i = 0; i < n; i++) {
    rq_t* pivot = &Q[i * n + i];

    if (rz_is_zero(&pivot->num)) {
      err = RABIN_ERR_DIV_BY_ZERO;
      goto cleanup;
    }
    if (pivot->num.is_neg) {
      err = RABIN_ERR_INVALID_ARG;
      goto cleanup;
    }

    for (u64 j = i + 1; j < n; j++) {
      // q_{j,i} <- q_{i,j}, then q_{i,j} <- q_{i,j} / q_{i,i}
      TRY(rq_copy(&Q[j * n + i], &Q[i * n + j]));
      TRY(rq_div(&Q[i * n + j], &Q[i * n + j], pivot));
    }

    // 3. [Main loop]
    for (u64 k = i + 1; k < n; k++) {
      if (rz_is_zero(&Q[k * n + i].num)) continue;

      for (u64 l = k; l < n; l++) {
        if (rz_is_zero(&Q[i * n + l].num)) continue;

        // q_{k,l} <- q_{k,l} - q_{k,i} * q_{i,l}
        TRY(rq_mul(&tmp, &Q[k * n + i], &Q[i * n + l]));
        TRY(rq_sub(&Q[k * n + l], &Q[k * n + l], &tmp));
      }
    }
  }

  // 4. [Find matrix R]
  for (u64 i = 0; i < n; i++) {
    TRY(rq_sqrt_exact(&Q[i * n + i], &Q[i * n + i]));

    for (u64 j = i + 1; j < n; j++) {
      TRY(rq_mul(&Q[i * n + j], &Q[i * n + j], &Q[i * n + i]));
    }
  }

  // Flush R into C on a single common denominator; C is untouched on error.
  TRY(rz_set_u64(&D, 1));
  for (u64 i = 0; i < n; i++) {
    for (u64 j = i; j < n; j++) {
      TRY(rz_lcm(&D, &D, &Q[i * n + j].den));
    }
  }

  for (u64 i = 0; i < n; i++) {
    for (u64 j = 0; j < n; j++) {
      if (j < i) {
        TRY(rz_set_u64(RMAT_GET(&C->M, i, j), 0));
        continue;
      }

      TRY(rz_div_exact(&s, &D, &Q[i * n + j].den));
      TRY(rz_mul(&t, &Q[i * n + j].num, &s));
      TRY(rz_copy(RMAT_GET(&C->M, i, j), &t));
    }
  }

  TRY(rz_copy(&C->den, &D));
  TRY(rmatq_normalize(C));

cleanup:
  if (Q != NULL) {
    for (u64 k = 0; k < inited; k++) rq_clear(&Q[k]);
    free(Q);
  }
  rq_clear(&tmp);
  rz_clear_multi(&D, &s, &t, NULL);
  return err;
}
