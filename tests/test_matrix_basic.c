#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "../include/rmat.h"
#include "../include/rmatq.h"
#include "../include/rvec.h"
#include "../include/rz.h"

static void fill(rmat_t* M, const long long vals[])
{
  for (u64 i = 0; i < M->rows; i++) {
    for (u64 j = 0; j < M->cols; j++) {
      assert(rz_set_i64(RMAT_GET(M, i, j), vals[i * M->cols + j]) == 0);
    }
  }
}

static bool eq_entry(const rmat_t* M, u64 i, u64 j, long long v)
{
  return rz_is_eq_i64(RMAT_GET(M, i, j), v);
}

int main(void)
{
  rmat_t A = {0}, B = {0}, R = {0}, E = {0};

  /* ---- rmat_zero / rmat_ones ---- */
  rmat_init(&A, 2, 3);
  {
    long long v[6] = {1, 2, 3, 4, 5, 6};
    fill(&A, v);
  }
  assert(rmat_zero(&A) == RABIN_SUCCESS);
  for (u64 i = 0; i < 6; i++) assert(rz_is_zero(&A.data[i]));
  assert(rmat_ones(&A) == RABIN_SUCCESS);
  for (u64 i = 0; i < 6; i++) assert(rz_is_eq_i64(&A.data[i], 1));
  assert(rmat_zero(NULL) == RABIN_ERR_NULL_PTR);
  assert(rmat_ones(NULL) == RABIN_ERR_NULL_PTR);
  printf("[PASS] rmat_zero / rmat_ones\n");

  /* ---- rmat_id ---- */
  rmat_init(&E, 3, 3);
  assert(rmat_id(&E, 3) == RABIN_SUCCESS);
  assert(rmat_is_identity(&E));
  printf("[PASS] rmat_id\n");

  /* ---- rmat_diag ---- */
  {
    rvec_t v = {0};
    rmat_init(&R, 3, 3);
    assert(rvec_init(&v, 3) == RABIN_SUCCESS);
    assert(rz_set_i64(&v.data[0], 2) == 0);
    assert(rz_set_i64(&v.data[1], -5) == 0);
    assert(rz_set_i64(&v.data[2], 7) == 0);
    assert(rmat_diag(&R, &v) == RABIN_SUCCESS);
    assert(eq_entry(&R, 0, 0, 2));
    assert(eq_entry(&R, 1, 1, -5));
    assert(eq_entry(&R, 2, 2, 7));
    assert(eq_entry(&R, 0, 1, 0));
    assert(eq_entry(&R, 2, 0, 0));
    /* dimension mismatch: result is not |v| x |v| */
    {
      rmat_t Bad = {0};
      rmat_init(&Bad, 2, 3);
      assert(rmat_diag(&Bad, &v) == RABIN_ERR_MATRIX_DIM);
      rmat_clear(&Bad);
    }
    assert(rmat_diag(NULL, &v) == RABIN_ERR_NULL_PTR);
    assert(rmat_diag(&R, NULL) == RABIN_ERR_NULL_PTR);
    rvec_clear(&v);
  }
  printf("[PASS] rmat_diag + dim errors\n");

  /* ---- rmat_swap_rows / rmat_swap_cols ---- */
  {
    long long m[6] = {1, 2, 3, 4, 5, 6};
    fill(&A, m);
    assert(rmat_swap_rows(&A, 0, 1) == RABIN_SUCCESS);
    assert(eq_entry(&A, 0, 0, 4) && eq_entry(&A, 1, 0, 1));
    assert(rmat_swap_cols(&A, 0, 2) == RABIN_SUCCESS);
    assert(eq_entry(&A, 0, 0, 6) && eq_entry(&A, 0, 2, 4));
    assert(rmat_swap_rows(&A, 1, 1) == RABIN_SUCCESS);
    assert(rmat_swap_rows(&A, 0, 5) == RABIN_ERR_INVALID_ARG);
    assert(rmat_swap_cols(&A, 0, 9) == RABIN_ERR_INVALID_ARG);
    assert(rmat_swap_rows(NULL, 0, 1) == RABIN_ERR_NULL_PTR);
    /* swap twice = identity */
    assert(rmat_swap_rows(&A, 0, 1) == RABIN_SUCCESS);
    assert(rmat_swap_rows(&A, 0, 1) == RABIN_SUCCESS);
    assert(eq_entry(&A, 0, 0, 6));
  }
  printf("[PASS] rmat_swap_rows / rmat_swap_cols + index errors\n");

  /* ---- rmat_hstack / rmat_vstack ---- */
  {
    long long a[6] = {1, 2, 3, 4, 5, 6};
    long long b[6] = {7, 8, 9, 10, 11, 12};
    rmat_t Bm = {0};
    fill(&A, a);
    rmat_init(&Bm, 2, 3);
    fill(&Bm, b);

    rmat_init(&R, 2, 6);
    assert(rmat_hstack(&R, &A, &Bm) == RABIN_SUCCESS);
    assert(R.rows == 2 && R.cols == 6);
    assert(eq_entry(&R, 0, 2, 3));
    assert(eq_entry(&R, 0, 3, 7));
    assert(eq_entry(&R, 1, 5, 12));

    rmat_t Vt = {0};
    rmat_init(&Vt, 4, 3);
    assert(rmat_vstack(&Vt, &A, &Bm) == RABIN_SUCCESS);
    assert(Vt.rows == 4 && Vt.cols == 3);
    assert(eq_entry(&Vt, 2, 0, 7));
    assert(eq_entry(&Vt, 3, 2, 12));

    /* dimension mismatch errors */
    {
      rmat_t Wide = {0};
      rmat_init(&Wide, 3, 2);
      assert(rmat_hstack(&R, &A, &Wide) == RABIN_ERR_MATRIX_DIM);
      assert(rmat_vstack(&Vt, &A, &Wide) == RABIN_ERR_MATRIX_DIM);
      {
        rmat_t SmallH = {0};
        rmat_init(&SmallH, 2, 5);
        assert(rmat_hstack(&SmallH, &A, &Bm) == RABIN_ERR_MATRIX_DIM);
        rmat_clear(&SmallH);
      }
      rmat_clear(&Wide);
    }
    assert(rmat_hstack(NULL, &A, &Bm) == RABIN_ERR_NULL_PTR);
    rmat_clear(&Bm);
    rmat_clear(&Vt);
  }
  printf("[PASS] rmat_hstack / rmat_vstack + dim errors\n");

  /* ---- rmat_pow ---- */
  {
    long long q[4] = {1, 1, 1, 0};
    rmat_t Q = {0}, Qn = {0};
    rmat_init(&Q, 2, 2);
    rmat_init(&Qn, 2, 2);
    fill(&Q, q);
    assert(rmat_pow(&Qn, &Q, 10) == RABIN_SUCCESS);
    assert(eq_entry(&Qn, 0, 0, 89)); /* F(11) */
    assert(eq_entry(&Qn, 0, 1, 55)); /* F(10) */
    assert(eq_entry(&Qn, 1, 0, 55));
    assert(eq_entry(&Qn, 1, 1, 34)); /* F(9) */
    assert(rmat_pow(&Qn, &Q, 0) == RABIN_SUCCESS);
    assert(rmat_is_identity(&Qn));
    assert(rmat_pow(&Qn, &Q, 1) == RABIN_SUCCESS);
    assert(rmat_equal(&Qn, &Q));
    /* alias R == A */
    assert(rmat_pow(&Q, &Q, 3) == RABIN_SUCCESS);
    {
      rmat_t Q2 = {0}, Q3 = {0};
      rmat_init(&Q2, 2, 2);
      rmat_init(&Q3, 2, 2);
      long long q2[4] = {1, 1, 1, 0};
      fill(&Q2, q2);
      assert(rmat_mul(&Q3, &Q2, &Q2) == RABIN_SUCCESS);
      assert(rmat_mul(&Q2, &Q3, &Q2) == RABIN_SUCCESS);
      assert(rmat_equal(&Q, &Q2));
      rmat_clear(&Q2);
      rmat_clear(&Q3);
    }
    /* non-square base fails */
    {
      rmat_t Rect = {0};
      rmat_init(&Rect, 2, 3);
      assert(rmat_pow(&Rect, &Rect, 2) == RABIN_ERR_MATRIX_DIM);
      rmat_clear(&Rect);
    }
    rmat_clear(&Q);
    rmat_clear(&Qn);
  }
  printf("[PASS] rmat_pow\n");

  /* ---- rmat_is_singular ---- */
  {
    rmat_t S = {0};
    rz_t d = {0};
    long long s[4] = {2, 4, 1, 2};
    long long t[4] = {2, 4, 1, 3};
    rmat_init(&S, 2, 2);
    rz_init(&d);
    fill(&S, s);
    assert(rmat_is_singular(&S) == true);
    assert(rmat_det(&d, &S) == RABIN_SUCCESS && rz_is_zero(&d));
    fill(&S, t);
    assert(rmat_is_singular(&S) == false);
    assert(rmat_det(&d, &S) == RABIN_SUCCESS && !rz_is_zero(&d));
    assert(rmat_is_singular(&A) == false); /* 2x3 non-square: false */
    rmat_clear(&S);
    rz_clear(&d);
  }
  printf("[PASS] rmat_is_singular\n");

  rmat_clear(&A);
  rmat_clear(&B);
  rmat_clear(&R);
  rmat_clear(&E);

  printf("test_matrix_basic: ALL PASSED\n");
  return 0;
}
