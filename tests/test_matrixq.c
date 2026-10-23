#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "../include/rmat.h"
#include "../include/rmatq.h"
#include "../include/rq.h"
#include "../include/rvec.h"
#include "../include/rz.h"

static void setq(rmatq_t* M, u64 i, u64 j, long long num, u64 den)
{
  rq_t v = {0};
  rq_t a = {0};
  rq_t b = {0};
  rz_t n = {0};
  rz_t d = {0};
  rq_init(&v);
  rq_init(&a);
  rq_init(&b);
  rz_init(&n);
  rz_init(&d);
  rz_set_i64(&n, num);
  rz_set_u64(&d, den);
  assert(rq_set(&a, &n) == 0);
  assert(rq_set(&b, &d) == 0);
  assert(rq_div(&v, &a, &b) == 0);
  assert(rmatq_set(M, &v, i, j) == 0);
  rq_clear(&v);
  rq_clear(&a);
  rq_clear(&b);
  rz_clear(&n);
  rz_clear(&d);
}

static bool qeq(const rmatq_t* M, u64 i, u64 j, long long num, u64 den)
{
  rq_t v = {0};
  rq_t want = {0};
  rq_t a = {0};
  rq_t b = {0};
  rz_t n = {0};
  rz_t d = {0};
  bool r = false;
  rq_init(&v);
  rq_init(&want);
  rq_init(&a);
  rq_init(&b);
  rz_init(&n);
  rz_init(&d);
  rz_set_i64(&n, num);
  rz_set_u64(&d, den);
  if (rq_set(&a, &n) != 0) goto done;
  if (rq_set(&b, &d) != 0) goto done;
  if (rq_div(&want, &a, &b) != 0) goto done;
  if (rmatq_get(&v, M, i, j) != 0) goto done;
  r = (rq_cmp(&v, &want) == 0);
done:
  rq_clear(&v);
  rq_clear(&want);
  rq_clear(&a);
  rq_clear(&b);
  rz_clear(&n);
  rz_clear(&d);
  return r;
}

int main(void)
{
  /* ---- zero / ones / id / predicates ---- */
  {
    rmatq_t A = {0};
    rmatq_init(&A, 2, 3);
    setq(&A, 0, 0, 5, 7);
    assert(rmatq_zero(&A) == RABIN_SUCCESS);
    assert(rmatq_is_zero(&A));
    assert(rmatq_ones(&A) == RABIN_SUCCESS);
    assert(!rmatq_is_zero(&A));
    assert(qeq(&A, 1, 2, 1, 1));
    assert(!rmatq_is_square(&A));
    assert(rmatq_zero(NULL) == RABIN_ERR_NULL_PTR);
    rmatq_clear(&A);

    rmatq_init(&A, 3, 3);
    assert(rmatq_id(&A, 3) == RABIN_SUCCESS);
    assert(rmatq_is_identity(&A));
    assert(!rmatq_is_zero(&A));
    assert(rmatq_is_square(&A));
    setq(&A, 0, 1, 1, 3);
    assert(!rmatq_is_identity(&A));
    assert(rmatq_id(NULL, 3) == RABIN_ERR_NULL_PTR);
    rmatq_clear(&A);
  }
  printf("[PASS] rmatq_zero/ones/id + predicates\n");

  /* ---- diag ---- */
  {
    rvec_t v = {0};
    rmatq_t D = {0};
    rmatq_init(&D, 3, 3);
    assert(rvec_init(&v, 3) == 0);
    rz_set_i64(&v.data[0], 2);
    rz_set_i64(&v.data[1], -4);
    rz_set_i64(&v.data[2], 0);
    assert(rmatq_diag(&D, &v) == RABIN_SUCCESS);
    assert(qeq(&D, 0, 0, 2, 1));
    assert(qeq(&D, 1, 1, -4, 1));
    assert(qeq(&D, 0, 2, 0, 1));
    assert(rmatq_is_symmetric(&D));
    {
      rmatq_t Bad = {0};
      rmatq_init(&Bad, 2, 2);
      assert(rmatq_diag(&Bad, &v) == RABIN_ERR_MATRIX_DIM);
      rmatq_clear(&Bad);
    }
    rvec_clear(&v);
    rmatq_clear(&D);
  }
  printf("[PASS] rmatq_diag\n");

  /* ---- swaps / transpose ---- */
  {
    rmatq_t A = {0}, T = {0};
    rmatq_init(&A, 2, 3);
    rmatq_init(&T, 3, 2);
    setq(&A, 0, 0, 1, 2);
    setq(&A, 0, 1, 3, 5);
    setq(&A, 0, 2, -7, 4);
    setq(&A, 1, 0, 2, 3);
    setq(&A, 1, 1, 8, 9);
    setq(&A, 1, 2, 1, 1);
    assert(rmatq_swap_rows(&A, 0, 1) == RABIN_SUCCESS);
    assert(qeq(&A, 0, 0, 2, 3));
    assert(rmatq_swap_cols(&A, 0, 2) == RABIN_SUCCESS);
    assert(qeq(&A, 0, 0, 1, 1));
    assert(rmatq_swap_rows(&A, 0, 9) == RABIN_ERR_INVALID_ARG);
    assert(rmatq_transpose(&T, &A) == RABIN_SUCCESS);
    assert(qeq(&T, 0, 0, 1, 1));
    assert(qeq(&T, 2, 0, 2, 3));
    {
      rmatq_t Wrong = {0};
      rmatq_init(&Wrong, 2, 2);
      assert(rmatq_transpose(&Wrong, &A) == RABIN_ERR_MATRIX_DIM);
      rmatq_clear(&Wrong);
    }
    assert(rmatq_transpose(NULL, &A) == RABIN_ERR_NULL_PTR);
    /* square self transpose */
    assert(rmatq_transpose(&A, &A) == RABIN_ERR_MATRIX_DIM);
    rmatq_clear(&A);
    rmatq_clear(&T);

    rmatq_init(&A, 2, 2);
    setq(&A, 0, 1, 5, 3);
    assert(rmatq_transpose(&A, &A) == RABIN_SUCCESS);
    assert(qeq(&A, 1, 0, 5, 3));
    assert(qeq(&A, 0, 1, 0, 1));
    rmatq_clear(&A);
  }
  printf("[PASS] rmatq swaps + transpose\n");

  /* ---- hstack / vstack with different denominators ---- */
  {
    rmatq_t P = {0}, Q = {0}, H = {0}, V = {0};
    rmatq_init(&P, 2, 2);
    rmatq_init(&Q, 2, 1);
    setq(&P, 0, 0, 1, 2);
    setq(&P, 0, 1, 1, 3);
    setq(&P, 1, 0, 1, 4);
    setq(&P, 1, 1, 1, 5);
    setq(&Q, 0, 0, 1, 6);
    setq(&Q, 1, 0, 1, 7);
    rmatq_init(&H, 2, 3);
    assert(rmatq_hstack(&H, &P, &Q) == RABIN_SUCCESS);
    assert(qeq(&H, 0, 2, 1, 6));
    assert(qeq(&H, 1, 0, 1, 4));
    rmatq_init(&V, 4, 2);
    assert(rmatq_vstack(&V, &P, &P) == RABIN_SUCCESS);
    assert(qeq(&V, 3, 1, 1, 5));
    assert(qeq(&V, 0, 0, 1, 2));
    assert(rmatq_hstack(&H, &P, &P) == RABIN_ERR_MATRIX_DIM);
    rmatq_clear(&P);
    rmatq_clear(&Q);
    rmatq_clear(&H);
    rmatq_clear(&V);
  }
  printf("[PASS] rmatq_hstack / rmatq_vstack\n");

  /* ---- scalar + normalize safety ---- */
  {
    rmatq_t S = {0}, T = {0};
    rq_t s = {0}, half = {0}, n1 = {0}, d1 = {0};
    rz_t a6 = {0}, b1 = {0};
    rmatq_init(&S, 2, 2);
    rmatq_init(&T, 2, 2);
    setq(&S, 0, 0, 1, 2);
    setq(&S, 0, 1, 3, 4);
    setq(&S, 1, 0, -1, 6);
    setq(&S, 1, 1, 5, 8);
    rq_init(&s);
    rq_init(&half);
    rq_init(&n1);
    rq_init(&d1);
    rz_init(&a6);
    rz_init(&b1);
    rq_set_u64(&s, 2);
    assert(rmatq_scalar(&T, &S, &s) == RABIN_SUCCESS);
    assert(qeq(&T, 0, 0, 1, 1));
    assert(qeq(&T, 0, 1, 3, 2));
    assert(rmatq_normalize(&T) == RABIN_SUCCESS);
    assert(qeq(&T, 0, 1, 3, 2));
    /* scalar by 1/2 recovers S */
    rq_set_u64(&n1, 1);
    rq_set_u64(&d1, 2);
    assert(rq_div(&half, &n1, &d1) == 0);
    assert(rmatq_scalar(&T, &T, &half) == RABIN_SUCCESS);
    assert(rmatq_equal(&T, &S));
    (void)a6;
    (void)b1;
    assert(rmatq_scalar(NULL, &S, &s) == RABIN_ERR_NULL_PTR);
    rq_clear(&s);
    rq_clear(&half);
    rq_clear(&n1);
    rq_clear(&d1);
    rz_clear(&a6);
    rz_clear(&b1);
    rmatq_clear(&S);
    rmatq_clear(&T);
  }
  printf("[PASS] rmatq_scalar + normalize\n");

  /* ---- pow vs repeated mul ---- */
  {
    rmatq_t S = {0}, P2 = {0}, P3 = {0}, M2 = {0}, M3 = {0};
    rmatq_init(&S, 2, 2);
    rmatq_init(&P2, 2, 2);
    rmatq_init(&P3, 2, 2);
    rmatq_init(&M2, 2, 2);
    rmatq_init(&M3, 2, 2);
    setq(&S, 0, 0, 1, 2);
    setq(&S, 0, 1, 3, 4);
    setq(&S, 1, 0, -1, 6);
    setq(&S, 1, 1, 5, 8);
    assert(rmatq_pow(&P2, &S, 2) == RABIN_SUCCESS);
    assert(rmatq_mul(&M2, &S, &S) == RABIN_SUCCESS);
    assert(rmatq_equal(&P2, &M2));
    assert(rmatq_pow(&P3, &S, 3) == RABIN_SUCCESS);
    assert(rmatq_mul(&M3, &M2, &S) == RABIN_SUCCESS);
    assert(rmatq_equal(&P3, &M3));
    assert(rmatq_pow(&P2, &S, 0) == RABIN_SUCCESS);
    assert(rmatq_is_identity(&P2));
    /* aliasing R == A */
    assert(rmatq_pow(&S, &S, 3) == RABIN_SUCCESS);
    assert(rmatq_equal(&S, &P3));
    rmatq_clear(&S);
    rmatq_clear(&P2);
    rmatq_clear(&P3);
    rmatq_clear(&M2);
    rmatq_clear(&M3);
  }
  printf("[PASS] rmatq_pow\n");

  /* ---- equal normalize-safe ---- */
  {
    rmatq_t A = {0}, B = {0};
    rmatq_init(&A, 2, 2);
    rmatq_init(&B, 2, 2);
    setq(&A, 0, 0, 1, 2);
    setq(&B, 0, 0, 5, 10);
    setq(&A, 0, 1, 2, 3);
    setq(&B, 0, 1, 8, 12);
    setq(&A, 1, 0, -1, 4);
    setq(&B, 1, 0, -3, 12);
    setq(&A, 1, 1, 0, 5);
    setq(&B, 1, 1, 0, 7);
    assert(rmatq_equal(&A, &B));
    setq(&B, 1, 1, 1, 7);
    assert(!rmatq_equal(&A, &B));
    assert(!rmatq_equal(NULL, &B));
    rmatq_clear(&A);
    rmatq_clear(&B);
  }
  printf("[PASS] rmatq_equal normalize-safe\n");

  /* ---- trace / det ---- */
  {
    rmatq_t S = {0};
    rq_t tr = {0}, d = {0}, want = {0}, a = {0}, b = {0};
    rmatq_init(&S, 2, 2);
    setq(&S, 0, 0, 1, 2);
    setq(&S, 0, 1, 3, 4);
    setq(&S, 1, 0, -1, 6);
    setq(&S, 1, 1, 5, 8);
    rq_init(&tr);
    rq_init(&d);
    rq_init(&want);
    rq_init(&a);
    rq_init(&b);
    assert(rmatq_trace(&tr, &S) == RABIN_SUCCESS);
    rq_set_u64(&a, 9);
    rq_set_u64(&b, 8);
    assert(rq_div(&want, &a, &b) == 0);
    assert(rq_cmp(&tr, &want) == 0);
    /* det = 1/2*5/8 - 3/4*(-1/6) = 5/16 + 1/8 = 7/16 */
    assert(rmatq_det(&d, &S) == RABIN_SUCCESS);
    rq_set_u64(&a, 7);
    rq_set_u64(&b, 16);
    assert(rq_div(&want, &a, &b) == 0);
    assert(rq_cmp(&d, &want) == 0);
    /* 1x1 convention: det([3/5]) = 3/5 */
    {
      rmatq_t One = {0};
      rq_t dd = {0}, w35 = {0};
      rmatq_init(&One, 1, 1);
      rq_init(&dd);
      rq_init(&w35);
      setq(&One, 0, 0, 3, 5);
      assert(rmatq_det(&dd, &One) == RABIN_SUCCESS);
      rq_set_u64(&a, 3);
      rq_set_u64(&b, 5);
      rq_div(&w35, &a, &b);
      assert(rq_cmp(&dd, &w35) == 0);
      rq_clear(&dd);
      rq_clear(&w35);
      rmatq_clear(&One);
    }
    assert(rmatq_trace(NULL, &S) == RABIN_ERR_NULL_PTR);
    rq_clear(&tr);
    rq_clear(&d);
    rq_clear(&want);
    rq_clear(&a);
    rq_clear(&b);
    rmatq_clear(&S);
  }
  printf("[PASS] rmatq_trace / rmatq_det\n");

  /* ---- rank / rref / lu ---- */
  {
    rmatq_t S = {0}, Z = {0}, Rr = {0}, L = {0}, U = {0}, PA = {0}, LU = {0};
    u64 rk = 999, Pp[2];
    rmatq_init(&S, 2, 2);
    rmatq_init(&Z, 2, 2);
    rmatq_init(&Rr, 2, 2);
    rmatq_init(&L, 2, 2);
    rmatq_init(&U, 2, 2);
    rmatq_init(&PA, 2, 2);
    rmatq_init(&LU, 2, 2);
    setq(&S, 0, 0, 1, 2);
    setq(&S, 0, 1, 3, 4);
    setq(&S, 1, 0, -1, 6);
    setq(&S, 1, 1, 5, 8);
    assert(rmatq_rank(&rk, &S) == RABIN_SUCCESS && rk == 2);
    assert(rmatq_copy(&Rr, &S) == RABIN_SUCCESS);
    assert(rmatq_rref(&Rr, &rk) == RABIN_SUCCESS && rk == 2);
    assert(rmatq_rref(&Rr, NULL) == RABIN_SUCCESS);
    assert(rmatq_lu(&L, &U, Pp, &S) == RABIN_SUCCESS);
    /* reconstruct P*S and compare with L*U */
    for (u64 i = 0; i < 2; i++) {
      for (u64 j = 0; j < 2; j++) {
        rq_t v = {0};
        rq_init(&v);
        assert(rmatq_get(&v, &S, Pp[i], j) == 0);
        assert(rmatq_set(&PA, &v, i, j) == 0);
        rq_clear(&v);
      }
    }
    assert(rmatq_mul(&LU, &L, &U) == RABIN_SUCCESS);
    assert(rmatq_equal(&LU, &PA));

    /* singular rational matrix */
    setq(&Z, 0, 0, 1, 3);
    setq(&Z, 0, 1, 2, 3);
    setq(&Z, 1, 0, 2, 3);
    setq(&Z, 1, 1, 4, 3);
    assert(rmatq_rank(&rk, &Z) == RABIN_SUCCESS && rk == 1);
    assert(rmatq_copy(&Rr, &Z) == RABIN_SUCCESS);
    assert(rmatq_rref(&Rr, &rk) == RABIN_SUCCESS && rk == 1);
    {
      rq_t v = {0};
      rq_init(&v);
      assert(rmatq_get(&v, &Rr, 1, 1) == 0 && rq_is_zero(&v));
      rq_clear(&v);
    }
    assert(rmatq_lu(&L, &U, Pp, &Z) == RABIN_SUCCESS);
    for (u64 i = 0; i < 2; i++) {
      for (u64 j = 0; j < 2; j++) {
        rq_t v = {0};
        rq_init(&v);
        assert(rmatq_get(&v, &Z, Pp[i], j) == 0);
        assert(rmatq_set(&PA, &v, i, j) == 0);
        rq_clear(&v);
      }
    }
    assert(rmatq_mul(&LU, &L, &U) == RABIN_SUCCESS);
    assert(rmatq_equal(&LU, &PA));
    rmatq_clear(&S);
    rmatq_clear(&Z);
    rmatq_clear(&Rr);
    rmatq_clear(&L);
    rmatq_clear(&U);
    rmatq_clear(&PA);
    rmatq_clear(&LU);
  }
  printf("[PASS] rmatq_rank / rref / lu\n");

  /* ---- solve consistency with rmatq_inv ---- */
  {
    rmatq_t A = {0}, B = {0}, X = {0}, Ainv = {0}, Xi = {0};
    rmatq_init(&A, 2, 2);
    rmatq_init(&B, 2, 1);
    rmatq_init(&X, 2, 1);
    rmatq_init(&Ainv, 2, 2);
    rmatq_init(&Xi, 2, 1);
    setq(&A, 0, 0, 1, 2);
    setq(&A, 0, 1, 3, 4);
    setq(&A, 1, 0, -1, 6);
    setq(&A, 1, 1, 5, 8);
    setq(&B, 0, 0, 1, 1);
    setq(&B, 1, 0, 0, 1);
    assert(rmatq_solve(&X, &A, &B) == RABIN_SUCCESS);
    assert(rmatq_inv(&Ainv, &A) == RABIN_SUCCESS);
    assert(rmatq_mul(&Xi, &Ainv, &B) == RABIN_SUCCESS);
    assert(rmatq_equal(&X, &Xi));
    /* verify A * X == B directly */
    {
      rmatq_t AX = {0};
      rmatq_init(&AX, 2, 1);
      assert(rmatq_mul(&AX, &A, &X) == RABIN_SUCCESS);
      assert(rmatq_equal(&AX, &B));
      rmatq_clear(&AX);
    }
    /* inconsistent system */
    {
      rmatq_t Y = {0}, YB = {0};
      rmatq_init(&Y, 2, 2);
      rmatq_init(&YB, 2, 1);
      setq(&Y, 0, 0, 1, 1);
      setq(&Y, 0, 1, 1, 1);
      setq(&Y, 1, 0, 2, 1);
      setq(&Y, 1, 1, 2, 1);
      setq(&YB, 0, 0, 1, 1);
      setq(&YB, 1, 0, 3, 1);
      assert(rmatq_solve(&X, &Y, &YB) == RABIN_ERR_INCONSISTENT);
      rmatq_clear(&Y);
      rmatq_clear(&YB);
    }
    /* singular but consistent */
    {
      rmatq_t Y = {0}, YB = {0};
      rmatq_init(&Y, 2, 2);
      rmatq_init(&YB, 2, 1);
      setq(&Y, 0, 0, 1, 1);
      setq(&Y, 0, 1, 1, 1);
      setq(&Y, 1, 0, 2, 1);
      setq(&Y, 1, 1, 2, 1);
      setq(&YB, 0, 0, 1, 1);
      setq(&YB, 1, 0, 2, 1);
      assert(rmatq_solve(&X, &Y, &YB) == RABIN_ERR_SINGULAR);
      rmatq_clear(&Y);
      rmatq_clear(&YB);
    }
    assert(rmatq_solve(NULL, &A, &B) == RABIN_ERR_NULL_PTR);
    rmatq_clear(&A);
    rmatq_clear(&B);
    rmatq_clear(&X);
    rmatq_clear(&Ainv);
    rmatq_clear(&Xi);
  }
  printf("[PASS] rmatq_solve == inv(A)*B\n");

  printf("test_matrixq: ALL PASSED\n");
  return 0;
}
