#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "../include/rmat.h"
#include "../include/rmatq.h"
#include "../include/rz.h"

static void fill(rmat_t* M, const long long vals[])
{
  for (u64 i = 0; i < M->rows; i++) {
    for (u64 j = 0; j < M->cols; j++) {
      assert(rz_set_i64(RMAT_GET(M, i, j), vals[i * M->cols + j]) == 0);
    }
  }
}

static u64 rs = 0x123456789abcdefULL;
static long long rnd_small(long long mod)
{
  rs ^= rs << 13;
  rs ^= rs >> 7;
  rs ^= rs << 17;
  return (long long)(rs % (u64)(2 * mod + 1)) - mod;
}

static void fill_random(rmat_t* M, long long mod)
{
  for (u64 i = 0; i < M->rows; i++) {
    for (u64 j = 0; j < M->cols; j++) {
      assert(rz_set_i64(RMAT_GET(M, i, j), rnd_small(mod)) == 0);
    }
  }
}

/* verify L * U == P * A for the factors of rmat_lu */
static void check_lu(const rmat_t* A, const rmat_t* L, const rmat_t* U,
                     const u64* P, u64 n)
{
  rmat_t PA = {0}, LU = {0};
  bool seen[64];
  assert(n <= 64);
  for (u64 i = 0; i < n; i++) seen[i] = false;
  for (u64 i = 0; i < n; i++) {
    assert(P[i] < n);
    assert(!seen[P[i]]);
    seen[P[i]] = true;
    assert(rz_is_eq_i64(RMAT_GET(L, i, i), 1));
    for (u64 j = i + 1; j < n; j++) {
      assert(rz_is_zero(RMAT_GET(L, i, j)));
    }
    for (u64 j = 0; j < i; j++) {
      assert(rz_is_zero(RMAT_GET(U, i, j)));
    }
  }
  rmat_init(&PA, n, n);
  rmat_init(&LU, n, n);
  for (u64 i = 0; i < n; i++) {
    for (u64 j = 0; j < n; j++) {
      assert(rz_copy(RMAT_GET(&PA, i, j), RMAT_GET(A, P[i], j)) == 0);
    }
  }
  assert(rmat_mul(&LU, L, U) == RABIN_SUCCESS);
  assert(rmat_equal(&LU, &PA));
  rmat_clear(&PA);
  rmat_clear(&LU);
}

int main(void)
{
  /* matrices used across sections */
  const long long a[9] = {2, 1, 3, 4, 2, 6, 1, 0, 5};  /* rank 2 */
  const long long b[9] = {1, 2, 3, 4, 5, 6, 7, 8, 10}; /* rank 3 */
  const long long nst[9] = {2, 1, 1, 1, 1, 1, 1, 1, 2};

  /* ---- rank ---- */
  {
    rmat_t M = {0};
    u64 rk = 999;
    rmat_init(&M, 3, 3);
    fill(&M, b);
    assert(rmat_rank(&rk, &M) == RABIN_SUCCESS && rk == 3);
    fill(&M, a);
    assert(rmat_rank(&rk, &M) == RABIN_SUCCESS && rk == 2);
    /* rows that are multiples */
    {
      const long long mult[6] = {3, -6, 3, 9, -18, 9};
      rmat_t R2 = {0};
      rmat_init(&R2, 2, 3);
      fill(&R2, mult);
      assert(rmat_rank(&rk, &R2) == RABIN_SUCCESS && rk == 1);
      rmat_clear(&R2);
    }
    /* zero matrix */
    rmat_zero(&M);
    assert(rmat_rank(&rk, &M) == RABIN_SUCCESS && rk == 0);
    /* leading column zero: pivoting needed for rank */
    {
      const long long piv[4] = {0, 1, 0, 2};
      rmat_t P4 = {0};
      rmat_init(&P4, 2, 2);
      fill(&P4, piv);
      assert(rmat_rank(&rk, &P4) == RABIN_SUCCESS && rk == 1);
      rmat_clear(&P4);
    }
    assert(rmat_rank(NULL, &M) == RABIN_ERR_NULL_PTR);
    assert(rmat_rank(&rk, NULL) == RABIN_ERR_NULL_PTR);
    rmat_clear(&M);
  }
  printf("[PASS] rmat_rank\n");

  /* ---- rref ---- */
  {
    rmat_t M = {0};
    rmatq_t Rr = {0};
    u64 rk = 999;
    rmat_init(&M, 3, 3);
    rmatq_init(&Rr, 3, 3);

    fill(&M, b);
    assert(rmat_rref(&Rr, &rk, &M) == RABIN_SUCCESS && rk == 3);
    /* square full rank -> identity */
    {
      rq_t v = {0};
      rq_init(&v);
      for (u64 i = 0; i < 3; i++) {
        for (u64 j = 0; j < 3; j++) {
          assert(rmatq_get(&v, &Rr, i, j) == 0);
          if (i == j) {
            assert(rq_is_zero(&v) == false);
            assert(rz_is_one(&v.num) && rz_is_one(&v.den));
          } else {
            assert(rq_is_zero(&v));
          }
        }
      }
      rq_clear(&v);
    }

    fill(&M, a);
    assert(rmat_rref(&Rr, &rk, &M) == RABIN_SUCCESS && rk == 2);
    {
      rq_t v = {0};
      rq_init(&v);
      assert(rmatq_get(&v, &Rr, 2, 0) == 0 && rq_is_zero(&v));
      assert(rmatq_get(&v, &Rr, 0, 0) == 0);
      assert(rz_is_one(&v.num));
      rq_clear(&v);
    }

    /* NULL rank allowed */
    assert(rmat_rref(&Rr, NULL, &M) == RABIN_SUCCESS);
    assert(rmat_rref(NULL, &rk, &M) == RABIN_ERR_NULL_PTR);
    rmat_clear(&M);
    rmatq_clear(&Rr);
  }
  printf("[PASS] rmat_rref\n");

  /* ---- rmat_lu ---- */
  {
    rmat_t M = {0}, L = {0}, U = {0};
    u64 P[3];
    rmat_init(&M, 3, 3);
    rmat_init(&L, 3, 3);
    rmat_init(&U, 3, 3);
    fill(&M, b);
    assert(rmat_lu(&L, &U, P, &M) == RABIN_SUCCESS);
    check_lu(&M, &L, &U, P, 3);
    fill(&M, nst);
    assert(rmat_lu(&L, &U, P, &M) == RABIN_SUCCESS);
    check_lu(&M, &L, &U, P, 3);
    /* singular: still a valid factorization, zero on U diagonal */
    fill(&M, a);
    assert(rmat_lu(&L, &U, P, &M) == RABIN_SUCCESS);
    check_lu(&M, &L, &U, P, 3);
    {
      bool zero_diag = false;
      for (u64 i = 0; i < 3; i++) {
        if (rz_is_zero(RMAT_GET(&U, i, i))) zero_diag = true;
      }
      assert(zero_diag);
    }
    /* pivoting required: leading entry zero */
    {
      const long long piv[4] = {0, 1, 1, 0};
      rmat_t M2 = {0}, L2 = {0}, U2 = {0};
      u64 P2[2];
      rmat_init(&M2, 2, 2);
      rmat_init(&L2, 2, 2);
      rmat_init(&U2, 2, 2);
      fill(&M2, piv);
      assert(rmat_lu(&L2, &U2, P2, &M2) == RABIN_SUCCESS);
      check_lu(&M2, &L2, &U2, P2, 2);
      assert(P2[0] != 0); /* row swap happened */
      rmat_clear(&M2);
      rmat_clear(&L2);
      rmat_clear(&U2);
    }
    assert(rmat_lu(NULL, &U, P, &M) == RABIN_ERR_NULL_PTR);
    rmat_clear(&M);
    rmat_clear(&L);
    rmat_clear(&U);
  }
  printf("[PASS] rmat_lu + LU recombination\n");

  /* ---- rmat_solve ---- */
  {
    rmat_t A = {0}, B = {0};
    rmatq_t X = {0};
    rmat_init(&A, 3, 3);
    rmat_init(&B, 3, 1);
    rmatq_init(&X, 3, 1);
    fill(&A, b);
    {
      const long long bb[3] = {3, 6, 9};
      fill(&B, bb);
      assert(rmat_solve(&X, &A, &B) == RABIN_SUCCESS);
      /* verify A * X == B over Q using copyZ + rmatq_mul */
      rmatq_t Aq = {0}, Bq = {0}, AX = {0};
      rmatq_init(&Aq, 3, 3);
      rmatq_init(&Bq, 3, 1);
      rmatq_init(&AX, 3, 1);
      assert(rmatq_copyZ(&Aq, &A) == 0);
      assert(rmatq_copyZ(&Bq, &B) == 0);
      assert(rmatq_mul(&AX, &Aq, &X) == 0);
      assert(rmatq_equal(&AX, &Bq));
      rmatq_clear(&Aq);
      rmatq_clear(&Bq);
      rmatq_clear(&AX);
    }
    /* singular but consistent */
    fill(&A, a);
    {
      const long long bb[3] = {3, 6, 9};
      fill(&B, bb);
      assert(rmat_solve(&X, &A, &B) == RABIN_ERR_SINGULAR);
    }
    /* inconsistent */
    {
      const long long ones[9] = {1, 1, 1, 1, 1, 1, 1, 1, 1};
      const long long bb[3] = {1, 2, 1};
      fill(&A, ones);
      fill(&B, bb);
      assert(rmat_solve(&X, &A, &B) == RABIN_ERR_INCONSISTENT);
    }
    assert(rmat_solve(NULL, &A, &B) == RABIN_ERR_NULL_PTR);
    rmat_clear(&A);
    rmat_clear(&B);
    rmatq_clear(&X);
  }
  printf("[PASS] rmat_solve (unique / singular / inconsistent)\n");

  /* ---- cofactor / adj ---- */
  {
    rmat_t A = {0}, Cof = {0}, ADJ = {0}, CHK = {0}, DI = {0};
    rz_t d = {0};
    rmat_init(&A, 3, 3);
    rmat_init(&Cof, 3, 3);
    rmat_init(&ADJ, 3, 3);
    rmat_init(&CHK, 3, 3);
    rmat_init(&DI, 3, 3);
    rz_init(&d);
    const long long* mats[2] = {b, a};
    for (int t = 0; t < 2; t++) {
      fill(&A, mats[t]);
      assert(rmat_cofactor(&Cof, &A) == RABIN_SUCCESS);
      assert(rmat_adj(&ADJ, &A) == RABIN_SUCCESS);
      assert(rmat_det(&d, &A) == RABIN_SUCCESS);
      assert(rmat_mul(&CHK, &A, &ADJ) == RABIN_SUCCESS);
      assert(rmat_id(&DI, 3) == RABIN_SUCCESS);
      for (u64 i = 0; i < 3; i++) {
        for (u64 j = 0; j < 3; j++) {
          assert(rz_mul(RMAT_GET(&DI, i, j), RMAT_GET(&DI, i, j), &d) == 0);
        }
      }
      assert(rmat_equal(&CHK, &DI));
    }
    /* 2x2 hand case: adj([[a,b],[c,d]]) = [[d,-b],[-c,a]] */
    {
      const long long m[4] = {3, 5, -2, 7};
      const long long want[4] = {7, -5, 2, 3};
      rmat_t A2 = {0}, J2 = {0};
      rmat_init(&A2, 2, 2);
      rmat_init(&J2, 2, 2);
      fill(&A2, m);
      assert(rmat_adj(&J2, &A2) == RABIN_SUCCESS);
      {
        rmat_t W = {0};
        rmat_init(&W, 2, 2);
        fill(&W, want);
        assert(rmat_equal(&J2, &W));
        rmat_clear(&W);
      }
      rmat_clear(&A2);
      rmat_clear(&J2);
    }
    assert(rmat_adj(NULL, &A) == RABIN_ERR_NULL_PTR);
    rmat_clear(&A);
    rmat_clear(&Cof);
    rmat_clear(&ADJ);
    rmat_clear(&CHK);
    rmat_clear(&DI);
    rz_clear(&d);
  }
  printf("[PASS] rmat_cofactor / rmat_adj\n");

  /* ---- rmat_det vs rmat_det_bareiss on random 4x4..6x6 ---- */
  {
    for (u64 n = 4; n <= 6; n++) {
      for (int trial = 0; trial < 3; trial++) {
        rmat_t M = {0};
        rz_t d1 = {0}, d2 = {0};
        rmat_init(&M, n, n);
        rz_init(&d1);
        rz_init(&d2);
        fill_random(&M, 30);
        assert(rmat_det(&d1, &M) == RABIN_SUCCESS);
        assert(rmat_det_bareiss(&d2, &M) == RABIN_SUCCESS);
        assert(rz_cmp(&d1, &d2) == 0);
        rmat_clear(&M);
        rz_clear(&d1);
        rz_clear(&d2);
      }
    }
  }
  printf("[PASS] rmat_det == rmat_det_bareiss (random, fixed seed)\n");

  printf("test_matrix_elim: ALL PASSED\n");
  return 0;
}
