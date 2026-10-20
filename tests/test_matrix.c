#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "../include/rmat.h"
#include "../include/rmatq.h"
#include "../include/rns.h"
#include "../include/rns_primes.h"
#include "../include/rpol.h"
#include "../include/rq.h"
#include "../include/rzmath.h"
#include "../include/rzrand.h"
#include "../include/u64.h"
double get_time()
{
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return ts.tv_sec + ts.tv_nsec * 1e-9;
}

void run_det_benchmark(u64 size, u64 bits)
{
  rmat_t M = {0};
  rz_t det, val;

  rz_init(&det);
  rz_init(&val);
  rmat_init(&M, size, size);

  rns_ctx_t ctx = {0};

  // Populate with random data to prevent "easy" zeros
  for (u64 i = 0; i < size; i++) {
    for (u64 j = 0; j < size; j++) {
      rz_gen_random(&val, bits);
      rmat_set(&M, &val, i, j);
    }
  }

  printf("Benchmarking %llu x %llu (%llu-bit entries) Determinant... ", size,
         size, bits);
  fflush(stdout);

  // rmat_print_python(&M);

  u64 k = rns_estimate_det(&M);
  // printf("k: %llu \n", k);

  if ((k + 1) > 500) {
    if (rns_ctx_init(&ctx, RNS_PRIMES2, k + 1) != RABIN_SUCCESS) return;
  } else {
    if (rns_ctx_init(&ctx, RNS_PRIMES, k + 1) != RABIN_SUCCESS) return;
  }

  double start = get_time();
  rmat_det_rns(&det, &M, &ctx);
  double end = get_time();

  printf("Time: %f seconds\n", end - start);

  // printf("rns: ");
  // rz_println(&det);
  printf("size of result: %llu bits \n", rz_bit_length(&det));

  rmat_clear(&M);
  rz_clear(&det);
  rz_clear(&val);
  rns_ctx_clear(&ctx);
}

void run_hadamard_benchmark(u64 size, u64 bits)
{
  rmat_t M = {0};
  rz_t det, val;

  rz_init(&det);
  rz_init(&val);
  rmat_init(&M, size, size);

  // Populate with random data to prevent "easy" zeros
  for (u64 i = 0; i < size; i++) {
    for (u64 j = 0; j < size; j++) {
      rz_gen_random(&val, bits);
      rmat_set(&M, &val, i, j);
    }
  }

  printf("Benchmarking %llu x %llu (%llu-bit entries) hadamard bound... ", size,
         size, bits);
  fflush(stdout);

  double start = get_time();
  rmat_hadamard(&det, &M);
  double end = get_time();

  printf("Time: %f seconds\n", end - start);

  rmat_clear(&M);
  rz_clear(&det);
  rz_clear(&val);
}

void run_hermite_benchmark(u64 size, u64 bits)
{
  rmat_t M = {0}, H = {0};
  rz_t det, val;

  rz_init(&det);
  rz_init(&val);
  rmat_init(&M, size, size);
  rmat_init(&H, size, size);

  // Populate with random data to prevent "easy" zeros
  for (u64 i = 0; i < size; i++) {
    for (u64 j = 0; j < size; j++) {
      rz_gen_random(&val, bits);
      rmat_set(&M, &val, i, j);
    }
  }

  rmat_det(&det, &M);

  det.is_neg = false;
  printf(
      "Benchmarking %llu x %llu (%llu-bit entries) Hermite normal form...\n ",
      size, size, bits);

  double start = get_time();
  rmat_hermite_mod_d(&H, &M, &det);
  double end = get_time();

  // rmat_print_tail(&H);

  printf("Time: %f seconds\n", end - start);

  rmat_clear(&M);
  rmat_clear(&H);
  rz_clear(&det);
  rz_clear(&val);
}

void run_mul_benchmark(u64 size, u64 bits)
{
  rmat_t M = {0}, A = {0}, B = {0};
  rz_t det, val;

  rz_init(&det);
  rz_init(&val);
  rmat_init(&M, size, size);
  rmat_init(&A, size, size);
  rmat_init(&B, size, size);

  // Populate with random data to prevent "easy" zeros
  for (u64 i = 0; i < size; i++) {
    for (u64 j = 0; j < size; j++) {
      rz_gen_random(&val, bits);
      rmat_set(&A, &val, i, j);
      rz_gen_random(&val, bits);
      rmat_set(&B, &val, i, j);
    }
  }

  printf("Benchmarking %llu x %llu (%llu-bit entries) matrix multiply ... ",
         size, size, bits);
  fflush(stdout);

  double start = get_time();
  rmat_mul(&M, &A, &B);
  double end = get_time();

  printf("Time: %f seconds\n", end - start);

  rmat_clear(&M);
  rmat_clear(&A);
  rmat_clear(&B);
  rz_clear(&det);
  rz_clear(&val);
}

void run_hermite_mod_d_benchmark(u64 size, u64 bits, bool with_exact)
{
  rmat_t M = {0}, H = {0};
  rz_t D, val;

  rz_init(&D);
  rz_init(&val);
  rmat_init(&M, size, size);
  rmat_init(&H, size, size);

  // Populate with random data to prevent "easy" zeros
  for (u64 i = 0; i < size; i++) {
    for (u64 j = 0; j < size; j++) {
      rz_gen_random(&val, bits);
      rmat_set(&M, &val, i, j);
    }
  }

  // D must be a multiple of the module determinant; for square M the
  // module determinant is |det M|, so use that
  for (int rep = 0; rep < 20; rep++) {
    rmat_det(&D, &M);
    if (!rz_is_zero(&D)) break;
    for (u64 i = 0; i < size; i++) {
      rz_gen_random(&val, bits);
      rmat_set(&M, &val, i, i);
    }
  }
  D.is_neg = false;

  if (with_exact) {
    printf(
        "Benchmarking %llu x %llu (%llu-bit entries) Hermite normal form "
        "(exact 2.4.5)... ",
        size, size, bits);
    fflush(stdout);

    double start = get_time();
    rmat_hermite_gcd(&H, &M);
    double end = get_time();

    printf("Time: %f seconds\n", end - start);
  }

  printf(
      "Benchmarking %llu x %llu (%llu-bit entries) Hermite normal form "
      "(mod D 2.4.8)... ",
      size, size, bits);
  fflush(stdout);

  double start = get_time();
  rmat_hermite_mod_d(&H, &M, &D);
  double end = get_time();

  printf("Time: %f seconds\n", end - start);

  rmat_clear(&M);
  rmat_clear(&H);
  rz_clear(&D);
  rz_clear(&val);
}

static void test_set_i64(rmat_t* M, u64 r, u64 c, i64 v)
{
  rz_t tmp;
  rz_init(&tmp);
  rz_set_i64(&tmp, v);
  rmat_set(M, &tmp, r, c);
  rz_clear(&tmp);
}

void test_hermite_mod_d(void)
{
  rz_t D, val;
  rz_init_multi(&D, &val, NULL);

  // 1. 2x2 with known HNF: A = [[1,2],[3,4]] -> W = [[2,1],[0,1]]
  {
    rmat_t A = {0}, W = {0}, H = {0};
    rmat_init(&A, 2, 2);
    rmat_init(&W, 2, 2);
    rmat_init(&H, 2, 2);
    test_set_i64(&A, 0, 0, 1);
    test_set_i64(&A, 0, 1, 2);
    test_set_i64(&A, 1, 0, 3);
    test_set_i64(&A, 1, 1, 4);

    // D = 4 = 2*|det A|, a multiple of the module determinant
    rz_set_u64(&D, 4);
    rmat_hermite_mod_d(&W, &A, &D);
    test_set_i64(&H, 0, 0, 2);
    test_set_i64(&H, 0, 1, 1);
    test_set_i64(&H, 1, 0, 0);
    test_set_i64(&H, 1, 1, 1);
    assert(rmat_equal(&W, &H));
    assert(rmat_hnf_check_structure(&W));

    // cross-check against the exact implementation
    rmat_t V = {0};
    rmat_init(&V, 2, 2);
    rmat_hermite(&V, &A);
    assert(rmat_equal(&W, &V));
    rmat_clear(&V);

    rmat_clear(&A);
    rmat_clear(&W);
    rmat_clear(&H);
    printf("[PASS] hermite_mod_d 2x2 known example\n");
  }

  // 2. 2x3 with unimodular column module: gcd of the 2x2 minors is 1,
  //    so D = 1 is valid and the HNF must be the 2x2 identity
  {
    rmat_t A = {0}, W = {0}, H = {0};
    rmat_init(&A, 2, 3);
    rmat_init(&W, 2, 2);
    rmat_init(&H, 2, 2);
    test_set_i64(&A, 0, 0, 6);
    test_set_i64(&A, 0, 1, 4);
    test_set_i64(&A, 0, 2, 9);
    test_set_i64(&A, 1, 0, 8);
    test_set_i64(&A, 1, 1, 1);
    test_set_i64(&A, 1, 2, 6);

    rz_set_u64(&D, 1);
    rmat_hermite_mod_d(&W, &A, &D);
    rmat_id(&H, 2);
    assert(rmat_equal(&W, &H));
    assert(rmat_hnf_check_structure(&W));
    printf("[PASS] hermite_mod_d 2x3 unimodular (D = 1)\n");

    rmat_clear(&A);
    rmat_clear(&W);
    rmat_clear(&H);
  }

  // 3. Random cross-checks against the exact implementation (HNF is unique).
  //    D must be a multiple of the module determinant Delta (GCD of the
  //    m x m minors); build A = [B | B R] so that Delta = |det B|:
  {
    u64 shapes[][2] = {{2, 2}, {2, 3}, {3, 3}, {3, 4}, {4, 4}, {2, 5}};
    u64 n_shapes = sizeof(shapes) / sizeof(shapes[0]);
    rz_t det, three;
    rz_init(&det);
    rz_set_u64(&three, 3);
    for (u64 s = 0; s < n_shapes; s++) {
      u64 m = shapes[s][0];
      u64 n = shapes[s][1];
      for (u64 trial = 0; trial < 4; trial++) {
        rmat_t B = {0}, R = {0}, C = {0}, A = {0}, H = {0}, W1 = {0}, W2 = {0};
        rmat_init(&B, m, m);
        rmat_init(&R, m, n - m);
        rmat_init(&C, m, n - m);
        rmat_init(&A, m, n);
        rmat_init(&H, m, m);
        rmat_init(&W1, m, m);
        rmat_init(&W2, m, m);

        for (int rep = 0; rep < 10; rep++) {
          for (u64 i = 0; i < m; i++) {
            for (u64 j = 0; j < m; j++) {
              rz_gen_random(&val, 20);
              rmat_set(&B, &val, i, j);
            }
          }
          rmat_det(&det, &B);
          if (!rz_is_zero(&det)) break;
        }
        assert(!rz_is_zero(&det));
        det.is_neg = false;

        for (u64 i = 0; i < m; i++) {
          for (u64 j = 0; j < n - m; j++) {
            rz_gen_random(&val, 20);
            rmat_set(&R, &val, i, j);
          }
        }
        if (m < n) {
          rmat_mul(&C, &B, &R);
        }
        // A = [B | B R]
        for (u64 i = 0; i < m; i++) {
          for (u64 j = 0; j < m; j++) {
            rmat_set(&A, &B.data[i * m + j], i, j);
          }
          for (u64 j = 0; j < n - m; j++) {
            rmat_set(&A, &C.data[i * (n - m) + j], i, m + j);
          }
        }

        rmat_hermite(&H, &A);
        assert(H.rows == m && H.cols == m);

        // D = |det B| = the module determinant itself
        rmat_hermite_mod_d(&W1, &A, &det);
        assert(rmat_equal(&W1, &H));
        assert(rmat_hnf_check_structure(&W1));
        assert(rmat_hnf_check_structure(&H));

        // D = 3 * |det B|, another valid multiple
        rz_mul(&D, &three, &det);
        rmat_hermite_mod_d(&W2, &A, &D);
        assert(rmat_equal(&W2, &H));
        assert(rmat_hnf_check_structure(&W2));

        rmat_clear(&B);
        rmat_clear(&R);
        rmat_clear(&C);
        rmat_clear(&A);
        rmat_clear(&H);
        rmat_clear(&W1);
        rmat_clear(&W2);
      }
      printf("[PASS] hermite_mod_d random %llux%llu vs exact\n", m, n);
    }
    rz_clear(&det);
    rz_clear(&three);
  }

  // 4. m > n falls back to the exact implementation
  {
    u64 m = 4, n = 3;
    rmat_t A = {0}, W = {0}, H = {0};
    rmat_init(&A, m, n);
    rmat_init(&W, m, n);
    rmat_init(&H, m, n);
    for (u64 i = 0; i < m; i++) {
      for (u64 j = 0; j < n; j++) {
        rz_gen_random(&val, 20);
        rmat_set(&A, &val, i, j);
      }
    }
    rz_set_u64(&D, 1000003);
    rmat_hermite_mod_d(&W, &A, &D);
    rmat_hermite(&H, &A);
    assert(rmat_equal(&W, &H));
    printf("[PASS] hermite_mod_d m > n fallback\n");

    rmat_clear(&A);
    rmat_clear(&W);
    rmat_clear(&H);
  }

  rz_clear_multi(&D, &val, NULL);
}

void test_hnf_verify(void)
{
  rz_t val;
  rz_init(&val);

  // 1. Structure: valid 2x2 with one violation per rule
  {
    // [[3,1],[0,2]]: upper triangular, positive pivots, 0 <= H[0][1] = 1 < 2
    rmat_t H = {0};
    rmat_init(&H, 2, 2);
    test_set_i64(&H, 0, 0, 3);
    test_set_i64(&H, 0, 1, 1);
    test_set_i64(&H, 1, 1, 2);
    assert(rmat_hnf_check_structure(&H));

    // violation: entry below the diagonal
    test_set_i64(&H, 1, 0, 5);
    assert(!rmat_hnf_check_structure(&H));
    test_set_i64(&H, 1, 0, 0);

    // violation: negative pivot
    test_set_i64(&H, 1, 1, -2);
    assert(!rmat_hnf_check_structure(&H));
    test_set_i64(&H, 1, 1, 2);

    // violation: entry right of the pivot not reduced (3 is not < 3)
    test_set_i64(&H, 0, 1, 3);
    assert(!rmat_hnf_check_structure(&H));

    // violation: negative entry right of the pivot
    test_set_i64(&H, 0, 1, -1);
    assert(!rmat_hnf_check_structure(&H));
    test_set_i64(&H, 0, 1, 1);
    assert(rmat_hnf_check_structure(&H));
    rmat_clear(&H);
    printf("[PASS] hnf_check_structure (2x2 rules)\n");

    // zero-row order: trailing zeros ok, a nonzero row after them is not
    rmat_t Z = {0};
    rmat_init(&Z, 3, 3);
    test_set_i64(&Z, 0, 0, 1);
    assert(rmat_hnf_check_structure(&Z));
    test_set_i64(&Z, 2, 2, 3);
    assert(!rmat_hnf_check_structure(&Z));
    rmat_clear(&Z);

    // zero pivot with a nonzero entry right of it: not a zero row
    rmat_t B = {0};
    rmat_init(&B, 2, 2);
    test_set_i64(&B, 0, 1, 1);
    assert(!rmat_hnf_check_structure(&B));
    rmat_clear(&B);

    // valid 3x3: [[2,1,0],[0,3,1],[0,0,5]]
    rmat_t T = {0};
    rmat_init(&T, 3, 3);
    test_set_i64(&T, 0, 0, 2);
    test_set_i64(&T, 0, 1, 1);
    test_set_i64(&T, 1, 1, 3);
    test_set_i64(&T, 1, 2, 1);
    test_set_i64(&T, 2, 2, 5);
    assert(rmat_hnf_check_structure(&T));
    rmat_clear(&T);
    printf("[PASS] hnf_check_structure (zero rows, 3x3)\n");
  }

  // 2. Unimodularity: det 1, det -1, shear; det 2 and non-square fail
  {
    rmat_t U = {0};
    rmat_init(&U, 2, 2);

    // identity: det 1
    test_set_i64(&U, 0, 0, 1);
    test_set_i64(&U, 1, 1, 1);
    assert(rmat_hnf_check_unimodular(&U));

    // swap: det -1, |det| = 1
    test_set_i64(&U, 0, 0, 0);
    test_set_i64(&U, 0, 1, 1);
    test_set_i64(&U, 1, 0, 1);
    test_set_i64(&U, 1, 1, 0);
    assert(rmat_hnf_check_unimodular(&U));

    // shear: det 1
    test_set_i64(&U, 0, 0, 1);
    test_set_i64(&U, 0, 1, 1);
    test_set_i64(&U, 1, 0, 0);
    test_set_i64(&U, 1, 1, 1);
    assert(rmat_hnf_check_unimodular(&U));

    // det 2: not unimodular
    test_set_i64(&U, 0, 0, 2);
    test_set_i64(&U, 1, 1, 1);
    assert(!rmat_hnf_check_unimodular(&U));

    rmat_clear(&U);

    // non-square: rejected without a determinant
    rmat_t R = {0};
    rmat_init(&R, 2, 3);
    assert(!rmat_hnf_check_unimodular(&R));
    rmat_clear(&R);
    printf("[PASS] hnf_check_unimodular\n");
  }

  // 3. Transformation: the library's known example, A * U = H
  //    A = [[1,2],[3,4]], U = [[-4,-1],[3,1]] (det -1),
  //    A * U = [[2,1],[0,1]] = H
  rmat_t A = {0}, U = {0}, H = {0};
  rmat_init(&A, 2, 2);
  rmat_init(&U, 2, 2);
  rmat_init(&H, 2, 2);
  test_set_i64(&A, 0, 0, 1);
  test_set_i64(&A, 0, 1, 2);
  test_set_i64(&A, 1, 0, 3);
  test_set_i64(&A, 1, 1, 4);
  test_set_i64(&U, 0, 0, -4);
  test_set_i64(&U, 0, 1, -1);
  test_set_i64(&U, 1, 0, 3);
  test_set_i64(&U, 1, 1, 1);
  test_set_i64(&H, 0, 0, 2);
  test_set_i64(&H, 0, 1, 1);
  test_set_i64(&H, 1, 1, 1);

  // sanity: A * U really is H (independent of the checker)
  rmat_t P = {0};
  rmat_init(&P, 2, 2);
  rmat_mul(&P, &A, &U);
  assert(rmat_equal(&P, &H));
  rmat_clear(&P);

  assert(rmat_hnf_check_transformation(&A, &H, &U));

  // corrupt H: product mismatch
  rmat_t Hbad = {0};
  rmat_init(&Hbad, 2, 2);
  rmat_copy(&Hbad, &H);
  test_set_i64(&Hbad, 0, 1, 7);
  assert(!rmat_hnf_check_transformation(&A, &Hbad, &U));
  rmat_clear(&Hbad);

  // shape mismatch: U must be 2x2 for a 2x2 A
  rmat_t Ubad = {0};
  rmat_init(&Ubad, 1, 2);
  assert(!rmat_hnf_check_transformation(&A, &H, &Ubad));
  rmat_clear(&Ubad);
  printf("[PASS] hnf_check_transformation\n");

  // 4. Master check: the full triple passes
  assert(rmat_hnf_verify(&A, &H, &U));

  // structure violation is caught first (lower triangle entry)
  rmat_t Hs = {0};
  rmat_init(&Hs, 2, 2);
  rmat_copy(&Hs, &H);
  test_set_i64(&Hs, 1, 0, 5);
  assert(!rmat_hnf_verify(&A, &Hs, &U));
  rmat_clear(&Hs);

  // structure-valid H, but A * U != H: [[2,0],[0,1]] is in HNF shape
  rmat_t Ht = {0};
  rmat_init(&Ht, 2, 2);
  test_set_i64(&Ht, 0, 0, 2);
  test_set_i64(&Ht, 1, 1, 1);
  assert(rmat_hnf_check_structure(&Ht));
  assert(!rmat_hnf_verify(&A, &Ht, &U));
  rmat_clear(&Ht);

  // H = A * U' with U' = 2 * U (det -4): structure and transformation
  // hold, unimodularity must be the failing check
  rmat_t U2 = {0}, H2 = {0};
  rmat_init(&U2, 2, 2);
  rmat_init(&H2, 2, 2);
  test_set_i64(&U2, 0, 0, -8);
  test_set_i64(&U2, 0, 1, -2);
  test_set_i64(&U2, 1, 0, 6);
  test_set_i64(&U2, 1, 1, 2);
  rmat_mul(&H2, &A, &U2);
  assert(rmat_hnf_check_structure(&H2));
  assert(rmat_hnf_check_transformation(&A, &H2, &U2));
  assert(!rmat_hnf_check_unimodular(&U2));
  assert(!rmat_hnf_verify(&A, &H2, &U2));
  rmat_clear(&U2);
  rmat_clear(&H2);
  printf("[PASS] hnf_verify (master)\n");

  // 5. Working HNF test: rmat_hermite on random square matrices
  //    must produce matrices in HNF structure
  for (u64 n = 2; n <= 6; n++) {
    for (int rep = 0; rep < 3; rep++) {
      rmat_t R = {0}, W = {0};
      rmat_init(&R, n, n);
      rmat_init(&W, n, n);
      for (u64 i = 0; i < n; i++) {
        for (u64 j = 0; j < n; j++) {
          rz_gen_random(&val, 20);
          rmat_set(&R, &val, i, j);
        }
      }
      rmat_hermite(&W, &R);
      assert(rmat_hnf_check_structure(&W));
      rmat_clear(&R);
      rmat_clear(&W);
    }
  }
  printf("[PASS] hnf structure of rmat_hermite output\n");

  rmat_clear(&A);
  rmat_clear(&U);
  rmat_clear(&H);
  rz_clear(&val);
}
/* ----------------- LLL / MLLL tests ------------------------------ */

/* rq-vector helpers (plain arrays of rationals) */
static rq_t* lllv_new(u64 m)
{
  rq_t* v = calloc(m, sizeof(rq_t));
  assert(v != NULL);
  for (u64 i = 0; i < m; i++) {
    assert(rq_init(&v[i]) == RABIN_SUCCESS);
  }
  return v;
}

static void lllv_free(rq_t* v, u64 m)
{
  for (u64 i = 0; i < m; i++) {
    rq_clear(&v[i]);
  }
  free(v);
}

static void lllv_set_col(rq_t* v, const rmat_t* M, u64 c)
{
  for (u64 i = 0; i < M->rows; i++) {
    assert(rq_set(&v[i], RMAT_GET(M, i, c)) == RABIN_SUCCESS);
  }
}

static void lllv_dot(rq_t* r, const rq_t* a, const rq_t* b, u64 m)
{
  rq_t t;
  assert(rq_init(&t) == RABIN_SUCCESS);
  assert(rq_set_u64(r, 0) == RABIN_SUCCESS);
  for (u64 i = 0; i < m; i++) {
    assert(rq_mul(&t, &a[i], &b[i]) == RABIN_SUCCESS);
    assert(rq_add(r, r, &t) == RABIN_SUCCESS);
  }
  rq_clear(&t);
}

/* a - b >= 0 ? (rq values are normalized after every op) */
static bool lllv_cmp_ge(const rq_t* a, const rq_t* b)
{
  rq_t t;
  assert(rq_init(&t) == RABIN_SUCCESS);
  assert(rq_sub(&t, a, b) == RABIN_SUCCESS);
  bool r = !t.num.is_neg;
  rq_clear(&t);
  return r;
}

/*
 * Exact LLL check: the columns of B (m x p) form an LLL-reduced basis for
 * delta = dnum/dden. Runs exact rational Gram-Schmidt and verifies
 *   |mu[k][j]| <= 1/2                          (all j < k)
 *   ||b*_k||^2 >= (delta - mu[k][k-1]^2) * ||b*_{k-1}||^2
 */
static bool test_check_lll(const rmat_t* B, u64 p, i64 dnum, i64 dden)
{
  u64 m = B->rows;
  if (p <= 1) return true;

  rq_t** bs = malloc(p * sizeof(rq_t*));
  rq_t** bk = malloc(p * sizeof(rq_t*));
  assert(bs != NULL && bk != NULL);
  for (u64 k = 0; k < p; k++) {
    bs[k] = lllv_new(m);
    bk[k] = lllv_new(m);
  }
  rq_t* mu = lllv_new(p * p);  // flattened: mu[k * p + j]

  bool ok = true;
  rq_t num, two, oneq;
  assert(rq_init_multi(&num, &two, &oneq, NULL) == RABIN_SUCCESS);
  assert(rq_set_u64(&two, 2) == RABIN_SUCCESS);
  assert(rq_set_u64(&oneq, 1) == RABIN_SUCCESS);

  for (u64 k = 0; k < p && ok; k++) {
    lllv_set_col(bk[k], B, k);
    for (u64 j = k; j-- > 0;) {
      rq_t n2, d2;
      assert(rq_init(&n2) == RABIN_SUCCESS);
      assert(rq_init(&d2) == RABIN_SUCCESS);
      lllv_dot(&n2, bk[k], bs[j], m);
      lllv_dot(&d2, bs[j], bs[j], m);
      assert(rq_div(&mu[k * p + j], &n2, &d2) == RABIN_SUCCESS);
      for (u64 i = 0; i < m; i++) {
        assert(rq_mul(&n2, &mu[k * p + j], &bs[j][i]) == RABIN_SUCCESS);
        assert(rq_sub(&bk[k][i], &bk[k][i], &n2) == RABIN_SUCCESS);
      }
      rq_clear(&n2);
      rq_clear(&d2);
    }
    for (u64 i = 0; i < m; i++) {
      assert(rq_copy(&bs[k][i], &bk[k][i]) == RABIN_SUCCESS);
    }

    // size reduction: |2 * mu| <= 1
    for (u64 j = 0; j < k; j++) {
      assert(rq_mul(&num, &two, &mu[k * p + j]) == RABIN_SUCCESS);
      if (!lllv_cmp_ge(&oneq, &num)) ok = false;
    }
  }

  // Lovasz condition
  rq_t s1, s2, dq, dq2, rhs, dn, dd;
  rz_t dnumz, ddenz;
  assert(rq_init_multi(&s1, &s2, &dq, &dq2, &rhs, &dn, &dd, NULL) ==
         RABIN_SUCCESS);
  assert(rz_init_multi(&dnumz, &ddenz, NULL) == RABIN_SUCCESS);
  assert(rz_set_i64(&dnumz, dnum) == RABIN_SUCCESS);
  assert(rz_set_i64(&ddenz, dden) == RABIN_SUCCESS);
  assert(rq_set(&dn, &dnumz) == RABIN_SUCCESS);
  assert(rq_set(&dd, &ddenz) == RABIN_SUCCESS);
  assert(rq_div(&dq, &dn, &dd) == RABIN_SUCCESS);
  for (u64 k = 1; k < p && ok; k++) {
    lllv_dot(&s2, bs[k], bs[k], m);
    lllv_dot(&s1, bs[k - 1], bs[k - 1], m);
    assert(rq_sqr(&dq2, &mu[k * p + (k - 1)]) == RABIN_SUCCESS);
    assert(rq_sub(&dq, &dq, &dq2) == RABIN_SUCCESS);
    assert(rq_mul(&rhs, &dq, &s1) == RABIN_SUCCESS);
    if (!lllv_cmp_ge(&s2, &rhs)) ok = false;
  }

  rz_clear_multi(&dnumz, &ddenz, NULL);
  assert(rq_clear_multi(&s1, &s2, &dq, &dq2, &rhs, &dn, &dd, NULL) ==
         RABIN_SUCCESS);
  assert(rq_clear_multi(&num, &two, &oneq, NULL) == RABIN_SUCCESS);
  for (u64 i = 0; i < p * p; i++) {
    rq_clear(&mu[i]);
  }
  free(mu);
  for (u64 k = 0; k < p; k++) {
    lllv_free(bk[k], m);
    lllv_free(bs[k], m);
  }
  free(bk);
  free(bs);
  return ok;
}

/*
 * Is the integer vector b (length m) an integer combination of the columns
 * of W (m x p, full column rank)? Solves W x = b over Q via the normal
 * equations x = (W^T W)^{-1} W^T b (the unique rational solution) and
 * checks integrality with the adjugate.
 */
static bool test_in_col_Z(const rmat_t* W, const rmat_t* b, u64 p)
{
  u64 m = W->rows;
  assert(b->rows == m && b->cols == 1);
  if (p == 0) return false;

  rmat_t WT = {0}, G = {0}, N = {0}, Adj = {0}, Anum = {0};
  rpol_t cp = {0};
  rz_t detG, q, r;

  assert(rmat_init(&WT, p, m) == RABIN_SUCCESS);
  for (u64 i = 0; i < p; i++) {
    for (u64 j = 0; j < m; j++) {
      assert(rmat_set(&WT, RMAT_GET(W, j, i), i, j) == RABIN_SUCCESS);
    }
  }
  assert(rmat_init(&G, p, p) == RABIN_SUCCESS);
  assert(rmat_mul(&G, &WT, W) == RABIN_SUCCESS);
  assert(rmat_init(&N, p, 1) == RABIN_SUCCESS);
  for (u64 i = 0; i < p; i++) {
    rz_t acc;
    rz_init(&acc);
    rz_set_u64(&acc, 0);
    for (u64 j = 0; j < m; j++) {
      rz_t t;
      rz_init(&t);
      rz_mul(&t, RMAT_GET(W, j, i), RMAT_GET(b, j, 0));
      rz_add(&acc, &acc, &t);
      rz_clear(&t);
    }
    rmat_set(&N, &acc, i, 0);
    rz_clear(&acc);
  }
  assert(rmat_init(&Adj, p, p) == RABIN_SUCCESS);
  assert(rmat_init(&Anum, p, 1) == RABIN_SUCCESS);
  assert(rmat_charpoly_adj(&cp, &Adj, &G) == RABIN_SUCCESS);
  assert(rmat_mul(&Anum, &Adj, &N) == RABIN_SUCCESS);
  assert(rmat_det(&detG, &G) == RABIN_SUCCESS);
  rz_init_multi(&q, &r, NULL);

  bool ok = true;
  for (u64 i = 0; i < p && ok; i++) {
    assert(rz_divmod(&q, &r, RMAT_GET(&Anum, i, 0), &detG) == RABIN_SUCCESS);
    ok = rz_is_zero(&r);
  }

  rz_clear_multi(&q, &r, NULL);
  rz_clear(&detG);
  rpol_clear(&cp);
  rmat_clear(&Anum);
  rmat_clear(&Adj);
  rmat_clear(&N);
  rmat_clear(&G);
  rmat_clear(&WT);
  return ok;
}

/*
 * Full invariant check of an rmat_mlll result:
 *   - B == A0 * H
 *   - |det H| == 1
 *   - first n - p columns of B are zero
 *   - last p columns of B are LLL-reduced (delta = dnum/dden)
 *   - last p columns of B form a Z-basis of col_Z(A0)
 *     (rank via HNF, inclusion via the adjugate test)
 *   - if R_known != NULL (n x (n-p), a known Z-basis of the relation
 *     space): first n - p columns of H are a Z-basis of it (kernel test,
 *     inclusion, covolume equality)
 */
static bool test_mlll_check(const rmat_t* A0, const rmat_t* B, const rmat_t* H,
                            u64 p, const rmat_t* R_known, i64 dnum, i64 dden)
{
  u64 m = A0->rows;
  u64 n = A0->cols;
  bool ok = true;

  // 1. B == A0 * H
  rmat_t T = {0};
  assert(rmat_init(&T, m, n) == RABIN_SUCCESS);
  assert(rmat_mul(&T, A0, H) == RABIN_SUCCESS);
  ok = rmat_equal(&T, B);

  // 2. |det H| == 1
  if (n > 0) {
    rz_t det, one;
    rz_init_multi(&det, &one, NULL);
    rz_set_u64(&one, 1);
    assert(rmat_det(&det, H) == RABIN_SUCCESS);
    det.is_neg = false;
    ok = ok && rz_is_eq_i64(&det, 1);
    rz_clear_multi(&det, &one, NULL);
  }

  // 3. first n - p columns of B are zero; 3b. relation columns of H are in
  // the kernel (T = A0 * H already computed)
  for (u64 x = 0; x < n - p; x++) {
    for (u64 y = 0; y < m; y++) {
      ok = ok && rz_is_zero(RMAT_GET(B, y, x));
      ok = ok && rz_is_zero(RMAT_GET(&T, y, x));
    }
  }
  rmat_clear(&T);

  // 4. LLL-reduced basis block
  rmat_t Bp = {0};
  if (p > 0) {
    assert(rmat_init(&Bp, m, p) == RABIN_SUCCESS);
    for (u64 x = 0; x < p; x++) {
      for (u64 y = 0; y < m; y++) {
        assert(rmat_set(&Bp, RMAT_GET(B, y, n - p + x), y, x) == RABIN_SUCCESS);
      }
    }
    ok = ok && test_check_lll(&Bp, p, dnum, dden);
  }

  // 5. basis block is a Z-basis of col_Z(A0): W = HNF(A0) is one; check
  // W has p columns and every column of Bp lies in col_Z(W)
  rmat_t W = {0};
  assert(rmat_hermite(&W, A0) == RABIN_SUCCESS);
  ok = ok && (W.cols == p);
  if (p > 0) {
    for (u64 x = 0; x < p && ok; x++) {
      rmat_t col = {0};
      assert(rmat_init(&col, m, 1) == RABIN_SUCCESS);
      for (u64 y = 0; y < m; y++) {
        assert(rmat_set(&col, RMAT_GET(&Bp, y, x), y, 0) == RABIN_SUCCESS);
      }
      ok = test_in_col_Z(&W, &col, p);
      rmat_clear(&col);
    }
  }

  // 6. relation block is a Z-basis of the relation space (vs R_known)
  if (R_known != NULL) {
    u64 q = n - p;
    if (q > 0) {
      rmat_t Hrel = {0};
      assert(rmat_init(&Hrel, n, q) == RABIN_SUCCESS);
      for (u64 x = 0; x < q; x++) {
        for (u64 y = 0; y < n; y++) {
          assert(rmat_set(&Hrel, RMAT_GET(H, y, x), y, x) == RABIN_SUCCESS);
        }
      }
      for (u64 x = 0; x < q && ok; x++) {
        rmat_t col = {0};
        assert(rmat_init(&col, n, 1) == RABIN_SUCCESS);
        for (u64 y = 0; y < n; y++) {
          assert(rmat_set(&col, RMAT_GET(&Hrel, y, x), y, 0) == RABIN_SUCCESS);
        }
        ok = test_in_col_Z(R_known, &col, q);
        rmat_clear(&col);
      }
      // covolume equality: det(Gram(Hrel)) == det(Gram(R_known))
      rmat_t HrelT = {0}, RkT = {0}, GH = {0}, GR = {0};
      rz_t dH, dR;
      assert(rmat_init(&HrelT, q, n) == RABIN_SUCCESS);
      assert(rmat_init(&RkT, q, n) == RABIN_SUCCESS);
      assert(rmat_init(&GH, q, q) == RABIN_SUCCESS);
      assert(rmat_init(&GR, q, q) == RABIN_SUCCESS);
      assert(rz_init(&dH) == RABIN_SUCCESS);
      assert(rz_init(&dR) == RABIN_SUCCESS);
      for (u64 i = 0; i < q; i++) {
        for (u64 j = 0; j < n; j++) {
          assert(rmat_set(&HrelT, RMAT_GET(&Hrel, j, i), i, j) ==
                 RABIN_SUCCESS);
          assert(rmat_set(&RkT, RMAT_GET(R_known, j, i), i, j) ==
                 RABIN_SUCCESS);
        }
      }
      assert(rmat_mul(&GH, &HrelT, &Hrel) == RABIN_SUCCESS);
      assert(rmat_mul(&GR, &RkT, R_known) == RABIN_SUCCESS);
      assert(rmat_det(&dH, &GH) == RABIN_SUCCESS);
      assert(rmat_det(&dR, &GR) == RABIN_SUCCESS);
      dH.is_neg = false;
      dR.is_neg = false;
      ok = ok && rz_is_eq_i64(&dH, 0) == false && (rz_cmp(&dH, &dR) == 0);
      rmat_clear(&HrelT);
      rmat_clear(&RkT);
      rmat_clear(&GH);
      rmat_clear(&GR);
      rz_clear_multi(&dH, &dR, NULL);
      rmat_clear(&Hrel);
    }
  }

  rmat_clear(&W);
  rmat_clear(&Bp);
  return ok;
}

void test_lll(void)
{
  // 1. 2x2 full rank: A = [[1, 2], [3, 4]]
  {
    rmat_t A = {0}, B = {0}, H = {0};
    rmat_init(&A, 2, 2);
    test_set_i64(&A, 0, 0, 1);
    test_set_i64(&A, 0, 1, 2);
    test_set_i64(&A, 1, 0, 3);
    test_set_i64(&A, 1, 1, 4);
    rmat_init(&B, 2, 2);
    assert(rmat_copy(&B, &A) == RABIN_SUCCESS);
    assert(rmat_lll(&B, 2, 0.75, &H) == RABIN_SUCCESS);
    // invariants: B = A * H, |det H| = 1, LLL conditions
    rmat_t T = {0};
    rmat_init(&T, 2, 2);
    rmat_mul(&T, &A, &H);
    assert(rmat_equal(&T, &B));
    rz_t det;
    rz_init(&det);
    rmat_det(&det, &H);
    det.is_neg = false;
    assert(rz_is_eq_i64(&det, 1));
    rz_clear(&det);
    assert(test_check_lll(&B, 2, 3, 4));
    rmat_clear(&T);
    rmat_clear(&A);
    rmat_clear(&B);
    rmat_clear(&H);
    printf("[PASS] lll 2x2 known example\n");
  }

  // 2. rank deficient 1 x 3 [3, 5, 7]: LLL requires full rank -> rejected
  {
    rmat_t B = {0}, H = {0};
    rmat_init(&B, 1, 3);
    test_set_i64(&B, 0, 0, 3);
    test_set_i64(&B, 0, 1, 5);
    test_set_i64(&B, 0, 2, 7);
    assert(rmat_lll(&B, 3, 0.75, &H) == RABIN_ERR_INVALID_ARG);
    rmat_clear(&B);
    printf("[PASS] lll rank deficient rejected\n");
  }

  // 3. H passed as an empty (zero-initialized) matrix must work
  {
    rmat_t B = {0}, H = {0};
    rmat_init(&B, 2, 2);
    test_set_i64(&B, 0, 0, 4);
    test_set_i64(&B, 0, 1, 1);
    test_set_i64(&B, 1, 0, 1);
    test_set_i64(&B, 1, 1, 3);
    assert(rmat_lll(&B, 2, 0.75, &H) == RABIN_SUCCESS);
    assert(H.rows == 2 && H.cols == 2);
    assert(test_check_lll(&B, 2, 3, 4));
    rmat_clear(&B);
    rmat_clear(&H);
    printf("[PASS] lll empty H output\n");
  }

  // 4. dependent columns -> RABIN_ERR_INVALID_ARG
  {
    rmat_t B = {0}, H = {0};
    rmat_init(&B, 2, 2);
    test_set_i64(&B, 0, 0, 1);
    test_set_i64(&B, 0, 1, 2);
    test_set_i64(&B, 1, 0, 2);
    test_set_i64(&B, 1, 1, 4);
    assert(rmat_lll(&B, 2, 0.75, &H) == RABIN_ERR_INVALID_ARG);
    rmat_clear(&B);
    printf("[PASS] lll dependent columns rejected\n");
  }

  // 5. argument validation
  {
    rmat_t B = {0}, H = {0};
    rmat_init(&B, 2, 2);
    assert(rmat_lll(&B, 2, 0.0, &H) == RABIN_ERR_INVALID_ARG);
    assert(rmat_lll(&B, 2, 1.5, &H) == RABIN_ERR_INVALID_ARG);
    assert(rmat_lll(&B, 3, 0.75, &H) == RABIN_ERR_INVALID_ARG);
    assert(rmat_lll(NULL, 2, 0.75, &H) == RABIN_ERR_NULL_PTR);
    assert(rmat_lll(&B, 2, 0.75, NULL) == RABIN_ERR_NULL_PTR);
    rmat_clear(&B);
    printf("[PASS] lll argument validation\n");
  }

  // 6. random 4 x 3 (tall, full rank), delta = 0.5 (double branch)
  for (u64 trial = 0; trial < 5; trial++) {
    rmat_t A = {0}, B = {0}, H = {0}, T = {0};
    rz_t val;
    rz_init(&val);
    rmat_init(&A, 4, 3);
    rmat_init(&B, 4, 3);
    for (u64 i = 0; i < 4; i++) {
      for (u64 j = 0; j < 3; j++) {
        rz_gen_random(&val, 16);
        rmat_set(&A, &val, i, j);
      }
    }
    assert(rmat_copy(&B, &A) == RABIN_SUCCESS);
    if (rmat_lll(&B, 3, 0.5, &H) != RABIN_SUCCESS) continue;  // unlucky rank
    rmat_init(&T, 4, 3);
    rmat_mul(&T, &A, &H);
    assert(rmat_equal(&T, &B));
    assert(test_check_lll(&B, 3, 1, 2));
    rmat_clear(&T);
    rmat_clear(&B);
    rmat_clear(&A);
    rmat_clear(&H);
    rz_clear(&val);
  }
  printf("[PASS] lll random 4x3\n");
}

void test_mlll(void)
{
  // 1. [[2, 4]]: rank 1, lattice 2Z
  {
    rmat_t A0 = {0}, B = {0}, H = {0};
    rmat_init(&A0, 1, 2);
    test_set_i64(&A0, 0, 0, 2);
    test_set_i64(&A0, 0, 1, 4);
    rmat_init(&B, 1, 2);
    assert(rmat_copy(&B, &A0) == RABIN_SUCCESS);
    u64 rank = 99;
    assert(rmat_mlll(&B, 2, 0.75, &H, &rank) == RABIN_SUCCESS);
    assert(rank == 1);
    // W from the HNF is [[2]], so the output is exactly [0 | 2]
    assert(rz_is_zero(RMAT_GET(&B, 0, 0)));
    assert(rz_is_eq_i64(RMAT_GET(&B, 0, 1), 2));
    assert(test_mlll_check(&A0, &B, &H, 1, NULL, 3, 4));
    rmat_clear(&A0);
    rmat_clear(&B);
    rmat_clear(&H);
    printf("[PASS] mlll 1x2 known example\n");
  }

  // 2. 1 x 3 gcd case: [3, 5, 7] -> rank 1, reduced column +-gcd = +-1
  {
    rmat_t A0 = {0}, B = {0}, H = {0};
    rmat_init(&A0, 1, 3);
    test_set_i64(&A0, 0, 0, 3);
    test_set_i64(&A0, 0, 1, 5);
    test_set_i64(&A0, 0, 2, 7);
    rmat_init(&B, 1, 3);
    assert(rmat_copy(&B, &A0) == RABIN_SUCCESS);
    u64 rank = 99;
    assert(rmat_mlll(&B, 3, 0.75, &H, &rank) == RABIN_SUCCESS);
    assert(rank == 1);
    // the surviving column is +-gcd(3,5,7) = +-1
    assert(rz_is_eq_i64(RMAT_GET(&B, 0, 2), 1) ||
           rz_is_eq_i64(RMAT_GET(&B, 0, 2), -1));
    assert(test_mlll_check(&A0, &B, &H, 1, NULL, 3, 4));
    rmat_clear(&A0);
    rmat_clear(&B);
    rmat_clear(&H);
    printf("[PASS] mlll 1x3 gcd case\n");
  }

  // 3. dependent 2x3: [[2, 0, 1], [0, 2, 1]], relations = Z * (-1,-1,2)
  {
    rmat_t A0 = {0}, B = {0}, H = {0}, Rk = {0};
    rmat_init(&A0, 2, 3);
    test_set_i64(&A0, 0, 0, 2);
    test_set_i64(&A0, 1, 0, 0);
    test_set_i64(&A0, 0, 1, 0);
    test_set_i64(&A0, 1, 1, 2);
    test_set_i64(&A0, 0, 2, 1);
    test_set_i64(&A0, 1, 2, 1);
    rmat_init(&B, 2, 3);
    assert(rmat_copy(&B, &A0) == RABIN_SUCCESS);
    u64 rank = 99;
    assert(rmat_mlll(&B, 3, 0.75, &H, &rank) == RABIN_SUCCESS);
    assert(rank == 2);
    rmat_init(&Rk, 3, 1);
    test_set_i64(&Rk, 0, 0, -1);
    test_set_i64(&Rk, 1, 0, -1);
    test_set_i64(&Rk, 2, 0, 2);
    assert(test_mlll_check(&A0, &B, &H, 2, &Rk, 3, 4));
    rmat_clear(&Rk);
    rmat_clear(&A0);
    rmat_clear(&B);
    rmat_clear(&H);
    printf("[PASS] mlll 2x3 dependent with known relations\n");
  }

  // 3. all-zero matrix: rank 0, H = I
  {
    rmat_t B = {0}, H = {0}, I3 = {0};
    rmat_init(&B, 2, 3);
    u64 rank = 99;
    assert(rmat_mlll(&B, 3, 0.75, &H, &rank) == RABIN_SUCCESS);
    assert(rank == 0);
    rmat_init(&I3, 3, 3);
    rmat_id(&I3, 3);
    assert(rmat_equal(&H, &I3));
    assert(rmat_is_zero(&B));
    rmat_clear(&I3);
    rmat_clear(&B);
    rmat_clear(&H);
    printf("[PASS] mlll all-zero matrix\n");
  }

  // 4. single vectors: [[5]] and [[-5]]
  {
    rmat_t B = {0}, H = {0};
    rmat_init(&B, 1, 1);
    test_set_i64(&B, 0, 0, 5);
    u64 rank = 99;
    assert(rmat_mlll(&B, 1, 0.75, &H, &rank) == RABIN_SUCCESS);
    assert(rank == 1);
    assert(rz_is_eq_i64(RMAT_GET(&B, 0, 0), 5));
    assert(rz_is_eq_i64(RMAT_GET(&H, 0, 0), 1));

    test_set_i64(&B, 0, 0, -5);
    assert(rmat_mlll(&B, 1, 0.75, &H, &rank) == RABIN_SUCCESS);
    assert(rank == 1);
    assert(rz_is_eq_i64(RMAT_GET(&B, 0, 0), 5));
    assert(rz_is_eq_i64(RMAT_GET(&H, 0, 0), -1));
    rmat_clear(&B);
    rmat_clear(&H);
    printf("[PASS] mlll single vector (+-5)\n");
  }

  // 5. zero coordinates (0 x 3)
  {
    rmat_t B = {0}, H = {0}, I3 = {0};
    rmat_init(&B, 0, 3);
    u64 rank = 99;
    assert(rmat_mlll(&B, 3, 0.75, &H, &rank) == RABIN_SUCCESS);
    assert(rank == 0);
    rmat_init(&I3, 3, 3);
    rmat_id(&I3, 3);
    assert(rmat_equal(&H, &I3));
    rmat_clear(&I3);
    rmat_clear(&B);
    rmat_clear(&H);
    printf("[PASS] mlll zero rows\n");
  }

  // 6. argument validation
  {
    rmat_t B = {0}, H = {0};
    rmat_init(&B, 3, 2);  // tall: m > n -> INVALID_ARG
    u64 rank = 99;
    assert(rmat_mlll(&B, 2, 0.75, &H, &rank) == RABIN_ERR_INVALID_ARG);
    rmat_init(&B, 2, 2);
    assert(rmat_mlll(&B, 2, 0.0, &H, &rank) == RABIN_ERR_INVALID_ARG);
    assert(rmat_mlll(&B, 2, 1.5, &H, &rank) == RABIN_ERR_INVALID_ARG);
    assert(rmat_mlll(&B, 3, 0.75, &H, &rank) == RABIN_ERR_INVALID_ARG);
    assert(rmat_mlll(NULL, 2, 0.75, &H, &rank) == RABIN_ERR_NULL_PTR);
    assert(rmat_mlll(&B, 2, 0.75, NULL, &rank) == RABIN_ERR_NULL_PTR);
    assert(rmat_mlll(&B, 2, 0.75, &H, NULL) == RABIN_ERR_NULL_PTR);
    rmat_clear(&B);
    printf("[PASS] mlll argument validation\n");
  }

  // 7. structured random tests with known rank and known relations:
  //    A = V * [I | C_r] with V m x p lower triangular (a known Z-basis of
  //    the column lattice); the known relations are e_j - sum_h c_h e_h.
  {
    struct {
      u64 m, n, p;
    } shapes[] = {{2, 2, 2}, {2, 3, 1}, {2, 3, 2}, {2, 4, 1}, {2, 4, 2},
                  {3, 4, 3}, {3, 3, 3}, {3, 4, 2}, {3, 5, 2}, {3, 5, 3},
                  {3, 5, 1}, {4, 6, 3}, {4, 6, 4}};
    u64 nshapes = sizeof(shapes) / sizeof(shapes[0]);
    rz_t val;
    rz_init(&val);
    for (u64 s = 0; s < nshapes; s++) {
      u64 m = shapes[s].m, n = shapes[s].n, p = shapes[s].p;
      for (int delta_i = 0; delta_i < 2; delta_i++) {
        double delta = (delta_i == 0) ? 0.75 : 0.5;
        i64 dnum = (delta_i == 0) ? 3 : 1;
        i64 dden = (delta_i == 0) ? 4 : 2;

        rmat_t V = {0};   // m x p known basis
        rmat_t C = {0};   // p x (n-p) random coefficients
        rmat_t A0 = {0};  // m x n
        rmat_t Rk = {0};  // n x (n-p) known relations
        rmat_t B = {0};
        rmat_t H = {0};
        rmat_init(&V, m, p);
        rmat_init(&C, p, n - p);
        rmat_init(&A0, m, n);
        rmat_init(&Rk, n, n - p);
        for (u64 i = 0; i < m; i++) {
          for (u64 j = 0; j < p; j++) {
            if (i >= j) {
              rz_gen_random(&val, 12);
              // keep the diagonal nonzero and small
              if (i == j) {
                rz_t d;
                rz_init(&d);
                rz_gen_random(&d, 3);
                if (rz_is_zero(&d)) rz_set_u64(&d, 1);
                d.is_neg = (rand() % 2 == 0);
                rmat_set(&V, &d, i, j);
                rz_clear(&d);
              } else {
                rmat_set(&V, &val, i, j);
              }
            }
          }
        }
        for (u64 i = 0; i < p; i++) {
          for (u64 j = 0; j < n - p; j++) {
            rz_set_i64(&val, (i64)(rand() % 7) - 3);
            rmat_set(&C, &val, i, j);
          }
        }
        // A0[:, j] = V[:, j] for j < p;  A0[:, j] = V * C[:, j-p] else
        for (u64 j = 0; j < n; j++) {
          for (u64 i = 0; i < m; i++) {
            rz_t acc;
            rz_init(&acc);
            rz_set_u64(&acc, 0);
            for (u64 h = 0; h < p; h++) {
              rz_t cv, t;
              rz_init_multi(&cv, &t, NULL);
              if (j < p) {
                if (h == j) rz_set_u64(&cv, 1);
              } else {
                rz_copy(&cv, RMAT_GET(&C, h, j - p));
              }
              rz_mul(&t, &cv, RMAT_GET(&V, i, h));
              rz_add(&acc, &acc, &t);
              rz_clear_multi(&cv, &t, NULL);
            }
            rmat_set(&A0, &acc, i, j);
            rz_clear(&acc);
          }
        }
        // Rk[:, j] = e_{p+j} - sum_h C[h][j] e_h
        for (u64 j = 0; j < n - p; j++) {
          for (u64 i = 0; i < n; i++) {
            rz_t v2;
            rz_init(&v2);
            rz_set_u64(&v2, 0);
            if (i == p + j) {
              rz_set_u64(&v2, 1);
            } else if (i < p) {
              rz_neg(&v2, RMAT_GET(&C, i, j));
            }
            rmat_set(&Rk, &v2, i, j);
            rz_clear(&v2);
          }
        }
        rmat_init(&B, m, n);
        assert(rmat_copy(&B, &A0) == RABIN_SUCCESS);
        u64 rank = 99;
        assert(rmat_mlll(&B, n, delta, &H, &rank) == RABIN_SUCCESS);
        assert(rank == p);
        assert(test_mlll_check(&A0, &B, &H, p, (n - p > 0) ? &Rk : NULL, dnum,
                               dden));
        rmat_clear(&B);
        rmat_clear(&H);
        rmat_clear(&Rk);
        rmat_clear(&A0);
        rmat_clear(&C);
        rmat_clear(&V);
      }
    }
    rz_clear(&val);
    printf("[PASS] mlll structured random (%llu shapes x 2 deltas)\n",
           (unsigned long long)nshapes);
  }
}

void test_pascal_det(u64 size)
{
  rmat_t P = {0};
  rz_t det, val, a, b;  // Move a and b here
  rz_init(&det);
  rz_init(&val);
  rz_init(&a);
  rz_init(&b);

  rmat_init(&P, size, size);

  rns_ctx_t ctx = {0};
  // Ensure RNS_PRIMES has at least 50 elements!
  if (rns_ctx_init(&ctx, RNS_PRIMES, 10) != RABIN_SUCCESS) return;

  for (u64 i = 0; i < size; i++) {
    for (u64 j = 0; j < size; j++) {
      if (i == 0 || j == 0) {
        rz_set_u64(&val, 1);
      } else {
        // Reuse a and b instead of re-allocating
        rmat_get(&a, &P, i - 1, j);
        rmat_get(&b, &P, i, j - 1);
        rz_add(&val, &a, &b);
      }
      rmat_set(&P, &val, i, j);
    }
  }

  // Standard Det (for comparison)
  rmat_det(&det, &P);
  printf("Pascal %llu x %llu | Standard Det: ", size, size);
  rz_println(&det);

  // RNS Det
  rmat_det_rns(&det, &P, &ctx);
  printf("Pascal %llu x %llu | RNS Det:      ", size, size);
  rz_println(&det);

  // Verification
  if (rz_is_eq_i64(&det, 1)) {
    printf("RESULT: PASS\n");
  } else {
    printf("RESULT: FAIL (Expected 1)\n");
  }

  // CLEANUP
  rmat_clear(&P);
  rz_clear(&det);
  rz_clear(&val);
  rz_clear(&a);
  rz_clear(&b);

  rns_ctx_clear(&ctx);
}

// GCD of all m x m minors of A (m = A->rows <= n = A->cols): the
// determinant Delta of the Z-module generated by the columns of A.
static void test_minors_gcd(rz_t* g, const rmat_t* A, u64 m, u64 n, u64 start,
                            u64 left, u64* cols)
{
  if (left == 0) {
    rmat_t S = {0};
    rz_t det;
    rmat_init(&S, m, m);
    rz_init(&det);
    for (u64 i = 0; i < m; i++) {
      for (u64 j = 0; j < m; j++) {
        rmat_set(&S, &A->data[i * n + cols[j]], i, j);
      }
    }
    rmat_det(&det, &S);
    det.is_neg = false;
    rz_gcd(g, g, &det);
    rz_clear(&det);
    rmat_clear(&S);
    return;
  }

  for (u64 c = start; c <= n - left; c++) {
    cols[m - left] = c;
    test_minors_gcd(g, A, m, n, c + 1, left - 1, cols);
  }
}

static void test_module_determinant(rz_t* delta, const rmat_t* A)
{
  u64 m = A->rows;
  u64 n = A->cols;
  rz_t g;
  u64* cols = malloc(m * sizeof(u64));
  rz_init(&g);
  rz_set_u64(&g, 0);
  test_minors_gcd(&g, A, m, n, 0, m, cols);
  rz_copy(delta, &g);
  free(cols);
  rz_clear(&g);
}

void test_hnf_mod_d_multilimb(void)
{
  rz_t D, val, bit, delta, mult;
  rz_init_multi(&D, &val, &bit, &delta, &mult, NULL);

  // 1. Zero-pivot multi-limb input - regression for the out-of-bounds top
  //    limb read in rz_gcd_extended_lehmer when a size-0 rz_t is the
  //    first gcd argument: A = [[p, q], [r, 0]] (all entries 2-limb) has
  //    HNF [[q, p mod q], [0, r]], here p < q so p mod q = p
  {
    rmat_t A = {0}, W = {0}, H = {0};
    rz_t p, q, r, zero;
    rz_init_multi(&p, &q, &r, &zero, NULL);
    rmat_init(&A, 2, 2);
    rmat_init(&W, 2, 2);
    rmat_init(&H, 2, 2);

    // p = 2^70, q = 2^70 + 3, r = 2^65 + 7
    rz_set_u64(&p, 1);
    rz_lshift(&p, &p, 70);
    rz_set_u64(&q, 1);
    rz_lshift(&q, &q, 70);
    rz_add_u64(&q, &q, 3);
    rz_set_u64(&r, 1);
    rz_lshift(&r, &r, 65);
    rz_add_u64(&r, &r, 7);
    rz_set_u64(&zero, 0);

    rmat_set(&A, &p, 0, 0);
    rmat_set(&A, &q, 0, 1);
    rmat_set(&A, &r, 1, 0);
    rmat_set(&A, &zero, 1, 1);

    // D = 3 * q * r = 3 * |det A|, a multiple of the module determinant
    rz_mul(&D, &q, &r);
    rz_set_u64(&mult, 3);
    rz_mul(&D, &mult, &D);

    rmat_hermite_mod_d(&W, &A, &D);

    rmat_set(&H, &q, 0, 0);
    rmat_set(&H, &p, 0, 1);
    rmat_set(&H, &zero, 1, 0);
    rmat_set(&H, &r, 1, 1);
    assert(rmat_equal(&W, &H));
    assert(rmat_hnf_check_structure(&W));

    // cross-check against the exact implementation
    rmat_t V = {0};
    rmat_init(&V, 2, 2);
    rmat_hermite(&V, &A);
    assert(rmat_equal(&W, &V));

    rmat_clear(&V);
    rmat_clear(&A);
    rmat_clear(&W);
    rmat_clear(&H);
    rz_clear_multi(&p, &q, &r, &zero, NULL);
    printf("[PASS] hermite_mod_d 2x2 multi-limb zero pivot\n");
  }

  // 2. rz_gcd_extended_lehmer with a zero argument (multi-limb partner):
  //    gcd(0, b) -> d = |b|, u = 0, v = sign(b); gcd(a, 0) -> d = |a|,
  //    u = sign(a), v = 0
  {
    rz_t zero, a, b, bneg, u, v, d, n1;
    rz_init_multi(&zero, &a, &b, &bneg, &u, &v, &d, &n1, NULL);
    rz_set_i64(&n1, -1);
    rz_set_u64(&zero, 0);
    rz_set_u64(&a, 1);
    rz_lshift(&a, &a, 70);
    rz_add_u64(&a, &a, 5);  // 2^70 + 5
    rz_set_u64(&b, 1);
    rz_lshift(&b, &b, 66);  // 2^66

    rz_gcd_extended_lehmer(&u, &v, &d, &zero, &b);
    assert(rz_is_eq_i64(&u, 0));
    assert(rz_is_eq_i64(&v, 1));
    assert(rz_cmp(&d, &b) == 0);

    rz_gcd_extended_lehmer(&u, &v, &d, &a, &zero);
    assert(rz_is_eq_i64(&u, 1));
    assert(rz_is_eq_i64(&v, 0));
    assert(rz_cmp(&d, &a) == 0);

    rz_copy(&bneg, &b);
    bneg.is_neg = true;
    rz_gcd_extended_lehmer(&u, &v, &d, &zero, &bneg);
    assert(rz_is_eq_i64(&u, 0));
    assert(rz_cmp(&v, &n1) == 0);
    assert(rz_cmp(&d, &b) == 0);
    printf("[PASS] gcd_extended_lehmer zero arguments (multi-limb)\n");
    rz_clear_multi(&zero, &a, &b, &bneg, &u, &v, &d, &n1, NULL);
  }

  // 3. Random multi-limb cross-checks against the exact implementations
  //    (the HNF is unique). All entries are > 64 bits so the multi-limb
  //    Lehmer path is actually taken; in a third of the trials the
  //    bottom-right entry (the zero-pivot shape of test 1) is forced to 0
  {
    u64 shapes[][2] = {{2, 2}, {2, 3}, {3, 3}, {3, 4}, {4, 4}, {4, 6}, {5, 7}};
    u64 n_shapes = sizeof(shapes) / sizeof(shapes[0]);
    u64 bits_list[] = {66, 90, 130, 200};
    u64 n_bits = sizeof(bits_list) / sizeof(bits_list[0]);
    u64 mults[] = {1, 3, 7};
    u64 n_mults = sizeof(mults) / sizeof(mults[0]);

    for (u64 s = 0; s < n_shapes; s++) {
      u64 m = shapes[s][0];
      u64 n = shapes[s][1];
      bool square = (m == n);
      for (u64 bl = 0; bl < n_bits; bl++) {
        u64 bits = bits_list[bl];
        for (u64 trial = 0; trial < n_mults; trial++) {
          rmat_t A = {0}, H = {0}, W1 = {0}, W2 = {0};
          rmat_init(&A, m, n);
          rmat_init(&H, m, m);
          rmat_init(&W1, m, m);
          rmat_init(&W2, m, m);

          for (u64 i = 0; i < m; i++) {
            for (u64 j = 0; j < n; j++) {
              rz_gen_random(&val, bits);
              rz_gen_random(&bit, 1);
              if (!rz_is_zero(&bit)) {
                rz_neg(&val, &val);
              }
              rmat_set(&A, &val, i, j);
            }
          }

          // force the zero pivot (rectangular matrices only: on a square
          // matrix it could cost rank)
          bool force_zero = !square && trial % 2 == 1;
          if (force_zero) {
            rz_set_u64(&val, 0);
            rmat_set(&A, &val, m - 1, n - 1);
          }

          // Delta = GCD of all m x m minors; if it is 0 the matrix is rank
          // deficient, so draw a new one
          test_module_determinant(&delta, &A);
          if (rz_is_zero(&delta)) {
            for (u64 i = 0; i < m; i++) {
              for (u64 j = 0; j < n; j++) {
                rz_gen_random(&val, bits);
                rz_gen_random(&bit, 1);
                if (!rz_is_zero(&bit)) {
                  rz_neg(&val, &val);
                }
                rmat_set(&A, &val, i, j);
              }
            }
            if (force_zero) {
              rz_set_u64(&val, 0);
              rmat_set(&A, &val, m - 1, n - 1);
            }
            test_module_determinant(&delta, &A);
          }
          assert(!rz_is_zero(&delta));

          rmat_hermite(&H, &A);  // 2.4.4 oracle
          assert(H.rows == m && H.cols == m);
          assert(rmat_hnf_check_structure(&H));

          rz_set_u64(&mult, mults[trial]);
          rz_mul(&D, &mult, &delta);

          rmat_hermite_mod_d(&W1, &A, &D);  // 2.4.8 fast path
          assert(rmat_equal(&W1, &H));
          assert(rmat_hnf_check_structure(&W1));

          rmat_hermite_gcd(&W2, &A);  // 2.4.5 exact
          assert(rmat_equal(&W2, &H));
          assert(rmat_hnf_check_structure(&W2));

          rmat_clear(&A);
          rmat_clear(&H);
          rmat_clear(&W1);
          rmat_clear(&W2);
        }
      }
      printf("[PASS] hermite_mod_d random %llux%llu multi-limb\n", m, n);
    }

    // 4. Aggressive case: 10x10 with 300-bit entries (5 limbs)
    {
      u64 m = 10, n = 10;
      u64 bits = 300;
      u64 mults[] = {1, 7};
      for (u64 trial = 0; trial < 2; trial++) {
        rmat_t A = {0}, H = {0}, W1 = {0}, W2 = {0};
        rmat_init(&A, m, n);
        rmat_init(&H, m, m);
        rmat_init(&W1, m, m);
        rmat_init(&W2, m, m);

        for (u64 i = 0; i < m; i++) {
          for (u64 j = 0; j < n; j++) {
            rz_gen_random(&val, bits);
            rz_gen_random(&bit, 1);
            if (!rz_is_zero(&bit)) {
              rz_neg(&val, &val);
            }
            rmat_set(&A, &val, i, j);
          }
        }

        // square: Delta = |det A|
        rmat_det(&delta, &A);
        if (rz_is_zero(&delta)) {
          for (u64 i = 0; i < m; i++) {
            rz_gen_random(&val, bits);
            rmat_set(&A, &val, i, i);
          }
          rmat_det(&delta, &A);
        }
        assert(!rz_is_zero(&delta));
        delta.is_neg = false;

        rmat_hermite(&H, &A);
        assert(rmat_hnf_check_structure(&H));

        rz_set_u64(&mult, mults[trial]);
        rz_mul(&D, &mult, &delta);

        rmat_hermite_mod_d(&W1, &A, &D);
        assert(rmat_equal(&W1, &H));
        assert(rmat_hnf_check_structure(&W1));

        rmat_hermite_gcd(&W2, &A);
        assert(rmat_equal(&W2, &H));
        assert(rmat_hnf_check_structure(&W2));

        rmat_clear(&A);
        rmat_clear(&H);
        rmat_clear(&W1);
        rmat_clear(&W2);
      }
      printf("[PASS] hermite_mod_d 10x10 300-bit\n");
    }
  }

  rz_clear_multi(&D, &val, &bit, &delta, &mult, NULL);
}

static void test_rmatq_cholesky_case(rmatq_t* C, const rmatq_t* A)
{
  rq_t s, t, a, b;
  rq_init_multi(&s, &t, &a, &b, NULL);

  for (u64 i = 0; i < A->M.rows; i++) {
    for (u64 j = 0; j < A->M.cols; j++) {
      rq_set_u64(&s, 0);
      for (u64 k = 0; k < A->M.rows; k++) {
        rmatq_get(&a, C, k, i);
        if (rq_is_zero(&a)) continue;
        rmatq_get(&b, C, k, j);
        rq_mul(&t, &a, &b);
        rq_add(&s, &s, &t);
      }
      rmatq_get(&a, A, i, j);
      int eq = (rq_is_zero(&s) && rq_is_zero(&a)) ||
               (!rq_is_zero(&s) && !rq_is_zero(&a) && rq_cmp(&s, &a) == 0);
      assert(eq);
    }
  }

  rq_clear_multi(&s, &t, &a, &b, NULL);
}

static void test_rmatq_cholesky(void)
{
  rmatq_t A = {0}, C = {0};
  rmatq_init(&A, 3, 3);
  rmatq_init(&C, 3, 3);

  // A = [[4,2,0],[2,5,2],[0,2,5]] = C^t C, C = [[2,1,0],[0,2,1],[0,0,2]].
  test_set_i64(&A.M, 0, 0, 4);
  test_set_i64(&A.M, 0, 1, 2);
  test_set_i64(&A.M, 1, 0, 2);
  test_set_i64(&A.M, 1, 1, 5);
  test_set_i64(&A.M, 1, 2, 2);
  test_set_i64(&A.M, 2, 1, 2);
  test_set_i64(&A.M, 2, 2, 5);

  assert(rmatq_cholesky(&C, &A) == RABIN_SUCCESS);

  rq_t got, want;
  rq_init_multi(&got, &want, NULL);
  const i64 expected[3][3] = {{2, 1, 0}, {0, 2, 1}, {0, 0, 2}};
  for (u64 i = 0; i < 3; i++) {
    for (u64 j = 0; j < 3; j++) {
      rmatq_get(&got, &C, i, j);
      rq_set_u64(&want, (u64)expected[i][j]);
      assert(rq_cmp(&got, &want) == 0);
    }
  }
  rq_clear_multi(&got, &want, NULL);

  test_rmatq_cholesky_case(&C, &A);
  rmatq_clear(&A);
  rmatq_clear(&C);

  // Pivot 2 is not a rational square: no Cholesky decomposition over Q.
  rmatq_init(&A, 2, 2);
  rmatq_init(&C, 2, 2);
  test_set_i64(&A.M, 0, 0, 2);
  test_set_i64(&A.M, 0, 1, 1);
  test_set_i64(&A.M, 1, 0, 1);
  test_set_i64(&A.M, 1, 1, 2);
  assert(rmatq_cholesky(&C, &A) == RABIN_ERR_INVALID_ARG);
  rmatq_clear(&A);
  rmatq_clear(&C);
}

static void test_rmatq_cholesky_4x4(void)
{
  static const i64 vals[4][4] = {
      {4, 12, -16, 0}, {12, 37, -43, 0}, {-16, -43, 98, -9}, {0, 0, -9, 13}};

  rmatq_t A = {0}, C = {0};
  rmatq_init(&A, 4, 4);
  rmatq_init(&C, 4, 4);

  for (u64 i = 0; i < 4; i++) {
    for (u64 j = 0; j < 4; j++) {
      test_set_i64(&A.M, i, j, vals[i][j]);
    }
  }

  printf("\nInput A:\n");
  rmatq_print(&A);

  // Pivots are 4, 1, 9, 4: all rational squares, C = L^t is integral.
  assert(rmatq_cholesky(&C, &A) == RABIN_SUCCESS);

  printf("C:\n");
  rmatq_print(&C);

  const i64 expected[4][4] = {
      {2, 6, -8, 0}, {0, 1, 5, 0}, {0, 0, 3, -3}, {0, 0, 0, 2}};
  rz_t v;
  rq_t got, want;
  rz_init(&v);
  rq_init_multi(&got, &want, NULL);
  for (u64 i = 0; i < 4; i++) {
    for (u64 j = 0; j < 4; j++) {
      rmatq_get(&got, &C, i, j);
      rz_set_i64(&v, expected[i][j]);
      rq_set(&want, &v);
      int eq =
          (rq_is_zero(&got) && rq_is_zero(&want)) ||
          (!rq_is_zero(&got) && !rq_is_zero(&want) && rq_cmp(&got, &want) == 0);
      assert(eq);
    }
  }
  rq_clear_multi(&got, &want, NULL);
  rz_clear(&v);

  test_rmatq_cholesky_case(&C, &A);

  rmatq_clear(&A);
  rmatq_clear(&C);
}

int main()
{
  printf("--- Starting BigMatrix Test Suite ---\n");

  // 1. Test Initialization and Free
  rmat_t A = {0}, B = {0}, R = {0};
  rmat_init(&A, 2, 2);
  rmat_init(&B, 2, 2);
  rmat_init(&R, 2, 2);

  assert(A.data != NULL);
  assert(A.rows == 2 && A.cols == 2);
  printf("[PASS] Initialization\n");

  // 2. Setup Data: A = [[1, 2], [3, 4]]
  rz_t val;
  rz_init(&val);

  rz_set_u64(&val, 1);
  rmat_set(&A, &val, 0, 0);
  rz_set_u64(&val, 2);
  rmat_set(&A, &val, 0, 1);
  rz_set_u64(&val, 3);
  rmat_set(&A, &val, 1, 0);
  rz_set_u64(&val, 4);
  rmat_set(&A, &val, 1, 1);
  // 3. Test rmat_get
  rz_t check;
  rz_init(&check);
  rmat_get(&check, &A, 0, 1);  // Should be 2
  // Replace with your library's comparison function, e.g., rz_cmp_u64
  assert(rz_is_eq_i64(&check, 2));
  printf("[PASS] Get/Set\n");

  // 4. Test rmat_copy
  rmat_copy(&B, &A);
  rmat_get(&check, &B, 1, 1);
  assert(rz_is_eq_i64(&check, 4));
  printf("[PASS] Copy\n");

  // 5. Test rmat_add: R = A + B (Since B=A, R should be [[2, 4], [6, 8]])
  rmat_add(&R, &A, &B);
  rmat_get(&check, &R, 1, 0);
  assert(rz_is_eq_i64(&check, 6));
  printf("[PASS] Addition\n");

  // 6. Test rmat_mul: R = A * Identity
  // Setup B as Identity: [[1, 0], [0, 1]]
  rz_set_u64(&val, 1);
  rmat_set(&B, &val, 0, 0);
  rz_set_u64(&val, 0);
  rmat_set(&B, &val, 0, 1);
  rz_set_u64(&val, 0);
  rmat_set(&B, &val, 1, 0);
  rz_set_u64(&val, 1);
  rmat_set(&B, &val, 1, 1);

  rmat_mul(&R, &A, &B);  // Result should be A
  rmat_get(&check, &R, 1, 1);
  assert(rz_is_eq_i64(&check, 4));
  printf("[PASS] Multiplication (Identity)\n");

  // 7. Test rmat_mul: Rectangular Multiplication
  // A (2x3) * B (3x2) = R (2x2)
  rmat_t RectA = {0}, RectB = {0}, RectR = {0};
  rmat_init(&RectA, 2, 3);
  rmat_init(&RectB, 3, 2);
  rmat_init(&RectR, 2, 2);

  // Populate with 1s
  rz_set_u64(&val, 1);
  for (u64 i = 0; i < 6; i++) {
    rz_copy(&RectA.data[i], &val);
    rz_copy(&RectB.data[i], &val);
  }

  rmat_mul(&RectR, &RectA, &RectB);
  // Each cell in a 2x3 * 3x2 of all 1s should be 3
  rmat_get(&check, &RectR, 0, 0);
  assert(rz_is_eq_i64(&check, 3));
  printf("[PASS] Multiplication (Rectangular)\n");

  rmat_print(&R);

  rmat_det(&val, &A);
  rz_println(&val);

  rmat_hadamard(&val, &A);
  rz_println(&val);

  // 7b. Test the structural check predicates
  {
    // 2x2 zero matrix: square, zero, diagonal, both triangular,
    // symmetric; not the identity
    rmat_t Z = {0};
    rmat_init(&Z, 2, 2);
    assert(rmat_is_square(&Z));
    assert(rmat_is_zero(&Z));
    assert(!rmat_is_identity(&Z));
    assert(rmat_is_diagonal(&Z));
    assert(rmat_is_upper_triangular(&Z));
    assert(rmat_is_lower_triangular(&Z));
    assert(rmat_is_symmetric(&Z));
    rmat_clear(&Z);
    printf("[PASS] Checks (zero matrix)\n");

    // 3x3 identity
    rmat_t I3 = {0};
    rmat_init(&I3, 3, 3);
    test_set_i64(&I3, 0, 0, 1);
    test_set_i64(&I3, 1, 1, 1);
    test_set_i64(&I3, 2, 2, 1);
    assert(rmat_is_square(&I3));
    assert(rmat_is_identity(&I3));
    assert(!rmat_is_zero(&I3));
    assert(rmat_is_diagonal(&I3));
    assert(rmat_is_upper_triangular(&I3));
    assert(rmat_is_lower_triangular(&I3));
    assert(rmat_is_symmetric(&I3));
    rmat_clear(&I3);
    printf("[PASS] Checks (identity)\n");

    // -1 on the diagonal must not count as 1 (rz_is_eq_i64(x, 1) would
    // match it; rz_is_one must not)
    rmat_t N3 = {0};
    rmat_init(&N3, 3, 3);
    test_set_i64(&N3, 0, 0, -1);
    test_set_i64(&N3, 1, 1, 1);
    test_set_i64(&N3, 2, 2, 1);
    assert(!rmat_is_identity(&N3));
    assert(rmat_is_diagonal(&N3));
    assert(rmat_is_symmetric(&N3));
    rmat_clear(&N3);
    printf("[PASS] Checks (negative diagonal)\n");

    // 3x3 upper triangular: not diagonal, not lower, not symmetric
    rmat_t U3 = {0};
    rmat_init(&U3, 3, 3);
    test_set_i64(&U3, 0, 0, 1);
    test_set_i64(&U3, 0, 1, 2);
    test_set_i64(&U3, 0, 2, 3);
    test_set_i64(&U3, 1, 1, 4);
    test_set_i64(&U3, 1, 2, 5);
    test_set_i64(&U3, 2, 2, 6);
    assert(rmat_is_upper_triangular(&U3));
    assert(!rmat_is_lower_triangular(&U3));
    assert(!rmat_is_diagonal(&U3));
    assert(!rmat_is_symmetric(&U3));
    assert(!rmat_is_identity(&U3));
    rmat_clear(&U3);
    printf("[PASS] Checks (upper triangular)\n");

    // 3x3 lower triangular: the mirror case
    rmat_t L3 = {0};
    rmat_init(&L3, 3, 3);
    test_set_i64(&L3, 0, 0, 1);
    test_set_i64(&L3, 1, 0, 2);
    test_set_i64(&L3, 1, 1, 4);
    test_set_i64(&L3, 2, 0, 3);
    test_set_i64(&L3, 2, 1, 5);
    test_set_i64(&L3, 2, 2, 6);
    assert(rmat_is_lower_triangular(&L3));
    assert(!rmat_is_upper_triangular(&L3));
    assert(!rmat_is_diagonal(&L3));
    assert(!rmat_is_symmetric(&L3));
    rmat_clear(&L3);
    printf("[PASS] Checks (lower triangular)\n");

    // 3x3 symmetric, then break one mirror pair
    rmat_t S3 = {0};
    rmat_init(&S3, 3, 3);
    test_set_i64(&S3, 0, 0, 1);
    test_set_i64(&S3, 0, 1, 2);
    test_set_i64(&S3, 0, 2, 3);
    test_set_i64(&S3, 1, 0, 2);
    test_set_i64(&S3, 1, 1, 4);
    test_set_i64(&S3, 1, 2, 5);
    test_set_i64(&S3, 2, 0, 3);
    test_set_i64(&S3, 2, 1, 5);
    test_set_i64(&S3, 2, 2, 6);
    assert(rmat_is_symmetric(&S3));
    test_set_i64(&S3, 0, 2, 9);
    assert(!rmat_is_symmetric(&S3));
    rmat_clear(&S3);
    printf("[PASS] Checks (symmetric)\n");

    // 2x3 zero rectangle: not square / not identity / not symmetric,
    // the rest vacuously true
    rmat_t R23 = {0};
    rmat_init(&R23, 2, 3);
    assert(!rmat_is_square(&R23));
    assert(rmat_is_zero(&R23));
    assert(!rmat_is_identity(&R23));
    assert(rmat_is_diagonal(&R23));
    assert(rmat_is_upper_triangular(&R23));
    assert(rmat_is_lower_triangular(&R23));
    assert(!rmat_is_symmetric(&R23));
    // a single sub-diagonal entry breaks upper triangularity (and
    // diagonal/zero), but not lower triangularity
    test_set_i64(&R23, 1, 0, 7);
    assert(!rmat_is_zero(&R23));
    assert(!rmat_is_upper_triangular(&R23));
    assert(rmat_is_lower_triangular(&R23));
    assert(!rmat_is_diagonal(&R23));
    rmat_clear(&R23);
    printf("[PASS] Checks (rectangular)\n");

    // 3x4 identity-pattern: the dimension short-circuit must reject it
    // as identity/symmetric even though the scanned sub-regions are clean
    rmat_t Q34 = {0};
    rmat_init(&Q34, 3, 4);
    test_set_i64(&Q34, 0, 0, 1);
    test_set_i64(&Q34, 1, 1, 1);
    test_set_i64(&Q34, 2, 2, 1);
    assert(!rmat_is_square(&Q34));
    assert(!rmat_is_identity(&Q34));
    assert(!rmat_is_symmetric(&Q34));
    assert(rmat_is_diagonal(&Q34));
    assert(rmat_is_upper_triangular(&Q34));
    assert(rmat_is_lower_triangular(&Q34));
    rmat_clear(&Q34);
    printf("[PASS] Checks (non-square short-circuit)\n");

    // degenerate sizes: 0x0 is vacuously everything, 0x2 is not square
    // nor symmetric
    rmat_t E00 = {0}, E02 = {0};
    rmat_init(&E00, 0, 0);
    assert(rmat_is_square(&E00));
    assert(rmat_is_zero(&E00));
    assert(rmat_is_identity(&E00));
    assert(rmat_is_diagonal(&E00));
    assert(rmat_is_upper_triangular(&E00));
    assert(rmat_is_lower_triangular(&E00));
    assert(rmat_is_symmetric(&E00));
    rmat_clear(&E00);

    rmat_init(&E02, 0, 2);
    assert(!rmat_is_square(&E02));
    assert(rmat_is_zero(&E02));
    assert(!rmat_is_identity(&E02));
    assert(rmat_is_diagonal(&E02));
    assert(rmat_is_upper_triangular(&E02));
    assert(rmat_is_lower_triangular(&E02));
    assert(!rmat_is_symmetric(&E02));
    rmat_clear(&E02);
    printf("[PASS] Checks (degenerate sizes)\n");
  }

  // 8. Cleanup
  rz_clear(&val);
  rz_clear(&check);
  rmat_clear(&A);
  rmat_clear(&B);
  rmat_clear(&R);
  rmat_clear(&RectA);
  rmat_clear(&RectB);
  rmat_clear(&RectR);

  rmat_t D = {0}, ADJ = {0};
  rmat_init(&D, 20, 20);
  rmat_init(&ADJ, 20, 20);
  rz_t tmp;
  rz_init(&tmp);

  // set matrix D
  for (i64 i = 0; i < 20; i++) {
    for (i64 j = 0; j < 20; j++) {
      rz_set_i64(&tmp, (i + j) % 20);
      rmat_set(&D, &tmp, i, j);
    }
  }

  rmat_print(&D);
  rmat_det(&tmp, &D);
  rz_println(&tmp);

  rmat_trace(&tmp, &D);
  rz_println(&tmp);

  rpol_t p = {0};
  rpol_init(&p);

  rmat_charpoly_adj(&p, &ADJ, &D);

  printf("Char poly: ");
  rpol_print(&p);
  printf("\n");

  printf("Adjoint: \n ");
  rmat_print(&ADJ);
  printf("\n");

  printf("Hermite: \n ");
  rmat_hermite(&ADJ, &D);
  assert(rmat_hnf_check_structure(&ADJ));
  rmat_print(&ADJ);
  printf("\n");

  printf("Hermite GCD: \n ");
  rmat_hermite_gcd(&ADJ, &D);
  rmat_print(&ADJ);
  printf("\n");

  test_hermite_mod_d();
  test_hnf_verify();
  test_hnf_mod_d_multilimb();
  test_lll();
  test_mlll();

  rpol_clear(&p);

  rmat_clear(&ADJ);

  rns_ctx_t ctx = {0};
  if (rns_ctx_init(&ctx, RNS_PRIMES, 10) != RABIN_SUCCESS) return 1;

  rmat_det_rns(&tmp, &D, &ctx);

  printf("Det rns: ");
  rz_println(&tmp);

  rmat_det_bareiss(&tmp, &D);
  printf("Det bareiss: ");
  rz_println(&tmp);

  rns_ctx_clear(&ctx);

  rmat_clear(&D);

  rmat_t E = {0};
  rmat_init(&E, 3, 4);
  rvec_t v = {0}, res = {0};
  rvec_init(&v, 4);
  rvec_init(&res, 3);

  rmat_set(&E, &tmp, 0, 0);
  rz_add_u64(&tmp, &tmp, 1);
  rmat_set(&E, &tmp, 1, 0);
  rz_add_u64(&tmp, &tmp, 1);
  rmat_set(&E, &tmp, 1, 1);
  rz_add_u64(&tmp, &tmp, 1);
  rmat_set(&E, &tmp, 2, 0);
  rz_add_u64(&tmp, &tmp, 1);
  rmat_set(&E, &tmp, 2, 1);
  rz_add_u64(&tmp, &tmp, 1);
  rmat_set(&E, &tmp, 2, 2);

  rmat_print(&E);

  rvec_set(&v, &tmp, 0);

  rvec_set(&v, &tmp, 1);
  rvec_set(&v, &tmp, 2);
  rvec_set(&v, &tmp, 3);

  rvec_print(&v);
  printf("\n");

  rmat_mv(&res, &E, &v);

  printf("Matrix vector: ");
  rvec_print(&res);
  printf("\n");

  u64* testmatrix = malloc(sizeof(u64) * 9);
  u64 values[] = {18446744073709551610, 1, 1, 1, 18446744073709551600, 1, 1, 1,
                  18446744073709551590};
  memcpy(testmatrix, values, sizeof(u64) * 9);

  u64_mat_t A_test;
  A_test.data = testmatrix;
  A_test.cols = 3;
  A_test.rows = 3;
  A_test.modulus = 100000007;

  u64 det = u64_mat_det(&A_test);
  printf("u64 determinant: %llu \n", det);

  free(testmatrix);

  test_pascal_det(50);

  test_rmatq_cholesky();
  printf("[PASS] Rmatq Cholesky\n");

  test_rmatq_cholesky_4x4();
  printf("[PASS] Rmatq Cholesky (4x4)\n");

  printf("[PASS] Cleanup / Free\n");
  printf("--- All Tests Passed! ---\n");

  u64 sizes[] = {2, 4, 6, 8, 16, 32, 64, 128, 256, 300, 512, 700, 994, 1024};
  int num_tests = sizeof(sizes) / sizeof(sizes[0]);

  // for (int i = 0; i < num_tests; i++) {
  //   run_hadamard_benchmark(sizes[i], 32);
  // }

  for (int i = 0; i < 9; i++) {
    run_det_benchmark(sizes[i], 32);
  }

  for (int i = 0; i < 7; i++) {
    run_hermite_benchmark(sizes[i], 100);
  }

  // HNF mod D (2.4.8) vs exact (2.4.5): on the smaller sizes both are timed
  // on the same random matrix, on the larger ones only the mod D version
  u64 hnf_sizes[] = {4, 8, 16, 32, 64, 128};
  int n_hnf = sizeof(hnf_sizes) / sizeof(hnf_sizes[0]);
  for (int i = 0; i < 8; i++) {
    run_hermite_mod_d_benchmark(sizes[i], 10, i < 4);
  }

  for (int i = 0; i < 10; i++) {
    run_mul_benchmark(sizes[i], 64);
  }

  rz_clear(&tmp);

  rmat_clear(&E);
  rvec_clear(&v);
  rvec_clear(&res);

  return 0;
}
