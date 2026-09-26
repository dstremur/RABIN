#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include "../include/bigmath.h"
#include "../include/bigmatrix.h"
#include "../include/bigpoly.h"
#include "../include/bigrand.h"
#include "../include/bigrns.h"
#include "../include/primes.h"
#include "../include/u64.h"
double get_time()
{
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return ts.tv_sec + ts.tv_nsec * 1e-9;
}

void run_det_benchmark(u64 size, u64 bits)
{
  bigmatrix M;
  bignum det, val;

  bn_init(&det);
  bn_init(&val);
  bigmatrix_init(&M, size, size);

  ctx_rns ctx;

  // Populate with random data to prevent "easy" zeros
  for (u64 i = 0; i < size; i++) {
    for (u64 j = 0; j < size; j++) {
      bn_gen_random(&val, bits);
      bigmatrix_set(&M, &val, i, j);
    }
  }

  printf("Benchmarking %llu x %llu (%llu-bit entries) Determinant... ", size,
         size, bits);
  fflush(stdout);

  // bigmatrix_print_python(&M);

  u64 k = rns_estimate_determinant(&M);
  // printf("k: %llu \n", k);

  if ((k + 1) > 500) {
    rns_context_init(&ctx, RNS_PRIMES2, k + 1);
  } else {
    rns_context_init(&ctx, RNS_PRIMES, k + 1);
  }

  double start = get_time();
  bigmatrix_det_rns(&det, &M, &ctx);
  double end = get_time();

  printf("Time: %f seconds\n", end - start);

  // printf("rns: ");
  // bn_println(&det);

  bigmatrix_free(&M);
  bn_free(&det);
  bn_free(&val);
  rns_context_free(&ctx);
}

void run_hadamard_benchmark(u64 size, u64 bits)
{
  bigmatrix M;
  bignum det, val;

  bn_init(&det);
  bn_init(&val);
  bigmatrix_init(&M, size, size);

  // Populate with random data to prevent "easy" zeros
  for (u64 i = 0; i < size; i++) {
    for (u64 j = 0; j < size; j++) {
      bn_gen_random(&val, bits);
      bigmatrix_set(&M, &val, i, j);
    }
  }

  printf("Benchmarking %llu x %llu (%llu-bit entries) hadamard bound... ", size,
         size, bits);
  fflush(stdout);

  double start = get_time();
  bigmatrix_hadamard(&det, &M);
  double end = get_time();

  printf("Time: %f seconds\n", end - start);

  bigmatrix_free(&M);
  bn_free(&det);
  bn_free(&val);
}

void run_hermite_benchmark(u64 size, u64 bits)
{
  bigmatrix M, H;
  bignum det, val;

  bn_init(&det);
  bn_init(&val);
  bigmatrix_init(&M, size, size);
  bigmatrix_init(&H, size, size);

  // Populate with random data to prevent "easy" zeros
  for (u64 i = 0; i < size; i++) {
    for (u64 j = 0; j < size; j++) {
      bn_gen_random(&val, bits);
      bigmatrix_set(&M, &val, i, j);
    }
  }

  bigmatrix_det(&det, &M); 

  det.is_neg = false;
  printf(
      "Benchmarking %llu x %llu (%llu-bit entries) Hermite normal form...\n ",
      size, size, bits);

  double start = get_time();
bigmatrix_hermite_mod_d(&H, &M, &det);
  double end = get_time();

//bigmatrix_print_tail(&H);

  printf("Time: %f seconds\n", end - start);

  bigmatrix_free(&M);
  bigmatrix_free(&H);
  bn_free(&det);
  bn_free(&val);
}

void run_mul_benchmark(u64 size, u64 bits)
{
  bigmatrix M, A, B;
  bignum det, val;

  bn_init(&det);
  bn_init(&val);
  bigmatrix_init(&M, size, size);
  bigmatrix_init(&A, size, size);
  bigmatrix_init(&B, size, size);

  // Populate with random data to prevent "easy" zeros
  for (u64 i = 0; i < size; i++) {
    for (u64 j = 0; j < size; j++) {
      bn_gen_random(&val, bits);
      bigmatrix_set(&A, &val, i, j);
      bn_gen_random(&val, bits);
      bigmatrix_set(&B, &val, i, j);
    }
  }

  printf("Benchmarking %llu x %llu (%llu-bit entries) matrix multiply ... ",
         size, size, bits);
  fflush(stdout);

  double start = get_time();
  bigmatrix_mul(&M, &A, &B);
  double end = get_time();

  printf("Time: %f seconds\n", end - start);

  bigmatrix_free(&M);
  bigmatrix_free(&A);
  bigmatrix_free(&B);
  bn_free(&det);
  bn_free(&val);
}

void run_hermite_mod_d_benchmark(u64 size, u64 bits, bool with_exact)
{
  bigmatrix M, H;
  bignum D, val;

  bn_init(&D);
  bn_init(&val);
  bigmatrix_init(&M, size, size);
  bigmatrix_init(&H, size, size);

  // Populate with random data to prevent "easy" zeros
  for (u64 i = 0; i < size; i++) {
    for (u64 j = 0; j < size; j++) {
      bn_gen_random(&val, bits);
      bigmatrix_set(&M, &val, i, j);
    }
  }

  // D must be a multiple of the module determinant; for square M the
  // module determinant is |det M|, so use that
  for (int rep = 0; rep < 20; rep++) {
    bigmatrix_det(&D, &M);
    if (!bn_is_zero(&D)) break;
    for (u64 i = 0; i < size; i++) {
      bn_gen_random(&val, bits);
      bigmatrix_set(&M, &val, i, i);
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
    bigmatrix_hermite_gcd(&H, &M);
    double end = get_time();

    printf("Time: %f seconds\n", end - start);
  }

  printf(
      "Benchmarking %llu x %llu (%llu-bit entries) Hermite normal form "
      "(mod D 2.4.8)... ",
      size, size, bits);
  fflush(stdout);

  double start = get_time();
  bigmatrix_hermite_mod_d(&H, &M, &D);
  double end = get_time();

  printf("Time: %f seconds\n", end - start);

  bigmatrix_free(&M);
  bigmatrix_free(&H);
  bn_free(&D);
  bn_free(&val);
}

static void test_set_i64(bigmatrix* M, u64 r, u64 c, i64 v)
{
  bignum tmp;
  bn_init(&tmp);
  bn_set_i64(&tmp, v);
  bigmatrix_set(M, &tmp, r, c);
  bn_free(&tmp);
}

void test_hermite_mod_d(void)
{
  bignum D, val;
  bn_init_multi(&D, &val, NULL);

  // 1. 2x2 with known HNF: A = [[1,2],[3,4]] -> W = [[2,1],[0,1]]
  {
    bigmatrix A, W, H;
    bigmatrix_init(&A, 2, 2);
    bigmatrix_init(&W, 2, 2);
    bigmatrix_init(&H, 2, 2);
    test_set_i64(&A, 0, 0, 1);
    test_set_i64(&A, 0, 1, 2);
    test_set_i64(&A, 1, 0, 3);
    test_set_i64(&A, 1, 1, 4);

    // D = 4 = 2*|det A|, a multiple of the module determinant
    bn_set_u64(&D, 4);
    bigmatrix_hermite_mod_d(&W, &A, &D);
    test_set_i64(&H, 0, 0, 2);
    test_set_i64(&H, 0, 1, 1);
    test_set_i64(&H, 1, 0, 0);
    test_set_i64(&H, 1, 1, 1);
    assert(bigmatrix_equal(&W, &H));
    assert(bigmatrix_hnf_check_structure(&W));

    // cross-check against the exact implementation
    bigmatrix V;
    bigmatrix_init(&V, 2, 2);
    bigmatrix_hermite(&V, &A);
    assert(bigmatrix_equal(&W, &V));
    bigmatrix_free(&V);

    bigmatrix_free(&A);
    bigmatrix_free(&W);
    bigmatrix_free(&H);
    printf("[PASS] hermite_mod_d 2x2 known example\n");
  }

  // 2. 2x3 with unimodular column module: gcd of the 2x2 minors is 1,
  //    so D = 1 is valid and the HNF must be the 2x2 identity
  {
    bigmatrix A, W, H;
    bigmatrix_init(&A, 2, 3);
    bigmatrix_init(&W, 2, 2);
    bigmatrix_init(&H, 2, 2);
    test_set_i64(&A, 0, 0, 6);
    test_set_i64(&A, 0, 1, 4);
    test_set_i64(&A, 0, 2, 9);
    test_set_i64(&A, 1, 0, 8);
    test_set_i64(&A, 1, 1, 1);
    test_set_i64(&A, 1, 2, 6);

    bn_set_u64(&D, 1);
    bigmatrix_hermite_mod_d(&W, &A, &D);
    bigmatrix_id(&H, 2);
    assert(bigmatrix_equal(&W, &H));
    assert(bigmatrix_hnf_check_structure(&W));
    printf("[PASS] hermite_mod_d 2x3 unimodular (D = 1)\n");

    bigmatrix_free(&A);
    bigmatrix_free(&W);
    bigmatrix_free(&H);
  }

  // 3. Random cross-checks against the exact implementation (HNF is unique).
  //    D must be a multiple of the module determinant Delta (GCD of the
  //    m x m minors); build A = [B | B R] so that Delta = |det B|:
  {
    u64 shapes[][2] = {{2, 2}, {2, 3}, {3, 3}, {3, 4}, {4, 4}, {2, 5}};
    u64 n_shapes = sizeof(shapes) / sizeof(shapes[0]);
    bignum det, three;
    bn_init(&det);
    bn_set_u64(&three, 3);
    for (u64 s = 0; s < n_shapes; s++) {
      u64 m = shapes[s][0];
      u64 n = shapes[s][1];
      for (u64 trial = 0; trial < 4; trial++) {
        bigmatrix B, R, C, A, H, W1, W2;
        bigmatrix_init(&B, m, m);
        bigmatrix_init(&R, m, n - m);
        bigmatrix_init(&C, m, n - m);
        bigmatrix_init(&A, m, n);
        bigmatrix_init(&H, m, m);
        bigmatrix_init(&W1, m, m);
        bigmatrix_init(&W2, m, m);

        for (int rep = 0; rep < 10; rep++) {
          for (u64 i = 0; i < m; i++) {
            for (u64 j = 0; j < m; j++) {
              bn_gen_random(&val, 20);
              bigmatrix_set(&B, &val, i, j);
            }
          }
          bigmatrix_det(&det, &B);
          if (!bn_is_zero(&det)) break;
        }
        assert(!bn_is_zero(&det));
        det.is_neg = false;

        for (u64 i = 0; i < m; i++) {
          for (u64 j = 0; j < n - m; j++) {
            bn_gen_random(&val, 20);
            bigmatrix_set(&R, &val, i, j);
          }
        }
        if (m < n) {
          bigmatrix_mul(&C, &B, &R);
        }
        // A = [B | B R]
        for (u64 i = 0; i < m; i++) {
          for (u64 j = 0; j < m; j++) {
            bigmatrix_set(&A, &B.data[i * m + j], i, j);
          }
          for (u64 j = 0; j < n - m; j++) {
            bigmatrix_set(&A, &C.data[i * (n - m) + j], i, m + j);
          }
        }

        bigmatrix_hermite(&H, &A);
        assert(H.r_size == m && H.c_size == m);

        // D = |det B| = the module determinant itself
        bigmatrix_hermite_mod_d(&W1, &A, &det);
        assert(bigmatrix_equal(&W1, &H));
        assert(bigmatrix_hnf_check_structure(&W1));
        assert(bigmatrix_hnf_check_structure(&H));

        // D = 3 * |det B|, another valid multiple
        bn_mul(&D, &three, &det);
        bigmatrix_hermite_mod_d(&W2, &A, &D);
        assert(bigmatrix_equal(&W2, &H));
        assert(bigmatrix_hnf_check_structure(&W2));

        bigmatrix_free(&B);
        bigmatrix_free(&R);
        bigmatrix_free(&C);
        bigmatrix_free(&A);
        bigmatrix_free(&H);
        bigmatrix_free(&W1);
        bigmatrix_free(&W2);
      }
      printf("[PASS] hermite_mod_d random %llux%llu vs exact\n", m, n);
    }
    bn_free(&det);
    bn_free(&three);
  }

  // 4. m > n falls back to the exact implementation
  {
    u64 m = 4, n = 3;
    bigmatrix A, W, H;
    bigmatrix_init(&A, m, n);
    bigmatrix_init(&W, m, n);
    bigmatrix_init(&H, m, n);
    for (u64 i = 0; i < m; i++) {
      for (u64 j = 0; j < n; j++) {
        bn_gen_random(&val, 20);
        bigmatrix_set(&A, &val, i, j);
      }
    }
    bn_set_u64(&D, 1000003);
    bigmatrix_hermite_mod_d(&W, &A, &D);
    bigmatrix_hermite(&H, &A);
    assert(bigmatrix_equal(&W, &H));
    printf("[PASS] hermite_mod_d m > n fallback\n");

    bigmatrix_free(&A);
    bigmatrix_free(&W);
    bigmatrix_free(&H);
  }

  bn_free_multi(&D, &val, NULL);
}

void test_hnf_verify(void)
{
  bignum val;
  bn_init(&val);

  // 1. Structure: valid 2x2 with one violation per rule
  {
    // [[3,1],[0,2]]: upper triangular, positive pivots, 0 <= H[0][1] = 1 < 2
    bigmatrix H;
    bigmatrix_init(&H, 2, 2);
    test_set_i64(&H, 0, 0, 3);
    test_set_i64(&H, 0, 1, 1);
    test_set_i64(&H, 1, 1, 2);
    assert(bigmatrix_hnf_check_structure(&H));

    // violation: entry below the diagonal
    test_set_i64(&H, 1, 0, 5);
    assert(!bigmatrix_hnf_check_structure(&H));
    test_set_i64(&H, 1, 0, 0);

    // violation: negative pivot
    test_set_i64(&H, 1, 1, -2);
    assert(!bigmatrix_hnf_check_structure(&H));
    test_set_i64(&H, 1, 1, 2);

    // violation: entry right of the pivot not reduced (3 is not < 3)
    test_set_i64(&H, 0, 1, 3);
    assert(!bigmatrix_hnf_check_structure(&H));

    // violation: negative entry right of the pivot
    test_set_i64(&H, 0, 1, -1);
    assert(!bigmatrix_hnf_check_structure(&H));
    test_set_i64(&H, 0, 1, 1);
    assert(bigmatrix_hnf_check_structure(&H));
    bigmatrix_free(&H);
    printf("[PASS] hnf_check_structure (2x2 rules)\n");

    // zero-row order: trailing zeros ok, a nonzero row after them is not
    bigmatrix Z;
    bigmatrix_init(&Z, 3, 3);
    test_set_i64(&Z, 0, 0, 1);
    assert(bigmatrix_hnf_check_structure(&Z));
    test_set_i64(&Z, 2, 2, 3);
    assert(!bigmatrix_hnf_check_structure(&Z));
    bigmatrix_free(&Z);

    // zero pivot with a nonzero entry right of it: not a zero row
    bigmatrix B;
    bigmatrix_init(&B, 2, 2);
    test_set_i64(&B, 0, 1, 1);
    assert(!bigmatrix_hnf_check_structure(&B));
    bigmatrix_free(&B);

    // valid 3x3: [[2,1,0],[0,3,1],[0,0,5]]
    bigmatrix T;
    bigmatrix_init(&T, 3, 3);
    test_set_i64(&T, 0, 0, 2);
    test_set_i64(&T, 0, 1, 1);
    test_set_i64(&T, 1, 1, 3);
    test_set_i64(&T, 1, 2, 1);
    test_set_i64(&T, 2, 2, 5);
    assert(bigmatrix_hnf_check_structure(&T));
    bigmatrix_free(&T);
    printf("[PASS] hnf_check_structure (zero rows, 3x3)\n");
  }

  // 2. Unimodularity: det 1, det -1, shear; det 2 and non-square fail
  {
    bigmatrix U;
    bigmatrix_init(&U, 2, 2);

    // identity: det 1
    test_set_i64(&U, 0, 0, 1);
    test_set_i64(&U, 1, 1, 1);
    assert(bigmatrix_hnf_check_unimodular(&U));

    // swap: det -1, |det| = 1
    test_set_i64(&U, 0, 0, 0);
    test_set_i64(&U, 0, 1, 1);
    test_set_i64(&U, 1, 0, 1);
    test_set_i64(&U, 1, 1, 0);
    assert(bigmatrix_hnf_check_unimodular(&U));

    // shear: det 1
    test_set_i64(&U, 0, 0, 1);
    test_set_i64(&U, 0, 1, 1);
    test_set_i64(&U, 1, 0, 0);
    test_set_i64(&U, 1, 1, 1);
    assert(bigmatrix_hnf_check_unimodular(&U));

    // det 2: not unimodular
    test_set_i64(&U, 0, 0, 2);
    test_set_i64(&U, 1, 1, 1);
    assert(!bigmatrix_hnf_check_unimodular(&U));

    bigmatrix_free(&U);

    // non-square: rejected without a determinant
    bigmatrix R;
    bigmatrix_init(&R, 2, 3);
    assert(!bigmatrix_hnf_check_unimodular(&R));
    bigmatrix_free(&R);
    printf("[PASS] hnf_check_unimodular\n");
  }

  // 3. Transformation: the library's known example, A * U = H
  //    A = [[1,2],[3,4]], U = [[-4,-1],[3,1]] (det -1),
  //    A * U = [[2,1],[0,1]] = H
  bigmatrix A, U, H;
  bigmatrix_init(&A, 2, 2);
  bigmatrix_init(&U, 2, 2);
  bigmatrix_init(&H, 2, 2);
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
  bigmatrix P;
  bigmatrix_init(&P, 2, 2);
  bigmatrix_mul(&P, &A, &U);
  assert(bigmatrix_equal(&P, &H));
  bigmatrix_free(&P);

  assert(bigmatrix_hnf_check_transformation(&A, &H, &U));

  // corrupt H: product mismatch
  bigmatrix Hbad;
  bigmatrix_init(&Hbad, 2, 2);
  bigmatrix_copy(&Hbad, &H);
  test_set_i64(&Hbad, 0, 1, 7);
  assert(!bigmatrix_hnf_check_transformation(&A, &Hbad, &U));
  bigmatrix_free(&Hbad);

  // shape mismatch: U must be 2x2 for a 2x2 A
  bigmatrix Ubad;
  bigmatrix_init(&Ubad, 1, 2);
  assert(!bigmatrix_hnf_check_transformation(&A, &H, &Ubad));
  bigmatrix_free(&Ubad);
  printf("[PASS] hnf_check_transformation\n");

  // 4. Master check: the full triple passes
  assert(bigmatrix_hnf_verify(&A, &H, &U));

  // structure violation is caught first (lower triangle entry)
  bigmatrix Hs;
  bigmatrix_init(&Hs, 2, 2);
  bigmatrix_copy(&Hs, &H);
  test_set_i64(&Hs, 1, 0, 5);
  assert(!bigmatrix_hnf_verify(&A, &Hs, &U));
  bigmatrix_free(&Hs);

  // structure-valid H, but A * U != H: [[2,0],[0,1]] is in HNF shape
  bigmatrix Ht;
  bigmatrix_init(&Ht, 2, 2);
  test_set_i64(&Ht, 0, 0, 2);
  test_set_i64(&Ht, 1, 1, 1);
  assert(bigmatrix_hnf_check_structure(&Ht));
  assert(!bigmatrix_hnf_verify(&A, &Ht, &U));
  bigmatrix_free(&Ht);

  // H = A * U' with U' = 2 * U (det -4): structure and transformation
  // hold, unimodularity must be the failing check
  bigmatrix U2, H2;
  bigmatrix_init(&U2, 2, 2);
  bigmatrix_init(&H2, 2, 2);
  test_set_i64(&U2, 0, 0, -8);
  test_set_i64(&U2, 0, 1, -2);
  test_set_i64(&U2, 1, 0, 6);
  test_set_i64(&U2, 1, 1, 2);
  bigmatrix_mul(&H2, &A, &U2);
  assert(bigmatrix_hnf_check_structure(&H2));
  assert(bigmatrix_hnf_check_transformation(&A, &H2, &U2));
  assert(!bigmatrix_hnf_check_unimodular(&U2));
  assert(!bigmatrix_hnf_verify(&A, &H2, &U2));
  bigmatrix_free(&U2);
  bigmatrix_free(&H2);
  printf("[PASS] hnf_verify (master)\n");

  // 5. Working HNF test: bigmatrix_hermite on random square matrices
  //    must produce matrices in HNF structure
  for (u64 n = 2; n <= 6; n++) {
    for (int rep = 0; rep < 3; rep++) {
      bigmatrix R, W;
      bigmatrix_init(&R, n, n);
      bigmatrix_init(&W, n, n);
      for (u64 i = 0; i < n; i++) {
        for (u64 j = 0; j < n; j++) {
          bn_gen_random(&val, 20);
          bigmatrix_set(&R, &val, i, j);
        }
      }
      bigmatrix_hermite(&W, &R);
      assert(bigmatrix_hnf_check_structure(&W));
      bigmatrix_free(&R);
      bigmatrix_free(&W);
    }
  }
  printf("[PASS] hnf structure of bigmatrix_hermite output\n");

  bigmatrix_free(&A);
  bigmatrix_free(&U);
  bigmatrix_free(&H);
  bn_free(&val);
}

void test_pascal_det(u64 size)
{
  bigmatrix P;
  bignum det, val, a, b;  // Move a and b here
  bn_init(&det);
  bn_init(&val);
  bn_init(&a);
  bn_init(&b);

  bigmatrix_init(&P, size, size);

  ctx_rns ctx;
  // Ensure RNS_PRIMES has at least 50 elements!
  rns_context_init(&ctx, RNS_PRIMES, 10);

  for (u64 i = 0; i < size; i++) {
    for (u64 j = 0; j < size; j++) {
      if (i == 0 || j == 0) {
        bn_set_u64(&val, 1);
      } else {
        // Reuse a and b instead of re-allocating
        bigmatrix_get(&a, &P, i - 1, j);
        bigmatrix_get(&b, &P, i, j - 1);
        bn_add(&val, &a, &b);
      }
      bigmatrix_set(&P, &val, i, j);
    }
  }

  // Standard Det (for comparison)
  bigmatrix_det(&det, &P);
  printf("Pascal %llu x %llu | Standard Det: ", size, size);
  bn_println(&det);

  // RNS Det
  bigmatrix_det_rns(&det, &P, &ctx);
  printf("Pascal %llu x %llu | RNS Det:      ", size, size);
  bn_println(&det);

  // Verification
  if (bn_is_eq_i64(&det, 1)) {
    printf("RESULT: PASS\n");
  } else {
    printf("RESULT: FAIL (Expected 1)\n");
  }

  // CLEANUP
  bigmatrix_free(&P);
  bn_free(&det);
  bn_free(&val);
  bn_free(&a);
  bn_free(&b);

  rns_context_free(&ctx);
}

// GCD of all m x m minors of A (m = A->r_size <= n = A->c_size): the
// determinant Delta of the Z-module generated by the columns of A.
static void test_minors_gcd(bignum* g, const bigmatrix* A, u64 m, u64 n,
                            u64 start, u64 left, u64* cols)
{
  if (left == 0) {
    bigmatrix S;
    bignum det;
    bigmatrix_init(&S, m, m);
    bn_init(&det);
    for (u64 i = 0; i < m; i++) {
      for (u64 j = 0; j < m; j++) {
        bigmatrix_set(&S, &A->data[i * n + cols[j]], i, j);
      }
    }
    bigmatrix_det(&det, &S);
    det.is_neg = false;
    bn_gcd(g, g, &det);
    bn_free(&det);
    bigmatrix_free(&S);
    return;
  }

  for (u64 c = start; c <= n - left; c++) {
    cols[m - left] = c;
    test_minors_gcd(g, A, m, n, c + 1, left - 1, cols);
  }
}

static void test_module_determinant(bignum* delta, const bigmatrix* A)
{
  u64 m = A->r_size;
  u64 n = A->c_size;
  bignum g;
  u64* cols = malloc(m * sizeof(u64));
  bn_init(&g);
  bn_set_u64(&g, 0);
  test_minors_gcd(&g, A, m, n, 0, m, cols);
  bn_copy(delta, &g);
  free(cols);
  bn_free(&g);
}

void test_hnf_mod_d_multilimb(void)
{
  bignum D, val, bit, delta, mult;
  bn_init_multi(&D, &val, &bit, &delta, &mult, NULL);

  // 1. Zero-pivot multi-limb input - regression for the out-of-bounds top
  //    limb read in bn_gcd_extended_lehmer when a size-0 bignum is the
  //    first gcd argument: A = [[p, q], [r, 0]] (all entries 2-limb) has
  //    HNF [[q, p mod q], [0, r]], here p < q so p mod q = p
  {
    bigmatrix A, W, H;
    bignum p, q, r, zero;
    bn_init_multi(&p, &q, &r, &zero, NULL);
    bigmatrix_init(&A, 2, 2);
    bigmatrix_init(&W, 2, 2);
    bigmatrix_init(&H, 2, 2);

    // p = 2^70, q = 2^70 + 3, r = 2^65 + 7
    bn_set_u64(&p, 1);
    bn_lshift(&p, &p, 70);
    bn_set_u64(&q, 1);
    bn_lshift(&q, &q, 70);
    bn_add_u64(&q, &q, 3);
    bn_set_u64(&r, 1);
    bn_lshift(&r, &r, 65);
    bn_add_u64(&r, &r, 7);
    bn_set_u64(&zero, 0);

    bigmatrix_set(&A, &p, 0, 0);
    bigmatrix_set(&A, &q, 0, 1);
    bigmatrix_set(&A, &r, 1, 0);
    bigmatrix_set(&A, &zero, 1, 1);

    // D = 3 * q * r = 3 * |det A|, a multiple of the module determinant
    bn_mul(&D, &q, &r);
    bn_set_u64(&mult, 3);
    bn_mul(&D, &mult, &D);

    bigmatrix_hermite_mod_d(&W, &A, &D);

    bigmatrix_set(&H, &q, 0, 0);
    bigmatrix_set(&H, &p, 0, 1);
    bigmatrix_set(&H, &zero, 1, 0);
    bigmatrix_set(&H, &r, 1, 1);
    assert(bigmatrix_equal(&W, &H));
    assert(bigmatrix_hnf_check_structure(&W));

    // cross-check against the exact implementation
    bigmatrix V;
    bigmatrix_init(&V, 2, 2);
    bigmatrix_hermite(&V, &A);
    assert(bigmatrix_equal(&W, &V));

    bigmatrix_free(&V);
    bigmatrix_free(&A);
    bigmatrix_free(&W);
    bigmatrix_free(&H);
    bn_free_multi(&p, &q, &r, &zero, NULL);
    printf("[PASS] hermite_mod_d 2x2 multi-limb zero pivot\n");
  }

  // 2. bn_gcd_extended_lehmer with a zero argument (multi-limb partner):
  //    gcd(0, b) -> d = |b|, u = 0, v = sign(b); gcd(a, 0) -> d = |a|,
  //    u = sign(a), v = 0
  {
    bignum zero, a, b, bneg, u, v, d, n1;
    bn_init_multi(&zero, &a, &b, &bneg, &u, &v, &d, &n1, NULL);
    bn_set_i64(&n1, -1);
    bn_set_u64(&zero, 0);
    bn_set_u64(&a, 1);
    bn_lshift(&a, &a, 70);
    bn_add_u64(&a, &a, 5);  // 2^70 + 5
    bn_set_u64(&b, 1);
    bn_lshift(&b, &b, 66);  // 2^66

    bn_gcd_extended_lehmer(&u, &v, &d, &zero, &b);
    assert(bn_is_eq_i64(&u, 0));
    assert(bn_is_eq_i64(&v, 1));
    assert(bn_cmp(&d, &b) == 0);

    bn_gcd_extended_lehmer(&u, &v, &d, &a, &zero);
    assert(bn_is_eq_i64(&u, 1));
    assert(bn_is_eq_i64(&v, 0));
    assert(bn_cmp(&d, &a) == 0);

    bn_copy(&bneg, &b);
    bneg.is_neg = true;
    bn_gcd_extended_lehmer(&u, &v, &d, &zero, &bneg);
    assert(bn_is_eq_i64(&u, 0));
    assert(bn_cmp(&v, &n1) == 0);
    assert(bn_cmp(&d, &b) == 0);
    printf("[PASS] gcd_extended_lehmer zero arguments (multi-limb)\n");
    bn_free_multi(&zero, &a, &b, &bneg, &u, &v, &d, &n1, NULL);
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
          bigmatrix A, H, W1, W2;
          bigmatrix_init(&A, m, n);
          bigmatrix_init(&H, m, m);
          bigmatrix_init(&W1, m, m);
          bigmatrix_init(&W2, m, m);

          for (u64 i = 0; i < m; i++) {
            for (u64 j = 0; j < n; j++) {
              bn_gen_random(&val, bits);
              bn_gen_random(&bit, 1);
              if (!bn_is_zero(&bit)) {
                bn_neg(&val, &val);
              }
              bigmatrix_set(&A, &val, i, j);
            }
          }

          // force the zero pivot (rectangular matrices only: on a square
          // matrix it could cost rank)
          bool force_zero = !square && trial % 2 == 1;
          if (force_zero) {
            bn_set_u64(&val, 0);
            bigmatrix_set(&A, &val, m - 1, n - 1);
          }

          // Delta = GCD of all m x m minors; if it is 0 the matrix is rank
          // deficient, so draw a new one
          test_module_determinant(&delta, &A);
          if (bn_is_zero(&delta)) {
            for (u64 i = 0; i < m; i++) {
              for (u64 j = 0; j < n; j++) {
                bn_gen_random(&val, bits);
                bn_gen_random(&bit, 1);
                if (!bn_is_zero(&bit)) {
                  bn_neg(&val, &val);
                }
                bigmatrix_set(&A, &val, i, j);
              }
            }
            if (force_zero) {
              bn_set_u64(&val, 0);
              bigmatrix_set(&A, &val, m - 1, n - 1);
            }
            test_module_determinant(&delta, &A);
          }
          assert(!bn_is_zero(&delta));

          bigmatrix_hermite(&H, &A);  // 2.4.4 oracle
          assert(H.r_size == m && H.c_size == m);
          assert(bigmatrix_hnf_check_structure(&H));

          bn_set_u64(&mult, mults[trial]);
          bn_mul(&D, &mult, &delta);

          bigmatrix_hermite_mod_d(&W1, &A, &D);  // 2.4.8 fast path
          assert(bigmatrix_equal(&W1, &H));
          assert(bigmatrix_hnf_check_structure(&W1));

          bigmatrix_hermite_gcd(&W2, &A);  // 2.4.5 exact
          assert(bigmatrix_equal(&W2, &H));
          assert(bigmatrix_hnf_check_structure(&W2));

          bigmatrix_free(&A);
          bigmatrix_free(&H);
          bigmatrix_free(&W1);
          bigmatrix_free(&W2);
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
        bigmatrix A, H, W1, W2;
        bigmatrix_init(&A, m, n);
        bigmatrix_init(&H, m, m);
        bigmatrix_init(&W1, m, m);
        bigmatrix_init(&W2, m, m);

        for (u64 i = 0; i < m; i++) {
          for (u64 j = 0; j < n; j++) {
            bn_gen_random(&val, bits);
            bn_gen_random(&bit, 1);
            if (!bn_is_zero(&bit)) {
              bn_neg(&val, &val);
            }
            bigmatrix_set(&A, &val, i, j);
          }
        }

        // square: Delta = |det A|
        bigmatrix_det(&delta, &A);
        if (bn_is_zero(&delta)) {
          for (u64 i = 0; i < m; i++) {
            bn_gen_random(&val, bits);
            bigmatrix_set(&A, &val, i, i);
          }
          bigmatrix_det(&delta, &A);
        }
        assert(!bn_is_zero(&delta));
        delta.is_neg = false;

        bigmatrix_hermite(&H, &A);
        assert(bigmatrix_hnf_check_structure(&H));

        bn_set_u64(&mult, mults[trial]);
        bn_mul(&D, &mult, &delta);

        bigmatrix_hermite_mod_d(&W1, &A, &D);
        assert(bigmatrix_equal(&W1, &H));
        assert(bigmatrix_hnf_check_structure(&W1));

        bigmatrix_hermite_gcd(&W2, &A);
        assert(bigmatrix_equal(&W2, &H));
        assert(bigmatrix_hnf_check_structure(&W2));

        bigmatrix_free(&A);
        bigmatrix_free(&H);
        bigmatrix_free(&W1);
        bigmatrix_free(&W2);
      }
      printf("[PASS] hermite_mod_d 10x10 300-bit\n");
    }
  }

  bn_free_multi(&D, &val, &bit, &delta, &mult, NULL);
}

int main()
{
  printf("--- Starting BigMatrix Test Suite ---\n");

  // 1. Test Initialization and Free
  bigmatrix A, B, R;
  bigmatrix_init(&A, 2, 2);
  bigmatrix_init(&B, 2, 2);
  bigmatrix_init(&R, 2, 2);

  assert(A.data != NULL);
  assert(A.r_size == 2 && A.c_size == 2);
  printf("[PASS] Initialization\n");

  // 2. Setup Data: A = [[1, 2], [3, 4]]
  bignum val;
  bn_init(&val);

  bn_set_u64(&val, 1);
  bigmatrix_set(&A, &val, 0, 0);
  bn_set_u64(&val, 2);
  bigmatrix_set(&A, &val, 0, 1);
  bn_set_u64(&val, 3);
  bigmatrix_set(&A, &val, 1, 0);
  bn_set_u64(&val, 4);
  bigmatrix_set(&A, &val, 1, 1);
  // 3. Test bigmatrix_get
  bignum check;
  bn_init(&check);
  bigmatrix_get(&check, &A, 0, 1);  // Should be 2
  // Replace with your library's comparison function, e.g., bn_cmp_u64
  assert(bn_is_eq_i64(&check, 2));
  printf("[PASS] Get/Set\n");

  // 4. Test bigmatrix_copy
  bigmatrix_copy(&B, &A);
  bigmatrix_get(&check, &B, 1, 1);
  assert(bn_is_eq_i64(&check, 4));
  printf("[PASS] Copy\n");

  // 5. Test bigmatrix_add: R = A + B (Since B=A, R should be [[2, 4], [6, 8]])
  bigmatrix_add(&R, &A, &B);
  bigmatrix_get(&check, &R, 1, 0);
  assert(bn_is_eq_i64(&check, 6));
  printf("[PASS] Addition\n");

  // 6. Test bigmatrix_mul: R = A * Identity
  // Setup B as Identity: [[1, 0], [0, 1]]
  bn_set_u64(&val, 1);
  bigmatrix_set(&B, &val, 0, 0);
  bn_set_u64(&val, 0);
  bigmatrix_set(&B, &val, 0, 1);
  bn_set_u64(&val, 0);
  bigmatrix_set(&B, &val, 1, 0);
  bn_set_u64(&val, 1);
  bigmatrix_set(&B, &val, 1, 1);

  bigmatrix_mul(&R, &A, &B);  // Result should be A
  bigmatrix_get(&check, &R, 1, 1);
  assert(bn_is_eq_i64(&check, 4));
  printf("[PASS] Multiplication (Identity)\n");

  // 7. Test bigmatrix_mul: Rectangular Multiplication
  // A (2x3) * B (3x2) = R (2x2)
  bigmatrix RectA, RectB, RectR;
  bigmatrix_init(&RectA, 2, 3);
  bigmatrix_init(&RectB, 3, 2);
  bigmatrix_init(&RectR, 2, 2);

  // Populate with 1s
  bn_set_u64(&val, 1);
  for (u64 i = 0; i < 6; i++) {
    bn_copy(&RectA.data[i], &val);
    bn_copy(&RectB.data[i], &val);
  }

  bigmatrix_mul(&RectR, &RectA, &RectB);
  // Each cell in a 2x3 * 3x2 of all 1s should be 3
  bigmatrix_get(&check, &RectR, 0, 0);
  assert(bn_is_eq_i64(&check, 3));
  printf("[PASS] Multiplication (Rectangular)\n");

  bigmatrix_print(&R);

  bigmatrix_det(&val, &A);
  bn_println(&val);

  bigmatrix_hadamard(&val, &A);
  bn_println(&val);

  // 7b. Test the structural check predicates
  {
    // 2x2 zero matrix: square, zero, diagonal, both triangular,
    // symmetric; not the identity
    bigmatrix Z;
    bigmatrix_init(&Z, 2, 2);
    assert(bigmatrix_is_square(&Z));
    assert(bigmatrix_is_zero(&Z));
    assert(!bigmatrix_is_identity(&Z));
    assert(bigmatrix_is_diagonal(&Z));
    assert(bigmatrix_is_upper_triangular(&Z));
    assert(bigmatrix_is_lower_triangular(&Z));
    assert(bigmatrix_is_symmetric(&Z));
    bigmatrix_free(&Z);
    printf("[PASS] Checks (zero matrix)\n");

    // 3x3 identity
    bigmatrix I3;
    bigmatrix_init(&I3, 3, 3);
    test_set_i64(&I3, 0, 0, 1);
    test_set_i64(&I3, 1, 1, 1);
    test_set_i64(&I3, 2, 2, 1);
    assert(bigmatrix_is_square(&I3));
    assert(bigmatrix_is_identity(&I3));
    assert(!bigmatrix_is_zero(&I3));
    assert(bigmatrix_is_diagonal(&I3));
    assert(bigmatrix_is_upper_triangular(&I3));
    assert(bigmatrix_is_lower_triangular(&I3));
    assert(bigmatrix_is_symmetric(&I3));
    bigmatrix_free(&I3);
    printf("[PASS] Checks (identity)\n");

    // -1 on the diagonal must not count as 1 (bn_is_eq_i64(x, 1) would
    // match it; bn_is_one must not)
    bigmatrix N3;
    bigmatrix_init(&N3, 3, 3);
    test_set_i64(&N3, 0, 0, -1);
    test_set_i64(&N3, 1, 1, 1);
    test_set_i64(&N3, 2, 2, 1);
    assert(!bigmatrix_is_identity(&N3));
    assert(bigmatrix_is_diagonal(&N3));
    assert(bigmatrix_is_symmetric(&N3));
    bigmatrix_free(&N3);
    printf("[PASS] Checks (negative diagonal)\n");

    // 3x3 upper triangular: not diagonal, not lower, not symmetric
    bigmatrix U3;
    bigmatrix_init(&U3, 3, 3);
    test_set_i64(&U3, 0, 0, 1);
    test_set_i64(&U3, 0, 1, 2);
    test_set_i64(&U3, 0, 2, 3);
    test_set_i64(&U3, 1, 1, 4);
    test_set_i64(&U3, 1, 2, 5);
    test_set_i64(&U3, 2, 2, 6);
    assert(bigmatrix_is_upper_triangular(&U3));
    assert(!bigmatrix_is_lower_triangular(&U3));
    assert(!bigmatrix_is_diagonal(&U3));
    assert(!bigmatrix_is_symmetric(&U3));
    assert(!bigmatrix_is_identity(&U3));
    bigmatrix_free(&U3);
    printf("[PASS] Checks (upper triangular)\n");

    // 3x3 lower triangular: the mirror case
    bigmatrix L3;
    bigmatrix_init(&L3, 3, 3);
    test_set_i64(&L3, 0, 0, 1);
    test_set_i64(&L3, 1, 0, 2);
    test_set_i64(&L3, 1, 1, 4);
    test_set_i64(&L3, 2, 0, 3);
    test_set_i64(&L3, 2, 1, 5);
    test_set_i64(&L3, 2, 2, 6);
    assert(bigmatrix_is_lower_triangular(&L3));
    assert(!bigmatrix_is_upper_triangular(&L3));
    assert(!bigmatrix_is_diagonal(&L3));
    assert(!bigmatrix_is_symmetric(&L3));
    bigmatrix_free(&L3);
    printf("[PASS] Checks (lower triangular)\n");

    // 3x3 symmetric, then break one mirror pair
    bigmatrix S3;
    bigmatrix_init(&S3, 3, 3);
    test_set_i64(&S3, 0, 0, 1);
    test_set_i64(&S3, 0, 1, 2);
    test_set_i64(&S3, 0, 2, 3);
    test_set_i64(&S3, 1, 0, 2);
    test_set_i64(&S3, 1, 1, 4);
    test_set_i64(&S3, 1, 2, 5);
    test_set_i64(&S3, 2, 0, 3);
    test_set_i64(&S3, 2, 1, 5);
    test_set_i64(&S3, 2, 2, 6);
    assert(bigmatrix_is_symmetric(&S3));
    test_set_i64(&S3, 0, 2, 9);
    assert(!bigmatrix_is_symmetric(&S3));
    bigmatrix_free(&S3);
    printf("[PASS] Checks (symmetric)\n");

    // 2x3 zero rectangle: not square / not identity / not symmetric,
    // the rest vacuously true
    bigmatrix R23;
    bigmatrix_init(&R23, 2, 3);
    assert(!bigmatrix_is_square(&R23));
    assert(bigmatrix_is_zero(&R23));
    assert(!bigmatrix_is_identity(&R23));
    assert(bigmatrix_is_diagonal(&R23));
    assert(bigmatrix_is_upper_triangular(&R23));
    assert(bigmatrix_is_lower_triangular(&R23));
    assert(!bigmatrix_is_symmetric(&R23));
    // a single sub-diagonal entry breaks upper triangularity (and
    // diagonal/zero), but not lower triangularity
    test_set_i64(&R23, 1, 0, 7);
    assert(!bigmatrix_is_zero(&R23));
    assert(!bigmatrix_is_upper_triangular(&R23));
    assert(bigmatrix_is_lower_triangular(&R23));
    assert(!bigmatrix_is_diagonal(&R23));
    bigmatrix_free(&R23);
    printf("[PASS] Checks (rectangular)\n");

    // 3x4 identity-pattern: the dimension short-circuit must reject it
    // as identity/symmetric even though the scanned sub-regions are clean
    bigmatrix Q34;
    bigmatrix_init(&Q34, 3, 4);
    test_set_i64(&Q34, 0, 0, 1);
    test_set_i64(&Q34, 1, 1, 1);
    test_set_i64(&Q34, 2, 2, 1);
    assert(!bigmatrix_is_square(&Q34));
    assert(!bigmatrix_is_identity(&Q34));
    assert(!bigmatrix_is_symmetric(&Q34));
    assert(bigmatrix_is_diagonal(&Q34));
    assert(bigmatrix_is_upper_triangular(&Q34));
    assert(bigmatrix_is_lower_triangular(&Q34));
    bigmatrix_free(&Q34);
    printf("[PASS] Checks (non-square short-circuit)\n");

    // degenerate sizes: 0x0 is vacuously everything, 0x2 is not square
    // nor symmetric
    bigmatrix E00, E02;
    bigmatrix_init(&E00, 0, 0);
    assert(bigmatrix_is_square(&E00));
    assert(bigmatrix_is_zero(&E00));
    assert(bigmatrix_is_identity(&E00));
    assert(bigmatrix_is_diagonal(&E00));
    assert(bigmatrix_is_upper_triangular(&E00));
    assert(bigmatrix_is_lower_triangular(&E00));
    assert(bigmatrix_is_symmetric(&E00));
    bigmatrix_free(&E00);

    bigmatrix_init(&E02, 0, 2);
    assert(!bigmatrix_is_square(&E02));
    assert(bigmatrix_is_zero(&E02));
    assert(!bigmatrix_is_identity(&E02));
    assert(bigmatrix_is_diagonal(&E02));
    assert(bigmatrix_is_upper_triangular(&E02));
    assert(bigmatrix_is_lower_triangular(&E02));
    assert(!bigmatrix_is_symmetric(&E02));
    bigmatrix_free(&E02);
    printf("[PASS] Checks (degenerate sizes)\n");
  }

  // 8. Cleanup
  bn_free(&val);
  bn_free(&check);
  bigmatrix_free(&A);
  bigmatrix_free(&B);
  bigmatrix_free(&R);
  bigmatrix_free(&RectA);
  bigmatrix_free(&RectB);
  bigmatrix_free(&RectR);

  bigmatrix D, ADJ;
  bigmatrix_init(&D, 20, 20);
  bigmatrix_init(&ADJ, 20, 20);
  bignum tmp;
  bn_init(&tmp);

  // set matrix D
  for (i64 i = 0; i < 20; i++) {
    for (i64 j = 0; j < 20; j++) {
      bn_set_i64(&tmp, (i + j) % 20);
      bigmatrix_set(&D, &tmp, i, j);
    }
  }

  bigmatrix_print(&D);
  bigmatrix_det(&tmp, &D);
  bn_println(&tmp);

  bigmatrix_trace(&tmp, &D);
  bn_println(&tmp);

  bigpoly p;
  bigpoly_init(&p);

  bigmatrix_charpoly_adj(&p, &ADJ, &D);

  printf("Char poly: ");
  bigpoly_print(&p);
  printf("\n");

  printf("Adjoint: \n ");
  bigmatrix_print(&ADJ);
  printf("\n");

  printf("Hermite: \n ");
  bigmatrix_hermite(&ADJ, &D);
  assert(bigmatrix_hnf_check_structure(&ADJ));
  bigmatrix_print(&ADJ);
  printf("\n");

  printf("Hermite GCD: \n ");
  bigmatrix_hermite_gcd(&ADJ, &D);
  bigmatrix_print(&ADJ);
  printf("\n");

  test_hermite_mod_d();
  test_hnf_verify();
  test_hnf_mod_d_multilimb();

  bigpoly_free(&p);

  bigmatrix_free(&ADJ);

  ctx_rns ctx;
  rns_context_init(&ctx, RNS_PRIMES, 10);

  bigmatrix_det_rns(&tmp, &D, &ctx);

  bn_println(&tmp);

  rns_context_free(&ctx);

  bigmatrix_free(&D);

  bigmatrix E;
  bigmatrix_init(&E, 3, 4);
  bigvector v, res;
  bigvector_init(&v, 4);
  bigvector_init(&res, 3);

  bigmatrix_set(&E, &tmp, 0, 0);
  bn_add_u64(&tmp, &tmp, 1);
  bigmatrix_set(&E, &tmp, 1, 0);
  bn_add_u64(&tmp, &tmp, 1);
  bigmatrix_set(&E, &tmp, 1, 1);
  bn_add_u64(&tmp, &tmp, 1);
  bigmatrix_set(&E, &tmp, 2, 0);
  bn_add_u64(&tmp, &tmp, 1);
  bigmatrix_set(&E, &tmp, 2, 1);
  bn_add_u64(&tmp, &tmp, 1);
  bigmatrix_set(&E, &tmp, 2, 2);

  bigmatrix_print(&E);

  bigvector_set(&v, &tmp, 0);

  bigvector_set(&v, &tmp, 1);
  bigvector_set(&v, &tmp, 2);
  bigvector_set(&v, &tmp, 3);

  bigvector_print(&v);
  printf("\n");

  bigmatrix_mv(&res, &E, &v);

  printf("Matrix vector: ");
  bigvector_print(&res);
  printf("\n");

  u64* testmatrix = malloc(sizeof(u64) * 9);
  u64 values[] = {18446744073709551610, 1, 1, 1, 18446744073709551600, 1, 1, 1,
                  18446744073709551590};
  memcpy(testmatrix, values, sizeof(u64) * 9);

  matrix_u64 A_test;
  A_test.data = testmatrix;
  A_test.c_size = 3;
  A_test.r_size = 3;
  A_test.modulus = 100000007;

  u64 det = matrix_u64_det(&A_test);
  printf("u64 determinant: %llu \n", det);

  free(testmatrix);

  test_pascal_det(50);

  printf("[PASS] Cleanup / Free\n");
  printf("--- All Tests Passed! ---\n");

  u64 sizes[] = {2, 4, 6, 8, 16, 32, 64, 128, 256, 300, 512, 700, 994, 1024};
  int num_tests = sizeof(sizes) / sizeof(sizes[0]);

  // for (int i = 0; i < num_tests; i++) {
  //   run_hadamard_benchmark(sizes[i], 32);
  // }

  // for (int i = 0; i < 11; i++) {
  //   run_det_benchmark(sizes[i], 32);
  // }

  for (int i = 0; i < 7; i++) {
    run_hermite_benchmark(sizes[i], 10);
  }

  // HNF mod D (2.4.8) vs exact (2.4.5): on the smaller sizes both are timed
  // on the same random matrix, on the larger ones only the mod D version
  u64 hnf_sizes[] = {4, 8, 16, 32, 64, 128};
  int n_hnf = sizeof(hnf_sizes) / sizeof(hnf_sizes[0]);
  for (int i = 0; i < 8; i++) {
    run_hermite_mod_d_benchmark(sizes[i], 10, i < 4);
  }

  for (int i = 0; i < 5; i++) {
    run_mul_benchmark(sizes[i], 64);
  }

  bn_free(&tmp);

  bigmatrix_free(&E);
  bigvector_free(&v);
  bigvector_free(&res);

  return 0;
}
