#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

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

  printf("Benchmarking %llu x %llu (%llu-bit entries) Hermite normal form... ",
         size, size, bits);
  fflush(stdout);

  double start = get_time();
  bigmatrix_smith(&H, &M);
  double end = get_time();

  bigmatrix_print_full(&H);

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

static bool test_matrix_equal(const bigmatrix* A, const bigmatrix* B)
{
  if (A->r_size != B->r_size || A->c_size != B->c_size) return false;

  for (u64 i = 0; i < A->r_size * A->c_size; i++) {
    if (bn_cmp(&A->data[i], &B->data[i]) != 0) return false;
  }
  return true;
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
    assert(test_matrix_equal(&W, &H));

    // cross-check against the exact implementation
    bigmatrix V;
    bigmatrix_init(&V, 2, 2);
    bigmatrix_hermite(&V, &A);
    assert(test_matrix_equal(&W, &V));
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
    assert(test_matrix_equal(&W, &H));
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
        assert(test_matrix_equal(&W1, &H));

        // D = 3 * |det B|, another valid multiple
        bn_mul(&D, &three, &det);
        bigmatrix_hermite_mod_d(&W2, &A, &D);
        assert(test_matrix_equal(&W2, &H));

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
    assert(test_matrix_equal(&W, &H));
    printf("[PASS] hermite_mod_d m > n fallback\n");

    bigmatrix_free(&A);
    bigmatrix_free(&W);
    bigmatrix_free(&H);
  }

  bn_free_multi(&D, &val, NULL);
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
  bigmatrix_print(&ADJ);
  printf("\n");

  printf("Hermite GCD: \n ");
  bigmatrix_smith(&ADJ, &D);
  bigmatrix_print(&ADJ);
  printf("\n");

  test_hermite_mod_d();

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

  for (int i = 0; i < 10; i++) {
    run_hermite_benchmark(sizes[i], 10);
  }

  // HNF mod D (2.4.8) vs exact (2.4.5): on the smaller sizes both are timed
  // on the same random matrix, on the larger ones only the mod D version
  u64 hnf_sizes[] = {4, 8, 16, 32, 64, 128};
  int n_hnf = sizeof(hnf_sizes) / sizeof(hnf_sizes[0]);
  for (int i = 0; i < 10; i++) {
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
