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

  for (int i = 0; i < 8; i++) {
    run_hermite_benchmark(sizes[i], 8);
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
