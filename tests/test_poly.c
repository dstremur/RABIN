#include <assert.h>
#include <stdio.h>
#include <time.h>

#include "../include/rabin.h"
#include "../include/rpol.h"
#define ITERATIONS 5
// #define POLY_SIZE 512

static int failures = 0;

#define CHECK(cond, msg)                                       \
  do {                                                         \
    if (!(cond)) {                                             \
      printf("FAILED: %s (%s:%d)\n", msg, __FILE__, __LINE__); \
      failures++;                                              \
    }                                                          \
  } while (0)

// Self-test formerly shipped as rpol_test() in the library itself,
// now with real pass/fail checks.
static void test_basic_ops(void)
{
  rpol_t p = {0}, q = {0}, r = {0}, expect = {0};
  rpol_init(&p);
  rpol_init(&q);
  rpol_init(&r);
  rpol_init(&expect);

  // p(x) = 2x^2 + 3x + 1, q(x) = 4x + 5
  i64 p_vals[] = {1, 3, 2};
  i64 q_vals[] = {5, 4};
  CHECK(rpol_set_i64(&p, p_vals, 2) == RABIN_SUCCESS, "rpol_set_i64 p");
  CHECK(rpol_set_i64(&q, q_vals, 1) == RABIN_SUCCESS, "rpol_set_i64 q");
  CHECK(rpol_equal(&p, &p), "rpol_equal reflexive");

  // p + q = 2x^2 + 7x + 6
  i64 sum_vals[] = {6, 7, 2};
  CHECK(rpol_set_i64(&expect, sum_vals, 2) == RABIN_SUCCESS,
        "rpol_set_i64 sum");
  CHECK(rpol_add(&r, &p, &q) == RABIN_SUCCESS, "rpol_add");
  CHECK(rpol_equal(&r, &expect), "p + q == 2x^2 + 7x + 6");

  // p - q = 2x^2 - x - 4
  rpol_t minus = {0};
  i64 exp_sub[] = {-4, -1, 2};
  rpol_init(&minus);
  CHECK(rpol_set_i64(&minus, exp_sub, 2) == RABIN_SUCCESS, "rpol_set_i64 p-q");
  CHECK(rpol_sub(&r, &p, &q) == RABIN_SUCCESS, "rpol_sub");
  CHECK(rpol_equal(&r, &minus), "p - q == 2x^2 - x - 4");

  // p * q = 8x^3 + 22x^2 + 19x + 5 (schoolbook)
  i64 prod_vals[] = {5, 19, 22, 8};
  CHECK(rpol_mul_school(&r, &p, &q) == RABIN_SUCCESS, "rpol_mul_school");
  CHECK(rpol_set_i64(&expect, prod_vals, 3) == RABIN_SUCCESS,
        "rpol_set_i64 product");
  CHECK(rpol_equal(&r, &expect), "schoolbook p*q == 8x^3 + 22x^2 + 19x + 5");

  // NTT-u64 path must agree with schoolbook
  rpol_t r_ntt = {0};
  rpol_init(&r_ntt);
  CHECK(rpol_mul_ntt_u64(&r_ntt, &p, &q) == RABIN_SUCCESS, "rpol_mul_ntt_u64");
  CHECK(rpol_equal(&r_ntt, &expect), "NTT p*q matches schoolbook");

  rpol_clear(&minus);
  rpol_clear(&r_ntt);
  rpol_clear(&p);
  rpol_clear(&q);
  rpol_clear(&r);
  rpol_clear(&expect);
}

void benchmark_mul(u64 n)
{
  rpol_t a = {0}, b = {0}, r_std = {0}, r_ntt = {0}, r_ntt_big = {0};
  rpol_init(&a);
  rpol_init(&b);
  rpol_init(&r_std);
  rpol_init(&r_ntt);
  rpol_init(&r_ntt_big);

  rpol_alloc(&a, n);
  rpol_alloc(&b, n);
  a.deg = b.deg = n - 1;

  // Initialize coefficients
  rz_t* C = malloc(sizeof(rz_t) * n);
  for (u64 i = 0; i < n; i++) {
    rz_init(&C[i]);
    rz_set_u64(&C[i], 2 * i + 1);  // Avoid all zeros
  }
  rpol_set(&a, C, n - 1);
  rpol_set(&b, C, n - 1);

  printf("--- Benchmarking Polynomial Multiplication (N = %llu) ---\n", n);

  // 1. Standard Multiplication Timing
  clock_t start_std = clock();
  for (int i = 0; i < ITERATIONS; i++) {
    rpol_mul(&r_std, &a, &b);  // Standard O(N^2)
  }
  clock_t end_std = clock();
  double time_std = (double)(end_std - start_std) / CLOCKS_PER_SEC;

  // 2. NTT Multiplication Timing
  clock_t start_ntt = clock();
  for (int i = 0; i < ITERATIONS; i++) {
    rpol_mul_ntt_u64(&r_ntt, &a, &b);  // NTT O(N log N)
  }
  clock_t end_ntt = clock();
  double time_ntt = (double)(end_ntt - start_ntt) / CLOCKS_PER_SEC;

  clock_t start_ntt_big = clock();
  for (int i = 0; i < ITERATIONS; i++) {
    // rpol_mul_ntt(&r_ntt_big, &a, &b);  // NTT O(N log N)
  }
  clock_t end_ntt_big = clock();
  double time_ntt_big = (double)(end_ntt_big - start_ntt_big) / CLOCKS_PER_SEC;

  if (rpol_equal(&r_std, &r_ntt)) {
    printf("Verification: [PASSED] (NTT results match Standard)\n");
  } else {
    printf("Verification: [FAILED] (NTT results differ from Standard!)\n");
    // Optional: print first few coefficients to debug
    // rpol_print(&r_std);
    // rpol_print(&r_ntt);
  }

  // Results Output
  printf("Standard Mul: %.5f seconds (%d iterations)\n", time_std / ITERATIONS,
         ITERATIONS);
  printf("NTT Mul:      %.5f seconds (%d iterations)\n", time_ntt / ITERATIONS,
         ITERATIONS);
  printf("NTT big Mul:      %.5f seconds (%d iterations)\n",
         time_ntt_big / ITERATIONS, ITERATIONS);
  if ((2 * n * n * n) < 5 * (1ULL << 55)) {
    printf("STILL FITS\n");
  } else {
    printf("TOO LARGE\n");
  }
  if (time_ntt < time_std) {
    printf("Result: NTT is %.2fx faster than Standard.\n", time_std / time_ntt);
    printf("Result: NTT is %.2fx faster than NTT rz_t.\n",
           time_ntt_big / time_ntt);
    printf("Result: NTT big is %.2fx faster than Std.\n",
           time_std / time_ntt_big);
  } else {
    printf("Result: Standard is faster at this size (NTT overhead).\n");
  }

  // Cleanup
  for (u64 i = 0; i < n; i++) rz_clear(&C[i]);
  free(C);
  rpol_clear(&a);
  rpol_clear(&b);
  rpol_clear(&r_std);
  rpol_clear(&r_ntt);
  rpol_clear(&r_ntt_big);
}

int main()
{
  rz_init_constants();

  test_basic_ops();

  u64 n = 16;
  rpol_t a = {0}, b = {0}, r = {0};
  rpol_init(&a);
  rpol_init(&b);
  rpol_init(&r);
  rpol_alloc(&a, n);
  rpol_alloc(&b, n);
  a.deg = b.deg = n - 1;

  rz_t* C = malloc(sizeof(rz_t) * n);

  for (u64 i = 0; i < n; i++) {
    rz_init(&C[i]);
    rz_set_u64(&C[i], i);
  }

  rpol_set(&a, C, n - 1);
  rpol_set(&b, C, n - 1);

  rpol_print(&a);
  rpol_print(&b);

  rpol_mul_ntt(&r, &a, &b);
  rpol_print(&r);

  for (u64 i = 0; i < n; i++) {
    rz_clear(&C[i]);
  }
  free(C);
  rpol_clear(&a);
  rpol_clear(&b);
  rpol_clear(&r);

  benchmark_mul(15);
  benchmark_mul(64);
  printf("\n");
  benchmark_mul(384);
  benchmark_mul(512);
  benchmark_mul(1024);
  benchmark_mul(2048);
  benchmark_mul(3000);
  benchmark_mul(4000);
  benchmark_mul(4096);
  benchmark_mul(8192);
  benchmark_mul(10000);
  benchmark_mul(1 << 14);
  benchmark_mul(1 << 15);
  benchmark_mul(1 << 16);
  /*  benchmark_mul(1 << 17);
    benchmark_mul(1 << 18);
    benchmark_mul(1 << 19);
   */
  rz_clear_constants();

  if (failures == 0) {
    printf("test_poly: ALL BASIC CHECKS PASSED\n");
  } else {
    printf("test_poly: %d CHECK(S) FAILED\n", failures);
  }
  return failures == 0 ? 0 : 1;
}
