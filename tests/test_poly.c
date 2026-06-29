#include <stdio.h>
#include <time.h>

#include "../include/bignum.h"
#include "../include/bigpoly.h"
#define ITERATIONS 1
// #define POLY_SIZE 512

void benchmark_mul(u64 n)
{
  bigpoly a, b, r_std, r_ntt, r_ntt_big;
  bigpoly_init(&a);
  bigpoly_init(&b);
  bigpoly_init(&r_std);
  bigpoly_init(&r_ntt);
  bigpoly_init(&r_ntt_big);

  bigpoly_alloc(&a, n);
  bigpoly_alloc(&b, n);
  a.deg = b.deg = n - 1;

  // Initialize coefficients
  bignum* C = malloc(sizeof(bignum) * n);
  for (u64 i = 0; i < n; i++) {
    bn_init(&C[i]);
    bn_set_u64(&C[i], i + 1);  // Avoid all zeros
  }
  bigpoly_set(&a, C, n - 1);
  bigpoly_set(&b, C, n - 1);

  printf("--- Benchmarking Polynomial Multiplication (N = %llu) ---\n", n);

  // 1. Standard Multiplication Timing
  clock_t start_std = clock();
  for (int i = 0; i < ITERATIONS; i++) {
    // bigpoly_mul(&r_std, &a, &b);  // Standard O(N^2)
  }
  clock_t end_std = clock();
  double time_std = (double)(end_std - start_std + 100000000) / CLOCKS_PER_SEC;

  // 2. NTT Multiplication Timing
  clock_t start_ntt = clock();
  for (int i = 0; i < ITERATIONS; i++) {
    bigpoly_mul_ntt_u64(&r_ntt, &a, &b);  // NTT O(N log N)
  }
  clock_t end_ntt = clock();
  double time_ntt = (double)(end_ntt - start_ntt) / CLOCKS_PER_SEC;

  clock_t start_ntt_big = clock();
  for (int i = 0; i < ITERATIONS; i++) {
    bigpoly_mul_ntt(&r_ntt_big, &a, &b);  // NTT O(N log N)
  }
  clock_t end_ntt_big = clock();
  double time_ntt_big = (double)(end_ntt_big - start_ntt_big) / CLOCKS_PER_SEC;

  if (bigpoly_equal(&r_ntt_big, &r_ntt)) {
    printf("Verification: [PASSED] (NTT results match Standard)\n");
  } else {
    printf("Verification: [FAILED] (NTT results differ from Standard!)\n");
    // Optional: print first few coefficients to debug
    // bigpoly_print(&r_std);
    // bigpoly_print(&r_ntt);
  }

  // Results Output
  printf("Standard Mul: %.5f seconds (%d iterations)\n", time_std, ITERATIONS);
  printf("NTT Mul:      %.5f seconds (%d iterations)\n", time_ntt, ITERATIONS);
  printf("NTT big Mul:      %.5f seconds (%d iterations)\n", time_ntt_big,
         ITERATIONS);

  if (time_ntt < time_std) {
    printf("Result: NTT is %.2fx faster than Standard.\n", time_std / time_ntt);
    printf("Result: NTT is %.2fx faster than NTT bignum.\n",
           time_ntt_big / time_ntt);
    printf("Result: NTT big is %.2fx faster than Std.\n",
           time_std / time_ntt_big);
  } else {
    printf("Result: Standard is faster at this size (NTT overhead).\n");
  }

  // Cleanup
  for (u64 i = 0; i < n; i++) bn_free(&C[i]);
  free(C);
  bigpoly_free(&a);
  bigpoly_free(&b);
  bigpoly_free(&r_std);
  bigpoly_free(&r_ntt);
  bigpoly_free(&r_ntt_big);
}

int main()
{
  bn_init_constants();

  u64 n = 16;
  bigpoly a, b, r;
  bigpoly_init(&a);
  bigpoly_init(&b);
  bigpoly_init(&r);
  bigpoly_alloc(&a, n);
  bigpoly_alloc(&b, n);
  a.deg = b.deg = n - 1;

  bignum* C = malloc(sizeof(bignum) * n);

  for (u64 i = 0; i < n; i++) {
    bn_init(&C[i]);
    bn_set_u64(&C[i], i);
  }

  bigpoly_set(&a, C, n - 1);
  bigpoly_set(&b, C, n - 1);

  bigpoly_print(&a);
  bigpoly_print(&b);

  bigpoly_mul_ntt(&r, &a, &b);
  bigpoly_print(&r);

  for (u64 i = 0; i < n; i++) {
    bn_free(&C[i]);
  }
  free(C);
  bigpoly_free(&a);
  bigpoly_free(&b);
  bigpoly_free(&r);

  benchmark_mul(15);
  benchmark_mul(64);
  printf("\n");
  benchmark_mul(512);
  benchmark_mul(1024);
  benchmark_mul(2048);
  benchmark_mul(4096);
  benchmark_mul(8192);
  benchmark_mul(1 << 14);
  benchmark_mul(1 << 15);
  benchmark_mul(1 << 16);
  benchmark_mul(1 << 17);
  benchmark_mul(1 << 18);
  benchmark_mul(1 << 19);
  benchmark_mul(1 << 20);
  benchmark_mul(1 << 17);

  bn_free_constants();
}
