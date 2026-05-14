#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "../include/bignum.h"

// Helper to calculate elapsed time in seconds
double get_elapsed_time(struct timespec start, struct timespec end)
{
  return (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
}

// Helper to generate a random bignum/mpz of a certain bit length
void generate_random_inputs(bignum* bn_a, bignum* bn_b, bignum* bn_m,
                            mpz_t mpz_a, mpz_t mpz_b, mpz_t mpz_m, int bits)
{
  // Generate for GMP first
  gmp_randstate_t state;
  gmp_randinit_default(state);

  mpz_urandomb(mpz_a, state, bits);
  mpz_urandomb(mpz_b, state, bits);
  mpz_urandomb(mpz_m, state, bits);
  // Ensure modulus is odd for Montgomery
  mpz_setbit(mpz_m, 0);

  // Convert GMP strings to your bignum format
  char* s_a = mpz_get_str(NULL, 10, mpz_a);
  char* s_b = mpz_get_str(NULL, 10, mpz_b);
  char* s_m = mpz_get_str(NULL, 10, mpz_m);

  bn_init_val(bn_a, s_a);
  bn_init_val(bn_b, s_b);
  bn_init_val(bn_m, s_m);

  free(s_a);
  free(s_b);
  free(s_m);
  gmp_randclear(state);
}

void benchmark_mod_exp(int bits, int iterations)
{
  printf("\n--- ModExp Benchmark (%d bits, %d iterations) ---\n", bits,
         iterations);

  bignum bn_a, bn_b, bn_m, bn_res;
  mpz_t mpz_a, mpz_b, mpz_m, mpz_res;

  bn_init(&bn_a);
  bn_init(&bn_b);
  bn_init(&bn_m);
  bn_init(&bn_res);
  mpz_inits(mpz_a, mpz_b, mpz_m, mpz_res, NULL);

  generate_random_inputs(&bn_a, &bn_b, &bn_m, mpz_a, mpz_b, mpz_m, bits);

  struct timespec start, end;
  double total_custom = 0, total_gmp = 0;

  for (int i = 0; i < iterations; i++) {
    // Benchmark Custom
    clock_gettime(CLOCK_MONOTONIC, &start);
    bn_mod_exp(&bn_res, &bn_a, &bn_b, &bn_m);
    clock_gettime(CLOCK_MONOTONIC, &end);
    total_custom += get_elapsed_time(start, end);

    // Benchmark GMP
    clock_gettime(CLOCK_MONOTONIC, &start);
    mpz_powm(mpz_res, mpz_a, mpz_b, mpz_m);
    clock_gettime(CLOCK_MONOTONIC, &end);
    total_gmp += get_elapsed_time(start, end);
  }

  double avg_custom = total_custom / iterations;
  double avg_gmp = total_gmp / iterations;

  printf("Custom Library: %.6f s (avg)\n", avg_custom);
  printf("GMP Library:    %.6f s (avg)\n", avg_gmp);
  printf("Ratio:          %.2fx slower than GMP\n", avg_custom / avg_gmp);

  bn_free(&bn_a);
  bn_free(&bn_b);
  bn_free(&bn_m);
  bn_free(&bn_res);
  mpz_clears(mpz_a, mpz_b, mpz_m, mpz_res, NULL);
}

void benchmark_pow(int bits, unsigned long exp, int iterations)
{
  printf("\n--- Pow Benchmark (Base: %d bits, Exp: %lu, %d iterations) ---\n",
         bits, exp, iterations);

  bignum bn_a, bn_b, bn_res;
  mpz_t mpz_a, mpz_res;

  bn_init(&bn_a);
  bn_init(&bn_b);
  bn_init(&bn_res);
  mpz_init(mpz_a);
  mpz_init(mpz_res);

  // Setup base
  gmp_randstate_t state;
  gmp_randinit_default(state);
  mpz_urandomb(mpz_a, state, bits);
  char* s_a = mpz_get_str(NULL, 10, mpz_a);
  bn_init_val(&bn_a, s_a);
  free(s_a);

  // Setup exponent
  bn_set_u64(&bn_b, exp);

  struct timespec start, end;
  double total_custom = 0, total_gmp = 0;

  for (int i = 0; i < iterations; i++) {
    // Benchmark Custom
    clock_gettime(CLOCK_MONOTONIC, &start);
    bn_pow(&bn_res, &bn_a, &bn_b);
    clock_gettime(CLOCK_MONOTONIC, &end);
    total_custom += get_elapsed_time(start, end);

    // Benchmark GMP
    clock_gettime(CLOCK_MONOTONIC, &start);
    mpz_pow_ui(mpz_res, mpz_a, exp);
    clock_gettime(CLOCK_MONOTONIC, &end);
    total_gmp += get_elapsed_time(start, end);
  }

  double avg_custom = total_custom / iterations;
  double avg_gmp = total_gmp / iterations;

  printf("Custom Library: %.6f s (avg)\n", avg_custom);
  printf("GMP Library:    %.6f s (avg)\n", avg_gmp);
  printf("Ratio:          %.2fx slower than GMP\n", avg_custom / avg_gmp);

  bn_free(&bn_a);
  bn_free(&bn_b);
  bn_free(&bn_res);
  mpz_clears(mpz_a, mpz_res, NULL);
  gmp_randclear(state);
}

void benchmark_mul(int bits, unsigned long exp, int iterations)
{
  printf("\n--- Mul Benchmark (Base: %d bits, Exp: %lu, %d iterations) ---\n",
         bits, exp, iterations);

  bignum bn_a, bn_b, bn_res;
  mpz_t mpz_a, mpz_res;

  bn_init(&bn_a);
  bn_init(&bn_b);
  bn_init(&bn_res);
  mpz_init(mpz_a);
  mpz_init(mpz_res);

  // Setup base
  gmp_randstate_t state;
  gmp_randinit_default(state);
  mpz_urandomb(mpz_a, state, bits);
  char* s_a = mpz_get_str(NULL, 10, mpz_a);
  bn_init_val(&bn_a, s_a);
  free(s_a);

  // Setup exponent
  bn_set_u64(&bn_b, exp);

  struct timespec start, end;
  double total_custom = 0, total_gmp = 0;

  for (int i = 0; i < iterations; i++) {
    // Benchmark Custom
    clock_gettime(CLOCK_MONOTONIC, &start);
    bn_mul(&bn_res, &bn_a, &bn_a);
    clock_gettime(CLOCK_MONOTONIC, &end);
    total_custom += get_elapsed_time(start, end);

    // Benchmark GMP
    clock_gettime(CLOCK_MONOTONIC, &start);
    mpz_mul(mpz_res, mpz_a, mpz_a);
    clock_gettime(CLOCK_MONOTONIC, &end);
    total_gmp += get_elapsed_time(start, end);
  }

  double avg_custom = total_custom / iterations;
  double avg_gmp = total_gmp / iterations;

  printf("Custom Library: %.6f s (avg)\n", avg_custom);
  printf("GMP Library:    %.6f s (avg)\n", avg_gmp);
  printf("Ratio:          %.2fx slower than GMP\n", avg_custom / avg_gmp);

  bn_free(&bn_a);
  bn_free(&bn_b);
  bn_free(&bn_res);
  mpz_clears(mpz_a, mpz_res, NULL);
  gmp_randclear(state);
}

int main()
{
  int iterations = 10;

  // Test ModExp with increasing bit sizes
  int sizes[] = {512, 1024, 2048, 4096, 5000};
  for (int i = 0; i < 5; i++) {
    benchmark_mod_exp(sizes[i], iterations);
  }

  // Test Pow (Note: pow grows extremely fast, keep exponent reasonable)
  benchmark_pow(256, 1000, iterations);
  benchmark_pow(512, 1000, iterations);
  benchmark_pow(1024, 1000, iterations);
  benchmark_pow(2048, 1000, iterations);

  benchmark_mul(2048, 1000, iterations);
  benchmark_mul(20480, 5000, iterations);
  benchmark_mul(204800, 10000, iterations);
  benchmark_mul(2048000, 100000, iterations);

  return 0;
}
