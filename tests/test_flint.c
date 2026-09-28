#include <flint/flint.h>
#include <flint/fmpz.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "../include/rabin.h"

// Helper to calculate elapsed time in seconds
double get_elapsed_time(struct timespec start, struct timespec end)
{
  return (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
}

// Helper to generate random inputs
void generate_random_inputs(rz_t* rz_a, rz_t* rz_b, rz_t* rz_m, fmpz_t f_a,
                            fmpz_t f_b, fmpz_t f_m, int bits)
{
  flint_rand_t state;
  flint_rand_init(state);

  fmpz_randbits(f_a, state, bits);
  fmpz_randbits(f_b, state, bits);
  fmpz_randbits(f_m, state, bits);

  // Ensure modulus is odd
  fmpz_setbit(f_m, 0);

  // Convert to strings for your rz_t
  char* s_a = fmpz_get_str(NULL, 10, f_a);
  char* s_b = fmpz_get_str(NULL, 10, f_b);
  char* s_m = fmpz_get_str(NULL, 10, f_m);

  rz_init_val(rz_a, s_a);
  rz_init_val(rz_b, s_b);
  rz_init_val(rz_m, s_m);

  flint_free(s_a);
  flint_free(s_b);
  flint_free(s_m);

  flint_rand_clear(state);
}

void benchmark_mod_exp(int bits, int iterations)
{
  printf("\n--- ModExp Benchmark (%d bits, %d iterations) ---\n", bits,
         iterations);

  rz_t rz_a, rz_b, rz_m, rz_res;
  fmpz_t f_a, f_b, f_m, f_res;

  rz_init(&rz_a);
  rz_init(&rz_b);
  rz_init(&rz_m);
  rz_init(&rz_res);
  fmpz_init(f_a);
  fmpz_init(f_b);
  fmpz_init(f_m);
  fmpz_init(f_res);

  generate_random_inputs(&rz_a, &rz_b, &rz_m, f_a, f_b, f_m, bits);

  struct timespec start, end;
  double total_custom = 0, total_flint = 0;

  for (int i = 0; i < iterations; i++) {
    // Custom
    clock_gettime(CLOCK_MONOTONIC, &start);
    rz_mod_exp(&rz_res, &rz_a, &rz_b, &rz_m);
    clock_gettime(CLOCK_MONOTONIC, &end);
    total_custom += get_elapsed_time(start, end);

    // FLINT
    clock_gettime(CLOCK_MONOTONIC, &start);
    fmpz_powm(f_res, f_a, f_b, f_m);
    clock_gettime(CLOCK_MONOTONIC, &end);
    total_flint += get_elapsed_time(start, end);
  }

  printf("Custom Library: %.6f s (avg)\n", total_custom / iterations);
  printf("FLINT Library:  %.6f s (avg)\n", total_flint / iterations);
  printf("Ratio:          %.2fx slower than FLINT\n",
         (total_custom / iterations) / (total_flint / iterations));

  rz_clear(&rz_a);
  rz_clear(&rz_b);
  rz_clear(&rz_m);
  rz_clear(&rz_res);
  fmpz_clear(f_a);
  fmpz_clear(f_b);
  fmpz_clear(f_m);
  fmpz_clear(f_res);
}

void benchmark_pow(int bits, unsigned long exp, int iterations)
{
  printf("\n--- Pow Benchmark (Base: %d bits, Exp: %lu, %d iterations) ---\n",
         bits, exp, iterations);

  rz_t rz_a, rz_b, rz_res;
  fmpz_t f_a, f_res;

  rz_init(&rz_a);
  rz_init(&rz_b);
  rz_init(&rz_res);
  fmpz_init(f_a);
  fmpz_init(f_res);

  flint_rand_t state;
  flint_rand_init(state);

  fmpz_randbits(f_a, state, bits);

  char* s_a = fmpz_get_str(NULL, 10, f_a);
  rz_init_val(&rz_a, s_a);
  flint_free(s_a);

  rz_set_u64(&rz_b, exp);

  struct timespec start, end;
  double total_custom = 0, total_flint = 0;

  for (int i = 0; i < iterations; i++) {
    // Custom
    clock_gettime(CLOCK_MONOTONIC, &start);
    rz_pow(&rz_res, &rz_a, &rz_b);
    clock_gettime(CLOCK_MONOTONIC, &end);
    total_custom += get_elapsed_time(start, end);

    // FLINT
    clock_gettime(CLOCK_MONOTONIC, &start);
    fmpz_pow_ui(f_res, f_a, exp);
    clock_gettime(CLOCK_MONOTONIC, &end);
    total_flint += get_elapsed_time(start, end);
  }

  printf("Custom Library: %.6f s (avg)\n", total_custom / iterations);
  printf("FLINT Library:  %.6f s (avg)\n", total_flint / iterations);
  printf("Ratio:          %.2fx slower than FLINT\n",
         (total_custom / iterations) / (total_flint / iterations));

  rz_clear(&rz_a);
  rz_clear(&rz_b);
  rz_clear(&rz_res);
  fmpz_clear(f_a);
  fmpz_clear(f_res);
  flint_rand_clear(state);
}

void benchmark_mul(int bits, int iterations)
{
  printf("\n--- Mul Benchmark (%d bits, %d iterations) ---\n", bits,
         iterations);

  rz_t rz_a, rz_res;
  fmpz_t f_a, f_res;

  rz_init(&rz_a);
  rz_init(&rz_res);
  fmpz_init(f_a);
  fmpz_init(f_res);

  flint_rand_t state;
  flint_rand_init(state);

  fmpz_randbits(f_a, state, bits);

  char* s_a = fmpz_get_str(NULL, 10, f_a);
  rz_init_val(&rz_a, s_a);
  flint_free(s_a);

  struct timespec start, end;
  double total_custom = 0, total_flint = 0;

  for (int i = 0; i < iterations; i++) {
    // Custom
    clock_gettime(CLOCK_MONOTONIC, &start);
    rz_mul(&rz_res, &rz_a, &rz_a);
    clock_gettime(CLOCK_MONOTONIC, &end);
    total_custom += get_elapsed_time(start, end);

    // FLINT
    clock_gettime(CLOCK_MONOTONIC, &start);
    fmpz_mul(f_res, f_a, f_a);
    clock_gettime(CLOCK_MONOTONIC, &end);
    total_flint += get_elapsed_time(start, end);
  }

  printf("Custom Library: %.6f s (avg)\n", total_custom / iterations);
  printf("FLINT Library:  %.6f s (avg)\n", total_flint / iterations);
  printf("Ratio:          %.2fx slower than FLINT\n",
         (total_custom / iterations) / (total_flint / iterations));

  rz_clear(&rz_a);
  rz_clear(&rz_res);
  fmpz_clear(f_a);
  fmpz_clear(f_res);
  flint_rand_clear(state);
}

int main()
{
  int iterations = 10;

  int sizes[] = {512, 1024, 2048, 4096, 5000};
  for (int i = 0; i < 5; i++) {
    benchmark_mod_exp(sizes[i], iterations);
  }

  benchmark_pow(256, 1000, iterations);
  benchmark_pow(512, 1000, iterations);
  benchmark_pow(1024, 1000, iterations);
  benchmark_pow(2048, 1000, iterations);

  benchmark_mul(2048, iterations);
  benchmark_mul(20480, iterations);
  benchmark_mul(204800, iterations);
  benchmark_mul(2048000, iterations);

  return 0;
}
