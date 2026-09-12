/*
 * test_mul.c
 *
 * Unified GMP-verified test suite for bn_mul.
 *
 *   1. Edge cases: 0, 1, -1, 2^64-1, 2^64 boundaries, limb-boundary carries
 *      + pointer aliasing (r == a and r == b)
 *   2. 1000 randomized cases (1..4096 bits) vs mpz_mul
 *   3. Time-based benchmark vs GMP at 512/1024/2048/4096 bits
 *   4. (kept) schoolbook vs Karatsuba comparison per size
 *
 * Copyright (C) 2026 Diego Strebel
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "../include/bignum.h"

// =============================================================================
// HELPERS & VALIDATION
// =============================================================================

double get_elapsed_time(struct timespec start, struct timespec end)
{
  return (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
}

void print_table_header()
{
  printf(
      "\n+--------------------+--------------+---------------+---------------+-"
      "-----------+\n");
  printf(
      "| Operation          | Bits/Size    | Custom (s)    | GMP (s)       | "
      "Ratio      |\n");
  printf(
      "+--------------------+--------------+---------------+---------------+---"
      "---------+\n");
}

void print_table_footer()
{
  printf(
      "+--------------------+--------------+---------------+---------------+---"
      "---------+\n");
}

void print_table_row(const char* op, const char* size_info, double avg_custom,
                     double avg_gmp)
{
  double ratio = (avg_gmp > 0) ? (avg_custom / avg_gmp) : 0.0;
  printf("| %-18s | %-12s | %10.7f    | %10.7f    | %7.2fx   |\n", op,
         size_info, avg_custom, avg_gmp, ratio);
}

void assert_match(const char* op, const char* a_str, const char* b_str,
                  const bignum* bn, mpz_t mpz)
{
  char* gmp_str = mpz_get_str(NULL, 10, mpz);
  char* custom_str = bn_to_string(bn);

  if (strcmp(gmp_str, custom_str) != 0) {
    fprintf(stderr, "\n[FATAL ERROR] Correctness failure in %s!\n", op);
    fprintf(stderr, "Input A: %s\n", a_str ? a_str : "(n/a)");
    fprintf(stderr, "Input B: %s\n", b_str ? b_str : "(n/a)");
    fprintf(stderr, "GMP Result:    %s\n", gmp_str);
    fprintf(stderr, "Custom Result: %s\n", custom_str);
    free(gmp_str);
    free(custom_str);
    exit(EXIT_FAILURE);
  }

  free(gmp_str);
  free(custom_str);
}

char* random_mpz_str(mpz_t z, int bits, gmp_randstate_t state)
{
  mpz_urandomb(z, state, bits);
  if (rand() & 1) mpz_neg(z, z);
  return mpz_get_str(NULL, 10, z);
}

// =============================================================================
// PART 1: EDGE CASES (incl. pointer aliasing)
// =============================================================================

static void run_case(const char* name, const char* a_str, const char* b_str)
{
  bignum a, b, res;
  mpz_t za, zb, zr;
  char detail[256];

  bn_init_multi(&a, &b, &res, NULL);
  mpz_inits(za, zb, zr, NULL);
  bn_init_val(&a, a_str);
  bn_init_val(&b, b_str);
  mpz_set_str(za, a_str, 10);
  mpz_set_str(zb, b_str, 10);

  snprintf(detail, sizeof(detail), "bn_mul [%s] a=%s b=%s", name, a_str, b_str);

  // non-aliased
  bn_mul(&res, &a, &b);
  mpz_mul(zr, za, zb);
  assert_match(detail, a_str, b_str, &res, zr);

  // aliasing: r == a
  bn_init_val(&a, a_str);
  bn_mul(&a, &a, &b);
  snprintf(detail, sizeof(detail), "bn_mul [%s, r==a] a=%s b=%s", name, a_str,
           b_str);
  assert_match(detail, a_str, b_str, &a, zr);

  // aliasing: r == b
  bn_init_val(&a, a_str);
  bn_init_val(&b, b_str);
  bn_mul(&b, &a, &b);
  snprintf(detail, sizeof(detail), "bn_mul [%s, r==b] a=%s b=%s", name, a_str,
           b_str);
  assert_match(detail, a_str, b_str, &b, zr);

  bn_free_multi(&a, &b, &res, NULL);
  mpz_clears(za, zb, zr, NULL);
}

static void run_edge_cases()
{
  printf("\n--- bn_mul: edge cases ---\n");

  run_case("zeros", "0", "0");
  run_case("one * one", "1", "1");
  run_case("neg-one * neg-one", "-1", "-1");
  run_case("one * neg-one", "1", "-1");
  run_case("multiply by zero", "987654321987654321", "0");
  run_case("u64-max * u64-max (128-bit result)", "18446744073709551615",
           "18446744073709551615");
  run_case("u64-max * 2 (carry cascade)", "18446744073709551615", "2");
  run_case("2^64 * 2^64 (limb boundary)", "18446744073709551616",
           "18446744073709551616");
  run_case("2^64 * u64-max", "18446744073709551616", "18446744073709551615");
  run_case("negative * positive", "-18446744073709551615", "12345");
  run_case("asymmetric sizes", "123456789012345678901234567890", "2");

  printf("all edge cases passed\n");
}

// =============================================================================
// PART 2: RANDOMIZED FUZZING vs GMP
// =============================================================================

#define FUZZ_ITERATIONS 1000

static void run_random(gmp_randstate_t state)
{
  printf("\n--- bn_mul: %d randomized cases vs mpz_mul ---\n", FUZZ_ITERATIONS);

  for (int i = 0; i < FUZZ_ITERATIONS; i++) {
    int bits_a = 1 + (rand() % 4096);
    int bits_b = 1 + (rand() % 4096);

    bignum a, b, res;
    mpz_t za, zb, zr;
    char* sa;
    char* sb;
    char detail[64];

    bn_init_multi(&a, &b, &res, NULL);
    mpz_inits(za, zb, zr, NULL);

    sa = random_mpz_str(za, bits_a, state);
    sb = random_mpz_str(zb, bits_b, state);
    bn_init_val(&a, sa);
    bn_init_val(&b, sb);

    bn_mul(&res, &a, &b);
    mpz_mul(zr, za, zb);

    snprintf(detail, sizeof(detail), "case %d (a=%d bits, b=%d bits)", i,
             bits_a, bits_b);
    assert_match(detail, sa, sb, &res, zr);

    free(sa);
    free(sb);
    bn_free_multi(&a, &b, &res, NULL);
    mpz_clears(za, zb, zr, NULL);
  }

  printf("all %d random cases passed\n", FUZZ_ITERATIONS);
}

// =============================================================================
// PART 3: TIME-BASED BENCHMARKS vs GMP
// =============================================================================

static void benchmark_mul(int bits, double target_sec, gmp_randstate_t state)
{
  bignum bn_a, bn_b, bn_res;
  mpz_t mpz_a, mpz_b, mpz_res;

  bn_init_multi(&bn_a, &bn_b, &bn_res, NULL);
  mpz_inits(mpz_a, mpz_b, mpz_res, NULL);

  char* sa = random_mpz_str(mpz_a, bits, state);
  char* sb = random_mpz_str(mpz_b, bits, state);
  bn_init_val(&bn_a, sa);
  bn_init_val(&bn_b, sb);
  free(sa);
  free(sb);

  struct timespec start, end;
  int ops_custom = 0, ops_gmp = 0;
  double total_custom = 0, total_gmp = 0;

  clock_gettime(CLOCK_MONOTONIC, &start);
  do {
    bn_mul(&bn_res, &bn_a, &bn_b);
    ops_custom++;
    clock_gettime(CLOCK_MONOTONIC, &end);
    total_custom = get_elapsed_time(start, end);
  } while (total_custom < target_sec);

  clock_gettime(CLOCK_MONOTONIC, &start);
  do {
    mpz_mul(mpz_res, mpz_a, mpz_b);
    ops_gmp++;
    clock_gettime(CLOCK_MONOTONIC, &end);
    total_gmp = get_elapsed_time(start, end);
  } while (total_gmp < target_sec);

  // Validate correctness before reporting
  assert_match("bn_mul (benchmark)", "(bench input)", "(bench input)", &bn_res,
               mpz_res);

  char size_info[32];
  snprintf(size_info, sizeof(size_info), "%d bits", bits);
  print_table_row("bn_mul", size_info, total_custom / ops_custom,
                  total_gmp / ops_gmp);

  bn_free_multi(&bn_a, &bn_b, &bn_res, NULL);
  mpz_clears(mpz_a, mpz_b, mpz_res, NULL);
}

// =============================================================================
// PART 4 (kept): SCHOOLBOOK vs KARATSUBA COMPARISON
// =============================================================================

static void benchmark_mul_internal(int limbs)
{
  bignum a, b, res_school, res_karat;
  bn_init_multi(&a, &b, &res_school, &res_karat, NULL);

  // Generate random numbers of a specific bit length (limbs * 64)
  bn_gen_random(&a, limbs * 64);
  bn_gen_random(&b, limbs * 64);

  // 1. Time Schoolbook
  clock_t start = clock();
  bn_mul_school(&res_school, &a, &b);
  clock_t end = clock();
  double time_school = (double)(end - start) / CLOCKS_PER_SEC;

  // 2. Time Karatsuba
  start = clock();
  bn_mul(&res_karat, &a, &b);
  end = clock();
  double time_karat = (double)(end - start) / CLOCKS_PER_SEC;

  // 3. Verify Correctness
  if (bn_cmp(&res_school, &res_karat) != 0) {
    printf("[FAIL] Mismatch at %d limbs!\n", limbs);
  } else {
    double speedup = (time_karat > 0.0) ? (time_school / time_karat) : 0.0;
    printf("Limbs: %4d | School: %.6fs | Karatsuba: %.6fs | Speedup: %.2fx\n",
           limbs, time_school, time_karat, speedup);
  }

  bn_free_multi(&a, &b, &res_school, &res_karat, NULL);
}

static void run_internal_comparison()
{
  printf("\n--- bn_mul: schoolbook vs Karatsuba (internal) ---\n");
  int sizes[] = {16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192};
  for (size_t i = 0; i < sizeof(sizes) / sizeof(sizes[0]); i++) {
    benchmark_mul_internal(sizes[i]);
  }
}

// =============================================================================
// MAIN
// =============================================================================

// Time budget scales down with size: the largest sizes need only 1-2
// iterations to be representative, and this keeps the total run time of
// the extended size list reasonable.
static double bench_budget(int bits)
{
  if (bits <= 32768) return 1.0;
  if (bits <= 131072) return 0.5;
  return 0.25;
}

int main()
{
  bn_init_constants();

  gmp_randstate_t state;
  gmp_randinit_default(state);
  gmp_randseed_ui(state, (unsigned long)time(NULL));
  srand((unsigned int)time(NULL));

  run_edge_cases();
  run_random(state);

  print_table_header();
  benchmark_mul(2048, bench_budget(2048), state);
  benchmark_mul(8192, bench_budget(8192), state);
  benchmark_mul(32768, bench_budget(32768), state);
  benchmark_mul(65536, bench_budget(65536), state);
  print_table_footer();

  run_internal_comparison();

  gmp_randclear(state);
  bn_free_constants();
  return 0;
}
