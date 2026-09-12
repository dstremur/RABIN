/*
 * test_div.c
 *
 * Unified GMP-verified test suite for bn_div (truncated division).
 *
 *   1. Edge cases: 0, 1, -1, 2^64-1, 2^64 boundaries, truncation sign
 *      semantics (q truncates toward 0) + pointer aliasing (q == a and
 *      q == b)
 *   2. 1000 randomized cases (1..4096 bits) vs mpz_tdiv_q
 *   3. Time-based benchmark vs GMP at 512/1024/2048/4096 bits
 *   4. (kept) 2N/N schoolbook vs Newton-Raphson comparison per size
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
  bignum a, b, q;
  mpz_t za, zb, zq;
  char detail[256];

  bn_init_multi(&a, &b, &q, NULL);
  mpz_inits(za, zb, zq, NULL);
  bn_init_val(&a, a_str);
  bn_init_val(&b, b_str);
  mpz_set_str(za, a_str, 10);
  mpz_set_str(zb, b_str, 10);

  snprintf(detail, sizeof(detail), "bn_div [%s] a=%s b=%s", name, a_str, b_str);

  // non-aliased
  bn_div(&q, &a, &b);
  mpz_tdiv_q(zq, za, zb);
  assert_match(detail, a_str, b_str, &q, zq);

  // aliasing: q == a
  bn_init_val(&a, a_str);
  bn_div(&a, &a, &b);
  snprintf(detail, sizeof(detail), "bn_div [%s, q==a] a=%s b=%s", name, a_str,
           b_str);
  assert_match(detail, a_str, b_str, &a, zq);

  // aliasing: q == b
  bn_init_val(&a, a_str);
  bn_init_val(&b, b_str);
  bn_div(&b, &a, &b);
  snprintf(detail, sizeof(detail), "bn_div [%s, q==b] a=%s b=%s", name, a_str,
           b_str);
  assert_match(detail, a_str, b_str, &b, zq);

  bn_free_multi(&a, &b, &q, NULL);
  mpz_clears(za, zb, zq, NULL);
}

static void run_edge_cases()
{
  printf("\n--- bn_div: edge cases ---\n");

  run_case("zero / one", "0", "1");
  run_case("one / one", "1", "1");
  run_case("one / -one", "1", "-1");
  run_case("seven / two (truncates)", "7", "2");
  run_case("neg seven / two (truncates to -3)", "-7", "2");
  run_case("seven / neg two (truncates to -3)", "7", "-2");
  run_case("neg seven / neg two", "-7", "-2");
  run_case("a < b (zero quotient)", "1", "18446744073709551615");
  run_case("a == b", "18446744073709551615", "18446744073709551615");
  run_case("u64-max / 2", "18446744073709551615", "2");
  run_case("u64-max / 3", "18446744073709551615", "3");
  run_case("2^64 / 2^64 (limb boundary)", "18446744073709551616",
           "18446744073709551616");
  run_case("2^64 / u64-max", "18446744073709551616", "18446744073709551615");
  run_case("neg u64-max / 2", "-18446744073709551615", "2");
  run_case("huge / 2 (asymmetric)", "123456789012345678901234567890", "2");

  printf("all edge cases passed\n");
}

// =============================================================================
// PART 2: RANDOMIZED FUZZING vs GMP
// =============================================================================

#define FUZZ_ITERATIONS 1000

static void run_random(gmp_randstate_t state)
{
  printf("\n--- bn_div: %d randomized cases vs mpz_tdiv_q ---\n",
         FUZZ_ITERATIONS);

  for (int i = 0; i < FUZZ_ITERATIONS; i++) {
    int bits_a = 1 + (rand() % 4096);
    int bits_b = 1 + (rand() % 4096);

    bignum a, b, q;
    mpz_t za, zb, zq;
    char* sa;
    char* sb;
    char detail[64];

    bn_init_multi(&a, &b, &q, NULL);
    mpz_inits(za, zb, zq, NULL);

    sa = random_mpz_str(za, bits_a, state);
    sb = random_mpz_str(zb, bits_b, state);

    // divisor must be nonzero
    if (mpz_cmp_ui(zb, 0) == 0) mpz_set_ui(zb, 1);
    sb = mpz_get_str(NULL, 10, zb);

    bn_init_val(&a, sa);
    bn_init_val(&b, sb);

    bn_div(&q, &a, &b);
    mpz_tdiv_q(zq, za, zb);

    snprintf(detail, sizeof(detail), "case %d (a=%d bits, b=%d bits)", i,
             bits_a, bits_b);
    assert_match(detail, sa, sb, &q, zq);

    free(sa);
    free(sb);
    bn_free_multi(&a, &b, &q, NULL);
    mpz_clears(za, zb, zq, NULL);
  }

  printf("all %d random cases passed\n", FUZZ_ITERATIONS);
}

// =============================================================================
// PART 3: TIME-BASED BENCHMARKS vs GMP
// =============================================================================

static void benchmark_div(int bits, double target_sec, gmp_randstate_t state)
{
  bignum bn_a, bn_b, bn_q;
  mpz_t mpz_a, mpz_b, mpz_q;

  bn_init_multi(&bn_a, &bn_b, &bn_q, NULL);
  mpz_inits(mpz_a, mpz_b, mpz_q, NULL);

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
    bn_div(&bn_q, &bn_a, &bn_b);
    ops_custom++;
    clock_gettime(CLOCK_MONOTONIC, &end);
    total_custom = get_elapsed_time(start, end);
  } while (total_custom < target_sec);

  clock_gettime(CLOCK_MONOTONIC, &start);
  do {
    mpz_tdiv_q(mpz_q, mpz_a, mpz_b);
    ops_gmp++;
    clock_gettime(CLOCK_MONOTONIC, &end);
    total_gmp = get_elapsed_time(start, end);
  } while (total_gmp < target_sec);

  // Validate correctness before reporting
  assert_match("bn_div (benchmark)", "(bench input)", "(bench input)", &bn_q,
               mpz_q);

  char size_info[32];
  snprintf(size_info, sizeof(size_info), "%d bits", bits);
  print_table_row("bn_div", size_info, total_custom / ops_custom,
                  total_gmp / ops_gmp);

  bn_free_multi(&bn_a, &bn_b, &bn_q, NULL);
  mpz_clears(mpz_a, mpz_b, mpz_q, NULL);
}

// =============================================================================
// PART 4 (kept): 2N / N SCHOOLBOOK vs NEWTON COMPARISON
// =============================================================================

// Global variable to prevent compiler dead-code elimination
static volatile u64 dummy_sum = 0;

static void benchmark_div_real(int divisor_limbs)
{
  // 1. The 2N / N Rule
  // To properly test division, the dividend (a) must be twice as large as the
  // divisor (b).
  int dividend_limbs = divisor_limbs * 2;

  bignum a, b, res_school, res_karat;
  bn_init_multi(&a, &b, &res_school, &res_karat, NULL);

  // Generate random numbers
  bn_gen_random(&a, dividend_limbs * 64);
  bn_gen_random(&b, divisor_limbs * 64);

  // Ensure the highest bit of 'b' is set so it is truly 'divisor_limbs' long
  bn_set_bit(&b, (divisor_limbs * 64) - 1);

  // 2. Time Schoolbook
  clock_t start = clock();
  bn_div(&res_school, &a, &b);
  clock_t end = clock();
  double time_school = ((double)(end - start)) / CLOCKS_PER_SEC;

  // 3. Time Newton-Raphson
  start = clock();
  bn_newton_div(&res_karat, &a, &b);
  end = clock();
  double time_karat = ((double)(end - start)) / CLOCKS_PER_SEC;

  // 4. Verify Correctness
  if (bn_cmp(&res_school, &res_karat) != 0) {
    printf("[FAIL] Mismatch! Divisor Limbs: %d\n", divisor_limbs);
    printf("       School quotient limbs: %llu\n",
           (unsigned long long)res_school.size);
    printf("       Newton quotient limbs: %llu\n",
           (unsigned long long)res_karat.size);
  } else {
    double speedup = (time_karat > 0.0) ? (time_school / time_karat) : 0.0;
    printf(
        "Divisor Limbs: %5d | Dividend Limbs: %5d | School: %8.6fs | Newton: "
        "%8.6fs | Speedup: %8.2fx\n",
        divisor_limbs, dividend_limbs, time_school, time_karat, speedup);
  }

  // 5. Anti-Optimization Hack
  if (res_school.size > 0) dummy_sum += res_school.limbs[0];
  if (res_karat.size > 0) dummy_sum += res_karat.limbs[0];

  bn_free_multi(&a, &b, &res_school, &res_karat, NULL);
}

static void run_internal_comparison()
{
  printf("\n--- bn_div: 2N/N schoolbook vs Newton (internal) ---\n");
  int sizes[] = {16,   32,   64,   128,  256,   512,
                 1024, 2048, 4096, 8192, 16000, 32000};
  for (size_t i = 0; i < sizeof(sizes) / sizeof(sizes[0]); i++) {
    benchmark_div_real(sizes[i]);
  }
  printf("Benchmark complete. (Checksum: %llu)\n",
         (unsigned long long)dummy_sum);
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
  benchmark_div(2048, bench_budget(2048), state);
  benchmark_div(8192, bench_budget(8192), state);
  benchmark_div(32768, bench_budget(32768), state);
  benchmark_div(65536, bench_budget(65536), state);
  benchmark_div(65539, bench_budget(65539), state);
  print_table_footer();

  run_internal_comparison();

  gmp_randclear(state);
  bn_free_constants();
  return 0;
}
