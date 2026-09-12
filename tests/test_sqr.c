/*
 * test_sqr.c
 *
 * Unified GMP-verified test suite for bn_sqr.
 *
 *   1. Edge cases: 0, 1, -1, 2^64-1, 2^64 boundaries + pointer aliasing
 *   2. 1000 randomized cases (1..4096 bits) vs mpz_mul(z, x, x)
 *   3. Time-based benchmark vs GMP at 512/1024/2048/4096 bits
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

void assert_sqr_match(const char* op, const char* a_str, const bignum* r,
                      mpz_t a)
{
  mpz_t expected;
  mpz_init(expected);
  mpz_mul(expected, a, a);

  char* exp_str = mpz_get_str(NULL, 10, expected);
  char* custom_str = bn_to_string(r);

  if (strcmp(exp_str, custom_str) != 0) {
    fprintf(stderr, "\n[FATAL ERROR] Correctness failure in %s!\n", op);
    fprintf(stderr, "Input A: %s\n", a_str ? a_str : "(n/a)");
    fprintf(stderr, "GMP Result:    %s\n", exp_str);
    fprintf(stderr, "Custom Result: %s\n", custom_str);
    free(exp_str);
    free(custom_str);
    mpz_clear(expected);
    exit(EXIT_FAILURE);
  }

  free(exp_str);
  free(custom_str);
  mpz_clear(expected);
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

static void run_case(const char* name, const char* a_str)
{
  bignum a, res;
  mpz_t za;
  char detail[256];

  bn_init_multi(&a, &res, NULL);
  mpz_init(za);
  bn_init_val(&a, a_str);
  mpz_set_str(za, a_str, 10);

  snprintf(detail, sizeof(detail), "bn_sqr [%s] a=%s", name, a_str);

  // non-aliased
  bn_sqr(&res, &a);
  assert_sqr_match(detail, a_str, &res, za);

  // aliasing: result == a
  bn_init_val(&a, a_str);
  bn_sqr(&a, &a);
  assert_sqr_match("bn_sqr (r==a)", a_str, &a, za);

  bn_free_multi(&a, &res, NULL);
  mpz_clear(za);
}

static void run_edge_cases()
{
  printf("\n--- bn_sqr: edge cases ---\n");

  run_case("zero", "0");
  run_case("one", "1");
  run_case("neg-one (positive result)", "-1");
  run_case("two", "2");
  run_case("u64-max (128-bit result)", "18446744073709551615");
  run_case("neg u64-max", "-18446744073709551615");
  run_case("2^64 (exact 2-limb square)", "18446744073709551616");

  printf("all edge cases passed\n");
}

// =============================================================================
// PART 2: RANDOMIZED FUZZING vs GMP
// =============================================================================

#define FUZZ_ITERATIONS 1000

static void run_random(gmp_randstate_t state)
{
  printf("\n--- bn_sqr: %d randomized cases vs mpz_mul ---\n", FUZZ_ITERATIONS);

  for (int i = 0; i < FUZZ_ITERATIONS; i++) {
    int bits = 1 + (rand() % 4096);

    bignum a, res;
    mpz_t za;
    char* sa;

    bn_init_multi(&a, &res, NULL);
    mpz_init(za);

    sa = random_mpz_str(za, bits, state);
    bn_init_val(&a, sa);

    bn_sqr(&res, &a);
    assert_sqr_match("bn_sqr (random)", sa, &res, za);

    free(sa);
    bn_free_multi(&a, &res, NULL);
    mpz_clear(za);
  }

  printf("all %d random cases passed\n", FUZZ_ITERATIONS);
}

// =============================================================================
// PART 3: TIME-BASED BENCHMARKS vs GMP
// =============================================================================

static void benchmark_sqr(int bits, double target_sec, gmp_randstate_t state)
{
  bignum bn_a, bn_res;
  mpz_t mpz_a, mpz_res;

  bn_init_multi(&bn_a, &bn_res, NULL);
  mpz_inits(mpz_a, mpz_res, NULL);

  char* sa = random_mpz_str(mpz_a, bits, state);
  bn_init_val(&bn_a, sa);
  free(sa);

  struct timespec start, end;
  int ops_custom = 0, ops_gmp = 0;
  double total_custom = 0, total_gmp = 0;

  clock_gettime(CLOCK_MONOTONIC, &start);
  do {
    bn_sqr(&bn_res, &bn_a);
    ops_custom++;
    clock_gettime(CLOCK_MONOTONIC, &end);
    total_custom = get_elapsed_time(start, end);
  } while (total_custom < target_sec);

  clock_gettime(CLOCK_MONOTONIC, &start);
  do {
    mpz_mul(mpz_res, mpz_a, mpz_a);
    ops_gmp++;
    clock_gettime(CLOCK_MONOTONIC, &end);
    total_gmp = get_elapsed_time(start, end);
  } while (total_gmp < target_sec);

  // Validate correctness before reporting
  assert_sqr_match("bn_sqr (benchmark)", "(bench input)", &bn_res, mpz_a);

  char size_info[32];
  snprintf(size_info, sizeof(size_info), "%d bits", bits);
  print_table_row("bn_sqr", size_info, total_custom / ops_custom,
                  total_gmp / ops_gmp);

  bn_free_multi(&bn_a, &bn_res, NULL);
  mpz_clears(mpz_a, mpz_res, NULL);
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
  gmp_randstate_t state;
  gmp_randinit_default(state);
  gmp_randseed_ui(state, (unsigned long)time(NULL));
  srand((unsigned int)time(NULL));

  run_edge_cases();
  run_random(state);

  print_table_header();
  benchmark_sqr(2048, bench_budget(2048), state);
  benchmark_sqr(8192, bench_budget(8192), state);
  benchmark_sqr(32768, bench_budget(32768), state);
  benchmark_sqr(65536, bench_budget(65536), state);
  print_table_footer();

  gmp_randclear(state);
  return 0;
}
