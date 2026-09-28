/*
 * test_lshift.c
 *
 * Unified GMP-verified test suite for rz_lshift.
 *
 *   1. Edge cases: 0, 1, -1, 2^64-1, 2^64 boundaries, limb-crossing and
 *      huge shifts + pointer aliasing
 *   2. 1000 randomized cases (1..4096 bits) vs mpz_mul_2exp
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

#include "../include/rabin.h"

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

void assert_lshift_match(const char* op, const char* a_str, int shift,
                         const rz_t* r, mpz_t a)
{
  mpz_t expected;
  mpz_init(expected);
  mpz_mul_2exp(expected, a, (size_t)shift);

  char* exp_str = mpz_get_str(NULL, 10, expected);
  char* custom_str = rz_to_string(r);

  if (strcmp(exp_str, custom_str) != 0) {
    fprintf(stderr, "\n[FATAL ERROR] Correctness failure in %s!\n", op);
    fprintf(stderr, "Input A: %s (shift=%d)\n", a_str ? a_str : "(n/a)", shift);
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

static void run_case(const char* name, const char* a_str, int shift)
{
  rz_t a, res;
  mpz_t za;
  char detail[256];

  rz_init_multi(&a, &res, NULL);
  mpz_init(za);
  rz_init_val(&a, a_str);
  mpz_set_str(za, a_str, 10);

  snprintf(detail, sizeof(detail), "rz_lshift [%s] a=%s shift=%d", name, a_str,
           shift);

  // non-aliased
  rz_lshift(&res, &a, shift);
  assert_lshift_match(detail, a_str, shift, &res, za);

  // aliasing: result == a
  rz_init_val(&a, a_str);
  rz_lshift(&a, &a, shift);
  assert_lshift_match("rz_lshift (r==a)", a_str, shift, &a, za);

  rz_clear_multi(&a, &res, NULL);
  mpz_clear(za);
}

static void run_edge_cases()
{
  printf("\n--- rz_lshift: edge cases ---\n");

  run_case("zero", "0", 123);
  run_case("one shift 0", "1", 0);
  run_case("one shift 1", "1", 1);
  run_case("one shift 64 (new limb)", "1", 64);
  run_case("one shift 65", "1", 65);
  run_case("one huge shift", "1", 100000);
  run_case("neg-one shift 3", "-1", 3);
  run_case("u64-max shift 0", "18446744073709551615", 0);
  run_case("u64-max shift 1 (carry out of limb)", "18446744073709551615", 1);
  run_case("u64-max shift 63", "18446744073709551615", 63);
  run_case("u64-max shift 64 (exact limb)", "18446744073709551615", 64);
  run_case("u64-max shift 65", "18446744073709551615", 65);
  run_case("2^64 shift 128 (two limbs)", "18446744073709551616", 128);
  run_case("neg u64-max shift 2", "-18446744073709551615", 2);

  // shift beyond the bit length grows the number
  run_case("shift == bit length", "3", 2);

  printf("all edge cases passed\n");
}

// =============================================================================
// PART 2: RANDOMIZED FUZZING vs GMP
// =============================================================================

#define FUZZ_ITERATIONS 1000

static void run_random(gmp_randstate_t state)
{
  printf("\n--- rz_lshift: %d randomized cases vs mpz_mul_2exp ---\n",
         FUZZ_ITERATIONS);

  for (int i = 0; i < FUZZ_ITERATIONS; i++) {
    int bits = 1 + (rand() % 4096);
    int shift = rand() % 8192;

    rz_t a, res;
    mpz_t za;
    char* sa;

    rz_init_multi(&a, &res, NULL);
    mpz_init(za);

    sa = random_mpz_str(za, bits, state);
    rz_init_val(&a, sa);

    rz_lshift(&res, &a, shift);
    assert_lshift_match("rz_lshift (random)", sa, shift, &res, za);

    free(sa);
    rz_clear_multi(&a, &res, NULL);
    mpz_clear(za);
  }

  printf("all %d random cases passed\n", FUZZ_ITERATIONS);
}

// =============================================================================
// PART 3: TIME-BASED BENCHMARKS vs GMP
// =============================================================================

static void benchmark_lshift(int bits, double target_sec, gmp_randstate_t state)
{
  const int shift = 123;  // fixed limb-crossing shift
  rz_t rz_a, rz_res;
  mpz_t mpz_a, mpz_res;

  rz_init_multi(&rz_a, &rz_res, NULL);
  mpz_inits(mpz_a, mpz_res, NULL);

  char* sa = random_mpz_str(mpz_a, bits, state);
  rz_init_val(&rz_a, sa);
  free(sa);

  struct timespec start, end;
  int ops_custom = 0, ops_gmp = 0;
  double total_custom = 0, total_gmp = 0;

  clock_gettime(CLOCK_MONOTONIC, &start);
  do {
    rz_lshift(&rz_res, &rz_a, shift);
    ops_custom++;
    clock_gettime(CLOCK_MONOTONIC, &end);
    total_custom = get_elapsed_time(start, end);
  } while (total_custom < target_sec);

  clock_gettime(CLOCK_MONOTONIC, &start);
  do {
    mpz_mul_2exp(mpz_res, mpz_a, (size_t)shift);
    ops_gmp++;
    clock_gettime(CLOCK_MONOTONIC, &end);
    total_gmp = get_elapsed_time(start, end);
  } while (total_gmp < target_sec);

  // Validate correctness before reporting
  assert_lshift_match("rz_lshift (benchmark)", "(bench input)", shift, &rz_res,
                      mpz_a);

  char size_info[32];
  snprintf(size_info, sizeof(size_info), "%d bits", bits);
  print_table_row("rz_lshift", size_info, total_custom / ops_custom,
                  total_gmp / ops_gmp);

  rz_clear_multi(&rz_a, &rz_res, NULL);
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
  benchmark_lshift(4096, bench_budget(4096), state);
  benchmark_lshift(65536, bench_budget(65536), state);
  benchmark_lshift(1000000, bench_budget(1000000), state);
  print_table_footer();

  gmp_randclear(state);
  return 0;
}
