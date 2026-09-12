/*
 * test_log_2.c
 *
 * Unified GMP-verified test suite for bn_log_2.
 *
 *   1. Edge cases: 0, 1, 2, powers of two, 2^64-1 / 2^64 boundaries
 *   2. 1000 randomized cases (1..4096 bits) vs mpz_sizeinbase(a, 2) - 1
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

// floor(log2(a)) for a > 0 equals GMP's bit length minus one
void assert_log2_match(const char* op, const char* a_str, const bignum* r,
                       mpz_t a)
{
  mpz_t expected;
  mpz_init(expected);
  mpz_set_ui(expected, (size_t)mpz_sizeinbase(a, 2) - 1);

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

  snprintf(detail, sizeof(detail), "bn_log_2 [%s] a=%s", name, a_str);

  // non-aliased
  bn_log_2(&res, &a);
  assert_log2_match(detail, a_str, &res, za);

  // aliasing: result == a (safe: a is fully read before r is written)
  bn_init_val(&a, a_str);
  bn_log_2(&a, &a);
  assert_log2_match("bn_log_2 (r==a)", a_str, &a, za);

  bn_free_multi(&a, &res, NULL);
  mpz_clear(za);
}

static void run_edge_cases()
{
  printf("\n--- bn_log_2: edge cases ---\n");

  // a = 0: r must be left unchanged
  {
    bignum a, res;
    bn_init_multi(&a, &res, NULL);
    bn_init_val(&a, "0");
    bn_set_u64(&res, 12345);
    bn_log_2(&res, &a);
    if (!bn_is_eq_i64(&res, 12345)) {
      fprintf(stderr,
              "\n[FATAL ERROR] bn_log_2(0) must leave r unchanged, r=%s\n",
              bn_to_string(&res));
      bn_free_multi(&a, &res, NULL);
      exit(EXIT_FAILURE);
    }
    bn_free_multi(&a, &res, NULL);
  }

  run_case("one", "1");
  run_case("two", "2");
  run_case("three", "3");
  run_case("u64-max (2^64-1)", "18446744073709551615");
  run_case("2^64", "18446744073709551616");
  run_case("2^64 + 1", "18446744073709551617");

  // exact powers of two at limb boundaries
  mpz_t pow2;
  mpz_init(pow2);
  int exps[] = {64, 127, 128, 129, 255, 256, 1023, 1024, 1025, 4095, 4096};
  for (size_t i = 0; i < sizeof(exps) / sizeof(exps[0]); i++) {
    mpz_set_ui(pow2, 2);
    mpz_pow_ui(pow2, pow2, (unsigned long)exps[i]);
    char* s = mpz_get_str(NULL, 10, pow2);
    char name[64];
    snprintf(name, sizeof(name), "2^%d", exps[i]);
    run_case(name, s);
    free(s);
  }
  mpz_clear(pow2);

  printf("all edge cases passed\n");
}

// =============================================================================
// PART 2: RANDOMIZED FUZZING vs GMP
// =============================================================================

#define FUZZ_ITERATIONS 1000

static void run_random(gmp_randstate_t state)
{
  printf("\n--- bn_log_2: %d randomized cases vs GMP ---\n", FUZZ_ITERATIONS);

  for (int i = 0; i < FUZZ_ITERATIONS; i++) {
    int bits = 1 + (rand() % 4096);

    bignum a, res;
    mpz_t za;
    char* sa;

    bn_init_multi(&a, &res, NULL);
    mpz_init(za);

    sa = random_mpz_str(za, bits, state);
    bn_init_val(&a, sa);

    bn_log_2(&res, &a);
    assert_log2_match("bn_log_2 (random)", sa, &res, za);

    free(sa);
    bn_free_multi(&a, &res, NULL);
    mpz_clear(za);
  }

  printf("all %d random cases passed\n", FUZZ_ITERATIONS);
}

// =============================================================================
// PART 3: TIME-BASED BENCHMARKS vs GMP
// =============================================================================

static void benchmark_log2(int bits, double target_sec, gmp_randstate_t state)
{
  bignum bn_a, bn_res;
  mpz_t mpz_a;

  bn_init_multi(&bn_a, &bn_res, NULL);
  mpz_init(mpz_a);

  char* sa = random_mpz_str(mpz_a, bits, state);
  bn_init_val(&bn_a, sa);
  free(sa);

  struct timespec start, end;
  int ops_custom = 0, ops_gmp = 0;
  double total_custom = 0, total_gmp = 0;

  clock_gettime(CLOCK_MONOTONIC, &start);
  do {
    bn_log_2(&bn_res, &bn_a);
    ops_custom++;
    clock_gettime(CLOCK_MONOTONIC, &end);
    total_custom = get_elapsed_time(start, end);
  } while (total_custom < target_sec);

  clock_gettime(CLOCK_MONOTONIC, &start);
  do {
    (void)mpz_sizeinbase(mpz_a, 2);
    ops_gmp++;
    clock_gettime(CLOCK_MONOTONIC, &end);
    total_gmp = get_elapsed_time(start, end);
  } while (total_gmp < target_sec);

  // Validate correctness before reporting
  assert_log2_match("bn_log_2 (benchmark)", "(bench input)", &bn_res, mpz_a);

  char size_info[32];
  snprintf(size_info, sizeof(size_info), "%d bits", bits);
  print_table_row("bn_log_2", size_info, total_custom / ops_custom,
                  total_gmp / ops_gmp);

  bn_free_multi(&bn_a, &bn_res, NULL);
  mpz_clear(mpz_a);
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
  benchmark_log2(4096, bench_budget(4096), state);
  benchmark_log2(65536, bench_budget(65536), state);
  benchmark_log2(1000000, bench_budget(1000000), state);
  print_table_footer();

  gmp_randclear(state);
  return 0;
}
