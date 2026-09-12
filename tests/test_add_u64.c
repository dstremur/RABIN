/*
 * test_add_u64.c
 *
 * Unified GMP-verified test suite for bn_add_u64.
 *
 *   1. Edge cases: 0, 1, -1, 2^64-1, 2^64 boundaries, u64 max addend,
 *      sign flips + pointer aliasing (r == a)
 *   2. 1000 randomized cases (a: 1..4096 bits, b: full 64-bit) vs mpz_add_ui
 *   3. Time-based benchmark vs GMP at 512/1024/2048/4096 bits
 *
 * Copyright (C) 2026 Diego Strebel
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <gmp.h>
#include <stdint.h>
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

void assert_match(const char* op, const char* a_str, uint64_t b,
                  const bignum* bn, mpz_t expected)
{
  char* gmp_str = mpz_get_str(NULL, 10, expected);
  char* custom_str = bn_to_string(bn);

  if (strcmp(gmp_str, custom_str) != 0) {
    fprintf(stderr, "\n[FATAL ERROR] Correctness failure in %s!\n", op);
    fprintf(stderr, "Input A: %s\n", a_str ? a_str : "(n/a)");
    fprintf(stderr, "Input B: %llu\n", (unsigned long long)b);
    fprintf(stderr, "GMP Result:    %s\n", gmp_str);
    fprintf(stderr, "Custom Result: %s\n", custom_str);
    free(gmp_str);
    free(custom_str);
    exit(EXIT_FAILURE);
  }

  free(gmp_str);
  free(custom_str);
}

// random signed mpz a + its decimal string
static char* random_a(mpz_t za, int bits, gmp_randstate_t state)
{
  mpz_urandomb(za, state, bits);
  if (rand() & 1) mpz_neg(za, za);
  return mpz_get_str(NULL, 10, za);
}

static uint64_t random_u64(gmp_randstate_t state)
{
  mpz_t t;
  mpz_init(t);
  mpz_urandomb(t, state, 64);
  uint64_t v = mpz_get_ui(t);
  mpz_clear(t);
  return v;
}

// =============================================================================
// PART 1: EDGE CASES (incl. pointer aliasing)
// =============================================================================

static void run_case(const char* name, const char* a_str, uint64_t b)
{
  bignum a, res;
  mpz_t za, zr;
  char detail[256];

  bn_init_multi(&a, &res, NULL);
  mpz_inits(za, zr, NULL);
  bn_init_val(&a, a_str);
  mpz_set_str(za, a_str, 10);

  snprintf(detail, sizeof(detail), "bn_add_u64 [%s] a=%s b=%llu", name, a_str,
           (unsigned long long)b);

  // non-aliased
  bn_add_u64(&res, &a, b);
  mpz_add_ui(zr, za, (unsigned long)b);
  assert_match(detail, a_str, b, &res, zr);

  // aliasing: r == a
  bn_init_val(&a, a_str);
  bn_add_u64(&a, &a, b);
  assert_match("bn_add_u64 (r==a)", a_str, b, &a, zr);

  bn_free_multi(&a, &res, NULL);
  mpz_clears(za, zr, NULL);
}

static void run_edge_cases()
{
  printf("\n--- bn_add_u64: edge cases ---\n");

  run_case("zero + zero", "0", 0);
  run_case("zero + one", "0", 1);
  run_case("one + zero", "1", 0);
  run_case("neg-one + one (to +0)", "-1", 1);
  run_case("u64-max + 1 (carry into 2nd limb)", "18446744073709551615", 1);
  run_case("u64-max + u64-max", "18446744073709551615", UINT64_MAX);
  run_case("2^64 + u64-max", "18446744073709551616", UINT64_MAX);
  run_case("neg small + big u64 (sign flip)", "-5", 10);
  run_case("neg big + small u64 (stays neg)", "-18446744073709551615", 7);
  run_case("neg u64-max + u64-max (to +0)", "-18446744073709551615",
           UINT64_MAX);

  printf("all edge cases passed\n");
}

// =============================================================================
// PART 2: RANDOMIZED FUZZING vs GMP
// =============================================================================

#define FUZZ_ITERATIONS 1000

static void run_random(gmp_randstate_t state)
{
  printf("\n--- bn_add_u64: %d randomized cases vs mpz_add_ui ---\n",
         FUZZ_ITERATIONS);

  for (int i = 0; i < FUZZ_ITERATIONS; i++) {
    int bits = 1 + (rand() % 4096);
    uint64_t b = random_u64(state);

    bignum a, res;
    mpz_t za, zr;
    char* sa;
    char detail[64];

    bn_init_multi(&a, &res, NULL);
    mpz_inits(za, zr, NULL);

    sa = random_a(za, bits, state);
    bn_init_val(&a, sa);

    bn_add_u64(&res, &a, b);
    mpz_add_ui(zr, za, (unsigned long)b);

    snprintf(detail, sizeof(detail), "case %d (a=%d bits)", i, bits);
    assert_match(detail, sa, b, &res, zr);

    free(sa);
    bn_free_multi(&a, &res, NULL);
    mpz_clears(za, zr, NULL);
  }

  printf("all %d random cases passed\n", FUZZ_ITERATIONS);
}

// =============================================================================
// PART 3: TIME-BASED BENCHMARKS vs GMP
// =============================================================================

static void benchmark_add_u64(int bits, double target_sec,
                              gmp_randstate_t state)
{
  const uint64_t b = 0x9e3779b97f4a7c15ULL;
  bignum bn_a, bn_res;
  mpz_t mpz_a, mpz_res;

  bn_init_multi(&bn_a, &bn_res, NULL);
  mpz_inits(mpz_a, mpz_res, NULL);

  char* sa = random_a(mpz_a, bits, state);
  bn_init_val(&bn_a, sa);
  free(sa);

  struct timespec start, end;
  int ops_custom = 0, ops_gmp = 0;
  double total_custom = 0, total_gmp = 0;

  clock_gettime(CLOCK_MONOTONIC, &start);
  do {
    bn_add_u64(&bn_res, &bn_a, b);
    ops_custom++;
    clock_gettime(CLOCK_MONOTONIC, &end);
    total_custom = get_elapsed_time(start, end);
  } while (total_custom < target_sec);

  clock_gettime(CLOCK_MONOTONIC, &start);
  do {
    mpz_add_ui(mpz_res, mpz_a, (unsigned long)b);
    ops_gmp++;
    clock_gettime(CLOCK_MONOTONIC, &end);
    total_gmp = get_elapsed_time(start, end);
  } while (total_gmp < target_sec);

  // Validate correctness before reporting
  mpz_t expected;
  mpz_init(expected);
  mpz_add_ui(expected, mpz_a, (unsigned long)b);
  assert_match("bn_add_u64 (benchmark)", "(bench input)", b, &bn_res, expected);
  mpz_clear(expected);

  char size_info[32];
  snprintf(size_info, sizeof(size_info), "%d bits", bits);
  print_table_row("bn_add_u64", size_info, total_custom / ops_custom,
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
  benchmark_add_u64(4096, bench_budget(4096), state);
  benchmark_add_u64(65536, bench_budget(65536), state);
  benchmark_add_u64(1000000, bench_budget(1000000), state);
  print_table_footer();

  gmp_randclear(state);
  return 0;
}
