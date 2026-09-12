/*
 * test_mul_i64.c
 *
 * Unified GMP-verified test suite for bn_mul_i64.
 *
 * GMP's mpz_mul_si only takes an int, so the 64-bit operand is built as an
 * mpz via its decimal string for the reference computation.
 *
 *   1. Edge cases: 0, 1, -1, INT64_MAX/INT64_MIN, 2^64-1, 2^64 boundaries
 *      + pointer aliasing (r == a)
 *   2. 1000 randomized cases (a: 1..4096 bits, c: full 64-bit) vs mpz_mul
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

void mpz_set_i64(mpz_t z, int64_t c)
{
  char buf[32];
  snprintf(buf, sizeof(buf), "%lld", (long long)c);
  mpz_set_str(z, buf, 10);
}

void assert_match(const char* op, const char* a_str, int64_t c,
                  const bignum* bn, mpz_t expected)
{
  char* gmp_str = mpz_get_str(NULL, 10, expected);
  char* custom_str = bn_to_string(bn);

  if (strcmp(gmp_str, custom_str) != 0) {
    fprintf(stderr, "\n[FATAL ERROR] Correctness failure in %s!\n", op);
    fprintf(stderr, "Input A: %s\n", a_str ? a_str : "(n/a)");
    fprintf(stderr, "Input C: %lld\n", (long long)c);
    fprintf(stderr, "GMP Result:    %s\n", gmp_str);
    fprintf(stderr, "Custom Result: %s\n", custom_str);
    free(gmp_str);
    free(custom_str);
    exit(EXIT_FAILURE);
  }

  free(gmp_str);
  free(custom_str);
}

static char* random_a(mpz_t za, int bits, gmp_randstate_t state)
{
  mpz_urandomb(za, state, bits);
  if (rand() & 1) mpz_neg(za, za);
  return mpz_get_str(NULL, 10, za);
}

static int64_t random_i64(gmp_randstate_t state)
{
  mpz_t t;
  mpz_init(t);
  mpz_urandomb(t, state, 64);
  uint64_t v = (uint64_t)mpz_get_ui(t);
  mpz_clear(t);
  return (int64_t)v;
}

// =============================================================================
// PART 1: EDGE CASES (incl. pointer aliasing)
// =============================================================================

static void run_case(const char* name, const char* a_str, int64_t c)
{
  bignum a, res;
  mpz_t za, zc, zr;
  char detail[256];

  bn_init_multi(&a, &res, NULL);
  mpz_inits(za, zc, zr, NULL);
  bn_init_val(&a, a_str);
  mpz_set_str(za, a_str, 10);
  mpz_set_i64(zc, c);

  snprintf(detail, sizeof(detail), "bn_mul_i64 [%s] a=%s c=%lld", name, a_str,
           (long long)c);

  // non-aliased
  bn_mul_i64(&res, &a, c);
  mpz_mul(zr, za, zc);
  assert_match(detail, a_str, c, &res, zr);

  // aliasing: r == a
  bn_init_val(&a, a_str);
  bn_mul_i64(&a, &a, c);
  assert_match("bn_mul_i64 (r==a)", a_str, c, &a, zr);

  bn_free_multi(&a, &res, NULL);
  mpz_clears(za, zc, zr, NULL);
}

static void run_edge_cases()
{
  printf("\n--- bn_mul_i64: edge cases ---\n");

  run_case("zero * anything", "0", 12345);
  run_case("one * one", "1", 1);
  run_case("neg-one * neg-one", "-1", -1);
  run_case("u64-max * 2", "18446744073709551615", 2);
  run_case("u64-max * -1", "18446744073709551615", -1);
  run_case("2^64 * 3", "18446744073709551616", 3);
  run_case("a * 0", "12345678901234567890", 0);
  run_case("a * INT64_MAX", "12345678901234567890", INT64_MAX);
  run_case("a * INT64_MIN", "12345678901234567890", INT64_MIN);
  run_case("neg a * INT64_MIN", "-18446744073709551615", INT64_MIN);

  printf("all edge cases passed\n");
}

// =============================================================================
// PART 2: RANDOMIZED FUZZING vs GMP
// =============================================================================

#define FUZZ_ITERATIONS 1000

static void run_random(gmp_randstate_t state)
{
  printf("\n--- bn_mul_i64: %d randomized cases vs mpz_mul ---\n",
         FUZZ_ITERATIONS);

  for (int i = 0; i < FUZZ_ITERATIONS; i++) {
    int bits = 1 + (rand() % 4096);
    int64_t c = random_i64(state);

    bignum a, res;
    mpz_t za, zc, zr;
    char* sa;
    char detail[64];

    bn_init_multi(&a, &res, NULL);
    mpz_inits(za, zc, zr, NULL);

    sa = random_a(za, bits, state);
    bn_init_val(&a, sa);
    mpz_set_i64(zc, c);

    bn_mul_i64(&res, &a, c);
    mpz_mul(zr, za, zc);

    snprintf(detail, sizeof(detail), "case %d (a=%d bits)", i, bits);
    assert_match(detail, sa, c, &res, zr);

    free(sa);
    bn_free_multi(&a, &res, NULL);
    mpz_clears(za, zc, zr, NULL);
  }

  printf("all %d random cases passed\n", FUZZ_ITERATIONS);
}

// =============================================================================
// PART 3: TIME-BASED BENCHMARKS vs GMP
// =============================================================================

static void benchmark_mul_i64(int bits, double target_sec,
                              gmp_randstate_t state)
{
  const int64_t c = 0x123456789abcdef0LL;
  bignum bn_a, bn_res;
  mpz_t mpz_a, mpz_c, mpz_res;

  bn_init_multi(&bn_a, &bn_res, NULL);
  mpz_inits(mpz_a, mpz_c, mpz_res, NULL);

  char* sa = random_a(mpz_a, bits, state);
  bn_init_val(&bn_a, sa);
  free(sa);
  mpz_set_i64(mpz_c, c);

  struct timespec start, end;
  int ops_custom = 0, ops_gmp = 0;
  double total_custom = 0, total_gmp = 0;

  clock_gettime(CLOCK_MONOTONIC, &start);
  do {
    bn_mul_i64(&bn_res, &bn_a, c);
    ops_custom++;
    clock_gettime(CLOCK_MONOTONIC, &end);
    total_custom = get_elapsed_time(start, end);
  } while (total_custom < target_sec);

  clock_gettime(CLOCK_MONOTONIC, &start);
  do {
    mpz_mul(mpz_res, mpz_a, mpz_c);
    ops_gmp++;
    clock_gettime(CLOCK_MONOTONIC, &end);
    total_gmp = get_elapsed_time(start, end);
  } while (total_gmp < target_sec);

  // Validate correctness before reporting
  assert_match("bn_mul_i64 (benchmark)", "(bench input)", c, &bn_res, mpz_res);

  char size_info[32];
  snprintf(size_info, sizeof(size_info), "%d bits", bits);
  print_table_row("bn_mul_i64", size_info, total_custom / ops_custom,
                  total_gmp / ops_gmp);

  bn_free_multi(&bn_a, &bn_res, NULL);
  mpz_clears(mpz_a, mpz_c, mpz_res, NULL);
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
  benchmark_mul_i64(4096, bench_budget(4096), state);
  benchmark_mul_i64(65536, bench_budget(65536), state);
  benchmark_mul_i64(1000000, bench_budget(1000000), state);
  print_table_footer();

  gmp_randclear(state);
  return 0;
}
