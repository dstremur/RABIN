/*
 * test_cmp.c
 *
 * Unified GMP-verified test suite for bn_cmp.
 *
 *   1. Edge cases: 0, 1, -1, 2^64-1, 2^64 boundaries + pointer aliasing
 *   2. 1000 randomized cases (1..4096 bits) vs mpz_cmp
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

// Both comparators return 1/0/-1; signs must agree exactly.
void assert_cmp_match(const char* op, const char* a_str, const char* b_str,
                      int custom, int gmp)
{
  if ((custom < 0) != (gmp < 0) || (custom > 0) != (gmp > 0) ||
      (custom == 0) != (gmp == 0)) {
    fprintf(stderr, "\n[FATAL ERROR] Correctness failure in %s!\n", op);
    fprintf(stderr, "Input A: %s\n", a_str);
    fprintf(stderr, "Input B: %s\n", b_str);
    fprintf(stderr, "GMP Result:    %d\n", gmp);
    fprintf(stderr, "Custom Result: %d\n", custom);
    exit(EXIT_FAILURE);
  }
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
  bignum a, b;
  mpz_t za, zb;
  char detail[256];

  bn_init_multi(&a, &b, NULL);
  mpz_inits(za, zb, NULL);
  bn_init_val(&a, a_str);
  bn_init_val(&b, b_str);
  mpz_set_str(za, a_str, 10);
  mpz_set_str(zb, b_str, 10);

  snprintf(detail, sizeof(detail), "bn_cmp [%s] a=%s b=%s", name, a_str, b_str);

  assert_cmp_match(detail, a_str, b_str, bn_cmp(&a, &b), mpz_cmp(za, zb));

  bn_free_multi(&a, &b, NULL);
  mpz_clears(za, zb, NULL);
}

static void run_edge_cases()
{
  printf("\n--- bn_cmp: edge cases ---\n");
  run_case("zeros", "0", "0");
  run_case("one vs zero", "1", "0");
  run_case("zero vs one", "0", "1");
  run_case("neg-one vs one", "-1", "1");
  run_case("one vs neg-one", "1", "-1");
  run_case("neg-one vs neg-one", "-1", "-1");
  run_case("u64-max vs u64-max", "18446744073709551615",
           "18446744073709551615");
  run_case("u64-max vs u64-max-1", "18446744073709551615",
           "18446744073709551614");
  run_case("2^64 vs u64-max (cross limb)", "18446744073709551616",
           "18446744073709551615");
  run_case("equal magnitude, opposite sign", "18446744073709551615",
           "-18446744073709551615");
  run_case("same magnitude, both negative", "-12345678901234567890",
           "-12345678901234567890");
  run_case("diff magnitude, both negative", "-123456789012345678901",
           "-12345678901234567890");

  // aliasing: a == b must compare equal (0)
  {
    bignum a;
    mpz_t za;
    bn_init(&a);
    mpz_init(za);
    bn_init_val(&a, "18446744073709551615");
    mpz_set_str(za, "18446744073709551615", 10);
    if (bn_cmp(&a, &a) != 0) {
      fprintf(stderr, "\n[FATAL ERROR] bn_cmp(&a, &a) must be 0, got %d\n",
              bn_cmp(&a, &a));
      exit(EXIT_FAILURE);
    }
    bn_free(&a);
    mpz_clear(za);
    printf("aliasing (a==b) passed\n");
  }

  printf("all edge cases passed\n");
}

// =============================================================================
// PART 2: RANDOMIZED FUZZING vs GMP
// =============================================================================

#define FUZZ_ITERATIONS 1000

static void run_random(gmp_randstate_t state)
{
  printf("\n--- bn_cmp: %d randomized cases vs mpz_cmp ---\n", FUZZ_ITERATIONS);

  for (int i = 0; i < FUZZ_ITERATIONS; i++) {
    int bits_a = 1 + (rand() % 4096);
    int bits_b = 1 + (rand() % 4096);

    bignum a, b;
    mpz_t za, zb;
    char* sa;
    char* sb;
    char detail[128];

    bn_init_multi(&a, &b, NULL);
    mpz_inits(za, zb, NULL);

    sa = random_mpz_str(za, bits_a, state);
    sb = random_mpz_str(zb, bits_b, state);
    bn_init_val(&a, sa);
    bn_init_val(&b, sb);

    snprintf(detail, sizeof(detail), "case %d (a=%d bits, b=%d bits)", i,
             bits_a, bits_b);
    assert_cmp_match(detail, sa, sb, bn_cmp(&a, &b), mpz_cmp(za, zb));

    free(sa);
    free(sb);
    bn_free_multi(&a, &b, NULL);
    mpz_clears(za, zb, NULL);
  }

  printf("all %d random cases passed\n", FUZZ_ITERATIONS);
}

// =============================================================================
// PART 3: TIME-BASED BENCHMARKS vs GMP
// =============================================================================

static void benchmark_cmp(int bits, double target_sec, gmp_randstate_t state)
{
  bignum bn_a, bn_b;
  mpz_t mpz_a, mpz_b;

  bn_init_multi(&bn_a, &bn_b, NULL);
  mpz_inits(mpz_a, mpz_b, NULL);

  char* sa = random_mpz_str(mpz_a, bits, state);
  char* sb = random_mpz_str(mpz_b, bits, state);
  bn_init_val(&bn_a, sa);
  bn_init_val(&bn_b, sb);
  free(sa);
  free(sb);

  struct timespec start, end;
  int ops_custom = 0, ops_gmp = 0;
  double total_custom = 0, total_gmp = 0;

  volatile int last_custom = 0, last_gmp = 0;

  clock_gettime(CLOCK_MONOTONIC, &start);
  do {
    last_custom = bn_cmp(&bn_a, &bn_b);
    ops_custom++;
    clock_gettime(CLOCK_MONOTONIC, &end);
    total_custom = get_elapsed_time(start, end);
  } while (total_custom < target_sec);

  clock_gettime(CLOCK_MONOTONIC, &start);
  do {
    last_gmp = mpz_cmp(mpz_a, mpz_b);
    ops_gmp++;
    clock_gettime(CLOCK_MONOTONIC, &end);
    total_gmp = get_elapsed_time(start, end);
  } while (total_gmp < target_sec);

  // Validate correctness before reporting
  if ((last_custom < 0) != (last_gmp < 0) ||
      (last_custom > 0) != (last_gmp > 0)) {
    fprintf(stderr, "\n[FATAL ERROR] bn_cmp mismatch in benchmark!\n");
    exit(EXIT_FAILURE);
  }

  char size_info[32];
  snprintf(size_info, sizeof(size_info), "%d bits", bits);
  print_table_row("bn_cmp", size_info, total_custom / ops_custom,
                  total_gmp / ops_gmp);

  bn_free_multi(&bn_a, &bn_b, NULL);
  mpz_clears(mpz_a, mpz_b, NULL);
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
  benchmark_cmp(4096, bench_budget(4096), state);
  benchmark_cmp(65536, bench_budget(65536), state);
  benchmark_cmp(1000000, bench_budget(1000000), state);
  print_table_footer();

  gmp_randclear(state);
  return 0;
}
