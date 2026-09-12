/*
 * test_gcd_extended.c
 *
 * Unified GMP-verified test suite for bn_gcd_extended.
 *
 * Verified against mpz_gcdext:
 *   - d must equal the GMP gcd (non-negative)
 *   - the Bezout identity u*a + v*b == d must hold (computed with the
 *     custom library)
 *
 *   1. Edge cases: (0,0), (0,x), (x,0), 1, -1, 2^64-1, 2^64 boundaries
 *   2. 1000 randomized cases (1..4096 bits) vs mpz_gcdext
 *   3. Time-based benchmark vs GMP at 512/1024/2048/4096 bits
 *
 * (No pointer-aliasing cases: a and b are const inputs.)
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

static void fail(const char* op, const char* a_str, const char* b_str,
                 const char* what, const char* gmp_val, const char* custom_val)
{
  fprintf(stderr, "\n[FATAL ERROR] Correctness failure in %s!\n", op);
  fprintf(stderr, "Input A: %s\n", a_str ? a_str : "(n/a)");
  fprintf(stderr, "Input B: %s\n", b_str ? b_str : "(n/a)");
  fprintf(stderr, "%s GMP:    %s\n", what, gmp_val);
  fprintf(stderr, "%s Custom: %s\n", what, custom_val);
  exit(EXIT_FAILURE);
}

static void assert_bn_eq(const char* op, const char* a_str, const char* b_str,
                         const bignum* bn, const char* what, mpz_t expected)
{
  char* gmp_str = mpz_get_str(NULL, 10, expected);
  char* custom_str = bn_to_string(bn);
  if (strcmp(gmp_str, custom_str) != 0) {
    fail(op, a_str, b_str, what, gmp_str, custom_str);
  }
  free(gmp_str);
  free(custom_str);
}

// d vs gmp gcd; then Bezout identity u*a + v*b == d using custom ops
static void assert_gcd_ext(const char* op, const char* a_str, const char* b_str,
                           bignum* u, bignum* v, const bignum* d, bignum* a,
                           bignum* b, mpz_t gg)
{
  assert_bn_eq(op, a_str, b_str, d, "GCD", gg);

  bignum t, acc;
  bn_init_multi(&t, &acc, NULL);
  bn_mul(&t, u, a);
  bn_mul(&acc, v, b);
  bn_add(&t, &t, &acc);  // u*a + v*b
  assert_bn_eq(op, a_str, b_str, &t, "U*A+V*B", gg);
  bn_free_multi(&t, &acc, NULL);
}

static char* random_mpz_str(mpz_t z, int bits, gmp_randstate_t state)
{
  mpz_urandomb(z, state, bits);
  if (rand() & 1) mpz_neg(z, z);
  return mpz_get_str(NULL, 10, z);
}

// =============================================================================
// PART 1: EDGE CASES
// =============================================================================

static void run_case(const char* name, const char* a_str, const char* b_str)
{
  bignum a, b, u, v, d;
  mpz_t za, zb, gg, gx, gy;

  bn_init_multi(&a, &b, &u, &v, &d, NULL);
  bn_init_val(&a, a_str);
  bn_init_val(&b, b_str);

  mpz_inits(za, zb, gg, gx, gy, NULL);
  mpz_set_str(za, a_str, 10);
  mpz_set_str(zb, b_str, 10);
  mpz_gcdext(gg, gx, gy, za, zb);

  char detail[256];
  snprintf(detail, sizeof(detail), "bn_gcd_extended [%s]", name);

  bn_gcd_extended(&u, &v, &d, &a, &b);
  assert_gcd_ext(detail, a_str, b_str, &u, &v, &d, &a, &b, gg);

  bn_free_multi(&a, &b, &u, &v, &d, NULL);
  mpz_clears(za, zb, gg, gx, gy, NULL);
}

static void run_edge_cases()
{
  printf("\n--- bn_gcd_extended: edge cases ---\n");

  run_case("both zero (gcd = 0)", "0", "0");
  run_case("zero + one", "0", "1");
  run_case("one + zero", "1", "0");
  run_case("one + one", "1", "1");
  run_case("classic 1071/462", "1071", "462");
  run_case("negative pair", "-1071", "-462");
  run_case("mixed signs", "-1071", "462");
  run_case("u64-max pair", "18446744073709551615", "18446744073709551615");
  run_case("u64-max + 2^64", "18446744073709551615", "18446744073709551616");
  run_case("2^64 pair", "18446744073709551616", "18446744073709551616");
  run_case("2^64 + 2^63", "18446744073709551616", "9223372036854775808");

  printf("all edge cases passed\n");
}

// =============================================================================
// PART 2: RANDOMIZED FUZZING vs GMP
// =============================================================================

#define FUZZ_ITERATIONS 1000

static void run_random(gmp_randstate_t state)
{
  printf("\n--- bn_gcd_extended: %d randomized cases vs mpz_gcdext ---\n",
         FUZZ_ITERATIONS);

  for (int i = 0; i < FUZZ_ITERATIONS; i++) {
    int bits_a = 1 + (rand() % 4096);
    int bits_b = 1 + (rand() % 4096);

    bignum a, b, u, v, d;
    mpz_t za, zb, gg, gx, gy;
    char* sa;
    char* sb;
    char detail[64];

    bn_init_multi(&a, &b, &u, &v, &d, NULL);
    mpz_inits(za, zb, gg, gx, gy, NULL);

    sa = random_mpz_str(za, bits_a, state);
    sb = random_mpz_str(zb, bits_b, state);

    bn_init_val(&a, sa);
    bn_init_val(&b, sb);

    mpz_gcdext(gg, gx, gy, za, zb);

    bn_gcd_extended(&u, &v, &d, &a, &b);

    snprintf(detail, sizeof(detail), "case %d (a=%d bits, b=%d bits)", i,
             bits_a, bits_b);
    assert_gcd_ext(detail, sa, sb, &u, &v, &d, &a, &b, gg);

    free(sa);
    free(sb);
    bn_free_multi(&a, &b, &u, &v, &d, NULL);
    mpz_clears(za, zb, gg, gx, gy, NULL);
  }

  printf("all %d random cases passed\n", FUZZ_ITERATIONS);
}

// =============================================================================
// PART 3: TIME-BASED BENCHMARKS vs GMP
// =============================================================================

static void benchmark_gcd_extended(int bits, double target_sec,
                                   gmp_randstate_t state)
{
  bignum bn_a, bn_b, bn_u, bn_v, bn_d;
  mpz_t mpz_a, mpz_b, mpz_g, mpz_x, mpz_y;

  bn_init_multi(&bn_a, &bn_b, &bn_u, &bn_v, &bn_d, NULL);
  mpz_inits(mpz_a, mpz_b, mpz_g, mpz_x, mpz_y, NULL);

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
    bn_gcd_extended(&bn_u, &bn_v, &bn_d, &bn_a, &bn_b);
    ops_custom++;
    clock_gettime(CLOCK_MONOTONIC, &end);
    total_custom = get_elapsed_time(start, end);
  } while (total_custom < target_sec);

  clock_gettime(CLOCK_MONOTONIC, &start);
  do {
    mpz_gcdext(mpz_g, mpz_x, mpz_y, mpz_a, mpz_b);
    ops_gmp++;
    clock_gettime(CLOCK_MONOTONIC, &end);
    total_gmp = get_elapsed_time(start, end);
  } while (total_gmp < target_sec);

  // Validate correctness before reporting
  assert_gcd_ext("bn_gcd_extended (benchmark)", "(bench input)",
                 "(bench input)", &bn_u, &bn_v, &bn_d, &bn_a, &bn_b, mpz_g);

  char size_info[32];
  snprintf(size_info, sizeof(size_info), "%d bits", bits);
  print_table_row("bn_gcd_extended", size_info, total_custom / ops_custom,
                  total_gmp / ops_gmp);

  bn_free_multi(&bn_a, &bn_b, &bn_u, &bn_v, &bn_d, NULL);
  mpz_clears(mpz_a, mpz_b, mpz_g, mpz_x, mpz_y, NULL);
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
  benchmark_gcd_extended(2048, bench_budget(2048), state);
  benchmark_gcd_extended(8192, bench_budget(8192), state);
  benchmark_gcd_extended(32768, bench_budget(32768), state);
  benchmark_gcd_extended(65536, bench_budget(65536), state);
  benchmark_gcd_extended(65539, bench_budget(65539), state);
  benchmark_gcd_extended(100000, bench_budget(100000), state);
  benchmark_gcd_extended(200000, bench_budget(200000), state);
  benchmark_gcd_extended(500000, bench_budget(500000), state);
  print_table_footer();

  gmp_randclear(state);
  return 0;
}
