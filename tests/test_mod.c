/*
 * test_mod.c
 *
 * Unified GMP-verified test suite for rz_mod (truncated remainder, C
 * semantics: sign(r) = sign(a)).
 *
 * Verified against mpz_tdiv_r (NOT mpz_mod, which is floor remainder).
 *
 *   1. Edge cases: 0, 1, -1, 2^64-1, 2^64 boundaries, negative divisor +
 *      pointer aliasing (r == a and r == b)
 *   2. 1000 randomized cases (1..4096 bits) vs mpz_tdiv_r
 *   3. Time-based benchmark vs GMP at 512/1024/2048/4096 bits
 *
 * Copyright (C) 2026 Diego Strebel
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <gmp.h>
#include <stdbool.h>
#include <stdint.h>
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

void assert_match(const char* op, const char* a_str, const char* b_str,
                  const rz_t* bn, mpz_t expected)
{
  char* gmp_str = mpz_get_str(NULL, 10, expected);
  char* custom_str = rz_to_string(bn);

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

static char* random_mpz_str(mpz_t z, int bits, gmp_randstate_t state)
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
  rz_t a, b, r;
  mpz_t za, zb, zr;
  char detail[256];

  rz_init_multi(&a, &b, &r, NULL);
  mpz_inits(za, zb, zr, NULL);
  rz_init_val(&a, a_str);
  rz_init_val(&b, b_str);
  mpz_set_str(za, a_str, 10);
  mpz_set_str(zb, b_str, 10);

  snprintf(detail, sizeof(detail), "rz_mod [%s] a=%s b=%s", name, a_str, b_str);

  // non-aliased
  rz_mod(&r, &a, &b);
  mpz_tdiv_r(zr, za, zb);
  assert_match(detail, a_str, b_str, &r, zr);

  // aliasing: r == a
  rz_init_val(&a, a_str);
  rz_mod(&a, &a, &b);
  assert_match("rz_mod (r==a)", a_str, b_str, &a, zr);

  // aliasing: r == b
  rz_init_val(&a, a_str);
  rz_init_val(&b, b_str);
  rz_mod(&b, &a, &b);
  assert_match("rz_mod (r==b)", a_str, b_str, &b, zr);

  rz_clear_multi(&a, &b, &r, NULL);
  mpz_clears(za, zb, zr, NULL);
}

static void run_edge_cases()
{
  printf("\n--- rz_mod: edge cases ---\n");

  run_case("zero mod five", "0", "5");
  run_case("one mod one (to 0)", "1", "1");
  run_case("seven mod two", "7", "2");
  run_case("neg seven mod two (trunc: -1, not +1)", "-7", "2");
  run_case("seven mod neg two (trunc: +1)", "7", "-2");
  run_case("neg seven mod neg two (trunc: -1)", "-7", "-2");
  run_case("neg one mod two (trunc: -1)", "-1", "2");
  run_case("a < b (r = a)", "1", "18446744073709551615");
  run_case("a == b (r = 0)", "18446744073709551615", "18446744073709551615");
  run_case("neg a < |b| (r = a)", "-18446744073709551615",
           "18446744073709551615");
  run_case("u64-max mod 3", "18446744073709551615", "3");
  run_case("2^64 mod 2^64 (to 0)", "18446744073709551616",
           "18446744073709551616");
  run_case("huge mod 2 (asymmetric)", "123456789012345678901234567890", "2");

  printf("all edge cases passed\n");
}

// =============================================================================
// PART 2: RANDOMIZED FUZZING vs GMP
// =============================================================================

#define FUZZ_ITERATIONS 1000

static void run_random(gmp_randstate_t state)
{
  printf("\n--- rz_mod: %d randomized cases vs mpz_tdiv_r ---\n",
         FUZZ_ITERATIONS);

  for (int i = 0; i < FUZZ_ITERATIONS; i++) {
    int bits_a = 1 + (rand() % 4096);
    int bits_b = 1 + (rand() % 4096);

    rz_t a, b, r;
    mpz_t za, zb, zr;
    char* sa;
    char* sb;
    char detail[64];

    rz_init_multi(&a, &b, &r, NULL);
    mpz_inits(za, zb, zr, NULL);

    sa = random_mpz_str(za, bits_a, state);
    sb = random_mpz_str(zb, bits_b, state);

    // divisor must be nonzero
    if (mpz_cmp_ui(zb, 0) == 0) {
      mpz_set_ui(zb, 1);
      free(sb);
      sb = mpz_get_str(NULL, 10, zb);
    }

    rz_init_val(&a, sa);
    rz_init_val(&b, sb);

    rz_mod(&r, &a, &b);
    mpz_tdiv_r(zr, za, zb);

    snprintf(detail, sizeof(detail), "case %d (a=%d bits, b=%d bits)", i,
             bits_a, bits_b);
    assert_match(detail, sa, sb, &r, zr);

    free(sa);
    free(sb);
    rz_clear_multi(&a, &b, &r, NULL);
    mpz_clears(za, zb, zr, NULL);
  }

  printf("all %d random cases passed\n", FUZZ_ITERATIONS);
}

// =============================================================================
// PART 3: TIME-BASED BENCHMARKS vs GMP
// =============================================================================

static void benchmark_mod(int bits, double target_sec, gmp_randstate_t state)
{
  rz_t rz_a, rz_b, rz_r;
  mpz_t mpz_a, mpz_b, mpz_r;

  rz_init_multi(&rz_a, &rz_b, &rz_r, NULL);
  mpz_inits(mpz_a, mpz_b, mpz_r, NULL);

  char* sa = random_mpz_str(mpz_a, bits, state);
  char* sb = random_mpz_str(mpz_b, bits, state);
  rz_init_val(&rz_a, sa);
  rz_init_val(&rz_b, sb);
  free(sa);
  free(sb);

  struct timespec start, end;
  int ops_custom = 0, ops_gmp = 0;
  double total_custom = 0, total_gmp = 0;

  clock_gettime(CLOCK_MONOTONIC, &start);
  do {
    rz_mod(&rz_r, &rz_a, &rz_b);
    ops_custom++;
    clock_gettime(CLOCK_MONOTONIC, &end);
    total_custom = get_elapsed_time(start, end);
  } while (total_custom < target_sec);

  clock_gettime(CLOCK_MONOTONIC, &start);
  do {
    mpz_tdiv_r(mpz_r, mpz_a, mpz_b);
    ops_gmp++;
    clock_gettime(CLOCK_MONOTONIC, &end);
    total_gmp = get_elapsed_time(start, end);
  } while (total_gmp < target_sec);

  // Validate correctness before reporting
  assert_match("rz_mod (benchmark)", "(bench input)", "(bench input)", &rz_r,
               mpz_r);

  char size_info[32];
  snprintf(size_info, sizeof(size_info), "%d bits", bits);
  print_table_row("rz_mod", size_info, total_custom / ops_custom,
                  total_gmp / ops_gmp);

  rz_clear_multi(&rz_a, &rz_b, &rz_r, NULL);
  mpz_clears(mpz_a, mpz_b, mpz_r, NULL);
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
  benchmark_mod(2048, bench_budget(2048), state);
  benchmark_mod(8192, bench_budget(8192), state);
  benchmark_mod(32768, bench_budget(32768), state);
  benchmark_mod(65536, bench_budget(65536), state);
  benchmark_mod(65539, bench_budget(65539), state);
  print_table_footer();

  gmp_randclear(state);
  return 0;
}
