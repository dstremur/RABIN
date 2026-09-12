/*
 * test_mod_u64.c
 *
 * Unified GMP-verified test suite for bn_mod_u64.
 *
 * bn_mod_u64 returns the nonnegative residue in [0, d), so it is verified
 * against mpz_mod_ui (floor-style reduction).
 *
 *   1. Edge cases: 0, negative a, d = 1, max u64 divisor
 *   2. 1000 randomized cases (a: 1..4096 bits signed, d: full 64-bit)
 *   3. Time-based benchmark vs GMP at 512/1024/2048/4096 bits
 *
 * (No pointer-aliasing cases: the function has no output pointer.)
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

void assert_mod_u64(const char* op, const char* a_str, uint64_t d,
                    uint64_t custom, uint64_t expected)
{
  if (custom != expected) {
    fprintf(stderr, "\n[FATAL ERROR] Correctness failure in %s!\n", op);
    fprintf(stderr, "Input A: %s\n", a_str ? a_str : "(n/a)");
    fprintf(stderr, "Input D: %llu\n", (unsigned long long)d);
    fprintf(stderr, "GMP Result:    %llu\n", (unsigned long long)expected);
    fprintf(stderr, "Custom Result: %llu\n", (unsigned long long)custom);
    exit(EXIT_FAILURE);
  }
}

void set_mpz_u64(mpz_t z, uint64_t val)
{
  char buf[32];
  snprintf(buf, sizeof(buf), "%llu", (unsigned long long)val);
  mpz_set_str(z, buf, 10);
}

uint64_t mpz_to_u64(mpz_t z)
{
  char* s = mpz_get_str(NULL, 10, z);
  uint64_t val = strtoull(s, NULL, 10);
  free(s);
  return val;
}

static char* random_a(mpz_t za, int bits, gmp_randstate_t state)
{
  mpz_urandomb(za, state, bits);
  if (rand() & 1) mpz_neg(za, za);
  return mpz_get_str(NULL, 10, za);
}

static uint64_t random_u64_nonzero(gmp_randstate_t state)
{
  mpz_t t;
  mpz_init(t);
  mpz_urandomb(t, state, 63);
  mpz_add_ui(t, t, 1);  // d in [1, 2^63)
  uint64_t v = mpz_get_ui(t);
  mpz_clear(t);
  return v;
}

// =============================================================================
// PART 1: EDGE CASES
// =============================================================================

static void run_case(const char* name, const char* a_str, uint64_t d)
{
  bignum a;
  mpz_t za, zd, zr;

  bn_init(&a);
  bn_init_val(&a, a_str);

  mpz_inits(za, zd, zr, NULL);
  mpz_set_str(za, a_str, 10);
  set_mpz_u64(zd, d);

  // GMP nonnegative residue matches bn_mod_u64
  mpz_mod_ui(zr, za, (unsigned long)d);
  uint64_t gmp_res = mpz_to_u64(zr);
  uint64_t custom_res = bn_mod_u64(&a, d);

  char detail[256];
  snprintf(detail, sizeof(detail), "bn_mod_u64 [%s]", name);
  assert_mod_u64(detail, a_str, d, custom_res, gmp_res);

  bn_free(&a);
  mpz_clears(za, zd, zr, NULL);
}

static void run_edge_cases()
{
  printf("\n--- bn_mod_u64: edge cases ---\n");

  run_case("zero mod", "0", 12345);
  run_case("neg a", "-999999999999999", 7);
  run_case("small mod", "10", 3);
  run_case("d = 1 (always 0)", "12345678901234567890", 1);
  run_case("max u64 divisor", "340282366920938463463374607431768211455",
           18446744073709551615ULL);
  run_case("u64-max mod u64-max (to 0)", "18446744073709551615",
           18446744073709551615ULL);
  run_case("2^64 mod u64-max", "18446744073709551616", 18446744073709551615ULL);
  run_case("neg one mod two (residue 1)", "-1", 2);

  printf("all edge cases passed\n");
}

// =============================================================================
// PART 2: RANDOMIZED FUZZING vs GMP
// =============================================================================

#define FUZZ_ITERATIONS 1000

static void run_random(gmp_randstate_t state)
{
  printf("\n--- bn_mod_u64: %d randomized cases vs mpz_mod_ui ---\n",
         FUZZ_ITERATIONS);

  for (int i = 0; i < FUZZ_ITERATIONS; i++) {
    int bits = 1 + (rand() % 4096);
    uint64_t d = random_u64_nonzero(state);

    bignum a;
    mpz_t za, zd, zr;
    char* sa;
    char detail[64];

    bn_init(&a);
    mpz_inits(za, zd, zr, NULL);

    sa = random_a(za, bits, state);
    bn_init_val(&a, sa);
    set_mpz_u64(zd, d);

    mpz_mod_ui(zr, za, (unsigned long)d);
    uint64_t gmp_res = mpz_to_u64(zr);
    uint64_t custom_res = bn_mod_u64(&a, d);

    snprintf(detail, sizeof(detail), "case %d (a=%d bits)", i, bits);
    assert_mod_u64(detail, sa, d, custom_res, gmp_res);

    free(sa);
    bn_free(&a);
    mpz_clears(za, zd, zr, NULL);
  }

  printf("all %d random cases passed\n", FUZZ_ITERATIONS);
}

// =============================================================================
// PART 3: TIME-BASED BENCHMARKS vs GMP
// =============================================================================

static void benchmark_mod_u64(int bits, double target_sec,
                              gmp_randstate_t state)
{
  uint64_t d = random_u64_nonzero(state);

  bignum bn_a;
  mpz_t mpz_a, mpz_d, mpz_res;

  bn_init(&bn_a);
  mpz_inits(mpz_a, mpz_d, mpz_res, NULL);

  char* sa = random_a(mpz_a, bits, state);
  bn_init_val(&bn_a, sa);
  free(sa);
  set_mpz_u64(mpz_d, d);

  struct timespec start, end;
  int ops_custom = 0, ops_gmp = 0;
  double total_custom = 0, total_gmp = 0;
  volatile uint64_t v_rem;

  clock_gettime(CLOCK_MONOTONIC, &start);
  do {
    v_rem = bn_mod_u64(&bn_a, d);
    ops_custom++;
    clock_gettime(CLOCK_MONOTONIC, &end);
    total_custom = get_elapsed_time(start, end);
  } while (total_custom < target_sec);

  clock_gettime(CLOCK_MONOTONIC, &start);
  do {
    mpz_mod_ui(mpz_res, mpz_a, (unsigned long)d);
    ops_gmp++;
    clock_gettime(CLOCK_MONOTONIC, &end);
    total_gmp = get_elapsed_time(start, end);
  } while (total_gmp < target_sec);

  // Validate correctness before reporting
  assert_mod_u64("bn_mod_u64 (benchmark)", "(bench input)", d, v_rem,
                 mpz_to_u64(mpz_res));

  char size_info[32];
  snprintf(size_info, sizeof(size_info), "%d bits", bits);
  print_table_row("bn_mod_u64", size_info, total_custom / ops_custom,
                  total_gmp / ops_gmp);

  bn_free(&bn_a);
  mpz_clears(mpz_a, mpz_d, mpz_res, NULL);
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
  benchmark_mod_u64(4096, bench_budget(4096), state);
  benchmark_mod_u64(65536, bench_budget(65536), state);
  benchmark_mod_u64(1000000, bench_budget(1000000), state);
  print_table_footer();

  gmp_randclear(state);
  return 0;
}
