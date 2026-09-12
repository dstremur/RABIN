/*
 * test_divmod_u64.c
 *
 * Unified GMP-verified test suite for bn_divmod_u64.
 *
 * The quotient follows truncated (C) semantics like mpz_tdiv_qr_ui; the
 * returned remainder is the magnitude |a mod d|, so the GMP reference
 * remainder is absolute-valued before comparison.
 *
 *   1. Edge cases: 0, 1, -1, d = 1, max u64 divisor + pointer aliasing
 *      (q == a)
 *   2. 1000 randomized cases (a: 1..4096 bits signed, d: full 64-bit)
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

// q compared as string, remainder as magnitude |r|
static void assert_divmod_u64(const char* op, const char* a_str, uint64_t d,
                              const bignum* q, uint64_t custom_rem, mpz_t gq,
                              mpz_t gr)
{
  char* gmp_q = mpz_get_str(NULL, 10, gq);
  char* custom_q = bn_to_string(q);
  mpz_abs(gr, gr);
  uint64_t gmp_rem = mpz_to_u64(gr);

  int ok = (custom_rem == gmp_rem) && (strcmp(gmp_q, custom_q) == 0);
  if (!ok) {
    fprintf(stderr, "\n[FATAL ERROR] Correctness failure in %s!\n", op);
    fprintf(stderr, "Input A: %s\n", a_str ? a_str : "(n/a)");
    fprintf(stderr, "Input D: %llu\n", (unsigned long long)d);
    fprintf(stderr, "Expected Q: %s, Rem: %llu\n", gmp_q,
            (unsigned long long)gmp_rem);
    fprintf(stderr, "Got Q:      %s, Rem: %llu\n", custom_q,
            (unsigned long long)custom_rem);
    exit(EXIT_FAILURE);
  }

  free(gmp_q);
  free(custom_q);
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
// PART 1: EDGE CASES (incl. pointer aliasing)
// =============================================================================

static void run_case(const char* name, const char* a_str, uint64_t d)
{
  bignum a, q, alias;
  mpz_t za, zd, gq, gr;
  char detail[256];

  bn_init_multi(&a, &q, &alias, NULL);
  bn_init_val(&a, a_str);
  bn_init_val(&alias, a_str);

  mpz_inits(za, zd, gq, gr, NULL);
  mpz_set_str(za, a_str, 10);
  set_mpz_u64(zd, d);

  // GMP truncated division with a 64-bit divisor
  mpz_tdiv_qr_ui(gq, gr, za, (unsigned long)d);

  snprintf(detail, sizeof(detail), "bn_divmod_u64 [%s]", name);

  // Test 1: Non-aliased
  uint64_t custom_rem = bn_divmod_u64(&q, &a, d);
  assert_divmod_u64(detail, a_str, d, &q, custom_rem, gq, gr);

  // Test 2: Aliased (q == a)
  uint64_t alias_rem = bn_divmod_u64(&alias, &alias, d);
  assert_divmod_u64("bn_divmod_u64 (q==a)", a_str, d, &alias, alias_rem, gq,
                    gr);

  bn_free_multi(&a, &q, &alias, NULL);
  mpz_clears(za, zd, gq, gr, NULL);
}

static void run_edge_cases()
{
  printf("\n--- bn_divmod_u64: edge cases ---\n");

  run_case("zero div", "0", 12345);
  run_case("neg a", "-999999999999999", 7);
  run_case("divide by 1", "9876543210123456789", 1);
  run_case("neg one (q=0)", "-1", 2);
  run_case("max u64 divisor", "340282366920938463463374607431768211455",
           18446744073709551615ULL);
  run_case("u64-max div u64-max (q=1, r=0)", "18446744073709551615",
           18446744073709551615ULL);
  run_case("2^64 div u64-max (q=1, r=1)", "18446744073709551616",
           18446744073709551615ULL);
  run_case("u64-max div 2", "18446744073709551615", 2);

  printf("all edge cases passed\n");
}

// =============================================================================
// PART 2: RANDOMIZED FUZZING vs GMP
// =============================================================================

#define FUZZ_ITERATIONS 1000

static void run_random(gmp_randstate_t state)
{
  printf("\n--- bn_divmod_u64: %d randomized cases vs mpz_tdiv_qr_ui ---\n",
         FUZZ_ITERATIONS);

  for (int i = 0; i < FUZZ_ITERATIONS; i++) {
    int bits = 1 + (rand() % 4096);
    uint64_t d = random_u64_nonzero(state);

    bignum a, q, alias;
    mpz_t za, zd, gq, gr;
    char* sa;
    char detail[64];

    bn_init_multi(&a, &q, &alias, NULL);
    mpz_inits(za, zd, gq, gr, NULL);

    sa = random_a(za, bits, state);
    bn_init_val(&a, sa);
    bn_init_val(&alias, sa);
    set_mpz_u64(zd, d);

    mpz_tdiv_qr_ui(gq, gr, za, (unsigned long)d);

    // non-aliased
    uint64_t custom_rem = bn_divmod_u64(&q, &a, d);
    snprintf(detail, sizeof(detail), "case %d (a=%d bits)", i, bits);
    assert_divmod_u64(detail, sa, d, &q, custom_rem, gq, gr);

    // aliased (q == a) every 10th case
    if (i % 10 == 0) {
      uint64_t alias_rem = bn_divmod_u64(&alias, &alias, d);
      assert_divmod_u64("case (q==a)", sa, d, &alias, alias_rem, gq, gr);
    }

    free(sa);
    bn_free_multi(&a, &q, &alias, NULL);
    mpz_clears(za, zd, gq, gr, NULL);
  }

  printf("all %d random cases passed\n", FUZZ_ITERATIONS);
}

// =============================================================================
// PART 3: TIME-BASED BENCHMARKS vs GMP
// =============================================================================

static void benchmark_divmod_u64(int bits, double target_sec,
                                 gmp_randstate_t state)
{
  bignum bn_a, bn_q;
  mpz_t mpz_a, mpz_d, mpz_q, mpz_rem;

  bn_init_multi(&bn_a, &bn_q, NULL);
  mpz_inits(mpz_a, mpz_d, mpz_q, mpz_rem, NULL);

  char* sa = random_a(mpz_a, bits, state);
  bn_init_val(&bn_a, sa);
  free(sa);

  uint64_t d = random_u64_nonzero(state);
  set_mpz_u64(mpz_d, d);

  struct timespec start, end;
  int ops_custom = 0, ops_gmp = 0;
  double total_custom = 0, total_gmp = 0;
  volatile uint64_t v_rem;

  clock_gettime(CLOCK_MONOTONIC, &start);
  do {
    v_rem = bn_divmod_u64(&bn_q, &bn_a, d);
    ops_custom++;
    clock_gettime(CLOCK_MONOTONIC, &end);
    total_custom = get_elapsed_time(start, end);
  } while (total_custom < target_sec);

  clock_gettime(CLOCK_MONOTONIC, &start);
  do {
    mpz_tdiv_qr_ui(mpz_q, mpz_rem, mpz_a, (unsigned long)d);
    ops_gmp++;
    clock_gettime(CLOCK_MONOTONIC, &end);
    total_gmp = get_elapsed_time(start, end);
  } while (total_gmp < target_sec);

  // Validate correctness before reporting
  assert_divmod_u64("bn_divmod_u64 (benchmark)", "(bench input)", d, &bn_q,
                    v_rem, mpz_q, mpz_rem);

  char size_info[32];
  snprintf(size_info, sizeof(size_info), "%d bits", bits);
  print_table_row("bn_divmod_u64", size_info, total_custom / ops_custom,
                  total_gmp / ops_gmp);

  bn_free_multi(&bn_a, &bn_q, NULL);
  mpz_clears(mpz_a, mpz_d, mpz_q, mpz_rem, NULL);
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
  benchmark_divmod_u64(4096, bench_budget(4096), state);
  benchmark_divmod_u64(65536, bench_budget(65536), state);
  benchmark_divmod_u64(1000000, bench_budget(1000000), state);
  print_table_footer();

  gmp_randclear(state);
  return 0;
}
