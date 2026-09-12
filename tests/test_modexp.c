/*
 * test_modexp.c
 *
 * Unified GMP-verified test suite for bn_mod_exp (a^b mod m).
 *
 * Operands are non-negative: bn_mod_exp does not normalize negative bases
 * into [0, m), which is outside its contract (GMP's mpz_powm would).
 *
 *   1. Edge cases: 0^0, e = 0, e = 1, m = 1, m = 2, odd/even moduli,
 *      2^64-1 / 2^64 boundaries + pointer aliasing (r == a, r == b)
 *   2. 1000 randomized cases (a, m: 1..4096 bits, e: 1..512 bits) vs
 *      mpz_powm
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

void assert_match(const char* op, const char* a_str, const char* e_str,
                  const char* m_str, const bignum* bn, mpz_t expected)
{
  char* gmp_str = mpz_get_str(NULL, 10, expected);
  char* custom_str = bn_to_string(bn);

  if (strcmp(gmp_str, custom_str) != 0) {
    fprintf(stderr, "\n[FATAL ERROR] Correctness failure in %s!\n", op);
    fprintf(stderr, "Input A: %s\n", a_str ? a_str : "(n/a)");
    fprintf(stderr, "Input E: %s\n", e_str ? e_str : "(n/a)");
    fprintf(stderr, "Input M: %s\n", m_str ? m_str : "(n/a)");
    fprintf(stderr, "GMP Result:    %s\n", gmp_str);
    fprintf(stderr, "Custom Result: %s\n", custom_str);
    free(gmp_str);
    free(custom_str);
    exit(EXIT_FAILURE);
  }

  free(gmp_str);
  free(custom_str);
}

// random nonnegative mpz of given bits + decimal string
static char* random_mpz_str(mpz_t z, int bits, gmp_randstate_t state)
{
  mpz_urandomb(z, state, bits);
  return mpz_get_str(NULL, 10, z);
}

// =============================================================================
// PART 1: EDGE CASES (incl. pointer aliasing)
// =============================================================================

static void run_case(const char* name, const char* a_str, const char* e_str,
                     const char* m_str)
{
  bignum a, e, m, r;
  mpz_t za, ze, zm, zr;
  char detail[512];

  bn_init_multi(&a, &e, &m, &r, NULL);
  mpz_inits(za, ze, zm, zr, NULL);
  bn_init_val(&a, a_str);
  bn_init_val(&e, e_str);
  bn_init_val(&m, m_str);
  mpz_set_str(za, a_str, 10);
  mpz_set_str(ze, e_str, 10);
  mpz_set_str(zm, m_str, 10);

  snprintf(detail, sizeof(detail), "bn_mod_exp [%s] a=%s e=%s m=%s", name,
           a_str, e_str, m_str);

  // non-aliased
  bn_mod_exp(&r, &a, &e, &m);
  mpz_powm(zr, za, ze, zm);
  assert_match(detail, a_str, e_str, m_str, &r, zr);

  // aliasing: r == a
  bn_init_val(&a, a_str);
  bn_mod_exp(&a, &a, &e, &m);
  assert_match("bn_mod_exp (r==a)", a_str, e_str, m_str, &a, zr);

  // aliasing: r == b (exponent)
  bn_init_val(&a, a_str);
  bn_init_val(&e, e_str);
  bn_mod_exp(&e, &a, &e, &m);
  assert_match("bn_mod_exp (r==e)", a_str, e_str, m_str, &e, zr);

  bn_free_multi(&a, &e, &m, &r, NULL);
  mpz_clears(za, ze, zm, zr, NULL);
}

static void run_edge_cases()
{
  printf("\n--- bn_mod_exp: edge cases ---\n");

  run_case("one ^ zero", "1", "0", "5");
  run_case("zero ^ zero (= 1)", "0", "0", "5");
  run_case("zero ^ one (= 0)", "0", "1", "5");
  run_case("zero ^ big (= 0)", "0", "65537", "5");
  run_case("three ^ zero (= 1)", "3", "0", "2");
  run_case("three ^ five mod two (= 1)", "3", "5", "2");
  run_case("seven ^ one mod three (= 1)", "7", "1", "3");
  run_case("m = 1 (always 0 for e > 0)", "3", "1", "1");
  run_case("zero ^ zero mod one (= 0)", "0", "0", "1");
  run_case("a == m, e > 0 (= 0)", "7", "3", "7");
  run_case("u64-max base, prime modulus", "18446744073709551615", "65",
           "18446744073709551557");
  run_case("even modulus (slow path)", "123456789", "1000",
           "10000000000000000000");
  run_case("2^64 modulus, 2^64-1 base", "18446744073709551615", "2",
           "18446744073709551616");
  run_case("u64-max exponent", "3", "18446744073709551615", "7");

  printf("all edge cases passed\n");
}

// =============================================================================
// PART 2: RANDOMIZED FUZZING vs GMP
// =============================================================================

#define FUZZ_ITERATIONS 1000

static void run_random(gmp_randstate_t state)
{
  printf("\n--- bn_mod_exp: %d randomized cases vs mpz_powm ---\n",
         FUZZ_ITERATIONS);

  for (int i = 0; i < FUZZ_ITERATIONS; i++) {
    int bits_a = 1 + (rand() % 4096);
    int bits_m = 1 + (rand() % 4096);
    int bits_e = 1 + (rand() % 512);  // keep exponents tractable

    bignum a, e, m, r;
    mpz_t za, ze, zm, zr;
    char* sa;
    char* se;
    char* sm;
    char detail[64];

    bn_init_multi(&a, &e, &m, &r, NULL);
    mpz_inits(za, ze, zm, zr, NULL);

    sa = random_mpz_str(za, bits_a, state);
    se = random_mpz_str(ze, bits_e, state);
    sm = random_mpz_str(zm, bits_m, state);
    const char* se_use = se;

    if (i % 50 == 0) {
      // sprinkle in zero exponents: a^0 = 1 mod m
      mpz_set_ui(ze, 0);
      se_use = "0";
    }

    bn_init_val(&a, sa);
    bn_init_val(&e, se_use);
    bn_init_val(&m, sm);

    bn_mod_exp(&r, &a, &e, &m);
    mpz_powm(zr, za, ze, zm);

    snprintf(detail, sizeof(detail),
             "case %d (a=%d bits, e=%d bits, m=%d bits)", i, bits_a, bits_e,
             bits_m);
    assert_match(detail, sa, se_use, sm, &r, zr);

    free(sa);
    free(se);
    free(sm);
    bn_free_multi(&a, &e, &m, &r, NULL);
    mpz_clears(za, ze, zm, zr, NULL);
  }

  printf("all %d random cases passed\n", FUZZ_ITERATIONS);
}

// =============================================================================
// PART 3: TIME-BASED BENCHMARKS vs GMP
// =============================================================================

static void benchmark_modexp(int bits, double target_sec, gmp_randstate_t state)
{
  bignum bn_a, bn_e, bn_m, bn_r;
  mpz_t mpz_a, mpz_e, mpz_m, mpz_r;

  bn_init_multi(&bn_a, &bn_e, &bn_m, &bn_r, NULL);
  mpz_inits(mpz_a, mpz_e, mpz_m, mpz_r, NULL);

  char* sa = random_mpz_str(mpz_a, bits, state);
  char* se = random_mpz_str(mpz_e, bits, state);
  mpz_setbit(mpz_m, 0);  // Ensure odd modulus (Montgomery path)
  char* sm = random_mpz_str(mpz_m, bits, state);
  bn_init_val(&bn_a, sa);
  bn_init_val(&bn_e, se);
  bn_init_val(&bn_m, sm);
  free(sa);
  free(se);
  free(sm);

  struct timespec start, end;
  int ops_custom = 0, ops_gmp = 0;
  double total_custom = 0, total_gmp = 0;

  clock_gettime(CLOCK_MONOTONIC, &start);
  do {
    bn_mod_exp(&bn_r, &bn_a, &bn_e, &bn_m);
    ops_custom++;
    clock_gettime(CLOCK_MONOTONIC, &end);
    total_custom = get_elapsed_time(start, end);
  } while (total_custom < target_sec);

  clock_gettime(CLOCK_MONOTONIC, &start);
  do {
    mpz_powm(mpz_r, mpz_a, mpz_e, mpz_m);
    ops_gmp++;
    clock_gettime(CLOCK_MONOTONIC, &end);
    total_gmp = get_elapsed_time(start, end);
  } while (total_gmp < target_sec);

  // Validate correctness before reporting
  assert_match("bn_mod_exp (benchmark)", "(bench input)", "(bench input)",
               "(bench input)", &bn_r, mpz_r);

  char size_info[32];
  snprintf(size_info, sizeof(size_info), "%d bits", bits);
  print_table_row("bn_mod_exp", size_info, total_custom / ops_custom,
                  total_gmp / ops_gmp);

  bn_free_multi(&bn_a, &bn_e, &bn_m, &bn_r, NULL);
  mpz_clears(mpz_a, mpz_e, mpz_m, mpz_r, NULL);
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
  benchmark_modexp(1024, bench_budget(1024), state);
  benchmark_modexp(2048, bench_budget(2048), state);
  benchmark_modexp(4096, bench_budget(4096), state);
  benchmark_modexp(8192, bench_budget(8192), state);
  benchmark_modexp(10000, bench_budget(10000), state);
  benchmark_modexp(16384, bench_budget(16384), state);
  print_table_footer();

  gmp_randclear(state);
  return 0;
}
