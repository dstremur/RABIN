/*
 * test_mod_inverse.c
 *
 * Unified GMP-verified test suite for bn_mod_inverse.
 *
 * Verified against mpz_invert: both must agree on existence, and on the
 * inverse value (in [0, m)) when it exists.
 *
 *   1. Edge cases: no-inverse (non-coprime, a = 0, m = 1), m = 2, a > m,
 *      large prime modulus
 *   2. 1000 randomized cases (a: 1..4096 bits signed, m: 1..4096 bits,
 *      prime or random) vs mpz_invert
 *   3. Time-based benchmark vs GMP at 512/1024/2048/4096 bits
 *
 * (No pointer-aliasing cases: a and m are const inputs, res aliasing is
 * not documented as supported.)
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

static void assert_inverse(const char* op, const char* a_str, const char* m_str,
                           bool custom_has_inv, const bignum* res,
                           int gmp_has_inv, mpz_t gres)
{
  if ((int)custom_has_inv != gmp_has_inv) {
    fprintf(stderr, "\n[FATAL ERROR] Correctness failure in %s!\n", op);
    fprintf(stderr, "Input A: %s\n", a_str ? a_str : "(n/a)");
    fprintf(stderr, "Input M: %s\n", m_str ? m_str : "(n/a)");
    fprintf(stderr, "GMP inv exists: %d\n", gmp_has_inv);
    fprintf(stderr, "Custom inv exists: %d\n", (int)custom_has_inv);
    exit(EXIT_FAILURE);
  }

  if (custom_has_inv) {
    char* gmp_res_str = mpz_get_str(NULL, 10, gres);
    char* custom_res_str = bn_to_string(res);
    if (strcmp(gmp_res_str, custom_res_str) != 0) {
      fprintf(stderr, "\n[FATAL ERROR] Correctness failure in %s (value)!\n",
              op);
      fprintf(stderr, "Input A: %s\n", a_str ? a_str : "(n/a)");
      fprintf(stderr, "Input M: %s\n", m_str ? m_str : "(n/a)");
      fprintf(stderr, "GMP Result:    %s\n", gmp_res_str);
      fprintf(stderr, "Custom Result: %s\n", custom_res_str);
      free(gmp_res_str);
      free(custom_res_str);
      exit(EXIT_FAILURE);
    }
    free(gmp_res_str);
    free(custom_res_str);
  }
}

// random bignum of given bits, optionally forced prime, >= min_val
static void gen_mpz(mpz_t z, int bits, bool prime, int min_val,
                    gmp_randstate_t state)
{
  mpz_urandomb(z, state, bits);
  if (prime) {
    mpz_nextprime(z, z);
  } else {
    // ensure z >= min_val
    if (mpz_cmp_si(z, min_val) < 0) mpz_set_si(z, min_val);
  }
}

// =============================================================================
// PART 1: EDGE CASES
// =============================================================================

static void run_case(const char* name, const char* a_str, const char* m_str)
{
  bignum a, m, res;
  mpz_t za, zm, zr;

  bn_init_multi(&a, &m, &res, NULL);
  bn_init_val(&a, a_str);
  bn_init_val(&m, m_str);

  mpz_inits(za, zm, zr, NULL);
  mpz_set_str(za, a_str, 10);
  mpz_set_str(zm, m_str, 10);

  bool custom_has_inv = bn_mod_inverse(&res, &a, &m);
  int gmp_has_inv = mpz_invert(zr, za, zm);

  char detail[256];
  snprintf(detail, sizeof(detail), "bn_mod_inverse [%s]", name);
  assert_inverse(detail, a_str, m_str, custom_has_inv, &res, gmp_has_inv, zr);

  bn_free_multi(&a, &m, &res, NULL);
  mpz_clears(za, zm, zr, NULL);
}

static void run_edge_cases()
{
  printf("\n--- bn_mod_inverse: edge cases ---\n");

  run_case("no inv (a even, m even)", "2", "4");
  run_case("no inv (a = 0)", "0", "13");
  run_case("no inv (m = 1)", "5", "1");
  run_case("no inv (a = m)", "13", "13");
  run_case("basic inv (3^-1 mod 11 = 4)", "3", "11");
  run_case("a = 1 (inv = 1)", "1", "13");
  run_case("m = 2", "1", "2");
  run_case("a > m (25^-1 mod 13)", "25", "13");
  run_case("negative a (-3^-1 mod 11)", "-3", "11");
  run_case("large prime modulus", "123456789", "18446744073709551557");
  run_case("u64-max a, prime m", "18446744073709551615",
           "18446744073709551557");

  printf("all edge cases passed\n");
}

// =============================================================================
// PART 2: RANDOMIZED FUZZING vs GMP
// =============================================================================

#define FUZZ_ITERATIONS 1000

static void run_random(gmp_randstate_t state)
{
  printf("\n--- bn_mod_inverse: %d randomized cases vs mpz_invert ---\n",
         FUZZ_ITERATIONS);

  for (int i = 0; i < FUZZ_ITERATIONS; i++) {
    int bits_a = 1 + (rand() % 4096);
    int bits_m = 1 + (rand() % 4096);
    bool prime_m = (rand() & 1);  // half prime, half random

    bignum a, m, res;
    mpz_t za, zm, zr;
    char* sa;
    char* sm;
    char detail[64];

    bn_init_multi(&a, &m, &res, NULL);
    mpz_inits(za, zm, zr, NULL);

    // a: signed
    mpz_urandomb(za, state, bits_a);
    if (rand() & 1) mpz_neg(za, za);
    sa = mpz_get_str(NULL, 10, za);

    gen_mpz(zm, bits_m, prime_m, 2, state);  // m >= 2
    sm = mpz_get_str(NULL, 10, zm);

    bn_init_val(&a, sa);
    bn_init_val(&m, sm);

    bool custom_has_inv = bn_mod_inverse(&res, &a, &m);
    int gmp_has_inv = mpz_invert(zr, za, zm);

    snprintf(detail, sizeof(detail), "case %d (a=%d bits, m=%d bits, %s)", i,
             bits_a, bits_m, prime_m ? "prime m" : "random m");
    assert_inverse(detail, sa, sm, custom_has_inv, &res, gmp_has_inv, zr);

    free(sa);
    free(sm);
    bn_free_multi(&a, &m, &res, NULL);
    mpz_clears(za, zm, zr, NULL);
  }

  printf("all %d random cases passed\n", FUZZ_ITERATIONS);
}

// =============================================================================
// PART 3: TIME-BASED BENCHMARKS vs GMP
// =============================================================================

static void benchmark_mod_inverse(int bits, double target_sec,
                                  gmp_randstate_t state)
{
  bignum bn_a, bn_m, bn_res;
  mpz_t mpz_a, mpz_m, mpz_res;

  bn_init_multi(&bn_a, &bn_m, &bn_res, NULL);
  mpz_inits(mpz_a, mpz_m, mpz_res, NULL);

  char* sa;
  char* sm;

  mpz_urandomb(mpz_a, state, bits);
  sa = mpz_get_str(NULL, 10, mpz_a);
  gen_mpz(mpz_m, bits, true, 2, state);  // prime ensures inverse exists
  sm = mpz_get_str(NULL, 10, mpz_m);
  bn_init_val(&bn_a, sa);
  bn_init_val(&bn_m, sm);
  free(sa);
  free(sm);

  struct timespec start, end;
  int ops_custom = 0, ops_gmp = 0;
  double total_custom = 0, total_gmp = 0;

  volatile bool custom_ok;
  volatile int gmp_ok;

  clock_gettime(CLOCK_MONOTONIC, &start);
  do {
    custom_ok = bn_mod_inverse(&bn_res, &bn_a, &bn_m);
    ops_custom++;
    clock_gettime(CLOCK_MONOTONIC, &end);
    total_custom = get_elapsed_time(start, end);
  } while (total_custom < target_sec);

  clock_gettime(CLOCK_MONOTONIC, &start);
  do {
    gmp_ok = mpz_invert(mpz_res, mpz_a, mpz_m);
    ops_gmp++;
    clock_gettime(CLOCK_MONOTONIC, &end);
    total_gmp = get_elapsed_time(start, end);
  } while (total_gmp < target_sec);

  // Validate correctness before reporting
  assert_inverse("bn_mod_inverse (benchmark)", "(bench input)", "(bench input)",
                 custom_ok, &bn_res, gmp_ok, mpz_res);

  char size_info[32];
  snprintf(size_info, sizeof(size_info), "%d bits", bits);
  print_table_row("bn_mod_inverse", size_info, total_custom / ops_custom,
                  total_gmp / ops_gmp);

  bn_free_multi(&bn_a, &bn_m, &bn_res, NULL);
  mpz_clears(mpz_a, mpz_m, mpz_res, NULL);
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
  benchmark_mod_inverse(2048, bench_budget(2048), state);
  benchmark_mod_inverse(8192, bench_budget(8192), state);
  benchmark_mod_inverse(32768, bench_budget(32768), state);
  benchmark_mod_inverse(65536, bench_budget(65536), state);
  benchmark_mod_inverse(65539, bench_budget(65539), state);
  print_table_footer();

  gmp_randclear(state);
  return 0;
}
