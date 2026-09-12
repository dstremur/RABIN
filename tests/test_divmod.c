/*
 * test_divmod.c
 *
 * Unified GMP-verified test suite for bn_divmod (truncated division).
 *
 *   1. Edge cases: 0, 1, -1, 2^64-1, 2^64 boundaries, truncation sign
 *      semantics (q toward 0, sign(r) = sign(a)) + pointer aliasing
 *      (q == a, q == b, r == a, r == b) and the NULL-output variants
 *   2. 1000 randomized cases (1..4096 bits) vs mpz_tdiv_qr
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

void assert_bn_eq_mpz(const char* op, const char* a_str, const char* b_str,
                      const bignum* bn, const char* label, mpz_t expected)
{
  char* gmp_str = mpz_get_str(NULL, 10, expected);
  char* custom_str = bn_to_string(bn);

  if (strcmp(gmp_str, custom_str) != 0) {
    fprintf(stderr, "\n[FATAL ERROR] Correctness failure in %s!\n", op);
    fprintf(stderr, "Input A: %s\n", a_str ? a_str : "(n/a)");
    fprintf(stderr, "Input B: %s\n", b_str ? b_str : "(n/a)");
    fprintf(stderr, "%s GMP:    %s\n", label, gmp_str);
    fprintf(stderr, "%s Custom: %s\n", label, custom_str);
    free(gmp_str);
    free(custom_str);
    exit(EXIT_FAILURE);
  }

  free(gmp_str);
  free(custom_str);
}

static void assert_divmod_values(const char* op, const char* a_str,
                                 const char* b_str, const bignum* q,
                                 const bignum* r, mpz_t zq, mpz_t zr)
{
  assert_bn_eq_mpz(op, a_str, b_str, q, "Q", zq);
  assert_bn_eq_mpz(op, a_str, b_str, r, "R", zr);
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
  bignum a, b, q, r, t, u;
  mpz_t za, zb, zq, zr;
  char detail[256];

  bn_init_multi(&a, &b, &q, &r, &t, &u, NULL);
  mpz_inits(za, zb, zq, zr, NULL);
  bn_init_val(&a, a_str);
  bn_init_val(&b, b_str);
  mpz_set_str(za, a_str, 10);
  mpz_set_str(zb, b_str, 10);

  snprintf(detail, sizeof(detail), "bn_divmod [%s] a=%s b=%s", name, a_str,
           b_str);

  // non-aliased
  bn_divmod(&q, &r, &a, &b);
  mpz_tdiv_qr(zq, zr, za, zb);
  assert_divmod_values(detail, a_str, b_str, &q, &r, zq, zr);

  // internal invariant: q * b + r == a (computed with the custom lib)
  bn_mul(&t, &q, &b);
  bn_add(&u, &t, &r);
  assert_bn_eq_mpz("bn_divmod invariant q*b+r==a", a_str, b_str, &u, "Q*B+R",
                   za);

  // aliasing: q == a
  bn_init_val(&a, a_str);
  bn_init_val(&b, b_str);
  bn_divmod(&a, &r, &a, &b);
  assert_divmod_values("bn_divmod (q==a)", a_str, b_str, &a, &r, zq, zr);

  // aliasing: r == b
  bn_init_val(&a, a_str);
  bn_init_val(&b, b_str);
  bn_divmod(&q, &b, &a, &b);
  assert_divmod_values("bn_divmod (r==b)", a_str, b_str, &q, &b, zq, zr);

  // aliasing: q == b
  bn_init_val(&a, a_str);
  bn_init_val(&b, b_str);
  bn_divmod(&b, &r, &a, &b);
  assert_divmod_values("bn_divmod (q==b)", a_str, b_str, &b, &r, zq, zr);

  // aliasing: r == a
  bn_init_val(&a, a_str);
  bn_init_val(&b, b_str);
  bn_divmod(&q, &a, &a, &b);
  assert_divmod_values("bn_divmod (r==a)", a_str, b_str, &q, &a, zq, zr);

  // NULL quotient: remainder only
  bn_init_val(&a, a_str);
  bn_init_val(&b, b_str);
  bn_divmod(NULL, &r, &a, &b);
  assert_bn_eq_mpz("bn_divmod (q=NULL)", a_str, b_str, &r, "R", zr);

  // NULL remainder: quotient only
  bn_init_val(&a, a_str);
  bn_init_val(&b, b_str);
  bn_divmod(&q, NULL, &a, &b);
  assert_bn_eq_mpz("bn_divmod (r=NULL)", a_str, b_str, &q, "Q", zq);

  bn_free_multi(&a, &b, &q, &r, &t, &u, NULL);
  mpz_clears(za, zb, zq, zr, NULL);
}

static void run_edge_cases()
{
  printf("\n--- bn_divmod: edge cases ---\n");

  run_case("zero / one", "0", "1");
  run_case("one / one", "1", "1");
  run_case("seven / two", "7", "2");
  run_case("neg seven / two (q=-3, r=-1)", "-7", "2");
  run_case("seven / neg two (q=-3, r=1)", "7", "-2");
  run_case("neg seven / neg two (q=3, r=-1)", "-7", "-2");
  run_case("neg one / two (q=0, r=-1)", "-1", "2");
  run_case("a < b (q=0, r=a)", "1", "18446744073709551615");
  run_case("a == b (q=1, r=0)", "18446744073709551615", "18446744073709551615");
  run_case("neg a < |b| (q=0, r=a)", "-18446744073709551615",
           "18446744073709551615");
  run_case("u64-max / 2", "18446744073709551615", "2");
  run_case("2^64 / 2^64 (limb boundary)", "18446744073709551616",
           "18446744073709551616");
  run_case("huge / 2 (asymmetric)", "123456789012345678901234567890", "2");

  printf("all edge cases passed\n");
}

// =============================================================================
// PART 2: RANDOMIZED FUZZING vs GMP
// =============================================================================

#define FUZZ_ITERATIONS 1000

static void run_random(gmp_randstate_t state)
{
  printf("\n--- bn_divmod: %d randomized cases vs mpz_tdiv_qr ---\n",
         FUZZ_ITERATIONS);

  for (int i = 0; i < FUZZ_ITERATIONS; i++) {
    int bits_a = 1 + (rand() % 4096);
    int bits_b = 1 + (rand() % 4096);

    bignum a, b, q, r;
    mpz_t za, zb, zq, zr;
    char* sa;
    char* sb;
    char detail[64];

    bn_init_multi(&a, &b, &q, &r, NULL);
    mpz_inits(za, zb, zq, zr, NULL);

    sa = random_mpz_str(za, bits_a, state);
    sb = random_mpz_str(zb, bits_b, state);

    // divisor must be nonzero
    if (mpz_cmp_ui(zb, 0) == 0) {
      mpz_set_ui(zb, 1);
      free(sb);
      sb = mpz_get_str(NULL, 10, zb);
    }

    bn_init_val(&a, sa);
    bn_init_val(&b, sb);

    bn_divmod(&q, &r, &a, &b);
    mpz_tdiv_qr(zq, zr, za, zb);

    snprintf(detail, sizeof(detail), "case %d (a=%d bits, b=%d bits)", i,
             bits_a, bits_b);
    assert_divmod_values(detail, sa, sb, &q, &r, zq, zr);

    free(sa);
    free(sb);
    bn_free_multi(&a, &b, &q, &r, NULL);
    mpz_clears(za, zb, zq, zr, NULL);
  }

  printf("all %d random cases passed\n", FUZZ_ITERATIONS);
}

// =============================================================================
// PART 3: TIME-BASED BENCHMARKS vs GMP
// =============================================================================

static void benchmark_divmod(int bits, double target_sec, gmp_randstate_t state)
{
  bignum bn_a, bn_b, bn_q, bn_r;
  mpz_t mpz_a, mpz_b, mpz_q, mpz_r;

  bn_init_multi(&bn_a, &bn_b, &bn_q, &bn_r, NULL);
  mpz_inits(mpz_a, mpz_b, mpz_q, mpz_r, NULL);

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
    bn_divmod(&bn_q, &bn_r, &bn_a, &bn_b);
    ops_custom++;
    clock_gettime(CLOCK_MONOTONIC, &end);
    total_custom = get_elapsed_time(start, end);
  } while (total_custom < target_sec);

  clock_gettime(CLOCK_MONOTONIC, &start);
  do {
    mpz_tdiv_qr(mpz_q, mpz_r, mpz_a, mpz_b);
    ops_gmp++;
    clock_gettime(CLOCK_MONOTONIC, &end);
    total_gmp = get_elapsed_time(start, end);
  } while (total_gmp < target_sec);

  // Validate correctness before reporting
  assert_divmod_values("bn_divmod (benchmark)", "(bench input)",
                       "(bench input)", &bn_q, &bn_r, mpz_q, mpz_r);

  char size_info[32];
  snprintf(size_info, sizeof(size_info), "%d bits", bits);
  print_table_row("bn_divmod", size_info, total_custom / ops_custom,
                  total_gmp / ops_gmp);

  bn_free_multi(&bn_a, &bn_b, &bn_q, &bn_r, NULL);
  mpz_clears(mpz_a, mpz_b, mpz_q, mpz_r, NULL);
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
  benchmark_divmod(2048, bench_budget(2048), state);
  benchmark_divmod(8192, bench_budget(8192), state);
  benchmark_divmod(32768, bench_budget(32768), state);
  benchmark_divmod(65536, bench_budget(65536), state);
  benchmark_divmod(65539, bench_budget(65539), state);
  print_table_footer();

  gmp_randclear(state);
  return 0;
}
