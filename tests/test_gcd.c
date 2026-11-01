/*
 * test_gcd.c
 *
 * Unified GMP-verified test suite for rz_gcd.
 *
 *   1. Edge cases: 0, 1, -1, (0,0), 2^64-1, 2^64 boundaries + pointer
 *      aliasing (d == a and d == b)
 *   2. 1000 randomized cases (1..4096 bits, occasional zeros) vs mpz_gcd
 *   3. Deterministic rz_gcd_lehmer() stress cases (fixed seed): bit-63
 *      limb overflow patterns, fibonacci quotient chains, pairs whose
 *      truncated top limbs agree
 *   4. Time-based benchmark vs GMP at 512/1024/2048/4096 bits
 *   5. (kept) internal comparison of gcd variants (normal / extended /
 *      lehmer / lehmer-extended)
 *
 * Copyright (C) 2026 Diego Strebel
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <gmp.h>
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

static char* random_a(mpz_t za, int bits, gmp_randstate_t state)
{
  mpz_urandomb(za, state, bits);
  if (rand() & 1) mpz_neg(za, za);
  return mpz_get_str(NULL, 10, za);
}

// =============================================================================
// PART 1: EDGE CASES (incl. pointer aliasing)
// =============================================================================

static void run_case(const char* name, const char* a_str, const char* b_str)
{
  rz_t a, b, d;
  mpz_t za, zb, zd;
  char detail[256];

  rz_init_multi(&a, &b, &d, NULL);
  mpz_inits(za, zb, zd, NULL);
  rz_init_val(&a, a_str);
  rz_init_val(&b, b_str);
  mpz_set_str(za, a_str, 10);
  mpz_set_str(zb, b_str, 10);

  snprintf(detail, sizeof(detail), "rz_gcd [%s] a=%s b=%s", name, a_str, b_str);

  // non-aliased
  rz_gcd(&d, &a, &b);
  mpz_gcd(zd, za, zb);
  assert_match(detail, a_str, b_str, &d, zd);

  // aliasing: d == a
  rz_init_val(&a, a_str);
  rz_gcd(&a, &a, &b);
  assert_match("rz_gcd (d==a)", a_str, b_str, &a, zd);

  // aliasing: d == b
  rz_init_val(&a, a_str);
  rz_init_val(&b, b_str);
  rz_gcd(&b, &a, &b);
  assert_match("rz_gcd (d==b)", a_str, b_str, &b, zd);

  rz_clear_multi(&a, &b, &d, NULL);
  mpz_clears(za, zb, zd, NULL);
}

static void run_edge_cases()
{
  printf("\n--- rz_gcd: edge cases ---\n");

  run_case("both zero (gcd(0,0) = 0)", "0", "0");
  run_case("zero + one", "0", "1");
  run_case("one + zero", "1", "0");
  run_case("one + neg-one", "1", "-1");
  run_case("classic 1071/462 (=21)", "1071", "462");
  run_case("u64-max + u64-max", "18446744073709551615", "18446744073709551615");
  run_case("u64-max + 2^64 (gcd = 1)", "18446744073709551615",
           "18446744073709551616");
  run_case("2^64 + 2^64", "18446744073709551616", "18446744073709551616");
  run_case("2^64 + 2^63 (gcd = 2^63)", "18446744073709551616",
           "9223372036854775808");
  run_case("negative inputs (result non-negative)", "-1071", "-462");
  run_case("fibonacci worst case", "107011735028014877650776236975863908754",
           "66208420403979399419517229710735292609");

  printf("all edge cases passed\n");
}

// =============================================================================
// PART 2: RANDOMIZED FUZZING vs GMP
// =============================================================================

#define FUZZ_ITERATIONS 1000

static void run_random(gmp_randstate_t state)
{
  printf("\n--- rz_gcd: %d randomized cases vs mpz_gcd ---\n", FUZZ_ITERATIONS);

  for (int i = 0; i < FUZZ_ITERATIONS; i++) {
    int bits_a = 1 + (rand() % 4096);
    int bits_b = 1 + (rand() % 4096);

    rz_t a, b, d;
    mpz_t za, zb, zd;
    char* sa;
    char* sb;
    char detail[64];

    rz_init_multi(&a, &b, &d, NULL);
    mpz_inits(za, zb, zd, NULL);

    sa = random_a(za, bits_a, state);
    sb = random_a(zb, bits_b, state);
    const char* sa_use = sa;
    const char* sb_use = sb;

    // sprinkle in zeros: gcd(0, x) = |x|
    if (i % 100 == 0) {
      sa_use = "0";
      mpz_set_ui(za, 0);
    }
    if (i % 100 == 1) {
      sb_use = "0";
      mpz_set_ui(zb, 0);
    }

    rz_init_val(&a, sa_use);
    rz_init_val(&b, sb_use);

    rz_gcd(&d, &a, &b);
    mpz_gcd(zd, za, zb);

    snprintf(detail, sizeof(detail), "case %d (a=%d bits, b=%d bits)", i,
             bits_a, bits_b);
    assert_match(detail, sa_use, sb_use, &d, zd);

    free(sa);
    free(sb);
    rz_clear_multi(&a, &b, &d, NULL);
    mpz_clears(za, zb, zd, NULL);
  }

  printf("all %d random cases passed\n", FUZZ_ITERATIONS);
}

// =============================================================================
// PART 3: DETERMINISTIC LEHMER STRESS CASES
// =============================================================================

// Regression net for the Lehmer quotient loop: operands whose top limb
// carries bit 63 (bit lengths that are exact multiples of 64, where a
// naive limb-to-i64 conversion overflows), quotient-1 chains from
// consecutive Fibonacci numbers that grow the transformation matrix, and
// same-size pairs whose truncated top views agree before diverging.
// Every pair is also run through all four input sign combinations.
static void lehmer_stress_pair(const char* what, const mpz_t za, const mpz_t zb)
{
  rz_t a, b, d;
  mpz_t zg;
  rz_init_multi(&a, &b, &d, NULL);
  mpz_init(zg);
  mpz_gcd(zg, za, zb);

  char op[128];
  snprintf(op, sizeof(op), "rz_gcd_lehmer [%s]", what);

  static const int signs[4][2] = {{1, 1}, {-1, 1}, {1, -1}, {-1, -1}};
  for (int s = 0; s < 4; s++) {
    mpz_t xa, xb;
    mpz_inits(xa, xb, NULL);
    mpz_set(xa, za);
    mpz_set(xb, zb);
    if (signs[s][0] < 0) mpz_neg(xa, xa);
    if (signs[s][1] < 0) mpz_neg(xb, xb);

    char* sa = mpz_get_str(NULL, 10, xa);
    char* sb = mpz_get_str(NULL, 10, xb);
    rz_init_val(&a, sa);
    rz_init_val(&b, sb);

    rz_gcd_lehmer(&d, &a, &b);
    assert_match(op, sa, sb, &d, zg);

    free(sa);
    free(sb);
    mpz_clears(xa, xb, NULL);
  }

  rz_clear_multi(&a, &b, &d, NULL);
  mpz_clear(zg);
}

static void run_lehmer_stress(void)
{
  printf("\n--- rz_gcd_lehmer: deterministic stress cases ---\n");

  // fixed seed: a failing case must reproduce without the run's seed
  gmp_randstate_t state;
  gmp_randinit_default(state);
  gmp_randseed_ui(state, 0x5EED5EEDUL);

  mpz_t a, b, s;
  mpz_inits(a, b, s, NULL);

  // Consecutive Fibonacci pairs: every quotient is 1, the longest matrix
  // accumulation; the scaled variant leaves a multi-limb common factor
  const unsigned long fibs[] = {3, 10, 40, 80, 160, 320};
  for (size_t i = 0; i < sizeof(fibs) / sizeof(fibs[0]); i++) {
    mpz_fib_ui(a, fibs[i] + 1);
    mpz_fib_ui(b, fibs[i]);
    lehmer_stress_pair("fibonacci chain", a, b);

    mpz_set_ui(s, 1);
    mpz_setbit(s, 64);  // s = 2^64 + 1
    mpz_mul(a, a, s);
    mpz_mul(b, b, s);
    lehmer_stress_pair("scaled fibonacci chain", a, b);
  }

  // Sizes that are exact multiples of 64 bits: the top limb carries
  // bit 63, the exact overflow point of a raw limb-to-i64 view
  const unsigned limb_counts[] = {1, 2, 5, 16, 64};
  for (size_t i = 0; i < sizeof(limb_counts) / sizeof(limb_counts[0]); i++) {
    unsigned bits = 64 * limb_counts[i];

    mpz_urandomb(a, state, bits - 1);
    mpz_setbit(a, bits - 1);
    mpz_urandomb(b, state, bits - 1);
    mpz_setbit(b, bits - 1);
    lehmer_stress_pair("top limb bit 63", a, b);

    // b = 2^bits - 1: all-ones limb views
    mpz_ui_pow_ui(b, 2, bits);
    mpz_sub_ui(b, b, 1);
    lehmer_stress_pair("all-ones divisor", a, b);

    // b = a with one low bit toggled: identical truncated views that
    // only diverge in the bottom limbs
    mpz_set(b, a);
    mpz_combit(b, 7);
    lehmer_stress_pair("identical top views", a, b);
  }

  mpz_clears(a, b, s, NULL);
  gmp_randclear(state);
  printf("all lehmer stress cases passed\n");
}

// =============================================================================
// PART 4: TIME-BASED BENCHMARKS vs GMP
// =============================================================================

static void benchmark_gcd(int bits, double target_sec, gmp_randstate_t state)
{
  rz_t rz_a, rz_b, rz_d;
  mpz_t mpz_a, mpz_b, mpz_d;

  rz_init_multi(&rz_a, &rz_b, &rz_d, NULL);
  mpz_inits(mpz_a, mpz_b, mpz_d, NULL);

  char* sa = random_a(mpz_a, bits, state);
  char* sb = random_a(mpz_b, bits, state);
  rz_init_val(&rz_a, sa);
  rz_init_val(&rz_b, sb);
  free(sa);
  free(sb);

  struct timespec start, end;
  int ops_custom = 0, ops_gmp = 0;
  double total_custom = 0, total_gmp = 0;

  clock_gettime(CLOCK_MONOTONIC, &start);
  do {
    rz_gcd(&rz_d, &rz_a, &rz_b);
    ops_custom++;
    clock_gettime(CLOCK_MONOTONIC, &end);
    total_custom = get_elapsed_time(start, end);
  } while (total_custom < target_sec);

  clock_gettime(CLOCK_MONOTONIC, &start);
  do {
    mpz_gcd(mpz_d, mpz_a, mpz_b);
    ops_gmp++;
    clock_gettime(CLOCK_MONOTONIC, &end);
    total_gmp = get_elapsed_time(start, end);
  } while (total_gmp < target_sec);

  // Validate correctness before reporting
  assert_match("rz_gcd (benchmark)", "(bench input)", "(bench input)", &rz_d,
               mpz_d);

  char size_info[32];
  snprintf(size_info, sizeof(size_info), "%d bits", bits);
  print_table_row("rz_gcd", size_info, total_custom / ops_custom,
                  total_gmp / ops_gmp);

  rz_clear_multi(&rz_a, &rz_b, &rz_d, NULL);
  mpz_clears(mpz_a, mpz_b, mpz_d, NULL);
}

// =============================================================================
// PART 5 (kept): INTERNAL VARIANT COMPARISON
// =============================================================================

static void benchmark_gcd_internal(int limbs)
{
  rz_t a, b, u_1, v_1, d_1, u_2, v_2, d_2;
  rz_init_multi(&a, &b, &u_1, &v_1, &d_1, &u_2, &v_2, &d_2, NULL);

  // Generate random numbers of a specific bit length (limbs * 64)
  rz_gen_random(&a, limbs * 64);
  rz_gen_random(&b, limbs * 64);

  // 1. Time Normal GCD
  clock_t start = clock();
  rz_gcd(&d_1, &a, &b);
  clock_t end = clock();
  double time_normal = (double)(end - start) / CLOCKS_PER_SEC;

  // 2. Time Karatsuba
  start = clock();
  rz_gcd_extended(&u_2, &v_2, &d_2, &a, &b);
  end = clock();
  double time_ex = (double)(end - start) / CLOCKS_PER_SEC;

  start = clock();
  rz_gcd_lehmer(&d_2, &a, &b);
  end = clock();
  double time_l = (double)(end - start) / CLOCKS_PER_SEC;

  start = clock();
  rz_gcd_extended_lehmer(&u_2, &v_2, &d_2, &a, &b);
  end = clock();
  double time_lehmer = (double)(end - start) / CLOCKS_PER_SEC;

  // 3. Verify Correctness
  if (rz_cmp(&d_1, &d_2) != 0) {
    printf("[FAIL] Mismatch at %d limbs!\n", limbs);
  } else {
    printf(
        "Limbs: %4d | Normal: %fs | Extended: %fs | Lehmer: %fs | Lehmer "
        "Extended: "
        "%fs "
        "\n",
        limbs, time_normal, time_ex, time_l, time_lehmer);
  }

  rz_clear_multi(&a, &b, &u_1, &v_1, &d_1, &u_2, &v_2, &d_2, NULL);
}

static void run_internal_comparison()
{
  printf("\n--- rz_gcd: variant comparison (internal) ---\n");
  int sizes[] = {16,   32,   64,   128,  256,   512,
                 1024, 2048, 4096, 8192, 16000, 32000};
  for (size_t i = 0; i < 9; i++) {
    benchmark_gcd_internal(sizes[i]);
  }
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
  rz_init_constants();

  gmp_randstate_t state;
  gmp_randinit_default(state);
  gmp_randseed_ui(state, (unsigned long)time(NULL));
  srand((unsigned int)time(NULL));

  run_edge_cases();
  run_random(state);
  run_lehmer_stress();

  print_table_header();
  benchmark_gcd(2048, bench_budget(2048), state);
  benchmark_gcd(8192, bench_budget(8192), state);
  benchmark_gcd(32768, bench_budget(32768), state);
  benchmark_gcd(65536, bench_budget(65536), state);
  benchmark_gcd(65539, bench_budget(65539), state);
  benchmark_gcd(100000, bench_budget(100000), state);
  benchmark_gcd(200000, bench_budget(200000), state);
  benchmark_gcd(500000, bench_budget(500000), state);
  print_table_footer();

  run_internal_comparison();

  gmp_randclear(state);
  rz_clear_constants();
  return 0;
}
