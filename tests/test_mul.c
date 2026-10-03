/*
 * test_mul.c
 *
 * Unified GMP-verified test suite for rz_mul.
 *
 *   1. Edge cases: 0, 1, -1, 2^64-1, 2^64 boundaries, limb-boundary carries
 *      + pointer aliasing (r == a and r == b)
 *   2. 1000 randomized cases (1..4096 bits) vs mpz_mul
 *   3. Time-based benchmark vs GMP at 2048..262144 bits
 *   4. (kept) schoolbook vs rz_mul (Karatsuba/NTT) comparison per size
 *   5. u64-NTT path: threshold-straddling, asymmetric, carry, squaring vs GMP
 *   6. Karatsuba vs NTT timing for tuning RZ_NTT_LIMIT
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
                  const rz_t* bn, mpz_t mpz)
{
  char* gmp_str = mpz_get_str(NULL, 10, mpz);
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
  rz_t a, b, res;
  mpz_t za, zb, zr;
  char detail[256];

  rz_init_multi(&a, &b, &res, NULL);
  mpz_inits(za, zb, zr, NULL);
  rz_init_val(&a, a_str);
  rz_init_val(&b, b_str);
  mpz_set_str(za, a_str, 10);
  mpz_set_str(zb, b_str, 10);

  snprintf(detail, sizeof(detail), "rz_mul [%s] a=%s b=%s", name, a_str, b_str);

  // non-aliased
  rz_mul(&res, &a, &b);
  mpz_mul(zr, za, zb);
  assert_match(detail, a_str, b_str, &res, zr);

  // aliasing: r == a
  rz_init_val(&a, a_str);
  rz_mul(&a, &a, &b);
  snprintf(detail, sizeof(detail), "rz_mul [%s, r==a] a=%s b=%s", name, a_str,
           b_str);
  assert_match(detail, a_str, b_str, &a, zr);

  // aliasing: r == b
  rz_init_val(&a, a_str);
  rz_init_val(&b, b_str);
  rz_mul(&b, &a, &b);
  snprintf(detail, sizeof(detail), "rz_mul [%s, r==b] a=%s b=%s", name, a_str,
           b_str);
  assert_match(detail, a_str, b_str, &b, zr);

  rz_clear_multi(&a, &b, &res, NULL);
  mpz_clears(za, zb, zr, NULL);
}

static void run_edge_cases()
{
  printf("\n--- rz_mul: edge cases ---\n");

  run_case("zeros", "0", "0");
  run_case("one * one", "1", "1");
  run_case("neg-one * neg-one", "-1", "-1");
  run_case("one * neg-one", "1", "-1");
  run_case("multiply by zero", "987654321987654321", "0");
  run_case("u64-max * u64-max (128-bit result)", "18446744073709551615",
           "18446744073709551615");
  run_case("u64-max * 2 (carry cascade)", "18446744073709551615", "2");
  run_case("2^64 * 2^64 (limb boundary)", "18446744073709551616",
           "18446744073709551616");
  run_case("2^64 * u64-max", "18446744073709551616", "18446744073709551615");
  run_case("negative * positive", "-18446744073709551615", "12345");
  run_case("asymmetric sizes", "123456789012345678901234567890", "2");

  printf("all edge cases passed\n");
}

// =============================================================================
// PART 2: RANDOMIZED FUZZING vs GMP
// =============================================================================

#define FUZZ_ITERATIONS 1000

static void run_random(gmp_randstate_t state)
{
  printf("\n--- rz_mul: %d randomized cases vs mpz_mul ---\n", FUZZ_ITERATIONS);

  for (int i = 0; i < FUZZ_ITERATIONS; i++) {
    int bits_a = 1 + (rand() % 4096);
    int bits_b = 1 + (rand() % 4096);

    rz_t a, b, res;
    mpz_t za, zb, zr;
    char* sa;
    char* sb;
    char detail[64];

    rz_init_multi(&a, &b, &res, NULL);
    mpz_inits(za, zb, zr, NULL);

    sa = random_mpz_str(za, bits_a, state);
    sb = random_mpz_str(zb, bits_b, state);
    rz_init_val(&a, sa);
    rz_init_val(&b, sb);

    rz_mul(&res, &a, &b);
    mpz_mul(zr, za, zb);

    snprintf(detail, sizeof(detail), "case %d (a=%d bits, b=%d bits)", i,
             bits_a, bits_b);
    assert_match(detail, sa, sb, &res, zr);

    free(sa);
    free(sb);
    rz_clear_multi(&a, &b, &res, NULL);
    mpz_clears(za, zb, zr, NULL);
  }

  printf("all %d random cases passed\n", FUZZ_ITERATIONS);
}

// =============================================================================
// PART 3: TIME-BASED BENCHMARKS vs GMP
// =============================================================================

static void benchmark_mul(int bits, double target_sec, gmp_randstate_t state)
{
  rz_t rz_a, rz_b, rz_res;
  mpz_t mpz_a, mpz_b, mpz_res;

  rz_init_multi(&rz_a, &rz_b, &rz_res, NULL);
  mpz_inits(mpz_a, mpz_b, mpz_res, NULL);

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
    rz_mul(&rz_res, &rz_a, &rz_b);
    ops_custom++;
    clock_gettime(CLOCK_MONOTONIC, &end);
    total_custom = get_elapsed_time(start, end);
  } while (total_custom < target_sec);

  clock_gettime(CLOCK_MONOTONIC, &start);
  do {
    mpz_mul(mpz_res, mpz_a, mpz_b);
    ops_gmp++;
    clock_gettime(CLOCK_MONOTONIC, &end);
    total_gmp = get_elapsed_time(start, end);
  } while (total_gmp < target_sec);

  // Validate correctness before reporting
  assert_match("rz_mul (benchmark)", "(bench input)", "(bench input)", &rz_res,
               mpz_res);

  char size_info[32];
  snprintf(size_info, sizeof(size_info), "%d bits", bits);
  print_table_row("rz_mul", size_info, total_custom / ops_custom,
                  total_gmp / ops_gmp);

  rz_clear_multi(&rz_a, &rz_b, &rz_res, NULL);
  mpz_clears(mpz_a, mpz_b, mpz_res, NULL);
}

// =============================================================================
// PART 4 (kept): SCHOOLBOOK vs KARATSUBA COMPARISON
// =============================================================================

static void benchmark_mul_internal(int limbs)
{
  rz_t a, b, res_school, res_karat;
  rz_init_multi(&a, &b, &res_school, &res_karat, NULL);

  // Generate random numbers of a specific bit length (limbs * 64)
  rz_gen_random(&a, limbs * 64);
  rz_gen_random(&b, limbs * 64);

  // 1. Time Schoolbook
  clock_t start = clock();
  rz_mul_school(&res_school, &a, &b);
  clock_t end = clock();
  double time_school = (double)(end - start) / CLOCKS_PER_SEC;

  // 2. Time rz_mul (Karatsuba below RZ_NTT_LIMIT, NTT above it)
  start = clock();
  rz_mul(&res_karat, &a, &b);
  end = clock();
  double time_karat = (double)(end - start) / CLOCKS_PER_SEC;

  // 3. Verify Correctness
  if (rz_cmp(&res_school, &res_karat) != 0) {
    printf("[FAIL] Mismatch at %d limbs!\n", limbs);
  } else {
    double speedup = (time_karat > 0.0) ? (time_school / time_karat) : 0.0;
    printf("Limbs: %4d | School: %.6fs | rz_mul: %.6fs | Speedup: %.2fx\n",
           limbs, time_school, time_karat, speedup);
  }

  rz_clear_multi(&a, &b, &res_school, &res_karat, NULL);
}

static void run_internal_comparison()
{
  printf("\n--- rz_mul: schoolbook vs Karatsuba (internal) ---\n");
  int sizes[] = {16,   32,   64,   128,   256,   512,   1024,
                 2048, 4096, 8192, 10000, 20000, 40000, 80000};
  for (size_t i = 0; i < sizeof(sizes) / sizeof(sizes[0]); i++) {
    benchmark_mul_internal(sizes[i]);
  }
}

// =============================================================================
// PART 5: NTT PATH (threshold-straddling, asymmetry, carry, squaring)
// =============================================================================

// Fast GMP reference check: imports the limb arrays directly (O(n)),
// avoiding the quadratic decimal string round-trip at multi-Mbit sizes.
static void assert_match_gmp(const char* op, const rz_t* a, const rz_t* b,
                             const rz_t* res)
{
  mpz_t za, zb, zr, zref;
  mpz_inits(za, zb, zr, zref, NULL);

  mpz_import(za, a->size, 0, 8, 1, 64, a->limbs);
  mpz_import(zb, b->size, 0, 8, 1, 64, b->limbs);
  mpz_import(zr, res->size, 0, 8, 1, 64, res->limbs);
  if (a->is_neg) mpz_neg(za, za);
  if (b->is_neg) mpz_neg(zb, zb);
  if (res->is_neg) mpz_neg(zr, zr);

  mpz_mul(zref, za, zb);
  if (mpz_cmp(zr, zref) != 0) {
    fprintf(stderr, "\n[FATAL ERROR] Correctness failure in %s!\n", op);
    char* sa = rz_to_string(a);
    char* sb = rz_to_string(b);
    char* sr = rz_to_string(res);
    char* sg = mpz_get_str(NULL, 10, zref);
    fprintf(stderr, "Input A: %s\n", sa ? sa : "(null)");
    fprintf(stderr, "Input B: %s\n", sb ? sb : "(null)");
    fprintf(stderr, "GMP Result:    %s\n", sg);
    fprintf(stderr, "Custom Result: %s\n", sr ? sr : "(null)");
    free(sa);
    free(sb);
    free(sr);
    free(sg);
    mpz_clears(za, zb, zr, zref, NULL);
    exit(EXIT_FAILURE);
  }

  mpz_clears(za, zb, zr, zref, NULL);
}

static void run_ntt_case(int limbs_a, int limbs_b, const char* name)
{
  rz_t a, b, res;
  char detail[96];

  rz_init_multi(&a, &b, &res, NULL);

  rz_gen_random(&a, (u64)limbs_a * 64);
  rz_gen_random(&b, (u64)limbs_b * 64);

  snprintf(detail, sizeof(detail), "rz_mul [NTT %s]", name);
  rz_mul(&res, &a, &b);
  assert_match_gmp(detail, &a, &b, &res);

  rz_clear_multi(&a, &b, &res, NULL);
}

static void run_ntt_cases()
{
  printf("\n--- rz_mul: NTT threshold-straddling cases vs GMP ---\n");
  run_ntt_case(RZ_NTT_LIMIT / 2, RZ_NTT_LIMIT / 2, "limit/2");
  run_ntt_case(RZ_NTT_LIMIT / 2 + 1, RZ_NTT_LIMIT / 2 + 1, "limit/2+1");
  run_ntt_case(RZ_NTT_LIMIT - 1, RZ_NTT_LIMIT - 1, "limit-1");
  run_ntt_case(RZ_NTT_LIMIT, RZ_NTT_LIMIT, "limit");
  run_ntt_case(RZ_NTT_LIMIT + 1, RZ_NTT_LIMIT + 1, "limit+1");
  run_ntt_case(4 * RZ_NTT_LIMIT, 4 * RZ_NTT_LIMIT, "4*limit");
  run_ntt_case(16 * RZ_NTT_LIMIT, 16 * RZ_NTT_LIMIT, "16*limit");
  printf("all NTT threshold cases passed\n");
}

static void run_ntt_asymmetric()
{
  printf("\n--- rz_mul: NTT asymmetric pairs vs GMP ---\n");

  // 4:1 ratio straddles the asymmetry guard: dispatches to the NTT path
  run_ntt_case(4 * RZ_NTT_LIMIT, RZ_NTT_LIMIT, "4:1 via rz_mul");
  run_ntt_case(RZ_NTT_LIMIT, 4 * RZ_NTT_LIMIT, "1:4 via rz_mul");

  // Extreme ratio: rz_mul() falls back to schoolbook; rz_mul_fast()
  // exercises the explicit NTT kernel's zero-padding either way
  {
    rz_t a, b, res, res_fast;
    char detail[96];

    rz_init_multi(&a, &b, &res, &res_fast, NULL);

    rz_gen_random(&a, (u64)(4 * RZ_NTT_LIMIT) * 64);
    rz_set_u64(&b, 123456789);

    rz_mul(&res, &a, &b);
    snprintf(detail, sizeof(detail), "rz_mul [4*limit x 1, fallback]");
    assert_match_gmp(detail, &a, &b, &res);

    rz_mul_fast(&res_fast, &a, &b);
    snprintf(detail, sizeof(detail), "rz_mul_fast [4*limit x 1]");
    assert_match_gmp(detail, &a, &b, &res_fast);

    rz_mul(&res, &b, &a);
    snprintf(detail, sizeof(detail), "rz_mul [1 x 4*limit, fallback]");
    assert_match_gmp(detail, &b, &a, &res);

    rz_mul_fast(&res_fast, &b, &a);
    snprintf(detail, sizeof(detail), "rz_mul_fast [1 x 4*limit]");
    assert_match_gmp(detail, &b, &a, &res_fast);

    rz_clear_multi(&a, &b, &res, &res_fast, NULL);
  }
  printf("all NTT asymmetric cases passed\n");
}

static void run_ntt_carry()
{
  printf(
      "\n--- rz_mul: worst-case carry (all-0xFFFF limbs, self square) ---\n");
  int limbs = 2 * RZ_NTT_LIMIT;

  rz_t a, one, res;
  char detail[96];

  rz_init_multi(&a, &one, &res, NULL);

  // a = 2^(64*limbs) - 1, i.e. every limb is 0xFFFF
  rz_set_u64(&one, 1);
  rz_lshift(&a, &one, 64 * limbs);
  rz_sub(&a, &a, &one);

  rz_mul(&res, &a, &a);
  snprintf(detail, sizeof(detail), "rz_mul [all-0xFFFF %d limbs ^2]", limbs);
  assert_match_gmp(detail, &a, &a, &res);

  rz_clear_multi(&a, &one, &res, NULL);
  printf("worst-case carry passed\n");
}

static void run_ntt_sqr()
{
  printf("\n--- rz_mul: squaring dispatch (a == b) vs GMP ---\n");
  int limbs = 2 * RZ_NTT_LIMIT;

  rz_t a, res;
  char detail[96];

  rz_init_multi(&a, &res, NULL);

  rz_gen_random(&a, (u64)limbs * 64);

  rz_mul(&res, &a, &a); /* a == b -> rz_sqr() -> NTT path */
  snprintf(detail, sizeof(detail), "rz_sqr [NTT %d limbs]", limbs);
  assert_match_gmp(detail, &a, &a, &res);

  rz_clear_multi(&a, &res, NULL);
  printf("squaring dispatch passed\n");
}

// =============================================================================
// PART 6: NTT vs KARATSUBA TIMING (for tuning RZ_NTT_LIMIT)
// =============================================================================

static void benchmark_mul_ntt_vs_karat(int limbs)
{
  rz_t a, b, res_dispatch, res_karat, res_ntt;
  rz_init_multi(&a, &b, &res_dispatch, &res_karat, &res_ntt, NULL);

  rz_gen_random(&a, (u64)limbs * 64);
  rz_gen_random(&b, (u64)limbs * 64);

  // Establish a capacity-2n buffer via the dispatch path; keep a copy
  // as the correctness reference for the raw kernel below.
  rz_mul(&res_dispatch, &a, &b);
  rz_copy(&res_karat, &res_dispatch);

  // 1. Time the explicit NTT path
  clock_t start = clock();
  rz_mul_fast(&res_ntt, &a, &b);
  clock_t end = clock();
  double time_ntt = (double)(end - start) / CLOCKS_PER_SEC;

  // 2. Time Karatsuba (raw kernel, both operands zero-padded to `limbs`)
  u64* buf = rz_scratch_get(10 * (u64)limbs);
  u64* pad_a = buf;
  u64* pad_b = buf + (u64)limbs;
  u64* scratch = buf + 2 * (u64)limbs;
  memset(pad_a, 0, (size_t)limbs * sizeof(u64));
  memset(pad_b, 0, (size_t)limbs * sizeof(u64));
  memcpy(pad_a, a.limbs, a.size * sizeof(u64));
  memcpy(pad_b, b.limbs, b.size * sizeof(u64));

  start = clock();
  limbs_mul_karatsuba(res_dispatch.limbs, pad_a, pad_b, (u64)limbs, scratch);
  end = clock();
  double time_karat = (double)(end - start) / CLOCKS_PER_SEC;
  rz_scratch_release();

  res_dispatch.size = 2 * (u64)limbs;
  res_dispatch.is_neg = false;
  rz_trim(&res_dispatch);

  // 3. Verify all three results agree
  if (rz_cmp(&res_ntt, &res_karat) != 0 ||
      rz_cmp(&res_dispatch, &res_karat) != 0) {
    printf("[FAIL] Mismatch (NTT/Karatsuba) at %d limbs!\n", limbs);
  } else {
    double speedup = (time_ntt > 0.0) ? (time_karat / time_ntt) : 0.0;
    printf("Limbs: %5d | Karatsuba: %.6fs | NTT: %.6fs | Speedup: %.2fx\n",
           limbs, time_karat, time_ntt, speedup);
  }

  rz_clear_multi(&a, &b, &res_dispatch, &res_karat, &res_ntt, NULL);
}

static void run_ntt_benchmark()
{
  printf("\n--- rz_mul: Karatsuba vs NTT (internal) ---\n");
  int sizes[] = {8192,   16384,  32768,  65536, 100000,
                 200000, 256000, 300000, 500000};
  for (size_t i = 0; i < sizeof(sizes) / sizeof(sizes[0]); i++) {
    benchmark_mul_ntt_vs_karat(sizes[i]);
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

  print_table_header();
  benchmark_mul(2048, bench_budget(2048), state);
  benchmark_mul(8192, bench_budget(8192), state);
  benchmark_mul(32768, bench_budget(32768), state);
  benchmark_mul(65536, bench_budget(65536), state);
  benchmark_mul(262144, bench_budget(262144), state);
  benchmark_mul(362144, bench_budget(362144), state);
  benchmark_mul(462144, bench_budget(462144), state);
  benchmark_mul(562144, bench_budget(562144), state);
  benchmark_mul(662144, bench_budget(662144), state);
  print_table_footer();

  run_internal_comparison();
  run_ntt_cases();
  run_ntt_asymmetric();
  run_ntt_carry();
  run_ntt_sqr();
  run_ntt_benchmark();

  gmp_randclear(state);
  rz_clear_constants();
  return 0;
}
