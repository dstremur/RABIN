/*
 * test_logic.c
 *
 * Unified GMP-verified test suite for the logic functions (bn_and, bn_or,
 * bn_xor, bn_not).
 *
 * Verified against GMP:
 *   1. Edge cases: 0, 1, 2^64-1, 2^64/2^128/2^130 boundaries, size
 *      mismatch, limb misalignment + pointer aliasing
 *   2. 1000 randomized cases (a, m: 1..4096 bits) vs mpz_and / mpz_ior /
 *      mpz_xor / width-limited complement
 *   3. Time-based benchmark vs GMP at 2048/8192/32768/65536/65539/1000000
 *      bits
 *
 * Operands are non-negative: the library applies the logic operations to
 * the magnitude limbs (preserving a's sign), which differs from GMP's
 * two's-complement logic semantics for negative operands.
 *
 * bn_not complements within the current (trimmed) limb width, so its GMP
 * oracle is (2^(64*size) - 1) XOR a rather than mpz_com.
 *
 * Copyright (C) 2026 Diego Strebel
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "../include/biglogic.h"
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

static void assert_logic_match(const char* op, const char* a_str,
                               const char* m_str, const bignum* r,
                               mpz_t expected)
{
  char* exp_str = mpz_get_str(NULL, 10, expected);
  char* custom_str = bn_to_string(r);

  if (strcmp(exp_str, custom_str) != 0) {
    fprintf(stderr, "\n[FATAL ERROR] Correctness failure in %s!\n", op);
    fprintf(stderr, "Input A: %s\n", a_str ? a_str : "(n/a)");
    fprintf(stderr, "Input M: %s\n", m_str ? m_str : "(n/a)");
    fprintf(stderr, "GMP Result:    %s\n", exp_str);
    fprintf(stderr, "Custom Result: %s\n", custom_str);
    free(exp_str);
    free(custom_str);
    exit(EXIT_FAILURE);
  }

  free(exp_str);
  free(custom_str);
}

// GMP oracle for bn_not: complement within the current limb width,
// i.e. (2^(64*size) - 1) XOR a.
static void not_oracle(mpz_t ze, mpz_t za, size_t size)
{
  if (size == 0) {
    mpz_set_ui(ze, 0);
    return;
  }

  mpz_t mask;
  mpz_init(mask);
  mpz_set_ui(mask, 1);
  mpz_mul_2exp(mask, mask, 64 * size);
  mpz_sub_ui(mask, mask, 1);
  mpz_xor(ze, mask, za);
  mpz_clear(mask);
}

// =============================================================================
// PART 1: EDGE CASES (incl. pointer aliasing)
// =============================================================================

typedef enum { OP_AND, OP_OR, OP_XOR, OP_NOT } LogicOp;

static const char* op_name(LogicOp op)
{
  switch (op) {
    case OP_AND:
      return "bn_and";
    case OP_OR:
      return "bn_or";
    case OP_XOR:
      return "bn_xor";
    case OP_NOT:
      return "bn_not";
  }
  return "?";
}

static void run_case(LogicOp op, const char* a_str, const char* m_str)
{
  bignum a, m, res;
  mpz_t za, zm, ze;
  char detail[256];

  bn_init_multi(&a, &m, &res, NULL);
  mpz_inits(za, zm, ze, NULL);

  bn_init_val(&a, a_str);
  mpz_set_str(za, a_str, 10);
  if (m_str != NULL) {
    bn_init_val(&m, m_str);
    mpz_set_str(zm, m_str, 10);
  }

  snprintf(detail, sizeof(detail), "%s a=%s m=%s", op_name(op), a_str,
           m_str ? m_str : "-");

  // non-aliased
  switch (op) {
    case OP_AND:
      bn_and(&res, &a, &m);
      mpz_and(ze, za, zm);
      break;
    case OP_OR:
      bn_or(&res, &a, &m);
      mpz_ior(ze, za, zm);
      break;
    case OP_XOR:
      bn_xor(&res, &a, &m);
      mpz_xor(ze, za, zm);
      break;
    case OP_NOT:
      bn_not(&res, &a);
      not_oracle(ze, za, a.size);
      break;
  }
  assert_logic_match(detail, a_str, m_str, &res, ze);

  // aliasing: result == a
  bn_init_val(&a, a_str);
  snprintf(detail, sizeof(detail), "%s (r==a) a=%s m=%s", op_name(op), a_str,
           m_str ? m_str : "-");
  switch (op) {
    case OP_AND:
      bn_and(&a, &a, &m);
      break;
    case OP_OR:
      bn_or(&a, &a, &m);
      break;
    case OP_XOR:
      bn_xor(&a, &a, &m);
      break;
    case OP_NOT:
      bn_not(&a, &a);
      break;
  }
  assert_logic_match(detail, a_str, m_str, &a, ze);

  // aliasing: result == m
  if (m_str != NULL) {
    bn_init_val(&a, a_str);  // r==a case left a holding the result
    bn_init_val(&m, m_str);
    snprintf(detail, sizeof(detail), "%s (r==m) a=%s m=%s", op_name(op), a_str,
             m_str);
    switch (op) {
      case OP_AND:
        bn_and(&m, &a, &m);
        break;
      case OP_OR:
        bn_or(&m, &a, &m);
        break;
      case OP_XOR:
        bn_xor(&m, &a, &m);
        break;
      case OP_NOT:
        break;
    }
    assert_logic_match(detail, a_str, m_str, &m, ze);
  }

  bn_free_multi(&a, &m, &res, NULL);
  mpz_clears(za, zm, ze, NULL);
}

static void run_edge_cases()
{
  printf("\n--- logic functions: edge cases ---\n");

  const char* zero = "0";
  const char* one = "1";
  const char* five = "5";
  const char* u64max = "18446744073709551615";          // 2^64 - 1
  const char* two64 = "18446744073709551616";           // 2^64
  const char* two65 = "36893488147419103232";           // 2^65
  const char* two96 = "79228162514264337593543950336";  // 2^96
  const char* two128m1 =
      "340282366920938463463374607431768211455";                    // 2^128 - 1
  const char* two128 = "340282366920938463463374607431768211456";   // 2^128
  const char* two130 = "1361129467683753853853498429727072845824";  // 2^130

  // AND
  run_case(OP_AND, zero, zero);
  run_case(OP_AND, zero, one);
  run_case(OP_AND, one, one);
  run_case(OP_AND, u64max, u64max);
  run_case(OP_AND, u64max, one);
  run_case(OP_AND, u64max, zero);
  run_case(OP_AND, two64, one);     // cross-limb: 2^64 & 1 = 0
  run_case(OP_AND, two64, u64max);  // 2^64 & (2^64 - 1) = 0
  run_case(OP_AND, one, two128m1);  // size mismatch
  run_case(OP_AND, two130, two65);  // 3 limbs vs 2 limbs (odd alignment)
  run_case(OP_AND, two128, two96);  // limb-2 alignment

  // OR
  run_case(OP_OR, zero, zero);
  run_case(OP_OR, zero, one);
  run_case(OP_OR, one, one);
  run_case(OP_OR, u64max, zero);
  run_case(OP_OR, u64max, u64max);
  run_case(OP_OR, two64, one);     // 2^64 | 1 = 2^64 + 1
  run_case(OP_OR, one, two128m1);  // size mismatch
  run_case(OP_OR, two130, two65);  // 3 limbs vs 2 limbs (odd alignment)
  run_case(OP_OR, two128, two96);

  // XOR
  run_case(OP_XOR, one, one);  // = 0
  run_case(OP_XOR, zero, zero);
  run_case(OP_XOR, u64max, u64max);      // = 0
  run_case(OP_XOR, two128m1, two128m1);  // = 0
  run_case(OP_XOR, two64, u64max);       // = 2^128 - 1
  run_case(OP_XOR, two128, one);         // = 2^128 + 1
  run_case(OP_XOR, one, two128m1);       // size mismatch
  run_case(OP_XOR, two130, two65);       // 3 limbs vs 2 limbs (odd alignment)
  run_case(OP_XOR, two128, two96);

  // NOT
  run_case(OP_NOT, zero, NULL);    // 1 limb of zeros -> 2^64 - 1
  run_case(OP_NOT, one, NULL);     // 2^64 - 2
  run_case(OP_NOT, u64max, NULL);  // 0
  run_case(OP_NOT, two64, NULL);   // 2^128 - 1 - 2^64
  run_case(OP_NOT, two130, NULL);  // top limb trims
  run_case(OP_NOT, five, NULL);

  printf("all edge cases passed\n");
}

// =============================================================================
// PART 2: RANDOMIZED FUZZING vs GMP
// =============================================================================

#define FUZZ_ITERATIONS 1000

static void run_random(gmp_randstate_t state)
{
  printf(
      "\n--- logic functions: %d randomized cases vs GMP (non-negative) ---\n",
      FUZZ_ITERATIONS);

  for (int i = 0; i < FUZZ_ITERATIONS; i++) {
    int bits_a = 1 + (rand() % 4096);
    int bits_m = 1 + (rand() % 4096);

    bignum a, m, res;
    mpz_t za, zm, ze;
    char* sa;
    char* sm;
    char detail[64];

    bn_init_multi(&a, &m, &res, NULL);
    mpz_inits(za, zm, ze, NULL);

    mpz_urandomb(za, state, bits_a);  // non-negative
    mpz_urandomb(zm, state, bits_m);  // non-negative
    sa = mpz_get_str(NULL, 10, za);
    sm = mpz_get_str(NULL, 10, zm);
    bn_init_val(&a, sa);
    bn_init_val(&m, sm);

    // AND
    bn_and(&res, &a, &m);
    mpz_and(ze, za, zm);
    snprintf(detail, sizeof(detail), "case %d and (a=%d bits, m=%d bits)", i,
             bits_a, bits_m);
    assert_logic_match(detail, sa, sm, &res, ze);

    // OR
    bn_or(&res, &a, &m);
    mpz_ior(ze, za, zm);
    snprintf(detail, sizeof(detail), "case %d or (a=%d bits, m=%d bits)", i,
             bits_a, bits_m);
    assert_logic_match(detail, sa, sm, &res, ze);

    // XOR
    bn_xor(&res, &a, &m);
    mpz_xor(ze, za, zm);
    snprintf(detail, sizeof(detail), "case %d xor (a=%d bits, m=%d bits)", i,
             bits_a, bits_m);
    assert_logic_match(detail, sa, sm, &res, ze);

    // NOT
    bn_not(&res, &a);
    not_oracle(ze, za, a.size);
    snprintf(detail, sizeof(detail), "case %d not (a=%d bits)", i, bits_a);
    assert_logic_match(detail, sa, NULL, &res, ze);

    free(sa);
    free(sm);
    bn_free_multi(&a, &m, &res, NULL);
    mpz_clears(za, zm, ze, NULL);
  }

  printf("all %d random cases passed\n", FUZZ_ITERATIONS);
}

// =============================================================================
// PART 3: TIME-BASED BENCHMARKS vs GMP
// =============================================================================

static void benchmark_logic(int bits, double target_sec, gmp_randstate_t state)
{
  bignum bn_a, bn_m, bn_res;
  mpz_t mpz_a, mpz_m, mpz_mask, mpz_res;

  bn_init_multi(&bn_a, &bn_m, &bn_res, NULL);
  mpz_inits(mpz_a, mpz_m, mpz_mask, mpz_res, NULL);

  mpz_urandomb(mpz_a, state, bits);
  mpz_urandomb(mpz_m, state, bits);
  bn_init_val(&bn_a, mpz_get_str(NULL, 10, mpz_a));
  bn_init_val(&bn_m, mpz_get_str(NULL, 10, mpz_m));

  // Width-limited mask for the bn_not GMP baseline: 2^(64*size) - 1
  mpz_set_ui(mpz_mask, 1);
  mpz_mul_2exp(mpz_mask, mpz_mask, 64 * bn_a.size);
  mpz_sub_ui(mpz_mask, mpz_mask, 1);

  struct timespec start, end;
  int ops_custom = 0, ops_gmp = 0;
  double total_custom = 0, total_gmp = 0;

  // AND
  clock_gettime(CLOCK_MONOTONIC, &start);
  do {
    bn_and(&bn_res, &bn_a, &bn_m);
    ops_custom++;
    clock_gettime(CLOCK_MONOTONIC, &end);
    total_custom = get_elapsed_time(start, end);
  } while (total_custom < target_sec);

  clock_gettime(CLOCK_MONOTONIC, &start);
  do {
    mpz_and(mpz_res, mpz_a, mpz_m);
    ops_gmp++;
    clock_gettime(CLOCK_MONOTONIC, &end);
    total_gmp = get_elapsed_time(start, end);
  } while (total_gmp < target_sec);

  mpz_and(mpz_res, mpz_a, mpz_m);
  assert_logic_match("bn_and (benchmark)", "(bench input)", "(bench input)",
                     &bn_res, mpz_res);

  char size_info[32];
  snprintf(size_info, sizeof(size_info), "%d bits", bits);
  print_table_row("bn_and", size_info, total_custom / ops_custom,
                  total_gmp / ops_gmp);

  ops_custom = 0;
  ops_gmp = 0;
  total_custom = 0;
  total_gmp = 0;

  // OR
  clock_gettime(CLOCK_MONOTONIC, &start);
  do {
    bn_or(&bn_res, &bn_a, &bn_m);
    ops_custom++;
    clock_gettime(CLOCK_MONOTONIC, &end);
    total_custom = get_elapsed_time(start, end);
  } while (total_custom < target_sec);

  clock_gettime(CLOCK_MONOTONIC, &start);
  do {
    mpz_ior(mpz_res, mpz_a, mpz_m);
    ops_gmp++;
    clock_gettime(CLOCK_MONOTONIC, &end);
    total_gmp = get_elapsed_time(start, end);
  } while (total_gmp < target_sec);

  mpz_ior(mpz_res, mpz_a, mpz_m);
  assert_logic_match("bn_or (benchmark)", "(bench input)", "(bench input)",
                     &bn_res, mpz_res);
  print_table_row("bn_or", size_info, total_custom / ops_custom,
                  total_gmp / ops_gmp);

  ops_custom = 0;
  ops_gmp = 0;
  total_custom = 0;
  total_gmp = 0;

  // XOR
  clock_gettime(CLOCK_MONOTONIC, &start);
  do {
    bn_xor(&bn_res, &bn_a, &bn_m);
    ops_custom++;
    clock_gettime(CLOCK_MONOTONIC, &end);
    total_custom = get_elapsed_time(start, end);
  } while (total_custom < target_sec);

  clock_gettime(CLOCK_MONOTONIC, &start);
  do {
    mpz_xor(mpz_res, mpz_a, mpz_m);
    ops_gmp++;
    clock_gettime(CLOCK_MONOTONIC, &end);
    total_gmp = get_elapsed_time(start, end);
  } while (total_gmp < target_sec);

  mpz_xor(mpz_res, mpz_a, mpz_m);
  assert_logic_match("bn_xor (benchmark)", "(bench input)", "(bench input)",
                     &bn_res, mpz_res);
  print_table_row("bn_xor", size_info, total_custom / ops_custom,
                  total_gmp / ops_gmp);

  ops_custom = 0;
  ops_gmp = 0;
  total_custom = 0;
  total_gmp = 0;

  // NOT (GMP baseline: width-limited complement via XOR with the mask)
  clock_gettime(CLOCK_MONOTONIC, &start);
  do {
    bn_not(&bn_res, &bn_a);
    ops_custom++;
    clock_gettime(CLOCK_MONOTONIC, &end);
    total_custom = get_elapsed_time(start, end);
  } while (total_custom < target_sec);

  clock_gettime(CLOCK_MONOTONIC, &start);
  do {
    mpz_xor(mpz_res, mpz_mask, mpz_a);
    ops_gmp++;
    clock_gettime(CLOCK_MONOTONIC, &end);
    total_gmp = get_elapsed_time(start, end);
  } while (total_gmp < target_sec);

  mpz_xor(mpz_res, mpz_mask, mpz_a);
  assert_logic_match("bn_not (benchmark)", "(bench input)", "(n/a)", &bn_res,
                     mpz_res);
  print_table_row("bn_not", size_info, total_custom / ops_custom,
                  total_gmp / ops_gmp);

  bn_free_multi(&bn_a, &bn_m, &bn_res, NULL);
  mpz_clears(mpz_a, mpz_m, mpz_mask, mpz_res, NULL);
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
  if (bits <= 131072) return 2.0;
  return 4.0;
}

int main()
{
  gmp_randstate_t state;
  gmp_randinit_default(state);
  gmp_randseed_ui(state, (unsigned long)time(NULL));
  srand((unsigned int)time(NULL));

  printf("AVX-512 path: %s\n", bn_supports_avx512() ? "enabled" : "disabled");

  run_edge_cases();
  run_random(state);

  print_table_header();
  /*  benchmark_logic(2048, bench_budget(2048), state);
    benchmark_logic(8192, bench_budget(8192), state);
    benchmark_logic(32768, bench_budget(32768), state);
    benchmark_logic(65536, bench_budget(65536), state);
    benchmark_logic(65539, bench_budget(65539), state);
    */
  benchmark_logic(1000000, bench_budget(1000000), state);
  benchmark_logic(10000000, bench_budget(10000000), state);
  print_table_footer();

  gmp_randclear(state);
  return 0;
}
