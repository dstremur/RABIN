#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "../include/bignum.h"

// =============================================================================
// ADAPTER FOR YOUR LIBRARY
// =============================================================================

// IMPORTANT: Replace the inside of this function with your library's actual
// string export method. This is required for correctness checking!
char* bn_to_str(bignum* bn)
{
  // Example: return bn_get_string(bn, 10);
  // For now, returning a dummy string to prevent compilation failure if
  // missing. If your library does not have a string export, you will need to
  // write one.
  return bn_to_string(bn);
}

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

// Verifies that your bignum matches GMP's result exactly.
void assert_match(bignum* bn, mpz_t mpz, const char* context)
{
  char* gmp_str = mpz_get_str(NULL, 10, mpz);
  char* custom_str = bn_to_str(bn);

  if (strcmp(gmp_str, custom_str) != 0) {
    fprintf(stderr, "\n[FATAL ERROR] Correctness failure in %s!\n", context);
    fprintf(stderr, "GMP Result:    %s\n", gmp_str);
    fprintf(stderr, "Custom Result: %s\n", custom_str);
    exit(EXIT_FAILURE);
  }

  free(gmp_str);
  free(custom_str);
}

void generate_random_pair(bignum* bn_a, bignum* bn_b, mpz_t mpz_a, mpz_t mpz_b,
                          int bits_a, int bits_b, gmp_randstate_t state)
{
  mpz_urandomb(mpz_a, state, bits_a);
  mpz_urandomb(mpz_b, state, bits_b);

  // Ensure b != 0 for division
  if (mpz_cmp_ui(mpz_b, 0) == 0) mpz_set_ui(mpz_b, 1);

  // Ensure a > b for safe subtraction in unsigned libs
  if (mpz_cmp(mpz_a, mpz_b) < 0) mpz_swap(mpz_a, mpz_b);

  char* s_a = mpz_get_str(NULL, 10, mpz_a);
  char* s_b = mpz_get_str(NULL, 10, mpz_b);
  bn_init_val(bn_a, s_a);
  bn_init_val(bn_b, s_b);
  free(s_a);
  free(s_b);
}

// =============================================================================
// EDGE CASE SUITE
// =============================================================================

typedef enum { OP_ADD, OP_SUB, OP_MUL, OP_DIV } OpType;

void run_edge_case(const char* name, OpType op, const char* a_str,
                   const char* b_str)
{
  bignum bn_a, bn_b, bn_res;
  mpz_t mpz_a, mpz_b, mpz_res;

  bn_init(&bn_a);
  bn_init(&bn_b);
  bn_init(&bn_res);

  bn_init_val(&bn_a, a_str);
  bn_init_val(&bn_b, b_str);
  bn_init(&bn_res);

  mpz_init_set_str(mpz_a, a_str, 10);
  mpz_init_set_str(mpz_b, b_str, 10);
  mpz_init(mpz_res);

  switch (op) {
    case OP_ADD:
      bn_add(&bn_res, &bn_a, &bn_b);
      mpz_add(mpz_res, mpz_a, mpz_b);
      break;
    case OP_SUB:
      bn_sub(&bn_res, &bn_a, &bn_b);
      mpz_sub(mpz_res, mpz_a, mpz_b);
      break;
    case OP_MUL:
      bn_mul(&bn_res, &bn_a, &bn_b);
      mpz_mul(mpz_res, mpz_a, mpz_b);
      break;
    case OP_DIV:
      bn_div(&bn_res, &bn_a, &bn_b);
      mpz_tdiv_q(mpz_res, mpz_a, mpz_b);
      break;
  }

  assert_match(&bn_res, mpz_res, name);

  bn_free(&bn_a);
  bn_free(&bn_b);
  bn_free(&bn_res);
  mpz_clears(mpz_a, mpz_b, mpz_res, NULL);
}

void run_edge_case_suite()
{
  printf("\n--- Running Correctness & Edge Case Suite ---\n");

  // 1. Zeros
  run_edge_case("Zero Add", OP_ADD, "0", "0");
  run_edge_case("Multiply by Zero", OP_MUL, "987654321987654321", "0");
  run_edge_case("Divide by One", OP_DIV, "987654321987654321", "1");

  // 2. Carry Propagation (Power of 2 boundaries)
  // 2^64 - 1 + 1 (Tests word-boundary cascade)
  run_edge_case("Carry Propagation (64-bit)", OP_ADD, "18446744073709551615",
                "1");

  // 3. Asymmetric Math (Massive discrepancy in sizes)
  char massive_a[2000];
  memset(massive_a, '9', 1999);
  massive_a[1999] = '\0';
  run_edge_case("Asymmetric Multiply", OP_MUL, massive_a, "2");
  run_edge_case("Asymmetric Divide", OP_DIV, massive_a, "2");

  printf("All edge cases passed successfully!\n");
}

// =============================================================================
// BENCHMARKS (Time-Bound)
// =============================================================================

#define MAKE_BENCHMARK(func_name, op_name, custom_func, gmp_func)        \
  void func_name(int bits, double target_sec, gmp_randstate_t state)     \
  {                                                                      \
    bignum bn_a, bn_b, bn_res;                                           \
    mpz_t mpz_a, mpz_b, mpz_res;                                         \
    bn_init(&bn_a);                                                      \
    bn_init(&bn_b);                                                      \
    bn_init(&bn_res);                                                    \
    mpz_inits(mpz_a, mpz_b, mpz_res, NULL);                              \
    generate_random_pair(&bn_a, &bn_b, mpz_a, mpz_b, bits, bits, state); \
                                                                         \
    struct timespec start, end;                                          \
    int ops_custom = 0, ops_gmp = 0;                                     \
    double total_custom = 0, total_gmp = 0;                              \
                                                                         \
    /* Benchmark Custom */                                               \
    clock_gettime(CLOCK_MONOTONIC, &start);                              \
    do {                                                                 \
      custom_func(&bn_res, &bn_a, &bn_b);                                \
      ops_custom++;                                                      \
      clock_gettime(CLOCK_MONOTONIC, &end);                              \
      total_custom = get_elapsed_time(start, end);                       \
    } while (total_custom < target_sec);                                 \
                                                                         \
    /* Benchmark GMP */                                                  \
    clock_gettime(CLOCK_MONOTONIC, &start);                              \
    do {                                                                 \
      gmp_func(mpz_res, mpz_a, mpz_b);                                   \
      ops_gmp++;                                                         \
      clock_gettime(CLOCK_MONOTONIC, &end);                              \
      total_gmp = get_elapsed_time(start, end);                          \
    } while (total_gmp < target_sec);                                    \
                                                                         \
    /* Validate Correctness before reporting */                          \
    assert_match(&bn_res, mpz_res, op_name);                             \
                                                                         \
    char size_info[32];                                                  \
    snprintf(size_info, sizeof(size_info), "%d bits", bits);             \
    print_table_row(op_name, size_info, total_custom / ops_custom,       \
                    total_gmp / ops_gmp);                                \
                                                                         \
    bn_free(&bn_a);                                                      \
    bn_free(&bn_b);                                                      \
    bn_free(&bn_res);                                                    \
    mpz_clears(mpz_a, mpz_b, mpz_res, NULL);                             \
  }

MAKE_BENCHMARK(benchmark_add, "Addition", bn_add, mpz_add)
MAKE_BENCHMARK(benchmark_sub, "Subtraction", bn_sub, mpz_sub)
MAKE_BENCHMARK(benchmark_mul, "Multiply", bn_mul, mpz_mul)
MAKE_BENCHMARK(benchmark_div, "Division", bn_div, mpz_tdiv_q)

MAKE_BENCHMARK(benchmark_gcd, "GCD", bn_gcd_lehmer, mpz_gcd)

void benchmark_mod_exp(int bits, double target_sec, gmp_randstate_t state)
{
  bignum bn_a, bn_b, bn_m, bn_res;
  mpz_t mpz_a, mpz_b, mpz_m, mpz_res;

  bn_init(&bn_a);
  bn_init(&bn_b);
  bn_init(&bn_m);
  bn_init(&bn_res);
  mpz_inits(mpz_a, mpz_b, mpz_m, mpz_res, NULL);

  mpz_urandomb(mpz_a, state, bits);
  mpz_urandomb(mpz_b, state, bits);
  mpz_urandomb(mpz_m, state, bits);
  mpz_setbit(mpz_m, 0);  // Ensure odd modulus

  char* s_a = mpz_get_str(NULL, 10, mpz_a);
  char* s_b = mpz_get_str(NULL, 10, mpz_b);
  char* s_m = mpz_get_str(NULL, 10, mpz_m);
  bn_init_val(&bn_a, s_a);
  bn_init_val(&bn_b, s_b);
  bn_init_val(&bn_m, s_m);
  free(s_a);
  free(s_b);
  free(s_m);

  struct timespec start, end;
  int ops_custom = 0, ops_gmp = 0;
  double total_custom = 0, total_gmp = 0;

  clock_gettime(CLOCK_MONOTONIC, &start);
  do {
    bn_mod_exp(&bn_res, &bn_a, &bn_b, &bn_m);
    ops_custom++;
    clock_gettime(CLOCK_MONOTONIC, &end);
    total_custom = get_elapsed_time(start, end);
  } while (total_custom < target_sec);

  clock_gettime(CLOCK_MONOTONIC, &start);
  do {
    mpz_powm(mpz_res, mpz_a, mpz_b, mpz_m);
    ops_gmp++;
    clock_gettime(CLOCK_MONOTONIC, &end);
    total_gmp = get_elapsed_time(start, end);
  } while (total_gmp < target_sec);

  assert_match(&bn_res, mpz_res, "ModExp");

  char size_info[32];
  snprintf(size_info, sizeof(size_info), "%d bits", bits);
  print_table_row("ModExp", size_info, total_custom / ops_custom,
                  total_gmp / ops_gmp);

  bn_free(&bn_a);
  bn_free(&bn_b);
  bn_free(&bn_m);
  bn_free(&bn_res);
  mpz_clears(mpz_a, mpz_b, mpz_m, mpz_res, NULL);
}

void benchmark_bpsw(int bits, double target_sec, gmp_randstate_t state)
{
  bignum bn_n;
  mpz_t mpz_n;

  bn_init(&bn_n);
  mpz_init(mpz_n);

  // Generate a random probable prime or odd number of the given bit length
  mpz_urandomb(mpz_n, state, bits);
  mpz_setbit(mpz_n, 0);  // Ensure odd number

  char* s_n = mpz_get_str(NULL, 10, mpz_n);
  bn_init_val(&bn_n, s_n);
  free(s_n);

  struct timespec start, end;
  int ops_custom = 0, ops_gmp = 0;
  double total_custom = 0, total_gmp = 0;

  // Benchmark Custom BPSW test
  clock_gettime(CLOCK_MONOTONIC, &start);
  do {
    // bn_bpsw returns bool: true if probable prime, false if composite
    volatile bool custom_res = bn_bpsw(&bn_n);
    (void)custom_res;
    ops_custom++;
    clock_gettime(CLOCK_MONOTONIC, &end);
    total_custom = get_elapsed_time(start, end);
  } while (total_custom < target_sec);

  // Benchmark GMP probabilistic primality test (mpz_probab_prime_p with 25
  // rounds as a standard proxy)
  clock_gettime(CLOCK_MONOTONIC, &start);
  do {
    volatile int gmp_res = mpz_probab_prime_p(mpz_n, 25);
    (void)gmp_res;
    ops_gmp++;
    clock_gettime(CLOCK_MONOTONIC, &end);
    total_gmp = get_elapsed_time(start, end);
  } while (total_gmp < target_sec);

  // Correctness Cross-Check: Both should agree on primality status
  bool custom_prime = trialdiv(&bn_n, 10000) && bn_bpsw(&bn_n);
  int gmp_prime = mpz_probab_prime_p(mpz_n, 25);

  if ((custom_prime && gmp_prime == 0) || (!custom_prime && gmp_prime > 0)) {
    fprintf(stderr, "\n[FATAL ERROR] Correctness mismatch in BPSW test!\n");
    exit(EXIT_FAILURE);
  }

  char size_info[32];
  snprintf(size_info, sizeof(size_info), "%d bits", bits);
  print_table_row("BPSW Primality", size_info, total_custom / ops_custom,
                  total_gmp / ops_gmp);

  bn_free(&bn_n);
  mpz_clear(mpz_n);
}

// =============================================================================
// MAIN
// =============================================================================

int main()
{
  gmp_randstate_t state;
  gmp_randinit_default(state);
  gmp_randseed_ui(state, (unsigned long)time(NULL));

  // 1. Validate correctness on nasty edge cases first
  run_edge_case_suite();

  // Target time per test in seconds (0.5s is usually enough for a solid
  // average)
  double time_per_test = 1.0;

  print_table_header();

  // -- Massive Addition & Subtraction --
  // We go up to 1 million bits (approx 125 KB per number)
  benchmark_add(4096, time_per_test, state);
  benchmark_add(65536, time_per_test, state);
  benchmark_add(1000000, time_per_test, state);

  benchmark_sub(4096, time_per_test, state);
  benchmark_sub(65536, time_per_test, state);

  // -- Heavy Multiplication & Division --
  // If you are using O(n^2) schoolbook math, 65536 bits will be brutally slow.
  // The time-bound loop ensures it will just run 1 or 2 iterations and finish.
  benchmark_mul(2048, time_per_test, state);
  benchmark_mul(8192, time_per_test, state);
  benchmark_mul(32768, time_per_test, state);
  benchmark_mul(65536, time_per_test, state);

  benchmark_div(2048, time_per_test, state);
  benchmark_div(8192, time_per_test, state);
  benchmark_div(32768, time_per_test, state);
  benchmark_div(65536, time_per_test, state);
  benchmark_div(65539, time_per_test, state);

  // -- Mod Exp (Cryptographic sizes) --
  benchmark_mod_exp(1024, time_per_test, state);
  benchmark_mod_exp(2048, time_per_test, state);
  benchmark_mod_exp(4096, time_per_test, state);
  benchmark_mod_exp(8192, time_per_test, state);
  benchmark_mod_exp(10000, time_per_test, state);
  benchmark_mod_exp(16384, time_per_test, state);

  benchmark_bpsw(1024, time_per_test, state);
  benchmark_bpsw(2048, time_per_test, state);
  benchmark_bpsw(4096, time_per_test, state);
  benchmark_bpsw(8192, time_per_test, state);
  benchmark_bpsw(10000, time_per_test, state);
  benchmark_bpsw(16384, time_per_test, state);
  benchmark_bpsw(32768, time_per_test, state);
  benchmark_bpsw(65536, time_per_test, state);
  // benchmark_bpsw(250000, time_per_test, state);

  benchmark_gcd(2048, time_per_test, state);
  benchmark_gcd(8192, time_per_test, state);
  benchmark_gcd(32768, time_per_test, state);
  benchmark_gcd(65536, time_per_test, state);
  benchmark_gcd(65539, time_per_test, state);
  benchmark_gcd(100000, time_per_test, state);
  benchmark_gcd(200000, time_per_test, state);
  benchmark_gcd(500000, time_per_test, state);

  print_table_footer();

  gmp_randclear(state);
  return 0;
}
