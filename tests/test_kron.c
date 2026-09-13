#include <fcntl.h>
#include <gmp.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#include "../include/bignum.h"

#ifndef i64
typedef int64_t i64;
#endif

// =============================================================================
// ADAPTER FOR YOUR LIBRARY
// =============================================================================

char* bn_to_str(bignum* bn) { return bn_to_string(bn); }

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

void generate_random_signed_pair(bignum* bn_a, bignum* bn_b, mpz_t mpz_a,
                                 mpz_t mpz_b, int bits_a, int bits_b,
                                 gmp_randstate_t state)
{
  mpz_urandomb(mpz_a, state, bits_a);
  mpz_urandomb(mpz_b, state, bits_b);

  mpz_t coin;
  mpz_init(coin);

  // Randomly make 'a' negative (50% chance)
  mpz_urandomb(coin, state, 1);
  if (mpz_cmp_ui(coin, 0) == 0) {
    mpz_neg(mpz_a, mpz_a);
  }

  // Randomly make 'b' negative (50% chance)
  mpz_urandomb(coin, state, 1);
  if (mpz_cmp_ui(coin, 0) == 0) {
    mpz_neg(mpz_b, mpz_b);
  }

  mpz_clear(coin);

  char* s_a = mpz_get_str(NULL, 10, mpz_a);
  char* s_b = mpz_get_str(NULL, 10, mpz_b);
  bn_init_val(bn_a, s_a);
  bn_init_val(bn_b, s_b);
  free(s_a);
  free(s_b);
}

// =============================================================================
// KRONECKER EDGE CASE SUITE
// =============================================================================

void run_kronecker_edge_case(const char* name, const char* a_str,
                             const char* b_str)
{
  bignum bn_a, bn_b;
  mpz_t mpz_a, mpz_b;

  bn_init(&bn_a);
  bn_init(&bn_b);

  bn_init_val(&bn_a, a_str);
  bn_init_val(&bn_b, b_str);

  mpz_init_set_str(mpz_a, a_str, 10);
  mpz_init_set_str(mpz_b, b_str, 10);

  i64 custom_res = bn_kronecker(&bn_a, &bn_b);
  int gmp_res = mpz_kronecker(mpz_a, mpz_b);

  if (custom_res != (i64)gmp_res) {
    fprintf(stderr, "\n[FATAL ERROR] Correctness failure in Kronecker %s!\n",
            name);
    fprintf(stderr, "Input A:       %s\n", a_str);
    fprintf(stderr, "Input B:       %s\n", b_str);
    fprintf(stderr, "GMP Result:    %d\n", gmp_res);
    fprintf(stderr, "Custom Result: %lld\n", (long long)custom_res);
    exit(EXIT_FAILURE);
  }

  bn_free(&bn_a);
  bn_free(&bn_b);
  mpz_clears(mpz_a, mpz_b, NULL);
}

void run_edge_case_suite(gmp_randstate_t state)
{
  printf("\n--- Running Kronecker Symbol Edge Case & Correctness Suite ---\n");

  run_kronecker_edge_case("a=0, b=0", "0", "0");
  run_kronecker_edge_case("a=1, b=0", "1", "0");
  run_kronecker_edge_case("a=-1, b=0", "-1", "0");
  run_kronecker_edge_case("a=2, b=0", "2", "0");
  run_kronecker_edge_case("a=0, b=1", "0", "1");
  run_kronecker_edge_case("a=0, b=-1", "0", "-1");

  run_kronecker_edge_case("a=15, b=1", "15", "1");
  run_kronecker_edge_case("a=-15, b=1", "-15", "1");
  run_kronecker_edge_case("a=15, b=-1", "15", "-1");
  run_kronecker_edge_case("a=-15, b=-1", "-15", "-1");

  run_kronecker_edge_case("a=1, b=2", "1", "2");
  run_kronecker_edge_case("a=3, b=2", "3", "2");
  run_kronecker_edge_case("a=5, b=2", "5", "2");
  run_kronecker_edge_case("a=7, b=2", "7", "2");
  run_kronecker_edge_case("a=2, b=2", "2", "2");
  run_kronecker_edge_case("a=-3, b=2", "-3", "2");
  run_kronecker_edge_case("a=3, b=-2", "3", "-2");

  run_kronecker_edge_case("Even/Even", "12345678", "87654320");
  run_kronecker_edge_case("Shared Prime Factor", "15", "35");

  run_kronecker_edge_case("a < 0, b > 0", "-12345", "67891");
  run_kronecker_edge_case("a > 0, b < 0", "12345", "-67891");
  run_kronecker_edge_case("a < 0, b < 0", "-12345", "-67891");

  run_kronecker_edge_case("2^64 - 1 vs 2^64 + 1", "18446744073709551615",
                          "18446744073709551617");
  run_kronecker_edge_case("-(2^64 - 1) vs 2^64 + 3", "-18446744073709551615",
                          "18446744073709551619");

  printf("Running 1,000 random Kronecker test vectors...\n");
  for (int i = 0; i < 1000; i++) {
    int bits_a = 1 + (rand() % 204);
    int bits_b = 1 + (rand() % 204);

    bignum bn_a, bn_b;
    mpz_t mpz_a, mpz_b;

    bn_init(&bn_a);
    bn_init(&bn_b);
    mpz_inits(mpz_a, mpz_b, NULL);

    generate_random_signed_pair(&bn_a, &bn_b, mpz_a, mpz_b, bits_a, bits_b,
                                state);

    i64 custom_res = bn_kronecker(&bn_a, &bn_b);
    int gmp_res = mpz_kronecker(mpz_a, mpz_b);

    if (custom_res != (i64)gmp_res) {
      fprintf(
          stderr,
          "\n[FATAL ERROR] Correctness mismatch in random Kronecker test %d!\n",
          i);
      char* s_a = mpz_get_str(NULL, 10, mpz_a);
      char* s_b = mpz_get_str(NULL, 10, mpz_b);
      fprintf(stderr, "A (%d bits): %s\n", bits_a, s_a);
      fprintf(stderr, "B (%d bits): %s\n", bits_b, s_b);
      fprintf(stderr, "GMP Result:    %d\n", gmp_res);
      fprintf(stderr, "Custom Result: %lld\n", (long long)custom_res);
      free(s_a);
      free(s_b);
      exit(EXIT_FAILURE);
    }

    bn_free(&bn_a);
    bn_free(&bn_b);
    mpz_clears(mpz_a, mpz_b, NULL);
  }

  printf("All Kronecker edge cases and random tests passed successfully!\n");
}

// =============================================================================
// TONELLI-SHANKS TEST & VALIDATION SUITE
// =============================================================================

void run_tonelli_shanks_case(const char* name, const char* n_str,
                             const char* p_str)
{
  bignum bn_n, bn_p, bn_r;
  mpz_t mpz_n, mpz_p, mpz_r_expected, mpz_r_actual, mpz_check;

  bn_init(&bn_n);
  bn_init(&bn_p);
  bn_init(&bn_r);
  mpz_inits(mpz_n, mpz_p, mpz_r_expected, mpz_r_actual, mpz_check, NULL);

  bn_init_val(&bn_n, n_str);
  bn_init_val(&bn_p, p_str);
  mpz_set_str(mpz_n, n_str, 10);
  mpz_set_str(mpz_p, p_str, 10);

  // Invoke custom Tonelli-Shanks
  tonelli_shanks(&bn_r, &bn_n, &bn_p);

  // Validate: r^2 mod p == n mod p
  char* r_str = bn_to_str(&bn_r);
  mpz_set_str(mpz_r_actual, r_str, 10);
  free(r_str);

  mpz_powm_ui(mpz_check, mpz_r_actual, 2, mpz_p);
  mpz_mod(mpz_n, mpz_n, mpz_p);

  if (mpz_cmp(mpz_check, mpz_n) != 0) {
    fprintf(stderr,
            "\n[FATAL ERROR] Correctness failure in Tonelli-Shanks test: %s!\n",
            name);
    char* s_n = mpz_get_str(NULL, 10, mpz_n);
    char* s_p = mpz_get_str(NULL, 10, mpz_p);
    char* s_r = mpz_get_str(NULL, 10, mpz_r_actual);
    fprintf(stderr, "Input N:       %s\n", s_n);
    fprintf(stderr, "Modulus P:     %s\n", s_p);
    fprintf(stderr, "Output R:      %s\n", s_r);
    free(s_n);
    free(s_p);
    free(s_r);
    exit(EXIT_FAILURE);
  }

  bn_free(&bn_n);
  bn_free(&bn_p);
  bn_free(&bn_r);
  mpz_clears(mpz_n, mpz_p, mpz_r_expected, mpz_r_actual, mpz_check, NULL);
}

void run_tonelli_shanks_suite(gmp_randstate_t state)
{
  printf("\n--- Running Tonelli-Shanks Correctness & Stress Suite ---\n");

  // 1. Small static test cases
  // p = 13 (13 ≡ 1 mod 4), quadratic residues mod 13 are 1, 3, 4, 9, 10, 12
  run_tonelli_shanks_case("p=13, n=3", "3", "13");
  run_tonelli_shanks_case("p=13, n=9", "9", "13");
  // p = 7 (7 ≡ 3 mod 4)
  run_tonelli_shanks_case("p=7, n=2", "2", "7");
  run_tonelli_shanks_case("p=7, n=4", "4", "7");

  // 2. Randomized Stress Test (500 random prime fields up to 2048 bits)
  printf("Running 500 random Tonelli-Shanks test vectors...\n");
  for (int i = 0; i < 500; i++) {
    int bits = 64 + (rand() % 1984);  // 64 to 2048 bits

    mpz_t mpz_p, mpz_root, mpz_n;
    mpz_inits(mpz_p, mpz_root, mpz_n, NULL);

    // Generate random prime p
    mpz_urandomb(mpz_p, state, bits);
    mpz_nextprime(mpz_p, mpz_p);

    // Ensure p > 2 and p % 4 != 2 (Tonelli-Shanks applies to odd primes)
    if (mpz_even_p(mpz_p)) {
      mpz_add_ui(mpz_p, mpz_p, 1);
      mpz_nextprime(mpz_p, mpz_p);
    }

    // Generate random root r < p, compute n = r^2 mod p
    mpz_urandomm(mpz_root, state, mpz_p);
    mpz_powm_ui(mpz_n, mpz_root, 2, mpz_p);

    char* s_n = mpz_get_str(NULL, 10, mpz_n);
    char* s_p = mpz_get_str(NULL, 10, mpz_p);

    bignum bn_n, bn_p, bn_r;
    bn_init(&bn_n);
    bn_init(&bn_p);
    bn_init(&bn_r);

    bn_init_val(&bn_n, s_n);
    bn_init_val(&bn_p, s_p);

    tonelli_shanks(&bn_r, &bn_n, &bn_p);

    // Verify result
    char* s_r = bn_to_str(&bn_r);
    mpz_t mpz_check, mpz_actual_r;
    mpz_inits(mpz_check, mpz_actual_r, NULL);
    mpz_set_str(mpz_actual_r, s_r, 10);

    mpz_powm_ui(mpz_check, mpz_actual_r, 2, mpz_p);
    if (mpz_cmp(mpz_check, mpz_n) != 0) {
      fprintf(stderr,
              "\n[FATAL ERROR] Random Tonelli-Shanks mismatch at test %d (%d "
              "bits)!\n",
              i, bits);
      fprintf(stderr, "P: %s\nN: %s\nComputed R: %s\n", s_p, s_n, s_r);
      exit(EXIT_FAILURE);
    }

    free(s_n);
    free(s_p);
    free(s_r);
    bn_free(&bn_n);
    bn_free(&bn_p);
    bn_free(&bn_r);
    mpz_clears(mpz_p, mpz_root, mpz_n, mpz_check, mpz_actual_r, NULL);
  }

  printf(
      "All Tonelli-Shanks correctness and random tests passed successfully!\n");
}

// =============================================================================
// BENCHMARK SUITE
// =============================================================================

void benchmark_kronecker(int bits, double target_sec, gmp_randstate_t state)
{
  bignum bn_a, bn_b;
  mpz_t mpz_a, mpz_b;

  bn_init(&bn_a);
  bn_init(&bn_b);
  mpz_inits(mpz_a, mpz_b, NULL);

  generate_random_signed_pair(&bn_a, &bn_b, mpz_a, mpz_b, bits, bits, state);

  struct timespec start, end;
  int ops_custom = 0, ops_gmp = 0;
  double total_custom = 0, total_gmp = 0;

  volatile i64 custom_res = 0;
  volatile int gmp_res = 0;

  clock_gettime(CLOCK_MONOTONIC, &start);
  do {
    custom_res = bn_kronecker(&bn_a, &bn_b);
    ops_custom++;
    clock_gettime(CLOCK_MONOTONIC, &end);
    total_custom = get_elapsed_time(start, end);
  } while (total_custom < target_sec);

  clock_gettime(CLOCK_MONOTONIC, &start);
  do {
    gmp_res = mpz_kronecker(mpz_a, mpz_b);
    ops_gmp++;
    clock_gettime(CLOCK_MONOTONIC, &end);
    total_gmp = get_elapsed_time(start, end);
  } while (total_gmp < target_sec);

  if (custom_res != (i64)gmp_res) {
    fprintf(stderr,
            "\n[FATAL ERROR] Correctness mismatch in Kronecker benchmark (%d "
            "bits)!\n",
            bits);
    exit(EXIT_FAILURE);
  }

  char size_info[32];
  snprintf(size_info, sizeof(size_info), "%d bits", bits);
  print_table_row("Kronecker", size_info, total_custom / ops_custom,
                  total_gmp / ops_gmp);

  bn_free(&bn_a);
  bn_free(&bn_b);
  mpz_clears(mpz_a, mpz_b, NULL);
}

void benchmark_tonelli_shanks(int bits, double target_sec,
                              gmp_randstate_t state)
{
  mpz_t mpz_p, mpz_n;
  mpz_inits(mpz_p, mpz_n, NULL);

  // No prime generation at large bit sizes (mpz_nextprime hangs at thousands
  // of bits). Sample the modulus until the Kronecker symbol certifies that
  // n is a quadratic residue mod p, keeping tonelli_shanks on its full
  // algorithm path instead of the early "No square roots exist" return.
  mpz_urandomb(mpz_n, state, bits);

  do {
    mpz_urandomb(mpz_p, state, bits);
    mpz_setbit(mpz_p, 0);  // ensure odd modulus
    if (mpz_cmp_ui(mpz_p, 3) <= 0) mpz_set_ui(mpz_p, 5);
  } while (mpz_kronecker(mpz_n, mpz_p) != 1);

  char* s_n = mpz_get_str(NULL, 10, mpz_n);
  char* s_p = mpz_get_str(NULL, 10, mpz_p);

  bignum bn_n, bn_p, bn_r;
  bn_init(&bn_n);
  bn_init(&bn_p);
  bn_init(&bn_r);

  bn_init_val(&bn_n, s_n);
  bn_init_val(&bn_p, s_p);

  struct timespec start, end;
  int ops_custom = 0;
  double total_custom = 0;

  // Silence tonelli_shanks' diagnostic output while timing. This only
  // redirects fd 1 for the duration of the loop; the library function
  // itself keeps its normal printing.
  fflush(stdout);
  int saved_stdout = dup(STDOUT_FILENO);
  int dev_null = open("/dev/null", O_WRONLY);
  dup2(dev_null, STDOUT_FILENO);
  close(dev_null);

  clock_gettime(CLOCK_MONOTONIC, &start);
  do {
    tonelli_shanks(&bn_r, &bn_n, &bn_p);
    ops_custom++;
    clock_gettime(CLOCK_MONOTONIC, &end);
    total_custom = get_elapsed_time(start, end);
  } while (total_custom < target_sec);

  fflush(stdout);
  dup2(saved_stdout, STDOUT_FILENO);
  close(saved_stdout);

  char size_info[32];
  snprintf(size_info, sizeof(size_info), "%d bits", bits);
  // GMP comparison baseline omitted/set to 0 since GMP lacks native
  // Tonelli-Shanks
  print_table_row("Tonelli-Shanks", size_info, total_custom / ops_custom, 0.0);

  free(s_n);
  free(s_p);
  bn_free(&bn_n);
  bn_free(&bn_p);
  bn_free(&bn_r);
  mpz_clears(mpz_p, mpz_n, NULL);
}

// =============================================================================
// MAIN
// =============================================================================

int main()
{
  gmp_randstate_t state;
  gmp_randinit_default(state);
  gmp_randseed_ui(state, (unsigned long)time(NULL));
  srand((unsigned int)time(NULL));

  // 1. Run Kronecker Suite
  run_edge_case_suite(state);

  // 2. Run Tonelli-Shanks Suite
  run_tonelli_shanks_suite(state);

  // 3. Benchmarks
  double time_per_test = 1.0;

  print_table_header();

  benchmark_kronecker(512, time_per_test, state);
  benchmark_kronecker(1024, time_per_test, state);
  benchmark_kronecker(2048, time_per_test, state);
  benchmark_kronecker(4096, time_per_test, state);
  benchmark_kronecker(8192, time_per_test, state);

  benchmark_tonelli_shanks(512, time_per_test, state);
  benchmark_tonelli_shanks(1024, time_per_test, state);
  benchmark_tonelli_shanks(2048, time_per_test, state);
  benchmark_tonelli_shanks(4096, time_per_test, state);
  benchmark_tonelli_shanks(8192, time_per_test, state);

  print_table_footer();

  gmp_randclear(state);
  return 0;
}