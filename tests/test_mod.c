#include <gmp.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

/* Adjust based on your header layout */
#include "../include/bigcore.h"
#include "../include/bighelper.h"

#ifndef i64
typedef int64_t i64;
#endif

#ifndef u64
typedef uint64_t u64;
#endif

// =============================================================================
// HELPERS
// =============================================================================

char* bn_to_str(const bignum* bn)
{
  return bn_to_string(bn);
}  // Map to your library

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

void generate_random_bignum(bignum* bn_a, mpz_t mpz_a, int bits, bool allow_neg,
                            gmp_randstate_t state)
{
  mpz_urandomb(mpz_a, state, bits);
  if (allow_neg) {
    mpz_t coin;
    mpz_init(coin);
    mpz_urandomb(coin, state, 1);
    if (mpz_cmp_ui(coin, 0) == 0) {
      mpz_neg(mpz_a, mpz_a);
    }
    mpz_clear(coin);
  }
  char* s_a = mpz_get_str(NULL, 10, mpz_a);
  bn_init_val(bn_a, s_a);
  free(s_a);
}

// =============================================================================
// EDGE CASE SUITE
// =============================================================================

void test_mod_u64(const char* name, const char* a_str, uint64_t d)
{
  bignum bn_a;
  mpz_t mpz_a, mpz_d, mpz_rem;

  bn_init(&bn_a);
  bn_init_val(&bn_a, a_str);

  mpz_inits(mpz_a, mpz_d, mpz_rem, NULL);
  mpz_set_str(mpz_a, a_str, 10);
  set_mpz_u64(mpz_d, d);

  // GMP floor modulo matches bn_mod_u64 logic
  mpz_mod(mpz_rem, mpz_a, mpz_d);
  uint64_t gmp_res = mpz_to_u64(mpz_rem);
  uint64_t custom_res = bn_mod_u64(&bn_a, d);

  if (custom_res != gmp_res) {
    fprintf(stderr, "\n[FATAL ERROR] mod_u64 mismatch (%s)!\n", name);
    fprintf(stderr, "A: %s, D: %llu\nExpected: %llu, Got: %llu\n", a_str,
            (unsigned long long)d, (unsigned long long)gmp_res,
            (unsigned long long)custom_res);
    exit(EXIT_FAILURE);
  }

  bn_free(&bn_a);
  mpz_clears(mpz_a, mpz_d, mpz_rem, NULL);
}

void test_divmod_u64(const char* name, const char* a_str, uint64_t d)
{
  bignum bn_a, bn_q, bn_alias;
  mpz_t mpz_a, mpz_d, mpz_q, mpz_rem;

  bn_init_multi(&bn_a, &bn_q, &bn_alias, NULL);
  bn_init_val(&bn_a, a_str);
  bn_init_val(&bn_alias, a_str);

  mpz_inits(mpz_a, mpz_d, mpz_q, mpz_rem, NULL);
  mpz_set_str(mpz_a, a_str, 10);
  set_mpz_u64(mpz_d, d);

  // GMP truncated division matches bn_divmod_u64
  mpz_tdiv_qr(mpz_q, mpz_rem, mpz_a, mpz_d);
  mpz_abs(mpz_rem, mpz_rem);  // bn_divmod_u64 returns magnitude modulo

  uint64_t gmp_rem = mpz_to_u64(mpz_rem);
  char* gmp_q_str = mpz_get_str(NULL, 10, mpz_q);

  // Test 1: Non-Aliased
  uint64_t custom_rem = bn_divmod_u64(&bn_q, &bn_a, d);
  char* custom_q_str = bn_to_str(&bn_q);

  if (custom_rem != gmp_rem || strcmp(gmp_q_str, custom_q_str) != 0) {
    fprintf(stderr, "\n[FATAL ERROR] divmod_u64 non-aliased mismatch (%s)!\n",
            name);
    fprintf(stderr, "A: %s, D: %llu\n", a_str, (unsigned long long)d);
    fprintf(stderr, "Expected Q: %s, Rem: %llu\n", gmp_q_str,
            (unsigned long long)gmp_rem);
    fprintf(stderr, "Got Q:      %s, Rem: %llu\n", custom_q_str,
            (unsigned long long)custom_rem);
    exit(EXIT_FAILURE);
  }

  // Test 2: Aliased (q == a)
  uint64_t alias_rem = bn_divmod_u64(&bn_alias, &bn_alias, d);
  char* alias_q_str = bn_to_str(&bn_alias);

  if (alias_rem != gmp_rem || strcmp(gmp_q_str, alias_q_str) != 0) {
    fprintf(stderr, "\n[FATAL ERROR] divmod_u64 ALIASED mismatch (%s)!\n",
            name);
    fprintf(stderr, "A: %s, D: %llu\n", a_str, (unsigned long long)d);
    fprintf(stderr, "Expected Q: %s, Rem: %llu\n", gmp_q_str,
            (unsigned long long)gmp_rem);
    fprintf(stderr, "Got Q:      %s, Rem: %llu\n", alias_q_str,
            (unsigned long long)alias_rem);
    exit(EXIT_FAILURE);
  }

  free(gmp_q_str);
  free(custom_q_str);
  free(alias_q_str);
  bn_free_multi(&bn_a, &bn_q, &bn_alias, NULL);
  mpz_clears(mpz_a, mpz_d, mpz_q, mpz_rem, NULL);
}

void test_mod_inverse(const char* name, const char* a_str, const char* m_str)
{
  bignum bn_a, bn_m, bn_res;
  mpz_t mpz_a, mpz_m, mpz_res;

  bn_init(&bn_a);
  bn_init(&bn_m);
  bn_init(&bn_res);
  bn_init_val(&bn_a, a_str);
  bn_init_val(&bn_m, m_str);

  mpz_inits(mpz_a, mpz_m, mpz_res, NULL);
  mpz_set_str(mpz_a, a_str, 10);
  mpz_set_str(mpz_m, m_str, 10);

  bool custom_has_inv = bn_mod_inverse(&bn_res, &bn_a, &bn_m);
  int gmp_has_inv = mpz_invert(mpz_res, mpz_a, mpz_m);

  if ((custom_has_inv ? 1 : 0) != (gmp_has_inv != 0)) {
    fprintf(stderr, "\n[FATAL ERROR] mod_inverse mismatch (%s)!\n", name);
    fprintf(stderr, "A: %s, M: %s\nExpected inv_exists: %d, Got: %d\n", a_str,
            m_str, gmp_has_inv != 0, custom_has_inv);
    exit(EXIT_FAILURE);
  }

  if (custom_has_inv) {
    char* gmp_res_str = mpz_get_str(NULL, 10, mpz_res);
    char* custom_res_str = bn_to_str(&bn_res);
    if (strcmp(gmp_res_str, custom_res_str) != 0) {
      fprintf(stderr, "\n[FATAL ERROR] mod_inverse value mismatch (%s)!\n",
              name);
      fprintf(stderr, "Expected: %s, Got: %s\n", gmp_res_str, custom_res_str);
      exit(EXIT_FAILURE);
    }
    free(gmp_res_str);
    free(custom_res_str);
  }

  bn_free(&bn_a);
  bn_free(&bn_m);
  bn_free(&bn_res);
  mpz_clears(mpz_a, mpz_m, mpz_res, NULL);
}

void run_edge_cases(gmp_randstate_t state)
{
  printf("\n--- Running BigMod Edge Cases ---\n");

  test_mod_u64("Zero mod", "0", 12345);
  test_mod_u64("Neg A", "-999999999999999", 7);
  test_mod_u64("Small Mod", "10", 3);
  test_mod_u64("Max U64 Divisor", "340282366920938463463374607431768211455",
               18446744073709551615ULL);

  test_divmod_u64("Zero div", "0", 12345);
  test_divmod_u64("Neg A", "-999999999999999", 7);
  test_divmod_u64("Divide by 1", "9876543210123456789", 1);
  test_divmod_u64("Max U64 Divisor", "340282366920938463463374607431768211455",
                  18446744073709551615ULL);

  test_mod_inverse("No inv (even)", "2", "4");
  test_mod_inverse("No inv (0)", "0", "13");
  test_mod_inverse("Basic inv", "3", "11");
  test_mod_inverse("Large prime", "123456789", "18446744073709551557");

  printf("Running 1,000 randomized verification tests...\n");
  for (int i = 0; i < 1000; i++) {
    int bits = 1 + (rand() % 4096);

    bignum bn_a, bn_m;
    mpz_t mpz_a, mpz_m, mpz_rand_u64;

    bn_init(&bn_a);
    bn_init(&bn_m);
    mpz_inits(mpz_a, mpz_m, mpz_rand_u64, NULL);

    // Rand A and U64 D
    generate_random_bignum(&bn_a, mpz_a, bits, true, state);
    mpz_urandomb(mpz_rand_u64, state, 63);
    mpz_add_ui(mpz_rand_u64, mpz_rand_u64, 1);  // Avoid div by 0
    uint64_t d = mpz_to_u64(mpz_rand_u64);

    char* a_str = bn_to_str(&bn_a);
    test_mod_u64("Random", a_str, d);
    test_divmod_u64("Random", a_str, d);
    free(a_str);

    // Inverse
    generate_random_bignum(&bn_a, mpz_a, bits, false, state);
    generate_random_bignum(&bn_m, mpz_m, bits, false, state);
    mpz_add_ui(mpz_a, mpz_a, 1);
    mpz_add_ui(mpz_m, mpz_m, 2);  // m >= 2

    a_str = mpz_get_str(NULL, 10, mpz_a);
    char* m_str = mpz_get_str(NULL, 10, mpz_m);
    test_mod_inverse("Random", a_str, m_str);

    free(a_str);
    free(m_str);
    bn_free(&bn_a);
    bn_free(&bn_m);
    mpz_clears(mpz_a, mpz_m, mpz_rand_u64, NULL);
  }
  printf("All edge cases and random tests passed!\n");
}

// =============================================================================
// BENCHMARKS
// =============================================================================

void benchmark_divmod_u64(int bits, double target_sec, gmp_randstate_t state)
{
  bignum bn_a, bn_q;
  mpz_t mpz_a, mpz_d, mpz_q, mpz_rem;

  bn_init(&bn_a);
  bn_init(&bn_q);
  mpz_inits(mpz_a, mpz_d, mpz_q, mpz_rem, NULL);

  generate_random_bignum(&bn_a, mpz_a, bits, false, state);

  mpz_t rand_d;
  mpz_init(rand_d);
  mpz_urandomb(rand_d, state, 63);
  mpz_add_ui(rand_d, rand_d, 1);
  uint64_t d = mpz_to_u64(rand_d);
  set_mpz_u64(mpz_d, d);
  mpz_clear(rand_d);

  char* a_str = mpz_get_str(NULL, 10, mpz_a);
  bn_init_val(&bn_a, a_str);
  free(a_str);

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
    mpz_tdiv_qr(mpz_q, mpz_rem, mpz_a, mpz_d);
    ops_gmp++;
    clock_gettime(CLOCK_MONOTONIC, &end);
    total_gmp = get_elapsed_time(start, end);
  } while (total_gmp < target_sec);

  char size_info[32];
  snprintf(size_info, sizeof(size_info), "%d bits", bits);
  print_table_row("bn_divmod_u64", size_info, total_custom / ops_custom,
                  total_gmp / ops_gmp);

  bn_free(&bn_a);
  bn_free(&bn_q);
  mpz_clears(mpz_a, mpz_d, mpz_q, mpz_rem, NULL);
}

void benchmark_mod_inverse(int bits, double target_sec, gmp_randstate_t state)
{
  bignum bn_a, bn_m, bn_res;
  mpz_t mpz_a, mpz_m, mpz_res;

  bn_init(&bn_a);
  bn_init(&bn_m);
  bn_init(&bn_res);
  mpz_inits(mpz_a, mpz_m, mpz_res, NULL);

  generate_random_bignum(&bn_a, mpz_a, bits, false, state);
  generate_random_bignum(&bn_m, mpz_m, bits, false, state);
  mpz_nextprime(mpz_m, mpz_m);  // Prime ensures inverse exists

  char* a_str = mpz_get_str(NULL, 10, mpz_a);
  char* m_str = mpz_get_str(NULL, 10, mpz_m);
  bn_init_val(&bn_a, a_str);
  bn_init_val(&bn_m, m_str);
  free(a_str);
  free(m_str);

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

  char size_info[32];
  snprintf(size_info, sizeof(size_info), "%d bits", bits);
  print_table_row("bn_mod_inverse", size_info, total_custom / ops_custom,
                  total_gmp / ops_gmp);

  bn_free(&bn_a);
  bn_free(&bn_m);
  bn_free(&bn_res);
  mpz_clears(mpz_a, mpz_m, mpz_res, NULL);
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

  run_edge_cases(state);

  double time_per_test = 0.5;
  print_table_header();

  benchmark_divmod_u64(512, time_per_test, state);
  benchmark_divmod_u64(1024, time_per_test, state);
  benchmark_divmod_u64(2048, time_per_test, state);
  benchmark_divmod_u64(4096, time_per_test, state);

  // Separator for visual clarity
  printf(
      "+--------------------+--------------+---------------+---------------+---"
      "---------+\n");

  benchmark_mod_inverse(512, time_per_test, state);
  benchmark_mod_inverse(1024, time_per_test, state);
  benchmark_mod_inverse(2048, time_per_test, state);
  benchmark_mod_inverse(4096, time_per_test, state);

  print_table_footer();

  gmp_randclear(state);
  return 0;
}