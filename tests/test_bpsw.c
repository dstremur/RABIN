/*
 * test_bpsw.c
 *
 * Regression tests for rz_bpsw(), rz_stronglucas(), rz_rabin*() and the
 * prime generators, plus the original latency benchmark.
 *
 * Every composite in the pseudoprime battery is a strong Miller-Rabin
 * pseudoprime to base 2 (verified by an independent u128 reference and
 * `openssl prime`), so BPSW must reject it through the Lucas stage with
 * the given Selfridge parameters.
 *
 * Copyright (C) 2026 Diego Strebel
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "../include/rabin.h"
#include "test_util.h"

/* Oracle-verified vectors: n, Selfridge Q = (1 - D)/4 with (D|n) = -1. */
typedef struct {
  u64 n;
  long q;
  const char* factors;
} pseudoprime_t;

static const pseudoprime_t PSEUDOPRIMES[] = {
    {2047ULL, -1, "23 * 89"},
    {1373653ULL, -1, "829 * 1657"},
    {3215031751ULL, 3, "151 * 751 * 28351"},
    {2152302898747ULL, -1, "6763 * 10627 * 29947"},
    {3474749660383ULL, -1, "1303 * 16927 * 157543"},
    {3825123056546413051ULL, 2, "149491 * 747451 * 34233211"},
};
#define NUM_PSEUDOPRIMES (sizeof(PSEUDOPRIMES) / sizeof(PSEUDOPRIMES[0]))

#define SIEVE_MAX 200000
static unsigned char sieve_comp[SIEVE_MAX + 1];

static void build_sieve(void)
{
  sieve_comp[0] = sieve_comp[1] = 1;
  for (u64 p = 2; p * p <= SIEVE_MAX; p++) {
    if (!sieve_comp[p]) {
      for (u64 m = p * p; m <= SIEVE_MAX; m += p) sieve_comp[m] = 1;
    }
  }
}

/* Deterministic primality for u64 values > SIEVE_MAX: trial division by
 * odds up to sqrt(v); only used for generated primes up to 2^40. */
static int is_prime_det(u64 v)
{
  if (v <= SIEVE_MAX) return !sieve_comp[v];
  if (v < 2 || !(v & 1)) return v == 2;
  for (u64 d = 3; d * d <= v; d += 2)
    if (v % d == 0) return 0;
  return 1;
}

TEST(domain_boundaries)
{
  rz_t n;
  rz_init(&n);

  /* primes are positive integers > 1 */
  rz_set_i64(&n, -2);
  TEST_ASSERT(rz_bpsw(&n) == false, "bpsw(-2) must be false");
  rz_set_i64(&n, -3);
  TEST_ASSERT(rz_bpsw(&n) == false, "bpsw(-3) must be false");
  rz_set_i64(&n, -5);
  TEST_ASSERT(rz_bpsw(&n) == false, "bpsw(-5) must be false");
  rz_set_u64(&n, 0);
  TEST_ASSERT(rz_bpsw(&n) == false, "bpsw(0) must be false");
  rz_set_u64(&n, 1);
  TEST_ASSERT(rz_bpsw(&n) == false, "bpsw(1) must be false");

  rz_set_u64(&n, 2);
  TEST_ASSERT(rz_bpsw(&n) == true, "bpsw(2) must be true");
  rz_set_u64(&n, 3);
  TEST_ASSERT(rz_bpsw(&n) == true, "bpsw(3) must be true");
  rz_set_u64(&n, 5);
  TEST_ASSERT(rz_bpsw(&n) == true, "bpsw(5) must be true");

  rz_set_u64(&n, 4);
  TEST_ASSERT(rz_bpsw(&n) == false, "bpsw(4) must be false");
  rz_set_u64(&n, 100);
  TEST_ASSERT(rz_bpsw(&n) == false, "bpsw(100) must be false");

  rz_clear(&n);
}

TEST(exhaustive_vs_sieve)
{
  unsigned prime_count = 0;
  for (u64 v = 0; v <= SIEVE_MAX; v++) {
    rz_t n;
    rz_init(&n);
    rz_set_u64(&n, v);
    int want = !sieve_comp[v];
    if (want) prime_count++;
    if (rz_bpsw(&n) != (bool)want) {
      TEST_ASSERT(false, "bpsw disagrees with sieve in [0, 200000]");
      fprintf(stderr, "  first disagreement at n = %llu\n",
              (unsigned long long)v);
      rz_clear(&n);
      return;
    }
    rz_clear(&n);
  }
  TEST_ASSERT(prime_count == 17984,
              "sieve must contain 17984 primes <= 200000");
}

TEST(strong_pseudoprime_battery)
{
  rz_t n, P, Q, two;
  rz_init_multi(&n, &P, &Q, &two, NULL);
  rz_set_u64(&P, 1);
  rz_set_u64(&two, 2);

  for (unsigned i = 0; i < NUM_PSEUDOPRIMES; i++) {
    rz_set_u64(&n, PSEUDOPRIMES[i].n);
    rz_set_i64(&Q, PSEUDOPRIMES[i].q);

    /* all are strong MR pseudoprimes to base 2: MR must NOT reject */
    TEST_ASSERT(rz_rabin_mont(&n, &two) == true,
                "base-2 MR must accept (strong pseudoprime)");

    /* the Lucas stage must reject with the oracle Selfridge parameters */
    TEST_ASSERT(rz_stronglucas(&n, &P, &Q) == false,
                "strong Lucas must reject pseudoprime");

    /* and the full BPSW must reject */
    TEST_ASSERT(rz_bpsw(&n) == false, "BPSW must reject pseudoprime");
  }

  rz_clear_multi(&n, &P, &Q, &two, NULL);
}

TEST(perfect_squares)
{
  /* regression: squares used to be checked only at round 5 of the Selfridge
   * search, where (D|n) never reaches -1 and the loop spins forever */
  u64 sq[] = {9ULL,
              25ULL,
              121ULL,
              2047ULL * 2047ULL,
              10007ULL * 10007ULL,
              1000003ULL * 1000003ULL};
  for (unsigned i = 0; i < 6; i++) {
    rz_t n;
    rz_init(&n);
    rz_set_u64(&n, sq[i]);
    TEST_ASSERT(rz_bpsw(&n) == false, "perfect square must be rejected");
    rz_clear(&n);
  }
}

TEST(known_large_primes)
{
  u64 pr[] = {104729ULL, 2147483647ULL, /* 2^31 - 1 */
              2305843009213693951ULL};  /* 2^61 - 1 */
  for (unsigned i = 0; i < 3; i++) {
    rz_t n;
    rz_init(&n);
    rz_set_u64(&n, pr[i]);
    TEST_ASSERT(rz_bpsw(&n) == true, "known Mersenne prime must pass");
    rz_clear(&n);
  }
}

TEST(rabin_witnesses)
{
  rz_t n, a;
  rz_init_multi(&n, &a, NULL);
  rz_set_u64(&a, 2);

  /* domain guards */
  rz_set_i64(&n, -7);
  TEST_ASSERT(rz_rabin(&n, &a) == false, "rz_rabin(-7) must be false");
  rz_set_u64(&n, 1);
  TEST_ASSERT(rz_rabin(&n, &a) == false, "rz_rabin(1) must be false");
  rz_set_u64(&n, 2);
  TEST_ASSERT(rz_rabin(&n, &a) == true, "rz_rabin(2) must be true");

  /* 341 = 11 * 31 is a Fermat but NOT a strong base-2 pseudoprime */
  rz_set_u64(&n, 341);
  rz_set_u64(&a, 2);
  TEST_ASSERT(rz_rabin(&n, &a) == false, "341 must fail strong MR base 2");
  TEST_ASSERT(rz_rabin_mont(&n, &a) == false, "341 must fail MR(mont)");

  /* 2047 fails MR to base 3 (oracle-verified) */
  rz_set_u64(&n, 2047);
  rz_set_u64(&a, 3);
  TEST_ASSERT(rz_rabin(&n, &a) == false, "2047 must fail MR base 3");
  TEST_ASSERT(rz_rabin_mont(&n, &a) == false, "2047 must fail MR(mont) b3");

  /* 1373653 is a strong base-2 pseudoprime: both variants accept */
  rz_set_u64(&n, 1373653);
  rz_set_u64(&a, 2);
  TEST_ASSERT(rz_rabin(&n, &a) == true, "1373653 passes MR base 2");
  TEST_ASSERT(rz_rabin_mont(&n, &a) == true, "1373653 passes MR(mont) b2");

  rz_clear_multi(&n, &a, NULL);
}

TEST(rabin_mont_ctx_mismatch)
{
  rz_t n, a;
  rz_init_multi(&n, &a, NULL);
  rz_mont_ctx ctx = {0};

  rz_set_u64(&n, 10007);
  rz_set_u64(&a, 2);
  TEST_ASSERT(rz_mont_ctx_init(&ctx, &n) == RABIN_SUCCESS,
              "ctx init for 10007");

  /* modulus under test must match the context modulus */
  rz_set_u64(&n, 10009);
  TEST_ASSERT(rz_rabin_mont_ctx(&n, &a, &ctx) == false,
              "ctx modulus mismatch must be false");

  rz_set_u64(&n, 10007);
  TEST_ASSERT(rz_rabin_mont_ctx(&n, &a, &ctx) == true,
              "matching ctx must accept prime 10007");

  rz_mont_ctx_clear(&ctx);
  rz_clear_multi(&n, &a, NULL);
}

TEST(mont_negative_input)
{
  /* regression: rz_mont_in silently dropped the sign of negative inputs,
   * which made composites pass the Lucas stage of BPSW */
  rz_t n, a, in, out;
  rz_init_multi(&n, &a, &in, &out, NULL);
  rz_mont_ctx ctx = {0};

  rz_set_u64(&n, 10007);
  TEST_ASSERT(rz_mont_ctx_init(&ctx, &n) == RABIN_SUCCESS, "ctx init");

  /* -1 mod 10007 = 10006 */
  rz_set_i64(&a, -1);
  TEST_ASSERT(rz_mont_in(&in, &a, &ctx) == RABIN_SUCCESS, "mont_in(-1)");
  TEST_ASSERT(rz_mont_out(&out, &in, &ctx) == RABIN_SUCCESS, "mont_out");
  TEST_ASSERT(rz_is_eq_i64(&out, 10006), "mont_in(-1) must roundtrip to n-1");

  /* -12345 mod 10007 = 10007 - (12345 % 10007) = 10007 - 2338 = 7669 */
  rz_set_i64(&a, -12345);
  TEST_ASSERT(rz_mont_in(&in, &a, &ctx) == RABIN_SUCCESS, "mont_in(-12345)");
  TEST_ASSERT(rz_mont_out(&out, &in, &ctx) == RABIN_SUCCESS, "mont_out");
  TEST_ASSERT(rz_is_eq_i64(&out, 7669), "mont_in(-12345) -> 7669");

  rz_mont_ctx_clear(&ctx);
  rz_clear_multi(&n, &a, &in, &out, NULL);
}

TEST(gen_prime_sizes)
{
  /* regression: the old low-limb shortcut accepted any candidate whose
   * least-significant limb equaled a small prime, and tiny bit sizes could
   * retry forever. Generated primes are checked deterministically. */
  int sizes[] = {2, 3, 5, 7, 11, 13, 16, 20, 24, 32, 40};
  for (unsigned s = 0; s < 11; s++) {
    for (int it = 0; it < 10; it++) {
      rz_t p;
      rz_init(&p);
      TEST_ASSERT(rz_gen_prime(&p, sizes[s]) == RABIN_SUCCESS,
                  "gen_prime must succeed");
      TEST_ASSERT(rz_bit_length(&p) == sizes[s], "bit length must match");
      u64 v = p.limbs[0]; /* sizes[s] <= 40 fits in one 64-bit limb */
      TEST_ASSERT(is_prime_det(v) == 1, "generated number must be prime");
      rz_clear(&p);
    }
  }
}

TEST(gen_safe_prime_small)
{
  int sizes[] = {12, 16};
  for (unsigned s = 0; s < 2; s++) {
    rz_t p;
    rz_init(&p);
    TEST_ASSERT(rz_gen_safe_prime(&p, sizes[s]) == RABIN_SUCCESS,
                "gen_safe_prime must succeed");
    u64 v = p.limbs[0];
    TEST_ASSERT(is_prime_det(v) == 1, "safe prime p must be prime");
    TEST_ASSERT(is_prime_det((v - 1) / 2) == 1, "(p-1)/2 must be prime");
    rz_clear(&p);
  }
}

int main(void)
{
  build_sieve();

  RUN_TEST(domain_boundaries);
  RUN_TEST(exhaustive_vs_sieve);
  RUN_TEST(strong_pseudoprime_battery);
  RUN_TEST(perfect_squares);
  RUN_TEST(known_large_primes);
  RUN_TEST(rabin_witnesses);
  RUN_TEST(rabin_mont_ctx_mismatch);
  RUN_TEST(mont_negative_input);
  RUN_TEST(gen_prime_sizes);
  RUN_TEST(gen_safe_prime_small);

  int rc = TEST_REPORT();

  /* benchmark (kept from the original suite) */
  int iterations = 100;
  int bits = 2048;
  rz_t n;
  rz_init(&n);
  if (rz_gen_prime(&n, bits) == RABIN_SUCCESS) {
    printf("Benchmarking BPSW: %d iterations at %d bits...\n", iterations,
           bits);

    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC, &start);
    for (int i = 0; i < iterations; i++) rz_bpsw(&n);
    clock_gettime(CLOCK_MONOTONIC, &end);

    double total_time =
        (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
    printf("Total Time:  %.4f s\n", total_time);
    printf("Throughput:  %.2f tests/sec\n", iterations / total_time);
  }
  rz_clear(&n);

  return rc;
}
