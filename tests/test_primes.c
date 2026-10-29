/*
 * test_primes.c
 *
 * rz_gen_prime(): correctness gate (every generated prime is BPSW-verified,
 * every bit length must always succeed) plus generation timing benchmark.
 *
 * Copyright (C) 2026 Diego Strebel
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <stdio.h>
#include <time.h>

#include "../include/rabin.h"
#include "test_util.h"

#define ITERATIONS 1000

TEST(primes_known_answers)
{
  static const char* primes[] = {"2", "3", "104729", "2147483647",
                                 "170141183460469231731687303715884105727"};
  static const char* composites[] = {"0",
                                     "1",
                                     "4",
                                     "2047",
                                     "104730",
                                     "3215031751",
                                     "170141183460469231731687303715884105728"};

  rz_t n;
  rz_init(&n);

  for (size_t i = 0; i < sizeof(primes) / sizeof(primes[0]); i++) {
    TEST_ASSERT(rz_init_val(&n, primes[i]) == RABIN_SUCCESS, "init_val");
    TEST_ASSERT(rz_bpsw(&n), "known prime must pass BPSW");
  }
  for (size_t i = 0; i < sizeof(composites) / sizeof(composites[0]); i++) {
    TEST_ASSERT(rz_init_val(&n, composites[i]) == RABIN_SUCCESS, "init_val");
    TEST_ASSERT(!rz_bpsw(&n), "known composite must fail BPSW");
  }

  rz_clear(&n);
}

TEST(primes_generated_are_prime)
{
  int bits = 2048;
  rz_t p;
  rz_init(&p);

  struct timespec start, end;
  double total_time = 0;

  for (int i = 0; i < 3; i++) {
    clock_gettime(CLOCK_MONOTONIC, &start);
    rabin_err_t err = rz_gen_prime(&p, bits);
    clock_gettime(CLOCK_MONOTONIC, &end);

    TEST_ASSERT(err == RABIN_SUCCESS, "rz_gen_prime succeeds");
    TEST_ASSERT(rz_bpsw(&p), "generated 2048-bit candidate is prime");

    double elapsed =
        (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
    total_time += elapsed;
    printf("Prime #%d: Found in %.4f seconds\n", i + 1, elapsed);
  }

  printf("Total Time:   %.4f seconds\n", total_time);
  rz_clear(&p);
}

TEST(primes_small_bit_lengths)
{
  rz_t p;
  rz_init(&p);

  struct timespec start, end;

  for (int i = 5; i <= 100; i++) {
    int success_count = 0;

    bool composite_seen = false;
    clock_gettime(CLOCK_MONOTONIC, &start);
    for (int j = 0; j < ITERATIONS; j++) {
      if (rz_gen_prime(&p, i) == RABIN_SUCCESS) {
        success_count++;
        if (!rz_bpsw(&p)) {
          composite_seen = true;
          break;
        }
      }
    }
    clock_gettime(CLOCK_MONOTONIC, &end);

    TEST_ASSERT(success_count == ITERATIONS,
                "rz_gen_prime must succeed for every small bit length");
    TEST_ASSERT(!composite_seen, "every generated small prime passes BPSW");

    if (success_count == ITERATIONS) {
      double total =
          (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
      printf("%-10d | %-15.6f | %-15.9f\n", i, total, total / ITERATIONS);
    }
  }

  rz_clear(&p);
}

TEST(primes_invalid_args)
{
  rz_t p;
  rz_init(&p);

  TEST_ASSERT(rz_gen_prime(&p, 1) == RABIN_ERR_INVALID_ARG,
              "rz_gen_prime(bits < 2) is rejected");

  rz_clear(&p);
}

int main(void)
{
  rz_init_constants();

  RUN_TEST(primes_known_answers);
  RUN_TEST(primes_generated_are_prime);
  RUN_TEST(primes_small_bit_lengths);
  RUN_TEST(primes_invalid_args);

  rz_clear_constants();
  return TEST_REPORT();
}
