/*
 * test_primes_parallel.c
 *
 * OpenMP stress test for rz_gen_prime(): correctness gate (BPSW-verify every
 * prime generated across threads) plus parallel throughput benchmark.
 *
 * Copyright (C) 2026 Diego Strebel
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <omp.h>
#include <stdio.h>
#include <time.h>

#include "../include/rabin.h"
#include "test_util.h"

static int failures = 0;

TEST(primes_parallel_generated_are_prime)
{
  int num_to_generate = 200;
  int bits = 2048;

  printf("Generating %d primes at %d-bits (Parallel)...\n", num_to_generate,
         bits);
  printf("--------------------------------------------------\n");

  struct timespec global_start, global_end;
  clock_gettime(CLOCK_MONOTONIC, &global_start);

#pragma omp parallel for schedule(dynamic, 1)
  for (int i = 0; i < num_to_generate; i++) {
    rz_t p;
    rz_init(&p);

    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC, &start);
    rabin_err_t err = rz_gen_prime(&p, bits);
    clock_gettime(CLOCK_MONOTONIC, &end);

    bool prime = (err == RABIN_SUCCESS) && rz_bpsw(&p);

    if (!prime) {
#pragma omp atomic
      failures++;
    }

#pragma omp critical
    {
      if (prime) {
        double elapsed =
            (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
        printf("Thread %d | Prime #%d: Found in %.4f seconds\n",
               omp_get_thread_num(), i + 1, elapsed);
      } else {
        printf("Failed to generate/verify prime #%d\n", i + 1);
      }
    }

    rz_clear(&p);
  }

  clock_gettime(CLOCK_MONOTONIC, &global_end);
  double total_wall_time = (global_end.tv_sec - global_start.tv_sec) +
                           (global_end.tv_nsec - global_start.tv_nsec) / 1e9;

  printf("--------------------------------------------------\n");
  printf("Total Wall-Clock Time: %.4f seconds\n", total_wall_time);
  printf("Efficiency: %.4f primes/second\n", num_to_generate / total_wall_time);

  TEST_ASSERT(failures == 0, "every parallel-generated prime is BPSW-prime");
}

int main(void)
{
  rz_init_constants();

  RUN_TEST(primes_parallel_generated_are_prime);

  rz_clear_constants();
  return TEST_REPORT();
}
