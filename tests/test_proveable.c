/*
 * test_proveable.c
 *
 * rz_gen_provable_arithmetic(): correctness gate (returned value must be a
 * prime, certificate must be populated) plus generation timing benchmark.
 *
 * Copyright (C) 2026 Diego Strebel
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <stdio.h>
#include <time.h>

#include "../include/rabin.h"
#include "../include/rzcert.h"
#include "test_util.h"

TEST(provable_primes_are_prime)
{
  int num_to_generate = 3;
  int bits = 1024;

  rz_t p;
  rz_init(&p);

  struct timespec start, end;
  double total_time = 0;

  for (int i = 0; i < num_to_generate; i++) {
    rzcert_t* cert = NULL;

    clock_gettime(CLOCK_MONOTONIC, &start);
    rabin_err_t err = rz_gen_provable_arithmetic(&p, bits, &cert);
    clock_gettime(CLOCK_MONOTONIC, &end);

    TEST_ASSERT(err == RABIN_SUCCESS, "rz_gen_provable_arithmetic succeeds");
    TEST_ASSERT(cert != NULL, "certificate is emitted");
    TEST_ASSERT(rz_bpsw(&p), "provable prime passes BPSW");

    double elapsed =
        (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
    total_time += elapsed;

    printf("Prime #%d: Found in %.4f seconds\n", i + 1, elapsed);
    print_pocklington_cert(cert);

    // Requires Phase 3 (rzcert_verify implementation):
    // TEST_ASSERT(rzcert_verify(cert, 256), "certificate verifies");

    rzcert_clear(cert);
  }

  printf("--------------------------------------------------\n");
  printf("Total Time:   %.4f seconds\n", total_time);
  printf("Average Time: %.4f seconds per prime\n",
         total_time / num_to_generate);

  rz_clear(&p);
}

int main(void)
{
  rz_init_constants();

  RUN_TEST(provable_primes_are_prime);

  rz_clear_constants();
  return TEST_REPORT();
}
