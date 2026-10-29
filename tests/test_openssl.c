/*
 * test_openssl.c
 *
 * Cross-check rz_gen_prime() output against the `openssl prime` CLI. Skips
 * gracefully when the CLI is not installed.
 *
 * Copyright (C) 2026 Diego Strebel
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../include/rabin.h"
#include "test_util.h"

static int verify_with_openssl(const char* num_str)
{
  char command[5000]; /* Increased for 2048-bit strings */
  char result_buf[1024];
  int is_prime = 0;

  /* Redirect 2>&1 to catch warnings and ensure they don't block the pipe */
  snprintf(command, sizeof(command), "openssl prime %s 2>&1", num_str);

  FILE* fp = popen(command, "r");
  if (fp == NULL) return -1;

  while (fgets(result_buf, sizeof(result_buf), fp) != NULL) {
    if (strstr(result_buf, "is prime") != NULL) {
      is_prime = 1;
      break;
    }
  }

  pclose(fp);
  return is_prime;
}

static int openssl_cli_available(void)
{
  return system("command -v openssl >/dev/null 2>&1") == 0;
}

static void run_prime_test(int num_tests, int bits)
{
  printf("Starting %d tests for %d-bit primes...\n", num_tests, bits);

  for (int i = 0; i < num_tests; i++) {
    rz_t p;
    rz_init(&p);

    TEST_ASSERT(rz_gen_prime(&p, bits) == RABIN_SUCCESS,
                "rz_gen_prime succeeds");

    char* p_str = rz_to_string(&p);
    TEST_ASSERT(p_str != NULL, "rz_to_string succeeds");

    TEST_ASSERT(rz_bpsw(&p), "generated prime passes BPSW");

    int verified = verify_with_openssl(p_str);
    if (verified) {
      printf("[PASS] Test %d: %d-bit candidate verified prime by OpenSSL\n",
             i + 1, bits);
    } else {
      printf("[FAIL] Test %d: %s is NOT prime according to OpenSSL\n", i + 1,
             p_str);
    }
    TEST_ASSERT(verified == 1, "openssl prime CLI confirms primality");

    free(p_str);
    rz_clear(&p);
  }
}

TEST(openssl_selfcheck_vector)
{
  TEST_ASSERT(verify_with_openssl("104729") == 1,
              "openssl confirms known prime 104729");
  TEST_ASSERT(verify_with_openssl("104730") == 0,
              "openssl rejects known composite 104730");
}

TEST(openssl_generated_primes) { run_prime_test(20, 2048); }

int main(void)
{
  if (!openssl_cli_available()) {
    TEST_SKIP("openssl CLI not available");
    return 0;
  }

  rz_init_constants();

  RUN_TEST(openssl_selfcheck_vector);
  RUN_TEST(openssl_generated_primes);

  rz_clear_constants();
  return TEST_REPORT();
}
