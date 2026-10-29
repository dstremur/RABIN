/*
 * test_util.h
 *
 * Minimal assertion harness shared by the RABIN test suites.
 *
 * Usage pattern:
 *
 *   TEST(my_op) { ... TEST_ASSERT(cond, "msg"); ... }
 *   int main(void) { RUN_TEST(my_op); return TEST_REPORT(); }
 *
 * TEST_REPORT() returns the number of failed tests (plus a flag for stray
 * assertions outside RUN_TEST), so the process exit code is 0 only when
 * every check passed.
 *
 * Copyright (C) 2026 Diego Strebel
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef TEST_UTIL_H
#define TEST_UTIL_H

#include <stdio.h>

static unsigned test_total = 0;
static unsigned test_failed = 0;
static unsigned assert_failed = 0;

/* Records a failure but keeps running so one test shows all its problems. */
#define TEST_ASSERT(cond, msg)                                     \
  do                                                               \
  {                                                                \
    if (!(cond))                                                   \
    {                                                              \
      fprintf(stderr, "[FAIL] %s:%d: %s\n", __FILE__, __LINE__,    \
              (msg));                                              \
      assert_failed++;                                             \
    }                                                              \
  } while (0)

/* Defines a test case body: TEST(name) { ... } */
#define TEST(name) static void name(void)

/* Runs one test case and prints its verdict. */
#define RUN_TEST(fn)                                               \
  do                                                               \
  {                                                                \
    unsigned _before = assert_failed;                              \
    test_total++;                                                  \
    fn();                                                          \
    if (assert_failed > _before)                                   \
    {                                                              \
      test_failed++;                                               \
      printf("[FAIL] %s\n", #fn);                                  \
    }                                                              \
    else                                                           \
    {                                                              \
      printf("[PASS] %s\n", #fn);                                  \
    }                                                              \
    fflush(stdout);                                                \
  } while (0)

/* Prints the summary and yields the process exit code. */
#define TEST_REPORT()                                              \
  (printf("%u/%u tests passed, %u assertion failures\n",           \
          test_total - test_failed, test_total, assert_failed),    \
   fflush(stdout),                                                 \
   (test_failed != 0 || assert_failed != 0) ? 1 : 0)

/* Optional-dependency skip: prints the reason and exits 0. */
#define TEST_SKIP(dep)                                             \
  do                                                               \
  {                                                                \
    printf("SKIP: %s not available\n", (dep));                     \
    return 0;                                                      \
  } while (0)

#endif /* TEST_UTIL_H */
