/*
 * test_rat.c
 *
 * Assert-based checks for rq_t arithmetic invariants (normalized results,
 * positive denominators) and rmatq_t add/sub with known answers.
 *
 * Copyright (C) 2026 Diego Strebel
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include "../include/rmat.h"
#include "../include/rmatq.h"
#include "../include/rq.h"
#include "test_util.h"

static void rq_set_pair(rq_t* r, i64 num, i64 den)
{
  rz_set_i64(&r->num, num);
  rz_set_i64(&r->den, den);
}

static void rq_set_num_den_str(rq_t* r, i64 num, const char* den)
{
  rz_set_i64(&r->num, num);
  rz_init_val(&r->den, den);
}

/* Compare two rationals independent of their current reduction level. */
static bool rq_equal(const rq_t* a, const rq_t* b)
{
  rq_t x, y;
  rq_init_multi(&x, &y, NULL);
  rq_copy(&x, a);
  rq_copy(&y, b);
  rq_normalize(&x);
  rq_normalize(&y);
  bool eq = (rq_cmp(&x, &y) == 0);
  rq_clear_multi(&x, &y, NULL);
  return eq;
}

/* A normalized rational must have a positive denominator. */
static bool rq_den_positive(const rq_t* r)
{
  return !r->den.is_neg && !rz_is_zero(&r->den);
}

TEST(rat_arithmetic)
{
  rq_t r, q, p, expect;
  rq_init_multi(&r, &q, &p, &expect, NULL);

  rq_set_pair(&r, 3, 8905493208);
  rq_set_pair(&q, 3, 23534563456);

  rq_print(&r);
  TEST_ASSERT(rq_mul(&p, &r, &q) == RABIN_SUCCESS, "rq_mul");
  rq_normalize(&p);
  TEST_ASSERT(rq_den_positive(&p), "product has positive denominator");
  rq_set_num_den_str(&expect, 9, "209586895010653006848");
  TEST_ASSERT(rq_equal(&p, &expect), "3/8905493208 * 3/23534563456 correct");

  rq_print(&p);
  rq_print(&r);

  rq_set_pair(&r, -1, 4);
  rq_set_pair(&q, 1, 4);

  TEST_ASSERT(rq_div(&p, &r, &q) == RABIN_SUCCESS, "rq_div");
  rq_normalize(&p);
  rq_set_pair(&expect, -1, 1);
  TEST_ASSERT(rq_equal(&p, &expect), "(-1/4) / (1/4) == -1");

  TEST_ASSERT(rq_add(&p, &r, &q) == RABIN_SUCCESS, "rq_add");
  rq_normalize(&p);
  TEST_ASSERT(rq_is_zero(&p), "(-1/4) + (1/4) == 0");

  /* Roundtrip: (x*y)/y == x with non-trivial rationals. */
  rq_set_pair(&r, 5, 6);
  rq_set_pair(&q, -7, 15);
  TEST_ASSERT(rq_mul(&p, &r, &q) == RABIN_SUCCESS, "roundtrip mul");
  TEST_ASSERT(rq_div(&p, &p, &q) == RABIN_SUCCESS, "roundtrip div");
  TEST_ASSERT(rq_equal(&p, &r), "(x*y)/y == x");

  /* Negation twice is identity. */
  TEST_ASSERT(rq_neg(&p) == RABIN_SUCCESS, "rq_neg");
  TEST_ASSERT(rq_neg(&p) == RABIN_SUCCESS, "rq_neg again");
  TEST_ASSERT(rq_equal(&p, &r), "-(-x) == x");

  rq_clear_multi(&r, &q, &p, &expect, NULL);
}

TEST(rat_matrix_addsub)
{
  rmatq_t A = {0}, B = {0}, R_add = {0}, R_sub = {0};
  rq_t entry, expect;
  rq_init_multi(&entry, &expect, NULL);

  TEST_ASSERT(rmatq_init(&A, 5, 5) == RABIN_SUCCESS, "rmatq_init A");
  TEST_ASSERT(rmatq_init(&B, 5, 5) == RABIN_SUCCESS, "rmatq_init B");
  TEST_ASSERT(rmatq_init(&R_add, 5, 5) == RABIN_SUCCESS, "rmatq_init R_add");
  TEST_ASSERT(rmatq_init(&R_sub, 5, 5) == RABIN_SUCCESS, "rmatq_init R_sub");

  rz_set_u64(&A.den, 2);
  rz_set_u64(&B.den, 3);

  rz_t val;
  rz_init(&val);
  for (u64 i = 0; i < 5; i++) {
    for (u64 j = 0; j < 5; j++) {
      rz_set_u64(&val, i + j + 1);
      rmat_set(&A.M, &val, i, j);
      rz_set_u64(&val, i * j + 2);
      rmat_set(&B.M, &val, i, j);
    }
  }
  rz_clear(&val);

  TEST_ASSERT(rmatq_add(&R_add, &A, &B) == RABIN_SUCCESS, "rmatq_add");
  TEST_ASSERT(rmatq_normalize(&R_add) == RABIN_SUCCESS, "normalize R_add");
  TEST_ASSERT(rmatq_sub(&R_sub, &A, &B) == RABIN_SUCCESS, "rmatq_sub");
  TEST_ASSERT(rmatq_normalize(&R_sub) == RABIN_SUCCESS, "normalize R_sub");

  printf("Matrix A + B:\n");
  rmatq_print(&R_add);
  printf("Matrix A - B:\n");
  rmatq_print(&R_sub);

  for (u64 i = 0; i < 5; i++) {
    for (u64 j = 0; j < 5; j++) {
      /* A[i][j] + B[i][j] = (3(i+j+1) + 2(ij+2)) / 6 */
      TEST_ASSERT(rmatq_get(&entry, &R_add, i, j) == RABIN_SUCCESS, "get add");
      TEST_ASSERT(rq_den_positive(&entry), "add entry den > 0");
      rq_set_pair(&expect, (i64)(3 * (i + j + 1) + 2 * (i * j + 2)), 6);
      TEST_ASSERT(rq_equal(&entry, &expect), "add entry matches known value");

      /* A[i][j] - B[i][j] = (3(i+j+1) - 2(ij+2)) / 6 */
      TEST_ASSERT(rmatq_get(&entry, &R_sub, i, j) == RABIN_SUCCESS, "get sub");
      rq_set_pair(&expect, (i64)(3 * (i + j + 1) - 2 * (i * j + 2)), 6);
      TEST_ASSERT(rq_equal(&entry, &expect), "sub entry matches known value");
    }
  }

  rq_clear_multi(&entry, &expect, NULL);
  rmatq_clear(&A);
  rmatq_clear(&B);
  rmatq_clear(&R_add);
  rmatq_clear(&R_sub);
}

int main(void)
{
  rz_init_constants();

  RUN_TEST(rat_arithmetic);
  RUN_TEST(rat_matrix_addsub);

  rz_clear_constants();
  return TEST_REPORT();
}
