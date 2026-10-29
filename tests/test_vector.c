/*
 * test_vector.c
 *
 * Assert-based checks for dynamic rvec_t append/neg/dot/norm.
 *
 * Copyright (C) 2026 Diego Strebel
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include "../include/rvec.h"
#include "test_util.h"

TEST(vector_append_and_neg)
{
  rvec_t a = {0};

  TEST_ASSERT(rvec_init_dynamic(&a) == RABIN_SUCCESS, "rvec_init_dynamic");
  TEST_ASSERT(a.size == 0, "fresh dynamic vector has size 0");

  for (int i = 0; i < 8; i++) {
    TEST_ASSERT(rvec_append(&a, &RZ_ONE) == RABIN_SUCCESS, "rvec_append");
  }
  TEST_ASSERT(a.size == 8, "size is 8 after eight appends");

  printf("Vector a:\n");
  rvec_println(&a);

  for (u64 i = 0; i < a.size; i++) {
    TEST_ASSERT(rz_cmp(&a.data[i], &RZ_ONE) == 0, "entry equals 1");
  }

  TEST_ASSERT(rvec_neg(&a, &a) == RABIN_SUCCESS, "rvec_neg");
  for (u64 i = 0; i < a.size; i++) {
    TEST_ASSERT(rz_cmp(&a.data[i], &RZ_ONE) == -1, "entry equals -1");
  }

  printf("Vector -a:\n");
  rvec_println(&a);

  rvec_clear(&a);
}

TEST(vector_dot_and_norm)
{
  rvec_t v = {0};
  rz_t dot = {0}, norm = {0}, expected = {0}, tmp = {0};

  TEST_ASSERT(rvec_init(&v, 8) == RABIN_SUCCESS, "rvec_init");
  for (u64 i = 0; i < 8; i++) {
    TEST_ASSERT(rvec_set(&v, &RZ_ONE, i) == RABIN_SUCCESS, "rvec_set");
  }

  rz_init_multi(&dot, &norm, &expected, &tmp, NULL);

  /* All entries are +-1, so dot(v, v) must be exactly 8. */
  TEST_ASSERT(rvec_neg(&v, &v) == RABIN_SUCCESS, "rvec_neg");
  TEST_ASSERT(rvec_dot(&dot, &v, &v) == RABIN_SUCCESS, "rvec_dot");
  rz_set_u64(&expected, 8);
  TEST_ASSERT(rz_cmp(&dot, &expected) == 0, "dot(v,v) == 8");

  TEST_ASSERT(rvec_norm(&norm, &v) == RABIN_SUCCESS, "rvec_norm");
  rz_isqrt(&expected, &dot);
  TEST_ASSERT(rz_cmp(&norm, &expected) == 0, "norm equals isqrt(dot(v,v))");

  /* norm(a) == isqrt(dot(a, a)) after another sign flip as well. */
  TEST_ASSERT(rvec_neg(&v, &v) == RABIN_SUCCESS, "rvec_neg mixed");
  TEST_ASSERT(rvec_dot(&tmp, &v, &v) == RABIN_SUCCESS, "rvec_dot mixed");
  TEST_ASSERT(rvec_norm(&norm, &v) == RABIN_SUCCESS, "rvec_norm mixed");
  rz_isqrt(&tmp, &tmp);
  TEST_ASSERT(rz_cmp(&norm, &tmp) == 0, "norm == isqrt(dot) for mixed vector");

  rvec_clear(&v);
  rz_clear_multi(&dot, &norm, &expected, &tmp, NULL);
}

int main(void)
{
  rz_init_constants();

  RUN_TEST(vector_append_and_neg);
  RUN_TEST(vector_dot_and_norm);

  rz_clear_constants();
  return TEST_REPORT();
}
