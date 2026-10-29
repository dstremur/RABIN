/*
 * test_u64.c
 *
 * Known-answer and cross-consistency checks for u64 Barrett and Montgomery
 * arithmetic, plus a timing benchmark.
 *
 * Copyright (C) 2026 Diego Strebel
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <stdio.h>
#include <time.h>

#include "../include/rabin.h"
#include "../include/u64.h"
#include "test_util.h"

#define ITERATIONS 1000000

static u64 mul_mod(u64 a, u64 b, u64 p) { return (u64)(((u128)a * b) % p); }

TEST(u64_mont_context)
{
  u64 p = 180143985094819841ULL;
  u64_mont_ctx_t mctx;
  u64_mont_init(&mctx, p);

  /* r2_mod_p must be exactly 2^128 mod p. */
  u64 r2 = 1;
  for (int i = 0; i < 128; i++) {
    r2 = (u64)(((u128)r2 * 2) % p);
  }
  TEST_ASSERT(mctx.p == p, "context keeps p");
  TEST_ASSERT(mctx.r2_mod_p == r2, "r2_mod_p == 2^128 mod p");
}

TEST(u64_barrett_known_answers)
{
  u64 p = 180143985094819841ULL;
  u128 mu = u64_compute_mu(p);

  const u64 xs[] = {
      0, 1, p - 1, p, p + 1, 9023482357230203423ULL, 8123482357230203421ULL};
  for (size_t i = 0; i < sizeof(xs) / sizeof(xs[0]); i++) {
    for (size_t j = 0; j < sizeof(xs) / sizeof(xs[0]); j++) {
      u128 c = (u128)xs[i] * xs[j];
      u64 got = u64_barrett(c, p, mu);
      u64 want = mul_mod(xs[i], xs[j], p);
      if (got != want) {
        printf("barrett mismatch: %llu*%llu mod p: got %llu want %llu\n",
               (unsigned long long)xs[i], (unsigned long long)xs[j],
               (unsigned long long)got, (unsigned long long)want);
      }
      TEST_ASSERT(got == want, "u64_barrett == exact remainder");
    }
  }
}

TEST(u64_montgomery_roundtrip)
{
  u64 p = 180143985094819841ULL;
  u64_mont_ctx_t mctx;
  u64_mont_init(&mctx, p);

  const u64 xs[] = {
      0, 1, 42, p - 1, 9023482357230203423ULL, 8123482357230203421ULL};
  for (size_t i = 0; i < sizeof(xs) / sizeof(xs[0]); i++) {
    u64 x = xs[i] % p;
    TEST_ASSERT(u64_mont_out(u64_mont_in(xs[i], &mctx), &mctx) == x,
                "mont_in/out roundtrip");
    for (size_t j = 0; j < sizeof(xs) / sizeof(xs[0]); j++) {
      u64 a_hat = u64_mont_in(xs[i], &mctx);
      u64 b_hat = u64_mont_in(xs[j], &mctx);
      u64 got = u64_mont_out(u64_mont_mul(a_hat, b_hat, &mctx), &mctx);
      TEST_ASSERT(got == mul_mod(xs[i], xs[j], p),
                  "mont mul == exact remainder");
    }
  }
}

TEST(u64_barrett_montgomery_agree)
{
  u64 p = 180143985094819841ULL;
  u128 mu = u64_compute_mu(p);
  u64_mont_ctx_t mctx;
  u64_mont_init(&mctx, p);

  u64 a = 9023482357230203423ULL;
  u64 b = 8123482357230203421ULL;

  clock_t start = clock();
  u64 r_barrett = 0;
  for (u64 i = 0; i < ITERATIONS; i++) {
    u128 T = (u128)(a + i) * (b + i);
    r_barrett = u64_barrett(T, p, mu);
  }
  clock_t end = clock();
  double time_barrett = (double)(end - start) / CLOCKS_PER_SEC;

  u64 r_hat = 0;
  start = clock();
  for (u64 i = 0; i < ITERATIONS; i++) {
    u64 aa = u64_mont_in(a + i, &mctx);
    u64 bb = u64_mont_in(b + i, &mctx);
    r_hat = u64_mont_mul(aa, bb, &mctx);
  }
  end = clock();
  double time_mont = (double)(end - start) / CLOCKS_PER_SEC;

  u64 r_mont = u64_mont_out(r_hat, &mctx);

  printf("Barrett time    : %f s\n", time_barrett);
  printf("Montgomery time : %f s\n", time_mont);

  TEST_ASSERT(r_barrett == r_mont, "barrett and montgomery agree");
  TEST_ASSERT(r_barrett == mul_mod(a + ITERATIONS - 1, b + ITERATIONS - 1, p),
              "final loop result matches exact remainder");
}

int main(void)
{
  rz_init_constants();

  RUN_TEST(u64_mont_context);
  RUN_TEST(u64_barrett_known_answers);
  RUN_TEST(u64_montgomery_roundtrip);
  RUN_TEST(u64_barrett_montgomery_agree);

  rz_clear_constants();
  return TEST_REPORT();
}
