#include <stdio.h>
#include <time.h>

#include "../include/bignum.h"
#include "../include/u64.h"
#define ITERATIONS 1000000

void benchmark_mul(void)
{
  u64 p = 180143985094819841ULL;

  u64 a = 9023482357230203423ULL;
  u64 b = 8123482357230203421ULL;

  u64 mu = compute_mu(p);

  mont_ctx mctx;
  mont_init(&mctx, p);

  /* ---------- Barrett ---------- */

  u64 r_barrett = 0;

  clock_t start = clock();

  for (u64 i = 0; i < ITERATIONS; i++) {
    unsigned __int128 T = (unsigned __int128)(a + i) * (b + i);

    r_barrett = barrett_reduction(T, p, mu);
  }

  clock_t end = clock();

  double time_barrett = (double)(end - start) / CLOCKS_PER_SEC;

  /* ---------- Montgomery ---------- */

  u64 a_hat = mont_in(a, &mctx);
  u64 b_hat = mont_in(b, &mctx);

  u64 r_hat = 0;

  start = clock();

  for (u64 i = 0; i < ITERATIONS; i++) {
    u64 aa = mont_in(a + i, &mctx);
    u64 bb = mont_in(b + i, &mctx);

    r_hat = mont_mul(aa, bb, &mctx);
  }

  end = clock();

  double time_mont = (double)(end - start) / CLOCKS_PER_SEC;

  u64 r_mont = mont_out(r_hat, &mctx);

  printf("Barrett     : %llu\n", (unsigned long long)r_barrett);

  printf("Montgomery  : %llu\n", (unsigned long long)r_mont);

  printf("Barrett time    : %f s\n", time_barrett);
  printf("Montgomery time : %f s\n", time_mont);

  printf("R2 = %llu\n", (unsigned long long)mctx.r2_mod_p);

  u64 lhs = a;
  u64 rhs = b;

  unsigned __int128 product = (unsigned __int128)lhs * rhs;

  u64 barrett = barrett_reduction(product, p, mu);

  u64 lhs_hat = mont_in(lhs, &mctx);
  u64 rhs_hat = mont_in(rhs, &mctx);
  u64 mont = mont_out(mont_mul(lhs_hat, rhs_hat, &mctx), &mctx);

  printf("barrett = %llu\n", (unsigned long long)barrett);
  printf("mont    = %llu\n", (unsigned long long)mont);
}

int main()
{
  bn_init_constants();

  benchmark_mul();

  bn_free_constants();
  return 0;
}