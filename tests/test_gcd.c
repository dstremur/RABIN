#include <stdio.h>
#include <time.h>

#include "../include/bignum.h"

void benchmark_gcd(int limbs)
{
  bignum a, b, u_1, v_1, d_1, u_2, v_2, d_2;
  bn_init_multi(&a, &b, &u_1, &v_1, &d_1, &u_2, &v_2, &d_2, NULL);

  // Generate random numbers of a specific bit length (limbs * 64)
  bn_gen_random(&a, limbs * 64);
  bn_gen_random(&b, limbs * 64);

  // 1. Time Normal GCD
  clock_t start = clock();
  bn_gcd(&d_1, &a, &b);
  clock_t end = clock();
  double time_normal = (double)(end - start) / CLOCKS_PER_SEC;

  // 2. Time Karatsuba
  start = clock();
  bn_gcd_extended(&u_2, &v_2, &d_2, &a, &b);
  end = clock();
  double time_ex = (double)(end - start) / CLOCKS_PER_SEC;

  start = clock();
  bn_gcd_lehmer(&d_2, &a, &b);
  end = clock();
  double time_l = (double)(end - start) / CLOCKS_PER_SEC;

  start = clock();
  bn_gcd_extended_lehmer(&u_2, &v_2, &d_2, &a, &b);
  end = clock();
  double time_lehmer = (double)(end - start) / CLOCKS_PER_SEC;

  // 3. Verify Correctness
  if (bn_cmp(&d_1, &d_2) != 0) {
    printf("[FAIL] Mismatch at %d limbs!\n", limbs);
  } else {
    printf(
        "Limbs: %4d | Normal: %fs | Extended: %fs | Lehmer: %fs | Lehmer "
        "Extended: "
        "%fs "
        "\n",
        limbs, time_normal, time_ex, time_l, time_lehmer);
  }

  bn_free_multi(&a, &b, &u_1, &v_1, &d_1, &u_2, &v_2, &d_2, NULL);
}

int main()
{
  bn_init_constants();
  // Test from small to large to find the crossover point
  int sizes[] = {16,    32,     64,     128,    256,    512,
                 1024,  2048,   4096,   8192,   16000,  32000,
                 64000, 100000, 200000, 256000, 300000, 400000};
  for (int i = 0; i < 12; i++) {
    benchmark_gcd(sizes[i]);
  }

  bn_free_constants();
  return 0;
}
