#include <stdio.h>
#include <time.h>

#include "../include/bignum.h"

void benchmark_mul(int limbs)
{
  bignum a, b, res_school, res_karat;
  bn_init_multi(&a, &b, &res_school, &res_karat, NULL);

  // Generate random numbers of a specific bit length (limbs * 64)
  bn_gen_random(&a, limbs * 64);
  bn_gen_random(&b, limbs * 64);

  // 1. Time Schoolbook
  clock_t start = clock();
  bn_mul_school(&res_school, &a, &b);
  clock_t end = clock();
  double time_school = (double)(end - start) / CLOCKS_PER_SEC;

  // 2. Time Karatsuba
  start = clock();
  bn_mul(&res_karat, &a, &b);
  end = clock();
  double time_karat = (double)(end - start) / CLOCKS_PER_SEC;

  // 3. Verify Correctness
  if (bn_cmp(&res_school, &res_karat) != 0) {
    printf("[FAIL] Mismatch at %d limbs!\n", limbs);
  } else {
    printf("Limbs: %4d | School: %fs | Karatsuba: %fs | Speedup: %.2fx\n",
           limbs, time_school, time_karat, time_school / time_karat);
  }

  bn_free(&a);
  bn_free(&b);
  bn_free(&res_school);
  bn_free(&res_karat);
}

int main()
{
  // Test from small to large to find the crossover point
  int sizes[] = {16,   32,   64,   128,   256,   512,   1024,
                 2048, 4096, 8192, 16000, 32000, 64000, 100000};
  for (int i = 0; i < 14; i++) {
    benchmark_mul(sizes[i]);
  }
  return 0;
}
