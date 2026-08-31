#include <fcntl.h>
#include <stdio.h>
#include <time.h>
#include <unistd.h>

#include "../include/bignum.h"

#define ITERATIONS 1000

int main()
{
  int num_to_generate = 1;
  int bits = 2048;

  bignum p;
  bn_init(&p);

  printf("Benchmarking: Generating %d primes at %d-bits...\n", num_to_generate,
         bits);
  printf("--------------------------------------------------\n");

  struct timespec start, end;
  double total_time = 0;

  for (int i = 0; i < num_to_generate; i++) {
    clock_gettime(CLOCK_MONOTONIC, &start);

    if (bn_gen_prime(&p, bits)) {
      clock_gettime(CLOCK_MONOTONIC, &end);

      double elapsed =
          (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
      total_time += elapsed;

      printf("Prime #%d: Found in %.4f seconds\n", i + 1, elapsed);
    } else {
      printf("Failed to generate prime #%d\n", i + 1);
    }
  }

  printf("--------------------------------------------------\n");
  printf("Total Time:   %.4f seconds\n", total_time);
  printf("Average Time: %.4f seconds per prime\n",
         total_time / num_to_generate);

  printf("Benchmarking bn_gen_prime base cases...\n");
  printf("-------------------------------------------------\n");
  printf("%-10s | %-15s | %-15s\n", "Bit Len", "Total Time (s)",
         "Avg Time (s)");
  printf("-------------------------------------------------\n");

  for (int i = 5; i <= 100; i++) {
    int success_count = 0;

    // Start timer just before the batch execution
    clock_gettime(CLOCK_MONOTONIC, &start);

    for (int j = 0; j < ITERATIONS; j++) {
      if (bn_gen_prime(&p, i)) {
        success_count++;
      }
    }

    // Stop timer immediately after
    clock_gettime(CLOCK_MONOTONIC, &end);

    if (success_count == ITERATIONS) {
      double total_time =
          (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
      double avg_time = total_time / ITERATIONS;

      printf("%-10d | %-15.6f | %-15.9f\n", i, total_time, avg_time);
    } else {
      printf("Failed at bit length %d (%d/%d successes)\n", i, success_count,
             ITERATIONS);
    }
  }

  bn_free(&p);
  return 0;
}
