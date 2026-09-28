#include <stdio.h>
#include <time.h>

#include "../include/rabin.h"

int main()
{
  int iterations = 100;  // Adjust based on how long you want to wait
  int bits = 2048;
  rz_t n;
  rz_init(&n);

  // Get a high-quality prime to ensure the test runs the full BPSW path
  if (rz_gen_prime(&n, bits) != RABIN_SUCCESS) return 1;

  printf("Benchmarking BPSW: %d iterations at %d bits...\n", iterations, bits);

  struct timespec start, end;
  clock_gettime(CLOCK_MONOTONIC, &start);

  for (int i = 0; i < iterations; i++) {
    rz_bpsw(&n);
  }

  clock_gettime(CLOCK_MONOTONIC, &end);

  double total_time =
      (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
  double throughput = iterations / total_time;

  printf("------------------------------------------\n");
  printf("Total Time:  %.4f s\n", total_time);
  printf("Throughput:  %.2f tests/sec\n", throughput);
  printf("Avg Latency: %.2f ms\n", (total_time / iterations) * 1000.0);
  printf("------------------------------------------\n");

  rz_clear(&n);
  return 0;
}
