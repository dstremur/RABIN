#include <fcntl.h>
#include <stdio.h>
#include <time.h>
#include <unistd.h>

#include "../include/bignum.h"

int main()
{
  int num_to_generate = 5;
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

  bn_free(&p);
  return 0;
}
