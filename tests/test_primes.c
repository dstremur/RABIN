#include <fcntl.h>
#include <stdio.h>
#include <time.h>
#include <unistd.h>

#include "../include/bignum.h"

int main()
{
  int num_to_generate = 200;
  int bits = 2048;

    int bounds[] = {50, 100, 200, 500, 750, 1000, 1250, 1500, 1750, 2000, 3000, 4000, 5000, 6000};
    int num_bounds = sizeof(bounds) / sizeof(bounds[0]);

    FILE* csv = fopen("benchmark.csv", "w");
    fprintf(csv, "trial_bound,total_time,avg_time\n");

	srand(42);

  bignum p;
  bn_init(&p);


  for (int i = 0; i < num_bounds; i++) { 
  printf("Benchmarking: Generating %d primes at %d-bits... Trial: %d\n", num_to_generate,
         bits, bounds[i]);
  printf("--------------------------------------------------\n");

  u64 bound = bounds[i];
  struct timespec start, end;
  double total_time = 0;

  for (int i = 0; i < num_to_generate; i++) {
    clock_gettime(CLOCK_MONOTONIC, &start);

    if (bn_gen_prime(&p, bits, bound)) {
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

  // write CSV row
        fprintf(csv, "%d,%.6f,%.6f\n",
               (int) bound, total_time, total_time / num_to_generate);
  }
  bn_free(&p);

  fclose(csv);

    printf("\nCSV written to benchmark.csv\n");
  return 0;
}
