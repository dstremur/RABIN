#include <fcntl.h>
#include <stdio.h>
#include <time.h>
#include <unistd.h>

#include "../include/rabin.h"
#include "../include/rzcert.h"

int main()
{
  rz_init_constants();

  int num_to_generate = 20;
  int bits = 2048;

  rz_t p;
  rz_init(&p);

  printf("Benchmarking: Generating %d primes at %d-bits...\n", num_to_generate,
         bits);
  printf("--------------------------------------------------\n");

  struct timespec start, end;
  double total_time = 0;

  rzcert_t* cert = NULL;
  for (int i = 0; i < num_to_generate; i++) {
    clock_gettime(CLOCK_MONOTONIC, &start);

    rz_gen_provable_arithmetic(&p, bits, &cert);
    clock_gettime(CLOCK_MONOTONIC, &end);

    double elapsed =
        (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
    total_time += elapsed;

    rz_println(&p);
    printf("Prime #%d: Found in %.4f seconds\n", i + 1, elapsed);

    print_pocklington_cert(cert);
    rzcert_clear(cert);
  }

  printf("--------------------------------------------------\n");
  printf("Total Time:   %.4f seconds\n", total_time);
  printf("Average Time: %.4f seconds per prime\n",
         total_time / num_to_generate);

  rz_clear(&p);

  rz_clear_constants();
  return 0;
}
