#include <fcntl.h>
#include <stdio.h>
#include <time.h>
#include <unistd.h>

#include "../include/bigcert.h"
#include "../include/bignum.h"

int main()
{
  bn_init_constants();

  int num_to_generate = 10;
  int bits = 2048;

  bignum p;
  bn_init(&p);

  printf("Benchmarking: Generating %d primes at %d-bits...\n", num_to_generate,
         bits);
  printf("--------------------------------------------------\n");

  struct timespec start, end;
  double total_time = 0;

  pocklington_cert* cert = NULL;
  for (int i = 0; i < num_to_generate; i++) {
    clock_gettime(CLOCK_MONOTONIC, &start);

    gen_provable_primes_arithmetic(&p, bits, &cert);
    clock_gettime(CLOCK_MONOTONIC, &end);

    double elapsed =
        (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
    total_time += elapsed;

    bn_println(&p);
    printf("Prime #%d: Found in %.4f seconds\n", i + 1, elapsed);

    // print_pocklington_cert(cert);
  }

  printf("--------------------------------------------------\n");
  printf("Total Time:   %.4f seconds\n", total_time);
  printf("Average Time: %.4f seconds per prime\n",
         total_time / num_to_generate);

  bn_free(&p);

  bn_free_constants();
  return 0;
}
