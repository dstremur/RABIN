#include <fcntl.h>
#include <stdio.h>
#include <time.h>
#include <unistd.h>
#include <omp.h>

#include "../include/bignum.h"

int main()
{
    int num_to_generate = 200;
    int bits = 2048;

    printf("Benchmarking: Generating %d primes at %d-bits (Parallel)...\n", num_to_generate, bits);
    printf("--------------------------------------------------\n");

    struct timespec global_start, global_end;
    clock_gettime(CLOCK_MONOTONIC, &global_start);

    // Use OpenMP to distribute the work of the for loop across threads
    #pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i < num_to_generate; i++) {
        struct timespec start, end;
        
        // Each thread must have its own bignum to work with
        bignum p;
        bn_init(&p);
        
        clock_gettime(CLOCK_MONOTONIC, &start);

        if (bn_gen_prime(&p, bits)) {
            clock_gettime(CLOCK_MONOTONIC, &end);

            double elapsed = (end.tv_sec - start.tv_sec) + 
                             (end.tv_nsec - start.tv_nsec) / 1e9;

            // Use a critical section for printing to prevent text interleaving
            #pragma omp critical
            {
                printf("Thread %d | Prime #%d: Found in %.4f seconds\n", 
                        omp_get_thread_num(), i + 1, elapsed);
            }
        } else {
            #pragma omp critical
            {
                printf("Failed to generate prime #%d\n", i + 1);
            }
        }

        bn_free(&p);
    }

    clock_gettime(CLOCK_MONOTONIC, &global_end);
    double total_wall_time = (global_end.tv_sec - global_start.tv_sec) + 
                             (global_end.tv_nsec - global_start.tv_nsec) / 1e9;

    printf("--------------------------------------------------\n");
    printf("Total Wall-Clock Time: %.4f seconds\n", total_wall_time);
    printf("Efficiency: %.4f primes/second\n", num_to_generate / total_wall_time);

    return 0;
}
