#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "../include/bignum.h"
// [Assume your bignum headers and functions are included here]

void benchmark_div_real(int divisor_limbs)
{
    // 1. The 2N / N Rule
    // To properly test division, the dividend (a) must be twice as large as the divisor (b).
    int dividend_limbs = divisor_limbs * 2; 

    bignum a, b, res_school, res_karat;
    bn_init_multi(&a, &b, &res_school, &res_karat, NULL);

    // Generate random numbers
    bn_gen_random(&a, dividend_limbs * 64);
    bn_gen_random(&b, divisor_limbs * 64);

    // Ensure the highest bit of 'b' is set so it's truly 'divisor_limbs' long
    bn_set_bit(&b, (divisor_limbs * 64) - 1); 

    // 2. Time Schoolbook
    clock_t start = clock();
    bn_div(&res_school, &a, &b);
    clock_t end = clock();
    double time_school = ((double)(end - start)) / CLOCKS_PER_SEC;

    // 3. Time Newton-Raphson
    start = clock();
    bn_newton_div(&res_karat, &a, &b);
    end = clock();
    double time_karat = ((double)(end - start)) / CLOCKS_PER_SEC;

    // 4. Verify Correctness
    if (bn_cmp(&res_school, &res_karat) != 0) {
        printf("[FAIL] Mismatch! Divisor Limbs: %d\n", divisor_limbs);
        
        // Print the sizes of the generated quotients to see how they failed
        printf("       School quotient limbs: %llu\n", res_school.size);
        printf("       Newton quotient limbs: %llu\n", res_karat.size);
    } else {
        double speedup = (time_karat > 0.0) ? (time_school / time_karat) : 0.0;
        printf("Divisor Limbs: %5d | Dividend Limbs: %5d | School: %8.6fs | Newton: %8.6fs | Speedup: %8.2fx\n",
               divisor_limbs, dividend_limbs, time_school, time_karat, speedup);
    }

    // 5. Anti-Optimization Hack
    // We add the lowest limb of the result to a global volatile variable 
    // so the compiler is FORCED to actually compute the division.
    extern volatile u64 dummy_sum; 
    if (res_school.size > 0) dummy_sum += res_school.limbs[0];
    if (res_karat.size > 0) dummy_sum += res_karat.limbs[0];

    bn_free_multi(&a, &b, &res_school, &res_karat, NULL);
}

// Global variable to prevent compiler dead-code elimination
volatile u64 dummy_sum = 0;

int main()
{
    printf("Starting 2N / N Division Benchmark...\n");
    printf("------------------------------------------------------------------------------------------------\n");
    
    int sizes[] = {16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16000};
    int num_sizes = sizeof(sizes) / sizeof(sizes[0]);

    for (int i = 0; i < num_sizes; i++) {
        benchmark_div_real(sizes[i]);
    }

    printf("------------------------------------------------------------------------------------------------\n");
    printf("Benchmark complete. (Checksum: %llu)\n", dummy_sum);
    
    return 0;
}
