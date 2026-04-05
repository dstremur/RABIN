#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <../include/bignum.h>
// Assuming your bignum functions are available
// char* bn_to_string(bignum* n);
// void generate_random_prime(bignum* p, int bits);

int verify_with_openssl(const char* num_str) {
    char command[5000]; // Increased for 2048-bit strings
    char result_buf[1024];
    int is_prime = 0;
    
    // Redirect 2>&1 to catch warnings and ensure they don't block the pipe
    // Also, use 'grep' logic internally or just parse the result more carefully
    snprintf(command, sizeof(command), "openssl prime %s 2>&1", num_str);

    FILE* fp = popen(command, "r");
    if (fp == NULL) return -1;

    // Read the output line by line
    while (fgets(result_buf, sizeof(result_buf), fp) != NULL) {
        // We are looking for the word "is prime" anywhere in the output
        if (strstr(result_buf, "is prime") != NULL) {
            is_prime = 1;
            break; 
        }
    }

    pclose(fp);
    return is_prime;
}
void run_prime_test(int num_tests, int bits) {
    printf("Starting %d tests for %d-bit primes...\n", num_tests, bits);
    
    for (int i = 0; i < num_tests; i++) {
        bignum p;
        bn_init(&p);
        
        // 1. Generate your prime
        bn_gen_prime(&p, bits);
        
        // 2. Convert to string
        char* p_str = bn_to_string(&p);
        
        // 3. Verify with OpenSSL
        if (verify_with_openssl(p_str)) {
            printf("[PASS] Test %d: %s is prime\n", i + 1, p_str);
        } else {
            printf("[FAIL] Test %d: %s is NOT prime according to OpenSSL\n", i + 1, p_str);
        }

        free(p_str);
        bn_free(&p);
    }
}

int main() {
    // Test 5 random 128-bit primes
    run_prime_test(20, 2048);
    return 0;
}
