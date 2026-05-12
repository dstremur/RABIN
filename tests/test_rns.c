#include "../include/bigrns.h"

#include <stdio.h>

int main() {
	// 1. Setup Context
ctx_rns ctx;
u64 primes[] = {0xFFFFFFFFFFFFFFC5ULL, 0xFFFFFFFFFFFFFF4DULL}; // Example 64-bit primes
rns_context_init(&ctx, primes, 2);

// 2. Prepare Numbers
bignum n;
bn_init(&n);
bn_init_val(&n, "100112324252335235235252350");

rns_num r;
r.residues = NULL; // bignum_to_rns will malloc if NULL

// 3. Convert and Convert Back
bignum_to_rns(&r, &n, &ctx);



bignum reconstructed;
bn_init(&reconstructed);

rns_to_bignum(&reconstructed, &r, &ctx);

// Now: reconstructed == my_big

printf("original :");
bn_println(&n);

printf("rns :");
bn_println(&reconstructed);

bn_free_multi(&n, &reconstructed, NULL);
}
