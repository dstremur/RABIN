#ifndef BIGRNS_H
#define BIGRNS_H 

#include "bignum.h"
#include "bigmatrix.h"

typedef struct {
		u64* primes;
		u64 count;
		bignum* crt_weights;
		bignum prod;
} ctx_rns; 

typedef struct {
    u64* data;
    u64 r_size;
    u64 c_size;
    u64 modulus; // The prime for this specific slice
} matrix_u64;

/*
 3. Practical Guideline: The "High-Water Mark"

If you are writing a performance-critical application, follow the High-Water Mark strategy:

    Check: Does my current ctx_rns->prod have enough bits for this new calculation? (Use your Hadamard estimator).

    Reuse: If yes, just use the existing context. No extra cost.

    Expand: If no, free the old context and initialize a new one with more primes.


Cache: Never shrink the context unless you are severely low on memory.

precompute barret reduction or montgomery

*/ 





u64 mod_add(u64 a, u64 b, u64 p);
u64 mod_sub(u64 a, u64 b, u64 p);
u64 mod_mul(u64 a, u64 b, u64 p);
u64 mod_pow(u64 base, u64 exp, u64 p);
u64 mod_inverse_prime(u64 n, u64 p);

typedef struct {
	u64* residues;
	u64 size;
} rns_num; 

void rns_context_init(ctx_rns* ctx, const u64* primes, u64 count);
void rns_context_free(ctx_rns* ctx);

void bignum_to_rns(rns_num* r, const bignum* a, ctx_rns* ctx);
void rns_to_bignum(bignum* a, const rns_num* r, ctx_rns* ctx);

void rns_add(rns_num* r, const rns_num* a, const rns_num* b, const ctx_rns* ctx);
void rns_mul(rns_num* r, const rns_num* a, const rns_num* b, const ctx_rns* ctx);

u64 matrix_u64_det(matrix_u64* A);
void bigmatrix_det_rns(bignum* det, const bigmatrix* A, const ctx_rns* ctx); 
u64 rns_estimate_determinant(bigmatrix* A);
u64 bigrns_estimate_primes(const bignum* a);



#endif 
