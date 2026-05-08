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

u64 mod_add(u64 a, u64 b, u64 p);
u64 mod_sub(u64 a, u64 b, u64 p);
u64 mod_mul(u64 a, u64 b, u64 p);
u64 mod_pow(u64 base, u64 exp, u64 p);
u64 mod_inverse_prime(u64 n, u64 p);

typedef struct {
	u64* residues;
} rns_num; 

void rns_context_init(ctx_rns* ctx, const u64* primes, u64 count);
void rns_context_free(ctx_rns* ctx);

void bignum_to_rns(rns_num* r, const bignum* a, ctx_rns* ctx);
void rns_to_bignum(rns_num* r, const bignum* a, ctx_rns* ctx);

void rns_add(rns_num* r, const rns_num* a, const rns_num* b, const ctx_rns* ctx);
void rns_mul(rns_num* r, const rns_num* a, const rns_num* b, const ctx_rns* ctx);

u64 matrix_u64_det(matrix_u64* A);
void bigmatrix_det_rns(bignum* det, const bigmatrix* A); 

u64 bigrns_estimate_primes(const bignum* a);



#endif 
