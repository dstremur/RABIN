#ifndef BIGNTT_H
#define BIGNTT_H

#include "bignum.h"
#include "bigpoly.h"

typedef struct {
    u64 n;                  // length must be power of 2
    bignum q;               // prime modulus
    bignum n_inv;           // n^-1 mod q 
    
    // Twiddle factors 
    bignum* omega_powers;   // w^i mod q for i = 0 to n-1
    bignum* omega_inv_powers; // w^-i mod q for INTT
    
    // Optional: 2n-th root of unity for negacyclic convolution
    bignum* psi_powers;     // psi^i mod q 
    bignum* psi_inv_powers; // psi^-i mod q
    
    // Bit-reversal permutation array to avoid recomputing indices
    u64* bit_rev_indices;   
} ntt_ctx;

bool ntt_ctx_init(ntt_ctx* ctx, u64 n, const bignum* q, const bignum* omega, const bignum* psi);

bool bigntt_find_prime(bignum* q, u64 n, u64 bits);
bool bigntt_find_generator(bignum* g, const bignum* q); 
bool bigntt_compute_roots(bignum* omega, bignum* psi, u64 n, const bignum* q, const bignum* g);

void ntt_ctx_free(ntt_ctx* ctx); 

// forward NTT using Cooley-Tukey butterflies
// input in normal order, output in bit-reversed order
void bigntt_forward(bigpoly* out_hat, const bigpoly* in, const ntt_ctx* ctx);

// inverse NTT using Gentleman-Sande (GS) butterflies
// input in bit-reversed order, output in normal order
void bigntt_inverse(bigpoly* out, const bigpoly* in_hat, const ntt_ctx* ctx);

void bigntt_pointwise_mul(bigpoly* r_hat, const bigpoly* a_hat, const bigpoly* b_hat, const ntt_ctx* ctx);

// Performs a cyclic (positive-wrapped) convolution: r = INTT(NTT(a) ◦ NTT(b))
void bigpoly_conv_cyclic_ntt(bigpoly* r, const bigpoly* a, const bigpoly* b, const ntt_ctx* ctx);

// Performs a negacyclic (negative-wrapped) convolution using psi
void bigpoly_conv_negacyclic_ntt(bigpoly* r, const bigpoly* a, const bigpoly* b, const ntt_ctx* ctx);

// Standard fast polynomial multiplication.
void bigpoly_mul_ntt(bigpoly* r, const bigpoly* a, const bigpoly* b, const bignum* q);

#endif 
