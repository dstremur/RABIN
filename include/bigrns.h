#ifndef BIGRNS_H
#define BIGRNS_H

#include "bigmatrix.h"
#include "bignum.h"

typedef struct {
  u64* primes;
  u64 count;
  bignum* crt_weights;
  bignum prod;
} ctx_rns;



/*
 3. Practical Guideline: The "High-Water Mark"

If you are writing a performance-critical application, follow the High-Water
Mark strategy:

    Check: Does my current ctx_rns->prod have enough bits for this new
calculation? (Use your Hadamard estimator).

    Reuse: If yes, just use the existing context. No extra cost.

    Expand: If no, free the old context and initialize a new one with more
primes.


Cache: Never shrink the context unless you are severely low on memory.

precompute barret reduction or montgomery

*/



typedef struct {
  u64* residues;
  u64 size;
} rns_num;

void rns_context_init(ctx_rns* ctx, const u64* primes, u64 count);
void rns_context_free(ctx_rns* ctx);

void bignum_to_rns(rns_num* r, const bignum* a, ctx_rns* ctx);
void rns_to_bignum(bignum* a, const rns_num* r, ctx_rns* ctx);

void rns_add(rns_num* r, const rns_num* a, const rns_num* b,
             const ctx_rns* ctx);
void rns_mul(rns_num* r, const rns_num* a, const rns_num* b,
             const ctx_rns* ctx);


void bigmatrix_det_rns(bignum* det, const bigmatrix* A, const ctx_rns* ctx);
u64 rns_estimate_determinant(bigmatrix* A);
u64 bigrns_estimate_primes(const bignum* a);

#endif
