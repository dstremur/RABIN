#ifndef BIGU64_H
#define BIGU64_H

#include "bignum.h"

typedef struct {
  u64* data;
  u64 r_size;
  u64 c_size;
  u64 modulus;  // The prime for this specific slice
} matrix_u64;


typedef struct {
  u64 p;
  unsigned __int128 mu;
} barrett_ctx;


static inline unsigned __int128 mul128_high(unsigned __int128 a, unsigned __int128 b) {
    u64 a_lo = (u64)a, a_hi = (u64)(a >> 64);
    u64 b_lo = (u64)b, b_hi = (u64)(b >> 64);

    unsigned __int128 p00 = (unsigned __int128)a_lo * b_lo;
    unsigned __int128 p10 = (unsigned __int128)a_hi * b_lo;
    unsigned __int128 p01 = (unsigned __int128)a_lo * b_hi;
    unsigned __int128 p11 = (unsigned __int128)a_hi * b_hi;

    // Intermediate sum of middle products
    unsigned __int128 middle = p10 + p01 + (p00 >> 64);
    
    // The high 128 bits consist of p11 plus the carry-out from the middle
    return p11 + (middle >> 64);
}

void barrett_init(barrett_ctx* ctx, u64 p);

u64 mod_mul_barrett(u64 a, u64 b, const barrett_ctx* ctx);
u64 barrett_reduce(unsigned __int128 a, const barrett_ctx* ctx);
u64 mod_add(u64 a, u64 b, u64 p);
u64 mod_sub(u64 a, u64 b, u64 p);
u64 mod_mul(u64 a, u64 b, u64 p);
u64 mod_pow(u64 base, u64 exp, u64 p);
u64 mod_inverse_prime(u64 n, u64 p);
u64 mod_inverse_euclid(u64 a, u64 p);

u64 matrix_u64_det(matrix_u64* A);
u64 matrix_u64_det_optimized(u64* data, u64 n, const barrett_ctx* ctx);


#endif