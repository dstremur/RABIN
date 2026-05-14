#ifndef BIGU64_H
#define BIGU64_H

#include "bignum.h"

#define INLINE static inline __attribute__((always_inline))

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

typedef struct {
  u64 p;
  u64 p_inv;
  u64 r2_mod_p;  // 2^128 mod p
} mont_ctx;

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
u64 matrix_u64_det_optimized(u64* data, u64 n, const mont_ctx* ctx);

void mont_init(mont_ctx* ctx, u64 p);

#endif