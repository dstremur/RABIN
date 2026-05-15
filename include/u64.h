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
  u64 p_inv;
  u64 r2_mod_p;  // 2^128 mod p
} mont_ctx;

u64 mod_add(u64 a, u64 b, u64 p);
u64 mod_sub(u64 a, u64 b, u64 p);
u64 mod_mul(u64 a, u64 b, u64 p);
u64 mod_pow(u64 base, u64 exp, u64 p);
u64 mod_inverse_prime(u64 n, u64 p);
u64 mod_inverse_euclid(u64 a, u64 p);

u64 matrix_u64_det(matrix_u64* A);
u64 matrix_u64_det_optimized(u64* data, u64 n, const mont_ctx* ctx);

void mont_init(mont_ctx* ctx, u64 p);

static inline u64 mod_mul_mont(u64 a, u64 b, const mont_ctx* ctx)
{
  unsigned __int128 T = (unsigned __int128)a * b;
  u64 m = (u64)T * ctx->p_inv;
  // t = T + m * p 
  unsigned __int128 t = T + (unsigned __int128)m * ctx->p;

  u64 res = (u64)(t >> 64);
  if (res >= ctx->p) res -= ctx->p;
  return res;
}

// Convert standard number -> Montgomery form
static inline u64 to_mont(u64 x, const mont_ctx* ctx)
{
  return mod_mul_mont(x, ctx->r2_mod_p, ctx);
}

// Convert Montgomery form -> standard number
static inline u64 from_mont(u64 x, const mont_ctx* ctx)
{
  return mod_mul_mont(x, 1, ctx);
}


#endif