#ifndef BIGU64_H
#define BIGU64_H

#include "bignum.h"
#include "stdio.h"
#include "string.h"
typedef unsigned __int128 u128;

#define INLINE static inline __attribute__((always_inline))

typedef struct {
  u64* data;
  u64 r_size;
  u64 c_size;
  u64 modulus;  // The prime for this specific slice
} matrix_u64;

typedef struct {
  u64 n;
  u64 k;
  u64 q;
  u64 n_inv;              // in Montgomery form
  u64* omega_powers;      // all in Montgomery form
  u64* omega_inv_powers;  // all in Montgomery form
  u64* bit_rev_indices;
  mont_ctx mctx;
} ntt_ctx_u64;

u64 mod_add(u64 a, u64 b, u64 p);
u64 mod_sub(u64 a, u64 b, u64 p);
u64 mod_mul(u64 a, u64 b, u64 p);
u64 mod_pow(u64 base, u64 exp, u64 p);
u64 mod_inverse_prime(u64 n, u64 p);
u64 mod_inverse_euclid(u64 a, u64 p);
u64 mod_inverse(u64 n, u64 p);

u64 matrix_u64_det(matrix_u64* A);
u64 matrix_u64_det_optimized(u64* data, u64 n, const mont_ctx* ctx);

void mont_init(mont_ctx* ctx, u64 p);
u64 mont_inverse(u64 a_mont, const mont_ctx* ctx);
u64 mont_mul(u64 a, u64 b, mont_ctx* ctx);
u64 mont_in(u64 a, mont_ctx* ctx);
u64 mont_out(u64 a, mont_ctx* ctx);
u64 mont_redc(unsigned __int128 T, mont_ctx* ctx);

u64 barrett_reduction(u128 c, u64 q, u64 mu);
u64 compute_mu(u64 q);
u64 goldilock_red(u128 c);

bool ntt_ctx_u64_init_golden(ntt_ctx_u64* ctx, u64 k);
void ntt_u64_cyclic_forward(u64* a_hat, const u64* a, ntt_ctx_u64* ctx);
void ntt_u64_cyclic_inverse(u64* a_hat, const u64* a, ntt_ctx_u64* ctx);
void ntt_u64_cyclic_inverse_montgomery_in(u64* a_hat, const u64* a,
                                          ntt_ctx_u64* ctx);
void ntt_ctx_u64_free(ntt_ctx_u64* ctx);
bool ntt_ctx_u64_init(ntt_ctx_u64* ctx, u64 p, u64 k, u64 omega, u64 psi);
bool ntt_ctx_u64_init_golden(ntt_ctx_u64* ctx, u64 k);

#endif