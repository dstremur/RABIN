#ifndef BIGNTT_H
#define BIGNTT_H

#include "bignum.h"
#include "bigpoly.h"

typedef struct {
  u64 k;
  u64 n;         // n = 2^k
  bignum q;      // prime modulus
  bignum n_inv;  // n^-1 mod q

  // Twiddle factors
  bignum* omega_powers;      // w^i mod q for i = 0 to n-1
  bignum* omega_inv_powers;  // w^-i mod q for INTT

  // Optional: 2n-th root of unity for negacyclic convolution
  bignum* psi_powers;      // psi^i mod q
  bignum* psi_inv_powers;  // psi^-i mod q

  // Bit-reversal permutation array to avoid recomputing indices
  u64* bit_rev_indices;

  bn_mont_ctx mctx;  // montgomery context
} ntt_ctx;

bool bigntt_ctx_init_simple(ntt_ctx* ctx, u64 k, u64 c);
bool bigntt_ctx_init(ntt_ctx* ctx, bignum* p, bignum* g, bignum* omega,
                     bignum* psi, u64 k, u64 c);
bool bigntt_ctx_init_golden(ntt_ctx* ctx, u64 k);
bool bigntt_find_prime(bignum* q, u64 n, u64 bits);
bool bigntt_find_generator(bignum* g, const bignum* q);
bool bigntt_compute_roots(bignum* omega, bignum* psi, u64 n, const bignum* q);

void bigntt_ctx_free(ntt_ctx* ctx);

/**
 * @brief forward NTT using Cooley-Tukey butterflies.
 *
 * input in normal order, output in bit-reversed order.
 *
 * @param[out] out_hat Output polynomial in bit-reversed order.
 * @param[in]  in      Input polynomial in normal order.
 * @param[in]  ctx     NTT context.
 */
void bigntt_forward(bigpoly* out_hat, const bigpoly* in, const ntt_ctx* ctx);
void bigntt_cyclic_forward(bigpoly* a_hat, bigpoly* a, ntt_ctx* ctx);
void bigntt_cyclic_inverse(bigpoly* a_hat, bigpoly* a, ntt_ctx* ctx);
void bigntt_cyclic_inverse_mont_in(bigpoly* a_hat, bigpoly* a, ntt_ctx* ctx);

/**
 * @brief inverse NTT using Gentleman-Sande (GS) butterflies.
 *
 * input in bit-reversed order, output in normal order.
 *
 * @param[out] out    Output polynomial in normal order.
 * @param[in]  in_hat Input polynomial in bit-reversed order.
 * @param[in]  ctx    NTT context.
 */
void bigntt_inverse(bigpoly* out, const bigpoly* in_hat, const ntt_ctx* ctx);

void bigntt_pointwise_mul(bigpoly* r_hat, const bigpoly* a_hat,
                          const bigpoly* b_hat, const ntt_ctx* ctx);

/**
 * @brief Performs a cyclic (positive-wrapped) convolution:
 * r = INTT(NTT(a) ◦ NTT(b)).
 *
 * @param[out] r   Result polynomial.
 * @param[in]  a   First polynomial.
 * @param[in]  b   Second polynomial.
 * @param[in]  ctx NTT context.
 */
void bigpoly_conv_cyclic_ntt(bigpoly* r, const bigpoly* a, const bigpoly* b,
                             const ntt_ctx* ctx);

/**
 * @brief Performs a negacyclic (negative-wrapped) convolution using psi.
 *
 * @param[out] r   Result polynomial.
 * @param[in]  a   First polynomial.
 * @param[in]  b   Second polynomial.
 * @param[in]  ctx NTT context.
 */
void bigpoly_conv_negacyclic_ntt(bigpoly* r, const bigpoly* a, const bigpoly* b,
                                 const ntt_ctx* ctx);

// Standard fast polynomial multiplication.

bool bn_check_ntt_safety(u64 ntt_size, u64 bit_width, const bignum* p);
#endif
