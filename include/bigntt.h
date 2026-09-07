#ifndef BIGNTT_H
#define BIGNTT_H

/*===========================================================================
 *  bigntt.h
 *
 *  Number-theoretic transforms over NTT-friendly primes: context
 *  construction, Montgomerized forward/inverse transforms, and
 *  multiplicative-generator / Proth-prime construction.
 *
 *  Layout:
 *    - ntt_ctx and lifecycle     (bigntt_ctx_init_simple, bigntt_ctx_init,
 *                                bigntt_ctx_init_golden, bigntt_ctx_free)
 *    - safety / sizing           (bn_check_ntt_safety)
 *    - transforms                (bigntt_cyclic_forward, bigntt_cyclic_inverse,
 *                                bigntt_cyclic_inverse_mont_in)
 *    - generators (factors known)  (bn_find_gen, bn_find_gen_fp,
 *                                bn_find_gen_proth, bn_gen_proth_ntt)
 *===========================================================================*/

#include "bigcore.h"
#include "bigpoly.h"
#include "bigvector.h"

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

/**
 * @brief Find a generator g of the multiplicative group F_p^* given the
 * prime factorization of p - 1.
 *
 * Let n_l = p->size, measured in 64-bit limbs.
 *
 * Repeatedly picks a random a in [2, p-2] and checks that
 * a^((p-1)/r) != 1 mod p for every prime factor r of p - 1; the first
 * candidate that passes all checks is a primitive root, stored in g.
 *
 * Complexity:
 *   Time: O(pi(p-1) * n_l^2) expected - one modular exponentiation per
 *         factor per candidate, O(1) candidates expected
 *   Auxiliary memory: O(n_l) limbs for temporaries
 *   Output memory: O(n_l) limbs
 *
 * @param[out]      g       Result storing the generator (primitive root).
 * @param[in]       p       Prime modulus.
 * @param[in]       factors Vector of prime factors of p - 1.
 */
void bn_find_gen(bignum* g, bignum* p, bigvector* factors);

/**
 * @brief Find a generator g of the multiplicative group F_p^*.
 *
 * Let n_l = p->size, measured in 64-bit limbs.
 *
 * Factorizes p - 1 with bn_factorize() and delegates to
 * bn_find_gen().
 *
 * Complexity:
 *   Time: factorization of p - 1 (subexponential in practice) plus
 *         O(n_l^2) for the generator search
 *   Auxiliary memory: O(n_l) limbs
 *   Output memory: O(n_l) limbs
 *
 * @param[out] g Result storing the generator (primitive root).
 * @param[in]  p Prime modulus.
 */
void bn_find_gen_fp(bignum* g, bignum* p);

/**
 * @brief Find a generator g of F_p^* for a Proth prime p = c * 2^k + 1.
 *
 * Let n_l = p->size, measured in 64-bit limbs.
 *
 * The factorization of p - 1 = c * 2^k is known up to the factorization
 * of c: it is {2} union factorization(c). The generator search then
 * proceeds as in bn_find_gen().
 *
 * Complexity:
 *   Time: factorization of c plus O(n_l^2) for the generator search
 *   Auxiliary memory: O(n_l) limbs
 *   Output memory: O(n_l) limbs
 *
 * @param[out] g Result storing the generator (primitive root).
 * @param[in]  p Proth prime modulus.
 * @param[in]  c Odd multiplier of the Proth prime (p = c * 2^k + 1).
 */
void bn_find_gen_proth(bignum* g, bignum* p, bignum* c);

/**
 * @brief Generate a Proth prime p = c * 2^k + 1 (with c starting at the given
 * odd value) together with a generator g and the roots psi, omega.
 *
 * Let k = exponent of the power of two.
 *
 * Scans odd c upward until c * 2^k + 1 passes BPSW, then factorizes
 * p - 1 = c * 2^k (as {2} union factorization(c)), finds a generator
 * g, and sets:
 *
 *   psi   = g^c mod p      (primitive (k+1)-th root of unity)
 *   omega = psi^2 mod p    (primitive k-th root of unity)
 *
 * Complexity:
 *   Time: O(k^3) expected per candidate for the BPSW test, plus
 *         factorization of c and O(k^2) for the generator search
 *   Auxiliary memory: O(k/64) limbs
 *   Output memory: O(k/64) limbs per output
 *
 * @param[out]     g     Result storing the generator.
 * @param[out]     p     Result storing the Proth prime.
 * @param[out]     omega Result storing the primitive k-th root of unity.
 * @param[out]     psi   Result storing the primitive (k+1)-th root of unity.
 * @param[in]      k     Exponent of the power of two.
 * @param[in]      c     Starting value for the odd multiplier c.
 */
void bn_gen_proth_ntt(bignum* g, bignum* p, bignum* omega, bignum* psi, u64 k,
                      u64 c);

/**
 * @brief Initialize a bignum NTT context by generating a fresh Proth prime.
 *
 * Let k = log2 of the transform length (n = 2^k).
 *
 * Generates a Proth prime p = c' * 2^(k+1) + 1 (starting from the
 * given c) with its generator and roots via bn_gen_proth_ntt(), then
 * delegates to bigntt_ctx_init().
 *
 * Returns true on success, false on allocation failure.
 *
 * Complexity:
 *   Time: O(k^3) expected for the prime generation, plus O(n * k^2)
 *         for the context tables
 *   Auxiliary memory: O(k/64) limbs
 *   Output memory: O(n) bignums per table (4 tables)
 *
 * @param[out] ctx NTT context to initialize.
 * @param[in]    k log2 of the transform length (n = 2^k).
 * @param[in]    c Starting value for the odd multiplier of the Proth prime.
 *
 * @return true  On success.
 * @return false On allocation failure.
 */
bool bigntt_ctx_init_simple(ntt_ctx* ctx, u64 k, u64 c);

/**
 * @brief Initialize a bignum NTT context for the Goldilocks field
 * p = 5 * 2^55 + 1.
 *
 * Let k = log2 of the transform length (n = 2^k, k <= 54).
 *
 * Uses the known primitive root g = 3 and derives:
 *
 *   psi   = g^(5 * 2^(54-k)) mod p   (primitive (k+1)-th root)
 *   omega = psi^2 mod p              (primitive k-th root)
 *
 * then delegates to bigntt_ctx_init().
 *
 * Returns false (and leaves ctx unchanged) if k > 54 or on allocation
 * failure.
 *
 * Complexity:
 *   Time: O(n * k^2) for the context tables, plus O(k^2) for the root
 *         derivation
 *   Auxiliary memory: O(k/64) limbs
 *   Output memory: O(n) bignums per table (4 tables)
 *
 * @param[out] ctx NTT context to initialize.
 * @param[in]    k log2 of the transform length (n = 2^k, k <= 54).
 *
 * @return true  On success.
 * @return false If k > 54 or on allocation failure.
 */
bool bigntt_ctx_init_golden(ntt_ctx* ctx, u64 k);

/**
 * @brief Check that an NTT of the given size cannot wrap around modulo p.
 *
 * Let ntt_size = N (number of coefficients) and bit_width = W.
 *
 * The largest coefficient of the true (unreduced) convolution of two
 * degree < N polynomials with coefficients < 2^W is bounded by
 * N * (2^W - 1)^2. This returns true iff that bound is strictly less
 * than p, i.e. the cyclic convolution modulo p is exact (no
 * wrap-around).
 *
 * Complexity:
 *   Time: O(k^2) where k is the size of p in limbs
 *   Auxiliary memory: O(k) limbs for temporaries
 *   Output memory: O(1)
 *
 * @param[in]  ntt_size Number of coefficients (transform length N).
 * @param[in] bit_width Bit width W of the input coefficients.
 * @param[in]         p Modulus prime.
 *
 * @return true  If the convolution is exact (no wrap-around modulo p).
 * @return false If the convolution may wrap around modulo p.
 */
bool bn_check_ntt_safety(u64 ntt_size, u64 bit_width, const bignum* p);

/**
 * @brief Initialize a bignum NTT context for the prime p and transform length
 * n = 2^k.
 *
 * Precomputes:
 *
 *   - the Montgomery context for p,
 *   - the omega^i and omega^{-i} tables (i = 0 .. n-1),
 *   - the psi^i and psi^{-i} tables (i = 0 .. n-1),
 *   - n_inv = n^{-1} mod p (for the inverse transform scaling),
 *   - the bit-reversal index table.
 *
 * All tables are converted into the Montgomery domain at the end.
 * The caller must supply omega, a primitive k-th root of unity mod p
 * (and psi, a primitive (k+1)-th root; g and c are accepted for
 * interface compatibility).
 *
 * Returns true on success. On failure (allocation failure or a
 * non-invertible root) the context is freed and false is returned.
 *
 * Complexity:
 *   Time: O(n * k^2) - n modular multiplications per table (4 tables)
 *   Auxiliary memory: O(k) limbs for temporaries
 *   Output memory: O(n) bignums per table (4 tables) plus O(n) u64
 *                  indices
 *
 * @param[out]  ctx   NTT context to initialize.
 * @param[in]     p   Prime modulus.
 * @param[in]     g   Generator (accepted for interface compatibility).
 * @param[in] omega Primitive k-th root of unity mod p.
 * @param[in]   psi   Primitive (k+1)-th root of unity mod p.
 * @param[in]     k   log2 of the transform length (n = 2^k).
 * @param[in]     c   Odd multiplier (accepted for interface compatibility).
 *
 * @return true  On success.
 * @return false On allocation failure or a non-invertible root.
 */
bool bigntt_ctx_init(ntt_ctx* ctx, bignum* p, bignum* g, bignum* omega,
                     bignum* psi, u64 k, u64 c);

/**
 * @brief Free all storage held by a bignum NTT context.
 *
 * Frees the four power tables (each entry is a bignum), the
 * bit-reversal table, the modulus, n_inv, and the Montgomery context.
 * Safe to call on a partially initialized context.
 *
 * Complexity:
 *   Time: O(n)
 *   Auxiliary memory: O(1)
 *   Output memory: O(1)
 *
 * @param[in,out] ctx NTT context to free.
 */
void bigntt_ctx_free(ntt_ctx* ctx);

/**
 * @brief Forward cyclic NTT of a bignum polynomial: a_hat = NTT(a).
 *
 * Let n = ctx->n = 2^k and k_l = the size of ctx->q in limbs.
 *
 * Iterative Cooley-Tukey (decimation in time):
 *
 *   1. bit-reversal permutation of the coefficients, converting each
 *      into the Montgomery domain on the fly (coefficients beyond
 *      a->deg are zero),
 *   2. log2(n) stages of butterflies with precomputed twiddle factors
 *      (omega powers, in Montgomery form); additions and subtractions
 *      are reduced with a single conditional add/subtract of q.
 *
 * The output a_hat is in the Montgomery domain with a_hat->deg = n - 1.
 *
 * Complexity:
 *   Time: O(n log n * k_l^2) - each butterfly costs a few k_l-limb
 *         operations
 *   Auxiliary memory: O(n) bignums for the working copy
 *   Output memory: O(n) bignums
 *
 * @param[out] a_hat Result storing the forward NTT (Montgomery domain).
 * @param[in]      a Input polynomial.
 * @param[in]    ctx NTT context.
 */
void bigntt_cyclic_forward(bigpoly* a_hat, bigpoly* a, ntt_ctx* ctx);

/**
 * @brief Inverse cyclic NTT of a bignum polynomial: a_hat = NTT^{-1}(a).
 *
 * Let n = ctx->n = 2^k and k_l = the size of ctx->q in limbs.
 *
 * Thin wrapper around bigntt_cyclic_inverse_mont_in() for inputs in
 * the normal domain. The output contains standard (non-Montgomery)
 * integers with a_hat->deg = n - 1.
 *
 * Complexity:
 *   Time: O(n log n * k_l^2)
 *   Auxiliary memory: O(n) bignums for the working copy
 *   Output memory: O(n) bignums
 *
 * @param[out] a_hat Result storing the inverse NTT (normal domain).
 * @param[in]      a Input polynomial (normal domain).
 * @param[in]    ctx NTT context.
 */
void bigntt_cyclic_inverse(bigpoly* a_hat, bigpoly* a, ntt_ctx* ctx);

/**
 * @brief Inverse cyclic NTT with Montgomery-domain input:
 * a_hat = NTT^{-1}(a), where a is already in the Montgomery domain.
 *
 * Let n = ctx->n = 2^k and k_l = the size of ctx->q in limbs.
 *
 * Same structure as the forward transform, but with the inverse
 * twiddle factors (omega^{-i}), followed by scaling with n^{-1} mod q
 * and conversion out of the Montgomery domain. This variant skips the
 * mont_in of step 1, so the input must already be in the Montgomery
 * domain (e.g. the output of bigntt_cyclic_forward()).
 *
 * The output contains standard (non-Montgomery) integers with
 * a_hat->deg = n - 1.
 *
 * Complexity:
 *   Time: O(n log n * k_l^2)
 *   Auxiliary memory: O(n) bignums for the working copy
 *   Output memory: O(n) bignums
 *
 * @param[out] a_hat Result storing the inverse NTT (normal domain).
 * @param[in]      a Input polynomial (Montgomery domain).
 * @param[in]    ctx NTT context.
 */
void bigntt_cyclic_inverse_mont_in(bigpoly* a_hat, bigpoly* a, ntt_ctx* ctx);

#endif
