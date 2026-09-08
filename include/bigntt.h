#ifndef BIGNTT_H
#define BIGNTT_H

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
 * @brief Find a generator \f$g\f$ of the multiplicative group \f$F_p^*\f$ given
 * the prime factorization of \f$p - 1\f$.
 *
 * Let \f$n_l =\f$ p->size, measured in 64-bit limbs.
 *
 * Repeatedly picks a random \f$a\f$ in \f$[2, p - 2]\f$ and checks that
 * \f$a^{(p - 1)/r} \neq 1 \bmod p\f$ for every prime factor \f$r\f$ of \f$p -
 * 1\f$; the first candidate that passes all checks is a primitive root, stored
 * in \f$g\f$.
 *
 * Complexity:
 *   - Time: \f$O(\pi(p-1) \cdot n_l^2)\f$ expected - one modular exponentiation
 * per factor per candidate, \f$O(1)\f$ candidates expected
 *   - Auxiliary memory: \f$O(n_l)\f$ limbs for temporaries
 *   - Output memory: \f$O(n_l)\f$ limbs
 *
 * @param[out]      g       Result storing the generator (primitive root).
 * @param[in]       p       Prime modulus.
 * @param[in]       factors Vector of prime factors of p - 1.
 */
void bn_find_gen(bignum* g, bignum* p, bigvector* factors);

/**
 * @brief Find a generator \f$g\f$ of the multiplicative group \f$F_p^*\f$.
 *
 * Let \f$n_l =\f$ p->size, measured in 64-bit limbs.
 *
 * Factorizes \f$p - 1\f$ with bn_factorize() and delegates to
 * bn_find_gen().
 *
 * Complexity:
 *   - Time: factorization of p - 1 (subexponential in practice) plus
 *           \f$O(n_l^2)\f$ for the generator search
 *   - Auxiliary memory: \f$O(n_l)\f$ limbs
 *   - Output memory: \f$O(n_l)\f$ limbs
 *
 * @param[out] g Result storing the generator (primitive root).
 * @param[in]  p Prime modulus.
 */
void bn_find_gen_fp(bignum* g, bignum* p);

/**
 * @brief Find a generator \f$g\f$ of \f$F_p^*\f$ for a Proth prime \f$p = c
 * \cdot 2^k + 1\f$.
 *
 * Let \f$n_l =\f$ p->size, measured in 64-bit limbs.
 *
 * The factorization of \f$p - 1 = c \cdot 2^k\f$ is known up to the
 * factorization of \f$c\f$: it is \f$\{2\}\f$ union factorization(\f$c\f$). The
 * generator search then proceeds as in bn_find_gen().
 *
 * Complexity:
 *   - Time: factorization of c plus \f$O(n_l^2)\f$ for the generator search
 *   - Auxiliary memory: \f$O(n_l)\f$ limbs
 *   - Output memory: \f$O(n_l)\f$ limbs
 *
 * @param[out] g Result storing the generator (primitive root).
 * @param[in]  p Proth prime modulus.
 * @param[in]  c Odd multiplier of the Proth prime (\f$p = c \cdot 2^k + 1\f$).
 */
void bn_find_gen_proth(bignum* g, bignum* p, bignum* c);

/**
 * @brief Generate a Proth prime \f$p = c \cdot 2^k + 1\f$ (with \f$c\f$
 * starting at the given odd value) together with a generator \f$g\f$ and the
 * roots \f$\psi\f$, \f$\omega\f$.
 *
 * Let \f$k =\f$ exponent of the power of two.
 *
 * Scans odd \f$c\f$ upward until \f$c \cdot 2^k + 1\f$ passes BPSW, then
 * factorizes
 * \f$p - 1 = c \cdot 2^k\f$ (as \f$\{2\}\f$ union factorization(\f$c\f$)),
 * finds a generator
 * \f$g\f$, and sets:
 *
 *   \f$\psi = g^c \bmod p\f$      (primitive \f$(k + 1)\f$-th root of unity)
 *   \f$\omega = \psi^2 \bmod p\f$    (primitive \f$k\f$-th root of unity)
 *
 * Complexity:
 *   - Time: \f$O(k^3)\f$ expected per candidate for the BPSW test, plus
 *           factorization of c and \f$O(k^2)\f$ for the generator search
 *   - Auxiliary memory: \f$O(k/64)\f$ limbs
 *   - Output memory: \f$O(k/64)\f$ limbs per output
 *
 * @param[out]     g     Result storing the generator.
 * @param[out]     p     Result storing the Proth prime.
 * @param[out]     omega Result storing the primitive \f$k\f$-th root of unity.
 * @param[out]     psi   Result storing the primitive \f$(k + 1)\f$-th root of
 * unity.
 * @param[in]      k     Exponent of the power of two.
 * @param[in]      c     Starting value for the odd multiplier \f$c\f$.
 */
void bn_gen_proth_ntt(bignum* g, bignum* p, bignum* omega, bignum* psi, u64 k,
                      u64 c);

/**
 * @brief Initialize a bignum NTT context by generating a fresh Proth prime.
 *
 * Let \f$k =\f$ log2 of the transform length (\f$n = 2^k\f$).
 *
 * Generates a Proth prime \f$p = c' \cdot 2^{(k + 1)} + 1\f$ (starting from the
 * given \f$c\f$) with its generator and roots via bn_gen_proth_ntt(), then
 * delegates to bigntt_ctx_init().
 *
 * Returns true on success, false on allocation failure.
 *
 * Complexity:
 *   - Time: \f$O(k^3)\f$ expected for the prime generation, plus \f$O(n \cdot
 * k^2)\f$ for the context tables
 *   - Auxiliary memory: \f$O(k/64)\f$ limbs
 *   - Output memory: \f$O(n)\f$ bignums per table (4 tables)
 *
 * @param[out] ctx NTT context to initialize.
 * @param[in]    k log2 of the transform length (\f$n = 2^k\f$).
 * @param[in]    c Starting value for the odd multiplier of the Proth prime.
 *
 * @return true  On success.
 * @return false On allocation failure.
 */
bool bigntt_ctx_init_simple(ntt_ctx* ctx, u64 k, u64 c);

/**
 * @brief Initialize a bignum NTT context for the Goldilocks field
 * \f$p = 5 \cdot 2^{55} + 1\f$.
 *
 * Let \f$k =\f$ log2 of the transform length (\f$n = 2^k\f$, \f$k \le 54\f$).
 *
 * Uses the known primitive root \f$g = 3\f$ and derives:
 *
 *   \f$\psi = g^{5 \cdot 2^{54 - k}} \bmod p\f$   (primitive \f$(k + 1)\f$-th
 * root)
 *   \f$\omega = \psi^2 \bmod p\f$              (primitive \f$k\f$-th root)
 *
 * then delegates to bigntt_ctx_init().
 *
 * Returns false (and leaves ctx unchanged) if \f$k > 54\f$ or on allocation
 * failure.
 *
 * Complexity:
 *   - Time: \f$O(n \cdot k^2)\f$ for the context tables, plus \f$O(k^2)\f$ for
 * the root derivation
 *   - Auxiliary memory: \f$O(k/64)\f$ limbs
 *   - Output memory: \f$O(n)\f$ bignums per table (4 tables)
 *
 * @param[out] ctx NTT context to initialize.
 * @param[in]    k log2 of the transform length (\f$n = 2^k\f$, \f$k \le 54\f$).
 *
 * @return true  On success.
 * @return false If \f$k > 54\f$ or on allocation failure.
 */
bool bigntt_ctx_init_golden(ntt_ctx* ctx, u64 k);

/**
 * @brief Check that an NTT of the given size cannot wrap around modulo p.
 *
 * Let \f$N =\f$ ntt_size (number of coefficients) and \f$W =\f$ bit_width.
 *
 * The largest coefficient of the true (unreduced) convolution of two
 * \f$degree < N\f$ polynomials with coefficients \f$< 2^W\f$ is bounded by
 * \f$N \cdot (2^W - 1)^2\f$. This returns true iff that bound is strictly less
 * than \f$p\f$, i.e. the cyclic convolution modulo \f$p\f$ is exact (no
 * wrap-around).
 *
 * Complexity:
 *   - Time: \f$O(k^2)\f$ where \f$k =\f$ the size of \f$p\f$ in limbs
 *   - Auxiliary memory: \f$O(k)\f$ limbs for temporaries
 *   - Output memory: \f$O(1)\f$
 *
 * @param[in]  ntt_size Number of coefficients (transform length \f$N\f$).
 * @param[in] bit_width Bit width \f$W\f$ of the input coefficients.
 * @param[in]         p Modulus prime.
 *
 * @return true  If the convolution is exact (no wrap-around modulo p).
 * @return false If the convolution may wrap around modulo p.
 */
bool bn_check_ntt_safety(u64 ntt_size, u64 bit_width, const bignum* p);

/**
 * @brief Initialize a bignum NTT context for the prime \f$p\f$ and transform
 * length
 * \f$n = 2^k\f$.
 *
 * Precomputes:
 *
 *   - the Montgomery context for \f$p\f$,
 *   - the \f$\omega^i\f$ and \f$\omega^{-i}\f$ tables (\f$i = 0 \ldots n -
 * 1\f$),
 *   - the \f$\psi^i\f$ and \f$\psi^{-i}\f$ tables (\f$i = 0 \ldots n - 1\f$),
 *   - \f$n_{inv} = n^{-1} \bmod p\f$ (for the inverse transform scaling),
 *   - the bit-reversal index table.
 *
 * All tables are converted into the Montgomery domain at the end.
 * The caller must supply \f$\omega\f$, a primitive \f$k\f$-th root of unity
 * \f$\bmod p\f$ (and \f$\psi\f$, a primitive \f$(k + 1)\f$-th root; \f$g\f$ and
 * \f$c\f$ are accepted for interface compatibility).
 *
 * Returns true on success. On failure (allocation failure or a
 * non-invertible root) the context is freed and false is returned.
 *
 * Complexity:
 *   - Time: \f$O(n \cdot k^2)\f$ - \f$n\f$ modular multiplications per table (4
 * tables)
 *   - Auxiliary memory: \f$O(k)\f$ limbs for temporaries
 *   - Output memory: \f$O(n)\f$ bignums per table (4 tables) plus \f$O(n)\f$
 * u64 indices
 *
 * @param[out]  ctx   NTT context to initialize.
 * @param[in]     p   Prime modulus.
 * @param[in]     g   Generator (accepted for interface compatibility).
 * @param[in] omega Primitive \f$k\f$-th root of unity \f$\bmod p\f$.
 * @param[in]   psi   Primitive \f$(k + 1)\f$-th root of unity \f$\bmod p\f$.
 * @param[in]     k   log2 of the transform length (\f$n = 2^k\f$).
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
 * bit-reversal table, the modulus, \f$n_{inv}\f$, and the Montgomery context.
 * Safe to call on a partially initialized context.
 *
 * Complexity:
 *   - Time: \f$O(n)\f$
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(1)\f$
 *
 * @param[in,out] ctx NTT context to free.
 */
void bigntt_ctx_free(ntt_ctx* ctx);

/**
 * @brief Forward cyclic NTT of a bignum polynomial: \f$a_{hat} = NTT(a)\f$.
 *
 * Let \f$n =\f$ ctx->n = \f$2^k\f$ and \f$k_{l} =\f$ the size of ctx->q in
 * limbs.
 *
 * Iterative Cooley-Tukey (decimation in time):
 *
 *   1. bit-reversal permutation of the coefficients, converting each
 *      into the Montgomery domain on the fly (coefficients beyond
 *      a->deg are zero),
 *   2. \f$\log_2(n)\f$ stages of butterflies with precomputed twiddle factors
 *      (\f$\omega\f$ powers, in Montgomery form); additions and subtractions
 *      are reduced with a single conditional add/subtract of \f$q\f$.
 *
 * The output \f$a_{hat}\f$ is in the Montgomery domain with a_hat->deg = \f$n -
 * 1\f$.
 *
 * Complexity:
 *   - Time: \f$O(n \log n \cdot k_{l}^2)\f$ - each butterfly costs a few
 * \f$k_{l}\f$-limb operations
 *   - Auxiliary memory: \f$O(n)\f$ bignums for the working copy
 *   - Output memory: \f$O(n)\f$ bignums
 *
 * @param[out] a_hat Result storing the forward NTT (Montgomery domain).
 * @param[in]      a Input polynomial.
 * @param[in]    ctx NTT context.
 */
void bigntt_cyclic_forward(bigpoly* a_hat, bigpoly* a, ntt_ctx* ctx);

/**
 * @brief Inverse cyclic NTT of a bignum polynomial: \f$a_{hat} =
 * NTT^{-1}(a)\f$.
 *
 * Let \f$n =\f$ ctx->n = \f$2^k\f$ and \f$k_{l} =\f$ the size of ctx->q in
 * limbs.
 *
 * Thin wrapper around bigntt_cyclic_inverse_mont_in() for inputs in
 * the normal domain. The output contains standard (non-Montgomery)
 * integers with a_hat->deg = \f$n - 1\f$.
 *
 * Complexity:
 *   - Time: \f$O(n \log n \cdot k_{l}^2)\f$
 *   - Auxiliary memory: \f$O(n)\f$ bignums for the working copy
 *   - Output memory: \f$O(n)\f$ bignums
 *
 * @param[out] a_hat Result storing the inverse NTT (normal domain).
 * @param[in]      a Input polynomial (normal domain).
 * @param[in]    ctx NTT context.
 */
void bigntt_cyclic_inverse(bigpoly* a_hat, bigpoly* a, ntt_ctx* ctx);

/**
 * @brief Inverse cyclic NTT with Montgomery-domain input:
 * \f$a_{hat} = NTT^{-1}(a)\f$, where \f$a\f$ is already in the Montgomery
 * domain.
 *
 * Let \f$n =\f$ ctx->n = \f$2^k\f$ and \f$k_{l} =\f$ the size of ctx->q in
 * limbs.
 *
 * Same structure as the forward transform, but with the inverse
 * twiddle factors (\f$\omega^{-i}\f$), followed by scaling with \f$n^{-1} \bmod
 * q\f$ and conversion out of the Montgomery domain. This variant skips the
 * mont_in of step 1, so the input must already be in the Montgomery
 * domain (e.g. the output of bigntt_cyclic_forward()).
 *
 * The output contains standard (non-Montgomery) integers with
 * a_hat->deg = \f$n - 1\f$.
 *
 * Complexity:
 *   - Time: \f$O(n \log n \cdot k_{l}^2)\f$
 *   - Auxiliary memory: \f$O(n)\f$ bignums for the working copy
 *   - Output memory: \f$O(n)\f$ bignums
 *
 * @param[out] a_hat Result storing the inverse NTT (normal domain).
 * @param[in]      a Input polynomial (Montgomery domain).
 * @param[in]    ctx NTT context.
 */
void bigntt_cyclic_inverse_mont_in(bigpoly* a_hat, bigpoly* a, ntt_ctx* ctx);

#endif
