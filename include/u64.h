#ifndef U64_H
#define U64_H

#include <stdio.h>
#include <string.h>

#include "rz.h"

typedef unsigned __int128 u128;

typedef struct {
  u64 p;
  u64 p_inv;
  u64 r2_mod_p;  // 2^128 mod p, R = 2^64
} u64_mont_ctx_t;

typedef struct {
  u64* data;
  u64 rows;
  u64 cols;
  u64 modulus;  // The prime for this specific slice
} u64_mat_t;

/* NTT context. Golden-prime contexts (u64_ntt_ctx_init_golden(),
   u64_ntt_ctx_golden_cached()) store n_inv and all tables in the
   Montgomery domain. Goldilocks contexts (u64_ntt_ctx_init_goldilocks(),
   u64_ntt_ctx_goldilocks_cached()) store them canonically in [0, q)
   and leave bit_rev_indices NULL and mctx unused. */
typedef struct {
  u64 n;
  u64 k;
  u64 q;
  u64 n_inv;  // Montgomery form for golden prime, canonical for Goldilocks
  u64* omega_powers;      // all in the domain of n_inv
  u64* omega_inv_powers;  // all in the domain of n_inv
  u64* bit_rev_indices;   // NULL in Goldilocks contexts (unused)

  u64* twiddle_forward;
  u64* twiddle_inverse;
  u64_mont_ctx_t mctx;  // used by golden-prime contexts only
} u64_ntt_ctx_t;

// u64.c
/**
 * @brief Modular addition: \f$(a + b) \bmod p\f$, with \f$0 \le a, b < p\f$.
 *
 * Adds in u64 and conditionally subtracts p, using a branch-free
 * overflow/carry mask.
 *
 * Complexity:
 *   - Time: \f$O(1)\f$
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(1)\f$
 *
 * @param[in] a First operand.
 * @param[in] b Second operand.
 * @param[in] p Modulus.
 *
 * @return \f$(a + b) \bmod p\f$.
 */
u64 u64_mod_add(u64 a, u64 b, u64 p);

/**
 * @brief Modular subtraction: \f$(a - b) \bmod p\f$, with \f$0 \le a, b < p\f$.
 *
 * Subtracts in u64 and conditionally adds p when a borrow occurred,
 * using a branch-free mask.
 *
 * Complexity:
 *   - Time: \f$O(1)\f$
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(1)\f$
 *
 * @param[in] a First operand.
 * @param[in] b Second operand.
 * @param[in] p Modulus.
 *
 * @return \f$(a - b) \bmod p\f$.
 */
u64 u64_mod_sub(u64 a, u64 b, u64 p);

/**
 * @brief Modular multiplication: \f$(a \cdot b) \bmod p\f$, with \f$0 \le a, b
 * < p\f$.
 *
 * Multiplies in 128 bits and reduces with a hardware 128/64 division.
 *
 * Complexity:
 *   - Time: \f$O(1)\f$
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(1)\f$
 *
 * @param[in] a First operand.
 * @param[in] b Second operand.
 * @param[in] p Modulus.
 *
 * @return \f$(a \cdot b) \bmod p\f$.
 */
u64 u64_mod_mul(u64 a, u64 b, u64 p);

/**
 * @brief Modular exponentiation: \f$base^exp \bmod p\f$.
 *
 * Right-to-left binary exponentiation (square-and-multiply).
 *
 * Complexity:
 *   - Time: \f$O(\log exp)\f$ modular multiplications
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(1)\f$
 *
 * @param[in] base Base.
 * @param[in]  exp Exponent.
 * @param[in]    p Modulus.
 *
 * @return \f$base^exp \bmod p\f$.
 */
u64 u64_mod_pow(u64 base, u64 exp, u64 p);

/**
 * @brief Modular inverse via the extended Euclidean algorithm:
 * \f$a^{-1} \bmod p\f$.
 *
 * \f$p\f$ needs to be prime (more generally, \f$\gcd(a, p)\f$ must be \f$1\f$).
 * Returns \f$0\f$ if \f$a\f$ is \f$0\f$ (should not happen with primes).
 *
 * Complexity:
 *   - Time: \f$O(\log p)\f$
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(1)\f$
 *
 * @param[in] a Value to invert.
 * @param[in] p Modulus (prime).
 *
 * @return \f$a^{-1} \bmod p\f$, or \f$0\f$ if \f$a\f$ is \f$0\f$.
 */
u64 u64_mod_inverse_euclid(u64 a, u64 p);

/**
 * @brief Modular inverse via Fermat's little theorem: \f$a^{-1} = a^{p - 2}
 * \bmod p\f$.
 *
 * \f$p\f$ needs to be prime.
 *
 * Complexity:
 *   - Time: \f$O(\log p)\f$ modular multiplications
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(1)\f$
 *
 * @param[in] n Value to invert.
 * @param[in] p Modulus (prime).
 *
 * @return \f$n^{-1} \bmod p\f$.
 */
u64 u64_mod_inverse(u64 n, u64 p);

/**
 * @brief Compute the Barrett reduction constant \f$\mu = \lfloor (2^{128} - 1)
 * / q \rfloor\f$.
 *
 * Complexity:
 *   - Time: \f$O(1)\f$
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(1)\f$
 *
 * @param[in] q Modulus, \f$q \ge 2\f$.
 *
 * @return The Barrett reduction constant \f$\mu\f$.
 */
u128 u64_compute_mu(u64 q);

/**
 * @brief Barrett reduction: \f$c \bmod q\f$, where \f$\mu = \lfloor (2^{128} -
 * 1) / q \rfloor\f$.
 *
 * Estimates the quotient as \f$q_{est} = \lfloor c \cdot \mu / 2^{128}
 * \rfloor\f$ using the high 128 bits of the 256-bit product, computes
 * \f$r = c - q_{est} \cdot q\f$, and subtracts \f$q\f$ at most once to correct
 * the (rare) underestimate.
 *
 * Complexity:
 *   - Time: \f$O(1)\f$
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(1)\f$
 *
 * @param[in]  c  Value to reduce (128-bit).
 * @param[in]  q  Modulus, \f$q \ge 2\f$.
 * @param[in] mu Barrett reduction constant from u64_compute_mu().
 *
 * @return \f$c \bmod q\f$.
 */
u64 u64_barrett(u128 c, u64 q, u128 mu);

/**
 * @brief Goldilocks field reduction: \f$c \bmod p\f$ with \f$p = 2^{64} -
 * 2^{32}
 * + 1\f$.
 *
 * Folds the 128-bit value using the field relation \f$2^{64} = 2^{32} - 1
 * \bmod p\f$ (hence \f$2^{96} = -1 \bmod p\f$) down to one 64-bit word,
 * followed by a single conditional subtraction.
 *
 * Complexity:
 *   - Time: \f$O(1)\f$
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(1)\f$
 *
 * @param[in] c Value to reduce (128-bit).
 *
 * @return \f$c \bmod p\f$, where \f$p = 2^{64} - 2^{32} + 1\f$.
 * @see u64_mul_goldilocks()
 */
u64 u64_goldilock_red(u128 c);

/**
 * @brief Goldilocks field multiplication: \f$a \cdot b \bmod p\f$ with
 * \f$0 \le a, b < p\f$ and \f$p = 2^{64} - 2^{32} + 1\f$.
 *
 * Multiplies to 128 bits and folds with u64 add/subtract chains only (no
 * division): splitting \f$hi = hi_l + hi_h \cdot 2^{32}\f$, the relation
 * \f$2^{64} = 2^{32} - 1,\ 2^{96} = -1 \bmod p\f$ gives exactly
 * \f$a \cdot b \equiv lo - hi_h + hi_l (2^{32} - 1) \bmod p\f$. The tail
 * is branchless (masked fixups only).
 *
 * Complexity:
 *   - Time: \f$O(1)\f$
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(1)\f$
 *
 * @param[in] a First operand, \f$0 \le a < p\f$.
 * @param[in] b Second operand, \f$0 \le b < p\f$.
 *
 * @return \f$a \cdot b \bmod p\f$.
 * @see u64_goldilock_red(), u64_ntt_goldilocks_forward()
 */
u64 u64_mul_goldilocks(u64 a, u64 b);

// u64_mont.c
/**
 * @brief Initialize a single-limb Montgomery context for the prime \f$p\f$.
 *
 * Computes:
 *
 *   ctx->p        = \f$p\f$
 *   ctx->p_inv    = \f$-p^{-1} \bmod 2^{64}\f$
 *   ctx->r2_mod_p = \f$R^2 \bmod p\f$, with \f$R = 2^{64}\f$
 *
 * Precondition: \f$p\f$ is odd (in practice, prime).
 *
 * Complexity:
 *   - Time: \f$O(1)\f$
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(1)\f$
 *
 * @param[out] ctx Montgomery context to initialize.
 * @param[in]  p   Prime modulus (odd).
 *
 * @par Algorithm Reference:
 * P. L. Montgomery, "Modular Multiplication and Exponentiation on
 * General-Purpose Machines," Software: Practice and Experience, vol. 16,
 * no. 5, 1985.
 * @see u64_mont_redc(), u64_mont_mul()
 */
void u64_mont_init(u64_mont_ctx_t* ctx, u64 p);

/**
 * @brief Single-limb Montgomery reduction (REDC).
 *
 * Given \f$T\f$ with \f$0 \le T < p \cdot 2^{64}\f$, this computes:
 *
 *   \f$t = T \cdot R^{-1} \bmod p\f$,   \f$R = 2^{64}\f$
 *
 * by choosing \f$m = (T \bmod 2^{64}) \cdot (-p^{-1} \bmod 2^{64})\f$ so that
 * \f$T + m \cdot p\f$ is divisible by \f$2^{64}\f$, shifting right by 64 bits,
 * and applying one final conditional subtraction of \f$p\f$.
 *
 * Complexity:
 *   - Time: \f$O(1)\f$
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(1)\f$
 *
 * @param[in] T   Value to reduce (128-bit, \f$0 \le T < p \cdot 2^{64}\f$).
 * @param[in] ctx Montgomery context.
 *
 * @return \f$T \cdot R^{-1} \bmod p\f$.
 */
u64 u64_mont_redc(unsigned __int128 T, const u64_mont_ctx_t* ctx);

/**
 * @brief Montgomery multiplication: \f$(a \cdot b \cdot R^{-1}) \bmod p\f$.
 *
 * If \f$a\f$ and \f$b\f$ are in the Montgomery domain, the result is the
 * Montgomery form of their represented values' product.
 *
 * Complexity:
 *   - Time: \f$O(1)\f$
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(1)\f$
 *
 * @param[in] a   First operand.
 * @param[in] b   Second operand.
 * @param[in] ctx Montgomery context.
 *
 * @return \f$(a \cdot b \cdot R^{-1}) \bmod p\f$.
 */
u64 u64_mont_mul(u64 a, u64 b, const u64_mont_ctx_t* ctx);

/**
 * @brief Convert a value from the normal domain into the Montgomery domain:
 * \f$a_{hat} = a \cdot R \bmod p\f$.
 *
 * Implemented as one Montgomery multiplication by \f$R^2 \bmod p\f$.
 *
 * Complexity:
 *   - Time: \f$O(1)\f$
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(1)\f$
 *
 * @param[in] a   Value in the normal domain.
 * @param[in] ctx Montgomery context.
 *
 * @return \f$a \cdot R \bmod p\f$ (Montgomery form).
 */
u64 u64_mont_in(u64 a, const u64_mont_ctx_t* ctx);

/**
 * @brief Convert a value from the Montgomery domain back to the normal
 * domain: \f$a = a_{hat} \cdot R^{-1} \bmod p\f$.
 *
 * Implemented as one Montgomery multiplication by 1.
 *
 * Complexity:
 *   - Time: \f$O(1)\f$
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(1)\f$
 *
 * @param[in] a_hat Value \f$a_{hat}\f$ in the Montgomery domain.
 * @param[in]   ctx Montgomery context.
 *
 * @return \f$a_{hat} \cdot R^{-1} \bmod p\f$ (normal domain).
 */
u64 u64_mont_out(u64 a_hat, const u64_mont_ctx_t* ctx);

/**
 * @brief Modular inverse of a value given in the Montgomery domain.
 *
 * Converts out of the Montgomery domain, inverts with the extended
 * Euclidean algorithm, and converts the result back in.
 *
 * Complexity:
 *   - Time: \f$O(\log p)\f$
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(1)\f$
 *
 * @param[in] a_mont Value in the Montgomery domain to invert.
 * @param[in]    ctx Montgomery context.
 *
 * @return The modular inverse, in the Montgomery domain.
 */
u64 u64_mont_inverse(u64 a_mont, const u64_mont_ctx_t* ctx);

// u64_ntt.c
/**
 * @brief Forward cyclic NTT: \f$a_{hat} = NTT(a)\f$.
 *
 * Let \f$n =\f$ ctx->n = \f$2^k\f$.
 *
 * Iterative Cooley-Tukey (decimation in time):
 *
 *   1. bit-reversal permutation of the input, converting each entry
 *      into the Montgomery domain on the fly,
 *   2. \f$\log_2(n)\f$ stages of butterflies, each using precomputed twiddle
 *      factors (powers of \f$\omega\f$, already in Montgomery form) and
 *      Montgomery multiplication.
 *
 * The output \f$a_{hat}\f$ is in the Montgomery domain; pair it with
 * u64_ntt_cyclic_inverse() (or the montgomery_in variant) to recover
 * standard integers.
 *
 * \f$a_{hat}\f$ must not alias \f$a\f$.
 *
 * Complexity:
 *   - Time: \f$O(n \log n)\f$
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(n)\f$
 *
 * @param[out] a_hat Result storing the forward NTT (Montgomery domain).
 * @param[in]      a Input array.
 * @param[in]    ctx NTT context.
 *
 * @par Algorithm Reference:
 * C. M. Cooley and J. W. Tukey, "An Algorithm for the Machine
 * Calculation of Complex Fourier Series," Mathematics of Computation,
 * vol. 19, no. 90, 1965.
 * @see u64_ntt_cyclic_inverse(), u64_ntt_cyclic_inverse_montgomery_in()
 */
void u64_ntt_cyclic_forward(u64* a_hat, const u64* a, const u64_ntt_ctx_t* ctx);

/**
 * @brief Initialize an NTT context for the golden NTT prime \f$p = 5 \cdot
 * 2^{55}
 * + 1\f$.
 *
 * Let \f$k =\f$ the logarithm to base \f$2\f$ of the transform length
 * (\f$n = 2^k\f$, \f$k \le 54\f$).
 *
 * Derives the roots of unity from the primitive root \f$g = 3\f$:
 *
 *   \f$\psi   = g^{5 \cdot 2^{54 - k}} \bmod p\f$   (primitive \f$(k + 1)\f$-th
 * root)
 *   \f$\omega = \psi^2 \bmod p\f$              (primitive \f$k\f$-th root)
 *
 * and delegates to u64_ntt_ctx_init().
 *
 * Complexity:
 *   - Time: \f$O(n)\f$ for the precomputed tables, plus \f$O(\log p)\f$ for the
 * root derivation
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(n)\f$ per precomputed table (3 tables)
 *
 * @param[out] ctx NTT context to initialize.
 * @param[in]    k Logarithm to base 2 of the transform length (\f$n = 2^k\f$,
 *                 \f$k \le 54\f$).
 *
 * @return RABIN_SUCCESS on success, or RABIN_ERR_INVALID_ARG if \f$k > 54\f$.
 */
rabin_err_t u64_ntt_ctx_init_golden(u64_ntt_ctx_t* ctx, u64 k);

/**
 * @brief Return the shared golden-prime NTT context for transform length
 * \f$2^k\f$ (\f$k \le 54\f$), initializing it on first use.
 *
 * Returns NULL if \f$k > 54\f$ or initialization fails. The returned context
 * is owned by the library and must not be freed or modified by the
 * caller.
 *
 * Complexity:
 *   - Time: \f$O(1)\f$ after the first call for a given \f$k\f$, \f$O(2^k)\f$
 * the first time
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(2^k)\f$ u64s per table, once per \f$k\f$
 *
 * @param[in] k Logarithm to base 2 of the transform length (\f$n = 2^k\f$,
 *              \f$k \le 54\f$).
 *
 * @return The shared context, or NULL if \f$k > 54\f$ or initialization fails.
 */
u64_ntt_ctx_t* u64_ntt_ctx_golden_cached(u64 k);

/**
 * @brief Initialize an NTT context for the Goldilocks prime \f$p = 2^{64} -
 * 2^{32} + 1\f$ with canonical (non-Montgomery) tables.
 *
 * Let \f$k =\f$ the logarithm to base \f$2\f$ of the transform length
 * (\f$n = 2^k\f$, \f$k \le 31\f$; \f$p - 1\f$ has 2-adicity 32 and table
 * size caps \f$k\f$). The \f$n\f$-th root of unity is derived from the
 * primitive root \f$g = 7\f$ as \f$\omega = g^{(p-1)/n} \bmod p\f$.
 *
 * Montgomery reduction is unusable for this prime (\f$p > 2^{63}\f$
 * breaks the \f$2q\f$ butterfly bounds), so all tables (\f$n_{inv}\f$,
 * \f$\omega\f$ and \f$\omega^{-1}\f$ powers, and the pre-ordered forward
 * and inverse twiddles) are stored as canonical residues in \f$[0, p)\f$
 * and the transforms use u64_mul_goldilocks(). bit_rev_indices is left
 * NULL (the paired transforms need no explicit permutation).
 *
 * Complexity:
 *   - Time: \f$O(n)\f$ for the precomputed tables, plus \f$O(\log p)\f$ for the
 * root derivation
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(n)\f$ per precomputed table (4 tables)
 *
 * @param[out] ctx NTT context to initialize.
 * @param[in]    k Logarithm to base 2 of the transform length (\f$n = 2^k\f$,
 *                 \f$k \le 31\f$).
 *
 * @return RABIN_SUCCESS on success, RABIN_ERR_NULL_PTR if ctx is NULL,
 * RABIN_ERR_INVALID_ARG if \f$k > 31\f$, or RABIN_ERR_OUT_OF_MEMORY.
 * @see u64_ntt_goldilocks_forward(), u64_ntt_goldilocks_inverse()
 */
rabin_err_t u64_ntt_ctx_init_goldilocks(u64_ntt_ctx_t* ctx, u64 k);

/**
 * @brief Return the shared Goldilocks NTT context for transform length
 * \f$2^k\f$ (\f$k \le 31\f$), initializing it on first use.
 *
 * Returns NULL if \f$k > 31\f$ or initialization fails. The returned context
 * is owned by the library and must not be freed or modified by the
 * caller.
 *
 * Complexity:
 *   - Time: \f$O(1)\f$ after the first call for a given \f$k\f$, \f$O(2^k)\f$
 * the first time
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(2^k)\f$ u64s per table, once per \f$k\f$
 *
 * @param[in] k Logarithm to base 2 of the transform length (\f$n = 2^k\f$,
 *              \f$k \le 31\f$).
 *
 * @return The shared context, or NULL if \f$k > 31\f$ or initialization fails.
 */
u64_ntt_ctx_t* u64_ntt_ctx_goldilocks_cached(u64 k);

/**
 * @brief Forward cyclic NTT over the Goldilocks prime \f$p = 2^{64} -
 * 2^{32} + 1\f$: \f$a_{hat} = NTT(a)\f$, canonical domain.
 *
 * Let \f$n =\f$ ctx->n = \f$2^k\f$. Iterative Gentleman-Sande (decimation
 * in frequency) with pre-ordered twiddle tables; all values stay reduced
 * in \f$[0, p)\f$ using fold-based u64 arithmetic (no Montgomery domain).
 * The output is in the same digit-reversed order as the output of
 * u64_ntt_cyclic_forward(), so it pairs with u64_ntt_goldilocks_inverse()
 * without an explicit permutation.
 *
 * \f$a_{hat}\f$ must not alias \f$a\f$. Inputs must be reduced:
 * \f$0 \le a[i] < p\f$.
 *
 * Complexity:
 *   - Time: \f$O(n \log n)\f$
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(n)\f$
 *
 * @param[out] a_hat Result storing the forward NTT (canonical domain).
 * @param[in]      a Input array, entries in \f$[0, p)\f$.
 * @param[in]    ctx Goldilocks NTT context from
 * u64_ntt_ctx_goldilocks_cached().
 *
 * @par Algorithm Reference:
 * W. Gentleman and G. Sande, "Basic Optimizations of Fast Fourier
 * Transform Algorithms," IEEE Trans. Audio Electroacoust., vol. 16,
 * no. 3, 1968.
 * @see u64_ntt_goldilocks_inverse(), u64_ntt_ctx_goldilocks_cached()
 */
void u64_ntt_goldilocks_forward(u64* a_hat, const u64* a,
                                const u64_ntt_ctx_t* ctx);

/**
 * @brief Inverse cyclic NTT over the Goldilocks prime \f$p = 2^{64} -
 * 2^{32} + 1\f$: \f$a_{hat} = NTT^{-1}(a)\f$, canonical domain.
 *
 * Let \f$n =\f$ ctx->n = \f$2^k\f$. Iterative Cooley-Tukey (decimation in
 * time) with pre-ordered inverse twiddles; the final stage is fused with
 * the \f$n^{-1} \bmod p\f$ scaling. All values stay reduced in \f$[0, p)\f$
 * using fold-based u64 arithmetic (no Montgomery domain). Accepts the
 * output order produced by u64_ntt_goldilocks_forward().
 *
 * \f$a_{hat}\f$ must not alias \f$a\f$. Inputs must be reduced:
 * \f$0 \le a[i] < p\f$.
 *
 * Complexity:
 *   - Time: \f$O(n \log n)\f$
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(n)\f$
 *
 * @param[out] a_hat Result storing the inverse NTT (canonical domain).
 * @param[in]      a Input array, entries in \f$[0, p)\f$.
 * @param[in]    ctx Goldilocks NTT context from
 * u64_ntt_ctx_goldilocks_cached().
 *
 * @see u64_ntt_goldilocks_forward(), u64_ntt_ctx_goldilocks_cached()
 */
void u64_ntt_goldilocks_inverse(u64* a_hat, const u64* a,
                                const u64_ntt_ctx_t* ctx);

/**
 * @brief Initialize a u64 NTT context for a prime \f$p\f$ and transform
 * length \f$n = 2^k\f$.
 *
 * Precomputes:
 *
 *   - the Montgomery context for p,
 *   - \f$n_{inv} = n^{-1} \bmod p\f$ in Montgomery form (for the inverse
 *     transform scaling),
 *   - the powers of \f$\omega\f$ and \f$\omega^{-1}\f$ tables (in Montgomery
 * form),
 *   - the bit-reversal index table.
 *
 * The caller must supply \f$\omega\f$, a primitive \f$k\f$-th root of unity
 * \f$\bmod p\f$ (\f$\psi\f$ is accepted for interface compatibility but only
 * \f$\omega\f$ is used for the tables).
 *
 * Complexity:
 *   - Time: \f$O(n)\f$
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(n)\f$ per precomputed table (3 tables)
 *
 * @param[out]  ctx   NTT context to initialize.
 * @param[in]     p   Prime modulus.
 * @param[in]     k   Logarithm to base 2 of the transform length
 *                    (\f$n = 2^k\f$).
 * @param[in] omega Primitive \f$k\f$-th root of unity \f$\bmod p\f$.
 * @param[in]   psi   (Accepted for interface compatibility; unused for
 *                    tables.)
 *
 * @return RABIN_SUCCESS on success, or RABIN_ERR_OUT_OF_MEMORY on allocation
 * failure.
 */
rabin_err_t u64_ntt_ctx_init(u64_ntt_ctx_t* ctx, u64 p, u64 k, u64 omega,
                             u64 psi);

/**
 * @brief Free the precomputed tables of a u64 NTT context.
 *
 * Complexity:
 *   - Time: \f$O(1)\f$
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(1)\f$
 *
 * @param[in,out] ctx NTT context to free.
 */
void u64_ntt_ctx_clear(u64_ntt_ctx_t* ctx);

/**
 * @brief Inverse cyclic NTT: \f$a_{hat} = NTT^{-1}(a)\f$.
 *
 * Let \f$n =\f$ ctx->n = \f$2^k\f$.
 *
 * Same structure as the forward transform, but with the inverse
 * twiddle factors, followed by scaling with \f$n^{-1} \bmod p\f$ and
 * conversion out of the Montgomery domain. The output \f$a_{hat}\f$ contains
 * standard (non-Montgomery) integers.
 *
 * \f$a_{hat}\f$ must not alias \f$a\f$.
 *
 * Complexity:
 *   - Time: \f$O(n \log n)\f$
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(n)\f$
 *
 * @param[out] a_hat Result storing the inverse NTT (normal domain).
 * @param[in]      a Input array.
 * @param[in]    ctx NTT context.
 */
void u64_ntt_cyclic_inverse(u64* a_hat, const u64* a, const u64_ntt_ctx_t* ctx);

/**
 * @brief Inverse cyclic NTT with Montgomery-domain input:
 * \f$a_{hat} = NTT^{-1}(a)\f$, where \f$a\f$ is already in the Montgomery
 * domain.
 *
 * Let \f$n =\f$ ctx->n = \f$2^k\f$.
 *
 * Identical to u64_ntt_cyclic_inverse() except that step 1 performs
 * only the bit-reversal permutation (no u64_mont_in()), so the input must
 * already be in the Montgomery domain. This saves one Montgomery
 * multiplication per entry when chaining forward and inverse
 * transforms. The output contains standard (non-Montgomery) integers.
 *
 * \f$a_{hat}\f$ must not alias \f$a\f$.
 *
 * Complexity:
 *   - Time: \f$O(n \log n)\f$
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(n)\f$
 *
 * @param[out] a_hat Result storing the inverse NTT (normal domain).
 * @param[in]      a Input array (Montgomery domain).
 * @param[in]    ctx NTT context.
 */
void u64_ntt_cyclic_inverse_montgomery_in(u64* a_hat, const u64* a,
                                          const u64_ntt_ctx_t* ctx);

// u64_matrix.c
/**
 * @brief Determinant of a square u64 matrix modulo the prime M->modulus.
 *
 * Let \f$n =\f$ M->rows = M->cols.
 *
 * Performs Gaussian elimination with partial pivoting on a copy of the
 * matrix: for each column a nonzero pivot is searched for (rows are
 * swapped if needed, negating the determinant), the determinant is
 * accumulated as the product of the pivots, and the rows below are
 * eliminated. Returns \f$0\f$ if the matrix is singular.
 *
 * The matrix must be square; a message is printed otherwise (the
 * computation still proceeds with \f$n =\f$ rows).
 *
 * Complexity:
 *   - Time: \f$O(n^3)\f$
 *   - Auxiliary memory: \f$O(n^2)\f$ for the working copy
 *   - Output memory: \f$O(1)\f$
 *
 * @param[in] M Square matrix (entries modulo M->modulus).
 *
 * @return The determinant modulo M->modulus, or \f$0\f$ if singular.
 */
u64 u64_mat_det(u64_mat_t* M);

/**
 * @brief Determinant of a u64 matrix modulo ctx->p, in the Montgomery domain.
 *
 * Let \f$n =\f$ the matrix dimension.
 *
 * Same Gaussian elimination as u64_mat_det(), but all entries are
 * converted to the Montgomery domain up front so that the inner loop
 * uses fast Montgomery multiplication. The row updates are processed
 * in TILE_SIZE-wide blocks to improve cache behaviour.
 *
 * Works in place: the input matrix is destroyed, so pass a copy.
 *
 * Complexity:
 *   - Time: \f$O(n^3)\f$
 *   - Auxiliary memory: \f$O(1)\f$ (in place)
 *   - Output memory: \f$O(1)\f$
 *
 * @param[in,out] mat Matrix (destroyed in place).
 * @param[in]      n  Matrix dimension.
 * @param[in]    ctx  Montgomery context.
 *
 * @return The determinant modulo ctx->p, or \f$0\f$ if singular.
 */
u64 u64_mat_det_optimized(u64* mat, u64 n, const u64_mont_ctx_t* ctx);

#endif
