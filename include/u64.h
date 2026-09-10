#ifndef BIGU64_H
#define BIGU64_H

#include <stdio.h>
#include <string.h>

#include "bigcore.h"

typedef unsigned __int128 u128;

typedef struct {
  u64 p;
  u64 p_inv;
  u64 r2_mod_p;  // 2^128 mod p, R = 2^64
} mont_ctx;

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
u64 mod_add(u64 a, u64 b, u64 p);

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
u64 mod_sub(u64 a, u64 b, u64 p);

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
u64 mod_mul(u64 a, u64 b, u64 p);

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
u64 mod_pow(u64 base, u64 exp, u64 p);

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
u64 mod_inverse_euclid(u64 a, u64 p);

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
u64 mod_inverse(u64 n, u64 p);

/**
 * @brief Compute the Barrett reduction constant \f$\mu = \lfloor (2^{128} - 1)
 * / q \rfloor\f$.
 *
 * Complexity:
 *   - Time: \f$O(1)\f$
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(1)\f$
 *
 * @param[in] q Modulus.
 *
 * @return The Barrett reduction constant \f$\mu\f$.
 */
u64 compute_mu(u64 q);

/**
 * @brief Barrett reduction: \f$c \bmod q\f$, where \f$\mu = \lfloor (2^{128} -
 * 1) / q \rfloor\f$.
 *
 * Estimates the quotient as \f$q_{est} = \lfloor c \cdot \mu / 2^{128}
 * \rfloor\f$ using the high 128 bits of the 192-bit product, computes
 * \f$r = c - q_{est} \cdot q\f$, and subtracts \f$q\f$ at most twice to correct
 * the (rare) overestimate.
 *
 * Complexity:
 *   - Time: \f$O(1)\f$
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(1)\f$
 *
 * @param[in]  c  Value to reduce (128-bit).
 * @param[in]  q  Modulus.
 * @param[in] mu Barrett reduction constant.
 *
 * @return \f$c \bmod q\f$.
 */
u64 barrett_reduction(u128 c, u64 q, u64 mu);

/**
 * @brief Goldilocks field reduction: \f$c \bmod p\f$ with \f$p = 2^{64} -
 * 2^{32}
 * + 1\f$.
 *
 * Splits the 128-bit value \f$c\f$ into 32-bit chunks and uses the field
 * relation \f$2^{32} = 1 \bmod p\f$ to fold the upper chunks down, followed
 * by a single conditional subtraction.
 *
 * Complexity:
 *   - Time: \f$O(1)\f$
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(1)\f$
 *
 * @param[in] c Value to reduce (128-bit).
 *
 * @return \f$c \bmod p\f$, where \f$p = 2^{64} - 2^{32} + 1\f$.
 */
u64 goldilock_red(u128 c);

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
 */
void mont_init(mont_ctx* ctx, u64 p);

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
u64 mont_redc(unsigned __int128 T, const mont_ctx* ctx);

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
u64 mont_mul(u64 a, u64 b, const mont_ctx* ctx);

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
u64 mont_in(u64 a, const mont_ctx* ctx);

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
u64 mont_out(u64 a_hat, const mont_ctx* ctx);

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
u64 mont_inverse(u64 a_mont, const mont_ctx* ctx);

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
 * ntt_u64_cyclic_inverse() (or the montgomery_in variant) to recover
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
 */
void ntt_u64_cyclic_forward(u64* a_hat, const u64* a, const ntt_ctx_u64* ctx);

/**
 * @brief Initialize an NTT context for the Goldilocks field \f$p = 5 \cdot
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
 * and delegates to ntt_ctx_u64_init().
 *
 * Returns false (and leaves ctx unchanged) if \f$k > 54\f$.
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
 * @return true  On success.
 * @return false If \f$k > 54\f$.
 */
bool ntt_ctx_u64_init_golden(ntt_ctx_u64* ctx, u64 k);

/**
 * @brief Return the shared Goldilocks NTT context for transform length
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
ntt_ctx_u64* ntt_ctx_u64_golden_cached(u64 k);

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
 * @return true  On success.
 * @return false On allocation failure.
 */
bool ntt_ctx_u64_init(ntt_ctx_u64* ctx, u64 p, u64 k, u64 omega, u64 psi);

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
void ntt_ctx_u64_free(ntt_ctx_u64* ctx);

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
void ntt_u64_cyclic_inverse(u64* a_hat, const u64* a, const ntt_ctx_u64* ctx);

/**
 * @brief Inverse cyclic NTT with Montgomery-domain input:
 * \f$a_{hat} = NTT^{-1}(a)\f$, where \f$a\f$ is already in the Montgomery
 * domain.
 *
 * Let \f$n =\f$ ctx->n = \f$2^k\f$.
 *
 * Identical to ntt_u64_cyclic_inverse() except that step 1 performs
 * only the bit-reversal permutation (no mont_in), so the input must
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
void ntt_u64_cyclic_inverse_montgomery_in(u64* a_hat, const u64* a,
                                          const ntt_ctx_u64* ctx);

// u64_matrix.c
/**
 * @brief Determinant of a square u64 matrix modulo the prime M->modulus.
 *
 * Let \f$n =\f$ M->r_size = M->c_size.
 *
 * Performs Gaussian elimination with partial pivoting on a copy of the
 * matrix: for each column a nonzero pivot is searched for (rows are
 * swapped if needed, negating the determinant), the determinant is
 * accumulated as the product of the pivots, and the rows below are
 * eliminated. Returns \f$0\f$ if the matrix is singular.
 *
 * The matrix must be square; a message is printed otherwise (the
 * computation still proceeds with \f$n =\f$ r_size).
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
u64 matrix_u64_det(matrix_u64* M);

/**
 * @brief Determinant of a u64 matrix modulo ctx->p, in the Montgomery domain.
 *
 * Let \f$n =\f$ the matrix dimension.
 *
 * Same Gaussian elimination as matrix_u64_det(), but all entries are
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
u64 matrix_u64_det_optimized(u64* mat, u64 n, const mont_ctx* ctx);

#endif
