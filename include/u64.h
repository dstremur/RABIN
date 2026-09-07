#ifndef BIGU64_H
#define BIGU64_H

/*===========================================================================
 *  bigu64.h
 *
 *  64-bit single-limb number theory: modular arithmetic, Montgomery
 *  reduction, Barrett / Goldilocks reductions, the 64-bit NTT stack,
 *  and u64 matrix determinants.
 *
 *  Layout:
 *    - context structs           (mont_ctx, ntt_ctx_u64, matrix_u64)
 *    - u64.c                     (mod_add, mod_sub, mod_mul, mod_pow,
 *                                mod_inverse_euclid, mod_inverse,
 *                                compute_mu, barrett_reduction,
 *                                goldilock_red)
 *    - u64_mont.c                (mont_init, mont_redc, mont_mul,
 *                                mont_in, mont_out, mont_inverse)
 *    - u64_ntt.c                 (ntt_ctx_u64_init, ntt_ctx_u64_init_golden,
 *                                ntt_ctx_u64_golden_cached, ntt_ctx_u64_free,
 *                                ntt_u64_cyclic_forward,
 *                                ntt_u64_cyclic_inverse,
 *                                ntt_u64_cyclic_inverse_montgomery_in)
 *    - u64_matrix.c              (matrix_u64_det, matrix_u64_det_optimized)
 *===========================================================================*/

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
 * @brief Modular addition: (a + b) mod p, with 0 <= a, b < p.
 *
 * Adds in u64 and conditionally subtracts p, using a branch-free
 * overflow/carry mask.
 *
 * Complexity:
 *   Time: O(1)
 *   Auxiliary memory: O(1)
 *   Output memory: O(1)
 *
 * @param[in] a First operand.
 * @param[in] b Second operand.
 * @param[in] p Modulus.
 *
 * @return (a + b) mod p.
 */
u64 mod_add(u64 a, u64 b, u64 p);

/**
 * @brief Modular subtraction: (a - b) mod p, with 0 <= a, b < p.
 *
 * Subtracts in u64 and conditionally adds p when a borrow occurred,
 * using a branch-free mask.
 *
 * Complexity:
 *   Time: O(1)
 *   Auxiliary memory: O(1)
 *   Output memory: O(1)
 *
 * @param[in] a First operand.
 * @param[in] b Second operand.
 * @param[in] p Modulus.
 *
 * @return (a - b) mod p.
 */
u64 mod_sub(u64 a, u64 b, u64 p);

/**
 * @brief Modular multiplication: (a * b) mod p, with 0 <= a, b < p.
 *
 * Multiplies in 128 bits and reduces with a hardware 128/64 division.
 *
 * Complexity:
 *   Time: O(1)
 *   Auxiliary memory: O(1)
 *   Output memory: O(1)
 *
 * @param[in] a First operand.
 * @param[in] b Second operand.
 * @param[in] p Modulus.
 *
 * @return (a * b) mod p.
 */
u64 mod_mul(u64 a, u64 b, u64 p);

/**
 * @brief Modular exponentiation: base^exp mod p.
 *
 * Right-to-left binary exponentiation (square-and-multiply).
 *
 * Complexity:
 *   Time: O(log exp) modular multiplications
 *   Auxiliary memory: O(1)
 *   Output memory: O(1)
 *
 * @param[in] base Base.
 * @param[in]  exp Exponent.
 * @param[in]    p Modulus.
 *
 * @return base^exp mod p.
 */
u64 mod_pow(u64 base, u64 exp, u64 p);

/**
 * @brief Modular inverse via the extended Euclidean algorithm:
 * a^(-1) mod p.
 *
 * p needs to be prime (more generally, gcd(a, p) must be 1). Returns
 * 0 if a is 0 (should not happen with primes).
 *
 * Complexity:
 *   Time: O(log p)
 *   Auxiliary memory: O(1)
 *   Output memory: O(1)
 *
 * @param[in] a Value to invert.
 * @param[in] p Modulus (prime).
 *
 * @return a^(-1) mod p, or 0 if a is 0.
 */
u64 mod_inverse_euclid(u64 a, u64 p);

/**
 * @brief Modular inverse via Fermat's little theorem: a^(-1) = a^(p-2) mod p.
 *
 * p needs to be prime.
 *
 * Complexity:
 *   Time: O(log p) modular multiplications
 *   Auxiliary memory: O(1)
 *   Output memory: O(1)
 *
 * @param[in] n Value to invert.
 * @param[in] p Modulus (prime).
 *
 * @return n^(-1) mod p.
 */
u64 mod_inverse(u64 n, u64 p);

/**
 * @brief Compute the Barrett reduction constant mu = floor((2^128 - 1) / q).
 *
 * Complexity:
 *   Time: O(1)
 *   Auxiliary memory: O(1)
 *   Output memory: O(1)
 *
 * @param[in] q Modulus.
 *
 * @return The Barrett reduction constant mu.
 */
u64 compute_mu(u64 q);

/**
 * @brief Barrett reduction: c mod q, where mu = floor((2^128 - 1) / q).
 *
 * Estimates the quotient as q_est = floor(c * mu / 2^128) using the
 * high 128 bits of the 192-bit product, computes r = c - q_est * q,
 * and subtracts q at most twice to correct the (rare) overestimate.
 *
 * Complexity:
 *   Time: O(1)
 *   Auxiliary memory: O(1)
 *   Output memory: O(1)
 *
 * @param[in]  c  Value to reduce (128-bit).
 * @param[in]  q  Modulus.
 * @param[in] mu Barrett reduction constant.
 *
 * @return c mod q.
 */
u64 barrett_reduction(u128 c, u64 q, u64 mu);

/**
 * @brief Goldilocks field reduction: c mod p with p = 2^64 - 2^32 + 1.
 *
 * Splits the 128-bit value c into 32-bit chunks and uses the field
 * relation 2^32 = 1 (mod p) to fold the upper chunks down, followed by
 * a single conditional subtraction.
 *
 * Complexity:
 *   Time: O(1)
 *   Auxiliary memory: O(1)
 *   Output memory: O(1)
 *
 * @param[in] c Value to reduce (128-bit).
 *
 * @return c mod p, where p = 2^64 - 2^32 + 1.
 */
u64 goldilock_red(u128 c);

// u64_mont.c
/**
 * @brief Initialize a single-limb Montgomery context for the prime p.
 *
 * Computes:
 *
 *   ctx->p        = p
 *   ctx->p_inv    = -p^{-1} mod 2^64
 *   ctx->r2_mod_p = R^2 mod p, with R = 2^64
 *
 * Precondition: p is odd (in practice, prime).
 *
 * Complexity:
 *   Time: O(1)
 *   Auxiliary memory: O(1)
 *   Output memory: O(1)
 *
 * @param[out] ctx Montgomery context to initialize.
 * @param[in]  p   Prime modulus (odd).
 */
void mont_init(mont_ctx* ctx, u64 p);

/**
 * @brief Single-limb Montgomery reduction (REDC).
 *
 * Given T with 0 <= T < p * 2^64, this computes:
 *
 *   t = T * R^{-1} mod p,   R = 2^64
 *
 * by choosing m = (T mod 2^64) * (-p^{-1} mod 2^64) so that
 * T + m*p is divisible by 2^64, shifting right by 64 bits, and
 * applying one final conditional subtraction of p.
 *
 * Complexity:
 *   Time: O(1)
 *   Auxiliary memory: O(1)
 *   Output memory: O(1)
 *
 * @param[in] T   Value to reduce (128-bit, 0 <= T < p * 2^64).
 * @param[in] ctx Montgomery context.
 *
 * @return T * R^{-1} mod p.
 */
u64 mont_redc(unsigned __int128 T, mont_ctx* ctx);

/**
 * @brief Montgomery multiplication: (a * b * R^{-1}) mod p.
 *
 * If a and b are in the Montgomery domain, the result is the
 * Montgomery form of their represented values' product.
 *
 * Complexity:
 *   Time: O(1)
 *   Auxiliary memory: O(1)
 *   Output memory: O(1)
 *
 * @param[in] a   First operand.
 * @param[in] b   Second operand.
 * @param[in] ctx Montgomery context.
 *
 * @return (a * b * R^{-1}) mod p.
 */
u64 mont_mul(u64 a, u64 b, mont_ctx* ctx);

/**
 * @brief Convert a value from the normal domain into the Montgomery domain:
 * a_hat = a * R mod p.
 *
 * Implemented as one Montgomery multiplication by R^2 mod p.
 *
 * Complexity:
 *   Time: O(1)
 *   Auxiliary memory: O(1)
 *   Output memory: O(1)
 *
 * @param[in] a   Value in the normal domain.
 * @param[in] ctx Montgomery context.
 *
 * @return a * R mod p (Montgomery form).
 */
u64 mont_in(u64 a, mont_ctx* ctx);

/**
 * @brief Convert a value from the Montgomery domain back to the normal
 * domain: a = a_hat * R^{-1} mod p.
 *
 * Implemented as one Montgomery multiplication by 1.
 *
 * Complexity:
 *   Time: O(1)
 *   Auxiliary memory: O(1)
 *   Output memory: O(1)
 *
 * @param[in] a_hat Value in the Montgomery domain.
 * @param[in]   ctx Montgomery context.
 *
 * @return a_hat * R^{-1} mod p (normal domain).
 */
u64 mont_out(u64 a_hat, mont_ctx* ctx);

/**
 * @brief Modular inverse of a value given in the Montgomery domain.
 *
 * Converts out of the Montgomery domain, inverts with the extended
 * Euclidean algorithm, and converts the result back in.
 *
 * Complexity:
 *   Time: O(log p)
 *   Auxiliary memory: O(1)
 *   Output memory: O(1)
 *
 * @param[in] a_mont Value in the Montgomery domain to invert.
 * @param[in]    ctx Montgomery context.
 *
 * @return The modular inverse, in the Montgomery domain.
 */
u64 mont_inverse(u64 a_mont, const mont_ctx* ctx);

// u64_ntt.c
/**
 * @brief Forward cyclic NTT: a_hat = NTT(a).
 *
 * Let n = ctx->n = 2^k.
 *
 * Iterative Cooley-Tukey (decimation in time):
 *
 *   1. bit-reversal permutation of the input, converting each entry
 *      into the Montgomery domain on the fly,
 *   2. log2(n) stages of butterflies, each using precomputed twiddle
 *      factors (omega powers, already in Montgomery form) and
 *      Montgomery multiplication.
 *
 * The output a_hat is in the Montgomery domain; pair it with
 * ntt_u64_cyclic_inverse() (or the montgomery_in variant) to recover
 * standard integers.
 *
 * a_hat must not alias a.
 *
 * Complexity:
 *   Time: O(n log n)
 *   Auxiliary memory: O(1)
 *   Output memory: O(n)
 *
 * @param[out] a_hat Result storing the forward NTT (Montgomery domain).
 * @param[in]      a Input array.
 * @param[in]    ctx NTT context.
 */
void ntt_u64_cyclic_forward(u64* a_hat, const u64* a, ntt_ctx_u64* ctx);

/**
 * @brief Initialize an NTT context for the Goldilocks field p = 5 * 2^55 + 1.
 *
 * Let k = log2 of the transform length (n = 2^k, k <= 54).
 *
 * Derives the roots of unity from the primitive root g = 3:
 *
 *   psi   = g^(5 * 2^(54-k)) mod p   (primitive (k+1)-th root)
 *   omega = psi^2 mod p              (primitive k-th root)
 *
 * and delegates to ntt_ctx_u64_init().
 *
 * Returns false (and leaves ctx unchanged) if k > 54.
 *
 * Complexity:
 *   Time: O(n) for the precomputed tables, plus O(log p) for the root
 *         derivation
 *   Auxiliary memory: O(1)
 *   Output memory: O(n) per precomputed table (3 tables)
 *
 * @param[out] ctx NTT context to initialize.
 * @param[in]    k log2 of the transform length (n = 2^k, k <= 54).
 *
 * @return true  On success.
 * @return false If k > 54.
 */
bool ntt_ctx_u64_init_golden(ntt_ctx_u64* ctx, u64 k);

/**
 * @brief Return the shared Goldilocks NTT context for transform length 2^k
 * (k <= 54), initializing it on first use.
 *
 * Returns NULL if k > 54 or initialization fails. The returned context
 * is owned by the library and must not be freed or modified by the
 * caller.
 *
 * Complexity:
 *   Time: O(1) after the first call for a given k, O(2^k) the first
 *         time
 *   Auxiliary memory: O(1)
 *   Output memory: O(2^k) u64s per table, once per k
 *
 * @param[in] k log2 of the transform length (n = 2^k, k <= 54).
 *
 * @return The shared context, or NULL if k > 54 or initialization fails.
 */
ntt_ctx_u64* ntt_ctx_u64_golden_cached(u64 k);

/**
 * @brief Initialize a u64 NTT context for a prime p and transform length
 * n = 2^k.
 *
 * Precomputes:
 *
 *   - the Montgomery context for p,
 *   - n_inv = n^{-1} mod p in Montgomery form (for the inverse
 *     transform scaling),
 *   - the omega and omega^{-1} power tables (in Montgomery form),
 *   - the bit-reversal index table.
 *
 * The caller must supply omega, a primitive k-th root of unity mod p
 * (psi is accepted for interface compatibility but only omega is used
 * for the tables).
 *
 * Complexity:
 *   Time: O(n)
 *   Auxiliary memory: O(1)
 *   Output memory: O(n) per precomputed table (3 tables)
 *
 * @param[out]  ctx   NTT context to initialize.
 * @param[in]     p   Prime modulus.
 * @param[in]     k   log2 of the transform length (n = 2^k).
 * @param[in] omega Primitive k-th root of unity mod p.
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
 *   Time: O(1)
 *   Auxiliary memory: O(1)
 *   Output memory: O(1)
 *
 * @param[in,out] ctx NTT context to free.
 */
void ntt_ctx_u64_free(ntt_ctx_u64* ctx);

/**
 * @brief Inverse cyclic NTT: a_hat = NTT^{-1}(a).
 *
 * Let n = ctx->n = 2^k.
 *
 * Same structure as the forward transform, but with the inverse
 * twiddle factors, followed by scaling with n^{-1} mod p and
 * conversion out of the Montgomery domain. The output a_hat contains
 * standard (non-Montgomery) integers.
 *
 * a_hat must not alias a.
 *
 * Complexity:
 *   Time: O(n log n)
 *   Auxiliary memory: O(1)
 *   Output memory: O(n)
 *
 * @param[out] a_hat Result storing the inverse NTT (normal domain).
 * @param[in]      a Input array.
 * @param[in]    ctx NTT context.
 */
void ntt_u64_cyclic_inverse(u64* a_hat, const u64* a, ntt_ctx_u64* ctx);

/**
 * @brief Inverse cyclic NTT with Montgomery-domain input:
 * a_hat = NTT^{-1}(a), where a is already in the Montgomery domain.
 *
 * Let n = ctx->n = 2^k.
 *
 * Identical to ntt_u64_cyclic_inverse() except that step 1 performs
 * only the bit-reversal permutation (no mont_in), so the input must
 * already be in the Montgomery domain. This saves one Montgomery
 * multiplication per entry when chaining forward and inverse
 * transforms. The output contains standard (non-Montgomery) integers.
 *
 * a_hat must not alias a.
 *
 * Complexity:
 *   Time: O(n log n)
 *   Auxiliary memory: O(1)
 *   Output memory: O(n)
 *
 * @param[out] a_hat Result storing the inverse NTT (normal domain).
 * @param[in]      a Input array (Montgomery domain).
 * @param[in]    ctx NTT context.
 */
void ntt_u64_cyclic_inverse_montgomery_in(u64* a_hat, const u64* a,
                                          ntt_ctx_u64* ctx);

// u64_matrix.c
/**
 * @brief Determinant of a square u64 matrix modulo the prime M->modulus.
 *
 * Let n = M->r_size = M->c_size.
 *
 * Performs Gaussian elimination with partial pivoting on a copy of the
 * matrix: for each column a nonzero pivot is searched for (rows are
 * swapped if needed, negating the determinant), the determinant is
 * accumulated as the product of the pivots, and the rows below are
 * eliminated. Returns 0 if the matrix is singular.
 *
 * The matrix must be square; a message is printed otherwise (the
 * computation still proceeds with n = r_size).
 *
 * Complexity:
 *   Time: O(n^3)
 *   Auxiliary memory: O(n^2) for the working copy
 *   Output memory: O(1)
 *
 * @param[in] M Square matrix (entries modulo M->modulus).
 *
 * @return The determinant modulo M->modulus, or 0 if singular.
 */
u64 matrix_u64_det(matrix_u64* M);

/**
 * @brief Determinant of a u64 matrix modulo ctx->p, in the Montgomery domain.
 *
 * Let n = matrix dimension.
 *
 * Same Gaussian elimination as matrix_u64_det(), but all entries are
 * converted to the Montgomery domain up front so that the inner loop
 * uses fast Montgomery multiplication. The row updates are processed
 * in TILE_SIZE-wide blocks to improve cache behaviour.
 *
 * Works in place: the input matrix is destroyed, so pass a copy.
 *
 * Complexity:
 *   Time: O(n^3)
 *   Auxiliary memory: O(1) (in place)
 *   Output memory: O(1)
 *
 * @param[in,out] mat Matrix (destroyed in place).
 * @param[in]      n  Matrix dimension.
 * @param[in]    ctx  Montgomery context.
 *
 * @return The determinant modulo ctx->p, or 0 if singular.
 */
u64 matrix_u64_det_optimized(u64* mat, u64 n, const mont_ctx* ctx);

#endif
