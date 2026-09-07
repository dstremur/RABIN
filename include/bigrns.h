#ifndef BIGRNS_H
#define BIGRNS_H

#include "bigcore.h"
#include "bigmatrix.h"
#include "u64.h"

typedef struct {
  u64* primes;
  u64 count;
  bignum* crt_weights;
  u64* garner_weights;
  bignum prod;
  mont_ctx* m_ctxs;
} ctx_rns;

/*
  3. Practical Guideline: The "High-Water Mark"

If you are writing a performance-critical application, follow the High-Water
Mark strategy:

    Check: Does my current ctx_rns->prod have enough bits for this new
calculation? (Use your Hadamard estimator).

    Reuse: If yes, just use the existing context. No extra cost.

    Expand: If no, free the old context and initialize a new one with more
primes.

Cache: Never shrink the context unless you are severely low on memory.

precompute barret reduction or montgomery

*/

typedef struct {
  u64* residues;
  u64 size;
} rns_num;

/**
 * @brief Estimate how many 62-bit primes are needed to represent a bignum.
 *
 * Let k = bit length of a.
 *
 * Returns ceil(k / 62), i.e. the number of 62-bit primes whose product
 * has at least k bits. (The header declares this as
 * bigrns_estimate_primes().)
 *
 * Complexity:
 *   Time: O(n) where n is the size of a in limbs
 *   Auxiliary memory: O(1)
 *   Output memory: O(1)
 *
 * @param[in] a Bignum to estimate the prime count for.
 *
 * @return The number of 62-bit primes needed (ceil(k / 62)).
 */
u64 rns_estimate_primes(const bignum* a);

/**
 * @brief Initialize an RNS context from a list of primes.
 *
 * Let c = count.
 *
 * Stores a copy of the primes, computes their product M, the CRT
 * weights w_i = (M / p_i) * (M / p_i)^{-1} mod p_i (as bignums), one
 * Montgomery context per prime, and the Garner weights
 * g_i = (p_0 * ... * p_{i-1})^{-1} mod p_i.
 *
 * Complexity:
 *   Time: O(c^2 * n) where n is the size of M in limbs (bignum
 *         multiplications and divisions by small primes)
 *   Auxiliary memory: O(n) limbs for temporaries
 *   Output memory: O(c) bignums for the CRT weights, O(c) u64s for
 *         the primes and Garner weights, O(c) mont_ctxs
 *
 * @param[out]   ctx    RNS context to initialize.
 * @param[in] primes Array of primes.
 * @param[in]  count  Number of primes.
 */
void rns_context_init(ctx_rns* ctx, const u64* primes, u64 count);

/**
 * @brief Free all storage owned by an RNS context.
 *
 * Complexity:
 *   Time: O(count)
 *   Auxiliary memory: O(1)
 *   Output memory: O(1)
 *
 * @param[in,out] ctx RNS context to free.
 */
void rns_context_free(ctx_rns* ctx);

/**
 * @brief Convert a bignum to its RNS residue vector: r_i = a mod p_i.
 *
 * Let c = ctx->count.
 *
 * Allocates the residue array on first use, then reduces a modulo
 * each prime.
 *
 * Complexity:
 *   Time: O(c * n) where n is the size of a in limbs
 *   Auxiliary memory: O(1)
 *   Output memory: O(c) u64s
 *
 * @param[out] r   Result residue vector.
 * @param[in]  a   Bignum to convert.
 * @param[in]  ctx RNS context.
 */
void bignum_to_rns(rns_num* r, const bignum* a, ctx_rns* ctx);

/**
 * @brief Component-wise addition of two RNS residue vectors: r = a + b.
 *
 * Let c = r->size.
 *
 * Adds the residues modulo the corresponding prime for each
 * component.
 *
 * Complexity:
 *   Time: O(c)
 *   Auxiliary memory: O(1)
 *   Output memory: O(c) u64s
 *
 * @param[out] r   Result residue vector.
 * @param[in]  a   First residue vector.
 * @param[in]  b   Second residue vector.
 * @param[in]  ctx RNS context.
 */
void rns_add(rns_num* r, const rns_num* a, const rns_num* b,
             const ctx_rns* ctx);

/**
 * @brief Reconstruct a bignum from its RNS residue vector via the CRT.
 *
 * Let c = r->size and n = size of the product M in limbs.
 *
 * Computes a = sum_i r_i * w_i mod M using the precomputed CRT
 * weights. The result is the unique value in [0, M) congruent to the
 * residues.
 *
 * Note: the start/end prints are leftover debug output.
 *
 * Complexity:
 *   Time: O(c * n^2) for the bignum multiplications, plus O(n^2) for
 *         the final reduction
 *   Auxiliary memory: O(n) limbs for temporaries
 *   Output memory: O(n) limbs
 *
 * @param[out] a   Result bignum.
 * @param[in]  r   Residue vector.
 * @param[in]  ctx RNS context.
 */
void rns_to_bignum(bignum* a, const rns_num* v, ctx_rns* ctx);

/**
 * @brief Estimate the number of primes needed for an RNS matrix determinant.
 *
 * Let n = A->r_size.
 *
 * Uses the Hadamard bound on |det(A)|, doubles it to cover the range
 * [-M, M] (so the sign can be recovered from the CRT result), and
 * returns the number of 62-bit primes whose product exceeds that
 * bound.
 *
 * Complexity:
 *   Time: O(n^3 * k^2) where k is the size of the entries in limbs
 *         (Hadamard bound), plus O(n) for the estimate
 *   Auxiliary memory: O(k) limbs
 *   Output memory: O(1)
 *
 * @param[in] A Matrix to estimate the prime count for.
 *
 * @return The number of 62-bit primes needed.
 */
u64 rns_estimate_determinant(bigmatrix* A);

/**
 * @brief Compute the determinant of a bignum matrix via RNS.
 *
 * Let n = A->r_size and c = ctx->count.
 *
 * For each prime p_i, reduces the whole matrix mod p_i and computes
 * the determinant in the u64 Montgomery domain (parallelized over the
 * primes with OpenMP), then recombines the residue vector with the
 * CRT. If the result exceeds M / 2 it is interpreted as negative
 * (two's-complement style) and M is subtracted.
 *
 * The primes must be chosen so that 2 * |det(A)| < M (see
 * rns_estimate_determinant()).
 *
 * Complexity:
 *   Time: O(c * n^3) for the u64 determinants (parallel over c), plus
 *         O(c * n^2) for the reductions and O(c * n_b^2) for the CRT
 *         where n_b is the size of M in limbs
 *   Auxiliary memory: O(n^2) u64s per thread for the reduced matrix
 *   Output memory: O(n_b) limbs
 *
 * @param[out] det Result storing the determinant.
 * @param[in]  A   Matrix.
 * @param[in]  ctx RNS context.
 */
void bigmatrix_det_rns(bignum* det, const bigmatrix* A, const ctx_rns* ctx);

#endif
