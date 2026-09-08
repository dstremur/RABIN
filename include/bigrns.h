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
 * Let \f$k =\f$ bit length of \f$a\f$.
 *
 * Returns \f$\lceil k / 62 \rceil\f$, i.e. the number of 62-bit primes whose
 * product has at least \f$k\f$ bits. (The header declares this as
 * bigrns_estimate_primes().)
 *
 * Complexity:
 *   - Time: \f$O(n)\f$ where \f$n =\f$ the size of \f$a\f$ in limbs
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(1)\f$
 *
 * @param[in] a Bignum to estimate the prime count for.
 *
 * @return The number of 62-bit primes needed (\f$\lceil k / 62 \rceil\f$).
 */
u64 rns_estimate_primes(const bignum* a);

/**
 * @brief Initialize an RNS context from a list of primes.
 *
 * Let \f$c =\f$ count.
 *
 * Stores a copy of the primes, computes their product \f$M\f$, the CRT
 * weights \f$w_i = (M / p_i) \cdot (M / p_i)^{-1} \bmod p_i\f$ (as bignums),
 * one Montgomery context per prime, and the Garner weights
 * \f$g_i = (p_0 \cdot \cdots \cdot p_{i - 1})^{-1} \bmod p_i\f$.
 *
 * Complexity:
 *   - Time: \f$O(c^2 \cdot n)\f$ where \f$n =\f$ the size of \f$M\f$ in limbs
 * (bignum multiplications and divisions by small primes)
 *   - Auxiliary memory: \f$O(n)\f$ limbs for temporaries
 *   - Output memory: \f$O(c)\f$ bignums for the CRT weights, \f$O(c)\f$ u64s
 * for the primes and Garner weights, \f$O(c)\f$ mont_ctxs
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
 *   - Time: \f$O(count)\f$
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(1)\f$
 *
 * @param[in,out] ctx RNS context to free.
 */
void rns_context_free(ctx_rns* ctx);

/**
 * @brief Convert a bignum to its RNS residue vector: \f$r_i = a \bmod p_i\f$.
 *
 * Let \f$c =\f$ ctx->count.
 *
 * Allocates the residue array on first use, then reduces \f$a\f$ modulo
 * each prime.
 *
 * Complexity:
 *   - Time: \f$O(c \cdot n)\f$ where \f$n =\f$ the size of \f$a\f$ in limbs
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(c)\f$ u64s
 *
 * @param[out] r   Result residue vector.
 * @param[in]  a   Bignum to convert.
 * @param[in]  ctx RNS context.
 */
void bignum_to_rns(rns_num* r, const bignum* a, ctx_rns* ctx);

/**
 * @brief Component-wise addition of two RNS residue vectors: \f$r = a + b\f$.
 *
 * Let \f$c =\f$ r->size.
 *
 * Adds the residues modulo the corresponding prime for each
 * component.
 *
 * Complexity:
 *   - Time: \f$O(c)\f$
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(c)\f$ u64s
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
 * Let \f$c =\f$ r->size and \f$n =\f$ size of the product \f$M\f$ in limbs.
 *
 * Computes \f$a = \sum_i r_i \cdot w_i \bmod M\f$ using the precomputed CRT
 * weights. The result is the unique value in \f$[0, M)\f$ congruent to the
 * residues.
 *
 * Note: the start/end prints are leftover debug output.
 *
 * Complexity:
 *   - Time: \f$O(c \cdot n^2)\f$ for the bignum multiplications, plus
 * \f$O(n^2)\f$ for the final reduction
 *   - Auxiliary memory: \f$O(n)\f$ limbs for temporaries
 *   - Output memory: \f$O(n)\f$ limbs
 *
 * @param[out] a   Result bignum.
 * @param[in]  r   Residue vector.
 * @param[in]  ctx RNS context.
 */
void rns_to_bignum(bignum* a, const rns_num* v, ctx_rns* ctx);

/**
 * @brief Estimate the number of primes needed for an RNS matrix determinant.
 *
 * Let \f$n =\f$ A->r_size.
 *
 * Uses the Hadamard bound on \f$|\det(A)|\f$, doubles it to cover the range
 * \f$[-M, M]\f$ (so the sign can be recovered from the CRT result), and
 * returns the number of 62-bit primes whose product exceeds that
 * bound.
 *
 * Complexity:
 *   - Time: \f$O(n^3 \cdot k^2)\f$ where \f$k =\f$ the size of the entries in
 * limbs (Hadamard bound), plus \f$O(n)\f$ for the estimate
 *   - Auxiliary memory: \f$O(k)\f$ limbs
 *   - Output memory: \f$O(1)\f$
 *
 * @param[in] A Matrix to estimate the prime count for.
 *
 * @return The number of 62-bit primes needed.
 */
u64 rns_estimate_determinant(bigmatrix* A);

/**
 * @brief Compute the determinant of a bignum matrix via RNS.
 *
 * Let \f$n =\f$ A->r_size and \f$c =\f$ ctx->count.
 *
 * For each prime \f$p_i\f$, reduces the whole matrix \f$\bmod p_i\f$ and
 * computes the determinant in the u64 Montgomery domain (parallelized over the
 * primes with OpenMP), then recombines the residue vector with the
 * CRT. If the result exceeds \f$M / 2\f$ it is interpreted as negative
 * (two's-complement style) and \f$M\f$ is subtracted.
 *
 * The primes must be chosen so that \f$2 \cdot |\det(A)| < M\f$ (see
 * rns_estimate_determinant()).
 *
 * Complexity:
 *   - Time: \f$O(c \cdot n^3)\f$ for the u64 determinants (parallel over
 * \f$c\f$), plus
 *           \f$O(c \cdot n^2)\f$ for the reductions and \f$O(c \cdot
 * n_{b}^2)\f$ for the CRT where \f$n_{b}\f$ is the size of \f$M\f$ in limbs
 *   - Auxiliary memory: \f$O(n^2)\f$ u64s per thread for the reduced matrix
 *   - Output memory: \f$O(n_b)\f$ limbs
 *
 * @param[out] det Result storing the determinant.
 * @param[in]  A   Matrix.
 * @param[in]  ctx RNS context.
 */
void bigmatrix_det_rns(bignum* det, const bigmatrix* A, const ctx_rns* ctx);

#endif
