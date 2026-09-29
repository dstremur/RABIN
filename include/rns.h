#ifndef RNS_H
#define RNS_H

#include "rabin_errors.h"
#include "rmat.h"
#include "rz.h"
#include "u64.h"

typedef struct {
  u64* primes;
  u64 count;
  rz_t* crt_weights;
  u64* garner_weights;
  rz_t prod;
  u64_mont_ctx_t* m_ctxs;
} rns_ctx_t;

/*
  3. Practical Guideline: The "High-Water Mark"

If you are writing a performance-critical application, follow the High-Water
Mark strategy:

    Check: Does my current rns_ctx_t->prod have enough bits for this new
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
} rns_num_t;

/**
 * @brief Estimate how many primes are needed to represent a rz_t.
 *
 * Let \f$k =\f$ bit length of \f$a\f$.
 *
 * Returns \f$\lceil k / 59 \rceil\f$, i.e. the number of primes from
 * RNS_PRIMES / RNS_PRIMES2 whose product has at least \f$k\f$ bits. Every
 * prime in those tables is \f$\ge 2^{59}\f$, so the product of
 * \f$\lceil k / 59 \rceil\f$ of them strictly exceeds \f$2^k \ge a\f$;
 * budgeting more bits per prime (e.g. 62, while the leading primes are only
 * 60-bit) lets the product fall short of \f$a\f$ and makes CRT
 * reconstruction wrap silently.
 *
 * Complexity:
 *   - Time: \f$O(n)\f$ where \f$n =\f$ the size of \f$a\f$ in limbs
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(1)\f$
 *
 * @param[in] a Bignum to estimate the prime count for.
 *
 * @return The number of primes needed (\f$\lceil k / 59 \rceil\f$).
 */
u64 rns_estimate_primes(const rz_t* a);

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
 * (rz_t multiplications and divisions by small primes)
 *   - Auxiliary memory: \f$O(n)\f$ limbs for temporaries
 *   - Output memory: \f$O(c)\f$ bignums for the CRT weights, \f$O(c)\f$ u64s
 * for the primes and Garner weights, \f$O(c)\f$ mont_ctxs
 *
 * @param[out] ctx RNS context to initialize.
 * @param[in]  primes Array of primes.
 * @param[in]  count  Number of primes.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR,
 * RABIN_ERR_INVALID_ARG, RABIN_ERR_OVERFLOW, or RABIN_ERR_OUT_OF_MEMORY.
 */
rabin_err_t rns_ctx_init(rns_ctx_t* ctx, const u64* primes, u64 count);

/**
 * @brief Free all storage owned by an RNS context.
 *
 * Complexity:
 *   - Time: \f$O(count)\f$
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(1)\f$
 *
 * @param[in,out] ctx RNS context to free.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR.
 */
rabin_err_t rns_ctx_clear(rns_ctx_t* ctx);

/**
 * @brief Convert a rz_t to its RNS residue vector: \f$r_i = a \bmod p_i\f$.
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
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR, RABIN_ERR_OVERFLOW,
 * or RABIN_ERR_OUT_OF_MEMORY.
 */
rabin_err_t rns_import(rns_num_t* r, const rz_t* a, rns_ctx_t* ctx);

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
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR /
 * RABIN_ERR_MATRIX_DIM.
 */
rabin_err_t rns_add(rns_num_t* r, const rns_num_t* a, const rns_num_t* b,
                    const rns_ctx_t* ctx);

/**
 * @brief Reconstruct a rz_t from its RNS residue vector via the CRT.
 *
 * Let \f$c =\f$ ctx->count and \f$n =\f$ size of the product \f$M\f$ in limbs.
 *
 * Uses the iterative CRT (Garner): starts from \f$x = v_0\f$ and, for each
 * prime \f$p_i\f$, adds \f$u \cdot \prod_{j < i} p_j\f$ where
 * \f$u = (v_i - x) \cdot g_i \bmod p_i\f$ with the precomputed Garner
 * weights \f$g_i\f$. For valid inputs the result is the unique value in
 * \f$[0, M)\f$ congruent to the residues.
 *
 * Complexity:
 *   - Time: \f$O(c \cdot n)\f$ for the rz_t multiplications plus
 * \f$O(c \cdot n)\f$ for the \f$O(1)\f$-size modular reductions
 *   - Auxiliary memory: \f$O(n)\f$ limbs for temporaries
 *   - Output memory: \f$O(n)\f$ limbs
 *
 * @param[out] a   Result rz_t.
 * @param[in]  v   Residue vector.
 * @param[in]  ctx RNS context.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR /
 * RABIN_ERR_OUT_OF_MEMORY.
 *
 * @par Algorithm Reference:
 * D. E. Knuth, "The Art of Computer Programming, Vol. 2: Seminumerical
 * Algorithms," 3rd ed., Addison-Wesley, 1997, Section 4.5.3 (Chinese
 * Remaindering, Garner's algorithm).
 * @see rns_import()
 */
rabin_err_t rns_export(rz_t* a, rns_num_t* v, const rns_ctx_t* ctx);

/**
 * @brief Estimate the number of primes needed for an RNS matrix determinant.
 *
 * Let \f$n =\f$ A->rows.
 *
 * Uses the Hadamard bound on \f$|\det(A)|\f$, doubles it to cover the range
 * \f$[-M, M]\f$ (so the sign can be recovered from the CRT result), and
 * returns the number of primes whose product exceeds that bound.
 *
 * Complexity:
 *   - Time: \f$O(n^3 \cdot k^2)\f$ where \f$k =\f$ the size of the entries in
 * limbs (Hadamard bound), plus \f$O(n)\f$ for the estimate
 *   - Auxiliary memory: \f$O(k)\f$ limbs
 *   - Output memory: \f$O(1)\f$
 *
 * @param[in] A Matrix to estimate the prime count for.
 *
 * @return The number of primes needed.
 */
u64 rns_estimate_det(const rmat_t* A);

/**
 * @brief Compute the determinant of a rz_t matrix via RNS.
 *
 * Let \f$n =\f$ A->rows and \f$c =\f$ ctx->count.
 *
 * For each prime \f$p_i\f$, reduces the whole matrix \f$\bmod p_i\f$ and
 * computes the determinant in the u64 Montgomery domain (parallelized over the
 * primes with OpenMP), then recombines the residue vector with the
 * CRT. If the result exceeds \f$M / 2\f$ it is interpreted as negative
 * (two's-complement style) and \f$M\f$ is subtracted.
 *
 * The primes must be chosen so that \f$2 \cdot |\det(A)| < M\f$ (see
 * rns_estimate_det()).
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
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR,
 * RABIN_ERR_INVALID_ARG, RABIN_ERR_OVERFLOW, RABIN_ERR_OUT_OF_MEMORY, or
 * RABIN_ERR_MATRIX_DIM.
 */
rabin_err_t rmat_det_rns(rz_t* det, const rmat_t* A, const rns_ctx_t* ctx);

#endif
