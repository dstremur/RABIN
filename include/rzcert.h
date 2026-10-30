#ifndef RZCERT_H
#define RZCERT_H

#include <stdlib.h>

#include "rabin_errors.h"
#include "rz.h"

typedef struct rzcert_t rzcert_t;
typedef struct rzcert_elem_t rzcert_elem_t;

struct rzcert_elem_t {
  rz_t q;            // prime
  rz_t alpha_q;      // base for pocklington theorem
  rzcert_t* q_cert;  // certificate for q's primality, if NULL then q is leaf
};

struct rzcert_t {
  rz_t N;               // prime number
  rzcert_elem_t* data;  // pointer to start node of pocklington certificates
  u64 size;
  u64 capacity;
};

/**
 * @brief Print a Pocklington certificate, recursively.
 *
 * Prints N, then for every element its factor q and base alpha_q, and
 * recurses into each child certificate. A certificate with size 0 is
 * printed as a base case.
 *
 * Complexity:
 *   - Time: \f$O(m)\f$, \f$m =\f$ total number of elements over all nested
 *     certificates
 *   - Auxiliary memory: \f$O(d)\f$ recursion stack, \f$d =\f$ certificate
 *     depth
 *   - Output memory: \f$O(0)\f$
 *
 * @param[in] cert Certificate to print.
 *
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR.
 */
rabin_err_t print_pocklington_cert(rzcert_t* cert);

/**
 * @brief Free a Pocklington certificate and all nested child
 * certificates.
 *
 * Frees cert->N, every element's q and alpha_q, each child certificate
 * (recursively), the data array, and the certificate itself.
 *
 * Complexity:
 *   - Time: \f$O(m)\f$, \f$m =\f$ total number of elements over all nested
 *     certificates
 *   - Auxiliary memory: \f$O(d)\f$ recursion stack, \f$d =\f$ certificate
 *     depth
 *   - Output memory: \f$O(0)\f$
 *
 * @param[in] cert Certificate to free.
 *
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR.
 */
rabin_err_t rzcert_clear(rzcert_t* cert);

/**
 * @brief Verify that cert is a valid Pocklington certificate for N.
 *
 * Checks the Pocklington conditions for N:
 *
 *   1. \f$\alpha_q^{N - 1} \equiv 1 \pmod N\f$ for each element
 *   2. \f$q \mid N - 1\f$ and \f$q > \sqrt{N} - 1\f$, the size condition
 *      being checked exactly as \f$(q + 1)^2 > N\f$
 *   3. \f$\gcd\left(\alpha_q^{\frac{N - 1}{q}}, N\right) = 1\f$ for each
 *      element
 *
 * Leaf factors are accepted only if they are provable by a base case;
 * base_case_bound is the size threshold (in bits) up to which a prime
 * factor is accepted without a child certificate.
 *
 * @warning Leaves within base_case_bound are accepted via rz_bpsw(), a
 *          probabilistic test: verification is then a proof modulo the
 *          (conjectured but unproven) correctness of BPSW, not a pure
 *          proof. Set base_case_bound = 0 to require full certificates,
 *          though no generator currently emits leaf-free certificates.
 * @warning A certificate with \f$N \le 1\f$, a missing child for a factor
 *          above base_case_bound, or a child whose proved value differs
 *          from its factor is rejected.
 *
 * Complexity:
 *   - Time: \f$O(m \cdot N^2 \log N)\f$, where \f$m\f$ = number of elements,
 *     \f$N\f$ = size in 64-bit limbs
 *   - Auxiliary memory: \f$O(N)\f$
 *   - Output memory: \f$O(0)\f$
 *
 * @param[in] cert            Certificate to verify.
 * @param[in] base_case_bound Maximum size (in bits) of a prime factor
 *                            accepted without a child certificate.
 *
 * @return true  If all Pocklington conditions hold (N is proven prime).
 * @return false If any condition fails or a leaf factor exceeds
 *               base_case_bound.
 *
 * @par Algorithm Reference:
 * J. B. Pocklington, "On the Factorization of Large Integers,"
 * Proceedings of the Cambridge Philosophical Society, vol. 78, 1975.
 */
bool rzcert_verify(const rzcert_t* cert, u64 base_case_bound);

/**
 * @brief Append a Pocklington element \f$(q, \alpha_q)\f$ to a certificate.
 *
 * Copies q and alpha_q into a new element of cert->data (the array grows
 * by doubling). q_cert is stored without copying: on success its ownership
 * transfers to the certificate and rzcert_clear() will free it; on failure
 * the caller keeps ownership. Pass q_cert = NULL to mark q as a leaf factor.
 *
 * Complexity:
 *   - Time: \f$O(\log c)\f$ amortized per call for \f$c = \f$ capacity
 *   - Auxiliary memory: \f$O(1)\f$ plus at most one array reallocation
 *   - Output memory: \f$O(\log |q|)\f$ limbs for the stored copies
 *
 * @param[in,out] cert    Certificate to extend.
 * @param[in]     q       Prime factor of \f$N - 1\f$ to certify.
 * @param[in]     alpha_q Pocklington witness base for q.
 * @param[in]     q_cert  Certificate for q, or NULL for a leaf factor.
 *
 * @return true on success; false on NULL inputs or allocation failure.
 *
 * @see rzcert_verify(), rzcert_clear()
 */
bool rzcert_add_elem(rzcert_t* cert, const rz_t* q, const rz_t* alpha_q,
                     rzcert_t* q_cert);

/**
 * @brief Return whether cert carries no Pocklington elements (a leaf).
 *
 * A base case has size 0 and is accepted by rzcert_verify() only via a
 * BPSW test when its value fits within base_case_bound.
 *
 * @param[in] cert Certificate to inspect.
 *
 * @return true if cert is non-NULL and has no elements, false otherwise.
 */
bool rzcert_is_base_case(const rzcert_t* cert);

/**
 * @brief Return the depth of the certificate tree.
 *
 * A leaf certificate has depth 1; otherwise 1 plus the maximum depth of
 * its child certificates.
 *
 * @param[in] cert Certificate to measure.
 *
 * @return Depth (number of certificate nodes on the longest root-to-leaf
 *         path), or 0 for NULL.
 */
u64 rzcert_get_depth(const rzcert_t* cert);

/**
 * @brief Return the total number of certificate nodes in the tree.
 *
 * Counts cert itself plus every nested child certificate.
 *
 * @param[in] cert Certificate to count.
 *
 * @return Node count, or 0 for NULL.
 */
u64 rzcert_get_node_count(const rzcert_t* cert);

/**
 * @brief Return the heap footprint of the certificate tree in bytes.
 *
 * Sums the certificate structs, element arrays (at capacity), and all
 * allocated limb buffers of every stored integer and child certificate.
 *
 * @param[in] cert Certificate to measure.
 *
 * @return Heap bytes, or 0 for NULL.
 */
size_t rzcert_get_byte_size(const rzcert_t* cert);

#endif
