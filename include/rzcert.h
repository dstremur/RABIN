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
 *   - Time: \f$O(m), m =\f$ total number of elements over all nested
 *         certificates
 *   - Auxiliary memory: \f$O(d)\f$ recursion stack, \f$d =\f$ certificate depth
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
 *   - Time: \f$O(m), m =\f$ total number of elements over all nested
 *         certificates
 *   - Auxiliary memory: \f$O(d)\f$ recursion stack, \f$d =\f$ certificate depth
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
 *   2. \f$q \mid N - 1\f$ and \f$q > \sqrt{N} - 1\f$
 *   3. \f$\gcd\left(\alpha_q^{\frac{N - 1}{q}}, N\right) = 1\f$ for each
 * element
 *
 * Leaf factors are accepted only if they are provable by a base case;
 * base_case_bound is the size threshold up to which a prime factor is
 * accepted without a child certificate.
 *
 * @note TODO: Not yet implemented (placeholder declaration).
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
 */
bool rzcert_verify(const rzcert_t* cert, u64 base_case_bound);

// TODO: Not yet implemented (placeholder declarations).
bool rzcert_add_elem(rzcert_t* cert, const rz_t* q, const rz_t* alpha_q,
                     rzcert_t* q_cert);

bool rzcert_is_base_case(const rzcert_t* cert);
u64 rzcert_get_depth(const rzcert_t* cert);
u64 rzcert_get_node_count(const rzcert_t* cert);
size_t rzcert_get_byte_size(const rzcert_t* cert);

#endif
