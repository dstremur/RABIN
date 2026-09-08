#ifndef BIGCERT_H
#define BIGCERT_H

#include <stdlib.h>

#include "bigcore.h"

typedef struct pocklington_cert pocklington_cert;
typedef struct pocklington_cert_elem pocklington_cert_elem;

struct pocklington_cert_elem {
  bignum q;        // prime
  bignum alpha_q;  // base for pocklington theorem
  pocklington_cert*
      q_cert;  // certificate for q's primality, if NULL then q is leaf
};

struct pocklington_cert {
  bignum N;  // prime number
  pocklington_cert_elem*
      data;  // pointer to start node of pocklington certificates
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
 * @param[in] cert Certificate to print (may be NULL, in which case
 *                 nothing is printed).
 */
void print_pocklington_cert(pocklington_cert* cert);

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
 * @param[in] cert Certificate to free (may be NULL, in which case this
 *                 is a no-op).
 */
void pocklington_cert_free(pocklington_cert* cert);

/**
 * @brief Verify that cert is a valid Pocklington certificate for N.
 *
 * Checks the Pocklington conditions for N:
 *
 *   1. \f$\alpha_q^{N - 1} \equiv 1 \pmod N\f$ for each element
 *   2. \f$q \mid N - 1\f$ and \f$q > \sqrt{N} - 1\f$ 
 *   3. \f$\gcd\left(\alpha_q^{\frac{N - 1}{q}}, N\right) = 1\f$ for each element
 *
 * Leaf factors are accepted only if they are provable by a base case;
 * base_case_bound is the size threshold up to which a prime factor is
 * accepted without a child certificate.
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
bool pocklington_cert_verify(const pocklington_cert* cert, u64 base_case_bound);

// TODO: implement these

bool pocklington_cert_add_elem(pocklington_cert* cert, const bignum* q,
                               const bignum* alpha_q, pocklington_cert* q_cert);

bool pocklington_cert_is_base_case(const pocklington_cert* cert);
u64 pocklington_cert_get_depth(const pocklington_cert* cert);
u64 pocklington_cert_get_node_count(const pocklington_cert* cert);
size_t pocklington_cert_get_byte_size(const pocklington_cert* cert);

#endif
