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
 * printed as a base case (verified via BPSW).
 *
 * Complexity:
 *   Time: O(m), m = total number of elements over all nested
 *         certificates
 *   Auxiliary memory: O(d) recursion stack, d = certificate depth
 *   Output memory: O(0)
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
 *   Time: O(m), m = total number of elements over all nested
 *         certificates
 *   Auxiliary memory: O(d) recursion stack, d = certificate depth
 *   Output memory: O(0)
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
 *   1. alpha_q^(N - 1) = 1 (mod N) for each element
 *   2. q | N - 1 and q > sqrt(N) - 1 (via the accumulated factored part)
 *   3. gcd(alpha_q^((N - 1) / q), N) = 1 for each element
 *
 * Leaf factors are accepted only if they are provable by a base case;
 * base_case_bound is the size threshold up to which a prime factor is
 * accepted without a child certificate.
 *
 * Complexity:
 *   Time: O(m * N^2 log N), m = number of elements, N = N->size in
 *         64-bit limbs
 *   Auxiliary memory: O(N)
 *   Output memory: O(0)
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
