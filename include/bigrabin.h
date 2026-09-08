#ifndef BIGRABIN_H
#define BIGRABIN_H

#include "bigcore.h"

/**
 * @brief Miller-Rabin test of \f$n\f$ to base \f$a\f$ (plain modular
 * arithmetic).
 *
 * Let \f$n_l =\f$ n->size, measured in 64-bit limbs.
 *
 * Writes \f$n - 1 = d \cdot 2^s\f$ with \f$d\f$ odd, computes \f$x = a^d \bmod
 * n\f$, and checks:
 *
 *   - \f$x = 1\f$ or \f$x = n - 1\f$: probably prime,
 *   - otherwise square \f$x\f$ up to \f$s - 1\f$ times; if \f$x\f$ becomes \f$n
 * - 1\f$ it is probably prime, if it becomes \f$1\f$ (a nontrivial square root
 * of \f$1\f$) or never reaches \f$n - 1\f$ it is composite.
 *
 * Returns true if \f$n\f$ passes the test (probably prime), false if \f$n\f$ is
 * even, \f$n \le 1\f$, or a witness for compositeness was found. \f$n = 2\f$
 * and
 * \f$n = 3\f$ are reported prime.
 *
 * Complexity:
 *   - Time: \f$O(s \cdot n_l^2)\f$ - one modular exponentiation plus \f$O(s)\f$
 *           squarings
 *   - Auxiliary memory: \f$O(n_l)\f$ limbs for temporaries
 *   - Output memory: \f$O(1)\f$
 *
 * @param[in] n Number to test.
 * @param[in] a Base of the test.
 *
 * @return true  If \f$n\f$ passes the test (probably prime).
 * @return false If \f$n\f$ is even, \f$n \le 1\f$, or a witness for
 * compositeness was found.
 */
bool bn_rabin(const bignum* n, const bignum* a);

/**
 * @brief Miller-Rabin test of \f$n\f$ to base \f$a\f$ (Montgomery domain).
 *
 * Let \f$n_l =\f$ n->size, measured in 64-bit limbs.
 *
 * Same test as bn_rabin(), but all modular squarings are Montgomery
 * multiplications, which avoids the division in each reduction. The
 * comparisons are done against one_mont (the Montgomery form of \f$1\f$) and
 * the Montgomery form of \f$n - 1\f$.
 *
 * Returns true if \f$n\f$ passes the test (probably prime), false if \f$n\f$ is
 * even, \f$n \le 1\f$, or a witness for compositeness was found. Small
 * single-limb \f$n\f$ are handled directly (\f$2\f$ and \f$3\f$ are prime).
 *
 * Complexity:
 *   - Time: \f$O(n_l^3)\f$ for the context initialization (see
 *           bn_mont_ctx_init), then \f$O(s \cdot n_l^2)\f$ for the test itself
 *   - Auxiliary memory: \f$O(n_l)\f$ limbs for the context and temporaries
 *   - Output memory: \f$O(1)\f$
 *
 * @param[in] n Number to test.
 * @param[in] a Base of the test.
 *
 * @return true  If \f$n\f$ passes the test (probably prime).
 * @return false If \f$n\f$ is even, \f$n \le 1\f$, or a witness for
 * compositeness was found.
 */
bool bn_rabin_mont(const bignum* n, const bignum* a);

/**
 * @brief Miller-Rabin test to base \f$a\f$ using a caller-supplied Montgomery
 * context.
 *
 * Identical to bn_rabin_mont(), but all multiplications are performed in
 * the Montgomery domain using the ready-made context ctx, so no context
 * has to be built and destroyed per call. This is the version to use when
 * testing many bases against the same modulus \f$n\f$.
 *
 * Complexity:
 *   - Time: \f$O(n^2 \log n)\f$, \f$n =\f$ n->size in 64-bit limbs
 *   - Auxiliary memory: \f$O(n)\f$
 *   - Output memory: \f$O(0)\f$
 *
 * @param[in] n Candidate modulus (odd, \f$> 1\f$); ctx must have been built for
 * \f$n\f$.
 * @param[in] a Base of the test (\f$0 < a < n - 1\f$).
 * @param[in] ctx Initialized Montgomery context for modulus \f$n\f$.
 *
 * @return true  If \f$n\f$ passes the test (probably prime).
 * @return false If \f$n\f$ is even, \f$n \le 1\f$, or a witness for
 * compositeness was found.
 */
bool bn_rabin_mont_ctx(const bignum* n, const bignum* a,
                       const bn_mont_ctx* ctx);

#endif
