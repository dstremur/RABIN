#ifndef BIGRABIN_H
#define BIGRABIN_H

#include "bigcore.h"

/**
 * @brief Miller-Rabin test of n to base a (plain modular arithmetic).
 *
 * Let n_l = n->size, measured in 64-bit limbs.
 *
 * Writes n - 1 = d * 2^s with d odd, computes x = a^d mod n, and
 * checks:
 *
 *   - x == 1 or x == n - 1: probably prime,
 *   - otherwise square x up to s - 1 times; if x becomes n - 1 it is
 *     probably prime, if it becomes 1 (a nontrivial square root of 1)
 *     or never reaches n - 1 it is composite.
 *
 * Returns true if n passes the test (probably prime), false if n is
 * even, n <= 1, or a witness for compositeness was found. n == 2 and
 * n == 3 are reported prime.
 *
 * Complexity:
 *   Time: O(s * n_l^2) - one modular exponentiation plus O(s)
 *         squarings
 *   Auxiliary memory: O(n_l) limbs for temporaries
 *   Output memory: O(1)
 *
 * @param[in] n Number to test.
 * @param[in] a Base of the test.
 *
 * @return true  If n passes the test (probably prime).
 * @return false If n is even, n <= 1, or a witness for compositeness
 *               was found.
 */
bool bn_rabin(const bignum* n, const bignum* a);

/**
 * @brief Miller-Rabin test of n to base a (Montgomery domain).
 *
 * Let n_l = n->size, measured in 64-bit limbs.
 *
 * Same test as bn_rabin(), but all modular squarings are Montgomery
 * multiplications, which avoids the division in each reduction. The
 * comparisons are done against one_mont (the Montgomery form of 1) and
 * the Montgomery form of n - 1.
 *
 * Returns true if n passes the test (probably prime), false if n is
 * even, n <= 1, or a witness for compositeness was found. Small
 * single-limb n are handled directly (2 and 3 are prime).
 *
 * Complexity:
 *   Time: O(n_l^3) for the context initialization (see
 *         bn_mont_ctx_init), then O(s * n_l^2) for the test itself
 *   Auxiliary memory: O(n_l) limbs for the context and temporaries
 *   Output memory: O(1)
 *
 * @param[in] n Number to test.
 * @param[in] a Base of the test.
 *
 * @return true  If n passes the test (probably prime).
 * @return false If n is even, n <= 1, or a witness for compositeness
 *               was found.
 */
bool bn_rabin_mont(const bignum* n, const bignum* a);

/**
 * @brief Miller-Rabin test to base a using a caller-supplied Montgomery
 * context.
 *
 * Identical to bn_rabin_mont(), but all multiplications are performed in
 * the Montgomery domain using the ready-made context ctx, so no context
 * has to be built and destroyed per call. This is the version to use when
 * testing many bases against the same modulus n.
 *
 * Complexity:
 *   Time: O(n^2 log n), n = n->size in 64-bit limbs
 *   Auxiliary memory: O(n)
 *   Output memory: O(0)
 *
 * @param[in] n Candidate modulus (odd, > 1); ctx must have been built for n.
 * @param[in] a Base of the test (0 < a < n-1).
 * @param[in] ctx Initialized Montgomery context for modulus n.
 *
 * @return true  If n passes the test (probably prime).
 * @return false If n is even, n <= 1, or a witness for compositeness
 *               was found.
 */
bool bn_rabin_mont_ctx(const bignum* n, const bignum* a,
                       const bn_mont_ctx* ctx);

#endif
