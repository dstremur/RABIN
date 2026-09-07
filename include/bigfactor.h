#ifndef BIGFACTOR_H
#define BIGFACTOR_H

#include "bigcore.h"
#include "bigvector.h"

/**
 * @brief Trial-divide n against the prime table up to the bound g.
 *
 * Returns false if any prime p < g (up to 50000 table entries) divides
 * n, true if no such divisor was found.
 *
 * Complexity:
 *   Time: O(min(50000, pi(g)) * n) - one O(n) modular reduction per
 *         prime
 *   Auxiliary memory: O(1)
 *   Output memory: O(1)
 *
 * @param[in,out] n Number to trial-divide (reduced in place).
 * @param[in]     g Upper bound on the trial primes.
 *
 * @return true  If no prime p < g divides n.
 * @return false If a prime p < g divides n.
 */
bool trialdiv(bignum* n, u64 g);

/**
 * @brief Trial-divide n against the prime table up to the bound g, storing
 * the first divisor found in f.
 *
 * Returns false and sets f to the smallest prime p < g (up to 70000
 * table entries) that divides n. Returns true if no such divisor was
 * found (f is left unchanged).
 *
 * Complexity:
 *   Time: O(min(70000, pi(g)) * n)
 *   Auxiliary memory: O(1)
 *   Output memory: O(1) limbs for f
 *
 * @param[out]    f Result storing the first divisor found.
 * @param[in,out] n Number to trial-divide (reduced in place).
 * @param[in]     g Upper bound on the trial primes.
 *
 * @return true  If no prime p < g divides n (f left unchanged).
 * @return false If a prime p < g divides n (f set to it).
 */
bool trialdiv_factor(bignum* f, bignum* n, u64 g);

/**
 * @brief One attempt at Pollard's rho with Floyd's cycle detection.
 *
 * Let n = n->size, measured in 64-bit limbs.
 *
 * Iterates x -> x^2 + c (mod n) with two pointers (one step and two
 * steps per round) and computes gcd(|x - y|, n) each round. A
 * nontrivial gcd (1 < d < n) is a factor of n.
 *
 * Returns true and stores the factor in f on success. Returns false if
 * the cycle closes without a factor (d == n) or after 1000000
 * iterations.
 *
 * Complexity:
 *   Time: O(n^{1/4} * n^2) expected - about n^{1/4} iterations, each
 *         costing a few n-limb multiplications and a gcd
 *   Auxiliary memory: O(n) limbs for temporaries
 *   Output memory: O(n) limbs for f
 *
 * @param[out] f Result storing the factor found.
 * @param[in]  n Number to factor.
 *
 * @return true  On success (a nontrivial factor is stored in f).
 * @return false If the cycle closes without a factor or after 1000000
 *               iterations.
 */
bool bn_pollard_rho_inner(bignum* f, const bignum* n);

/**
 * @brief Pollard's rho factorization of n.
 *
 * Let n = n->size, measured in 64-bit limbs.
 *
 * Tries bn_pollard_rho_inner() up to 64 times with fresh random
 * constants. Returns true and stores a nontrivial factor in f on
 * success, false if all attempts failed.
 *
 * Complexity:
 *   Time: O(n^{1/4} * n^2) expected per attempt, up to 64 attempts
 *   Auxiliary memory: O(n) limbs
 *   Output memory: O(n) limbs for f
 *
 * @param[out] f Result storing the factor found.
 * @param[in]  n Number to factor.
 *
 * @return true  On success (a nontrivial factor is stored in f).
 * @return false If all 64 attempts failed.
 */
bool bn_pollard_rho(bignum* f, const bignum* n);

/**
 * @brief Completely factorize n, appending the factors to the vector v.
 *
 * Let n = n->size, measured in 64-bit limbs.
 *
 * Strategy, applied recursively:
 *
 *   1. stop at 0 or 1,
 *   2. append 2 or 3 directly,
 *   3. append n if BPSW says it is prime,
 *   4. peel off all factors of 2,
 *   5. trial division up to 50000 (peel off all powers of the factor),
 *   6. Pollard's rho (up to 16 attempts), then Pollard p - 1,
 *   7. if the found factor f is prime (trial division up to 10000),
 *      peel off all powers of f and recurse on the cofactor; otherwise
 *      recurse on both f and the cofactor.
 *
 * n is modified (it is a working copy); v receives the factors in no
 * particular order, with consecutive duplicates suppressed.
 *
 * Complexity:
 *   Time: subexponential in practice; dominated by the Pollard rho /
 *         p - 1 attempts on the hardest composite branch
 *   Auxiliary memory: O(n) limbs for temporaries
 *   Output memory: O(n) limbs for the factor list
 *
 * @param[out]    v Vector receiving the factors.
 * @param[in,out] n Number to factorize (modified in place).
 */
void bn_factorize(bigvector* v, bignum* n);

/**
 * @brief One stage of Pollard's p - 1 method.
 *
 * Let n = n->size, measured in 64-bit limbs.
 *
 * Repeatedly (up to `iterations` times) picks a random base a in
 * [2, n-1] and computes a^M mod n where M is the product of the
 * highest prime powers q^e <= B for all primes q <= B. If p - 1 is
 * B-smooth for some prime factor p of n, then a^M = 1 mod p, so
 * gcd(a^M - 1, n) reveals p.
 *
 * Returns true and stores the factor in f on success, false otherwise.
 *
 * Complexity:
 *   Time: O(iterations * pi(B) * n^2) - one modular exponentiation per
 *         prime <= B per iteration
 *   Auxiliary memory: O(n) limbs for temporaries
 *   Output memory: O(n) limbs for f
 *
 * @param[out]         f          Result storing the factor found.
 * @param[in]          n          Number to factor.
 * @param[in]          B          Smoothness bound.
 * @param[in]          iterations Number of random-base attempts.
 *
 * @return true  On success (a nontrivial factor is stored in f).
 * @return false If no factor was found.
 */
bool bn_pollard_p_minus_one_stage_1(bignum* f, bignum* n, u64 B,
                                    u64 iterations);

/**
 * @brief Pollard's p - 1 factorization of n, with three escalating stages.
 *
 * Let n = n->size, measured in 64-bit limbs.
 *
 * Runs bn_pollard_p_minus_one_stage_1() with increasing smoothness
 * bounds and iteration counts:
 *
 *   1. B = 10000, 1000 rounds
 *   2. B = 10000, 100 rounds
 *   3. B = 70000, 100 rounds
 *
 * Returns true and stores a nontrivial factor in f as soon as any
 * stage succeeds.
 *
 * Complexity:
 *   Time: sum of the three stages, O(pi(B) * n^2) per iteration
 *   Auxiliary memory: O(n) limbs
 *   Output memory: O(n) limbs for f
 *
 * @param[out] f Result storing the factor found.
 * @param[in]  n Number to factor.
 *
 * @return true  On success (a nontrivial factor is stored in f).
 * @return false If all stages failed.
 */
bool bn_pollard_p_minus_one(bignum* f, bignum* n);

#endif
