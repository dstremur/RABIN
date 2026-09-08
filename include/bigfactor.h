#ifndef BIGFACTOR_H
#define BIGFACTOR_H

#include "bigcore.h"
#include "bigvector.h"

/**
 * @brief Trial-divide n against the prime table up to the bound g.
 *
 * Returns false if any prime \f$p < g\f$ (up to 50000 table entries)
 * divides \f$n\f$, true if no such divisor was found.
 *
 * Complexity:
 *   - Time: \f$O(\min(50000, \pi(g)) \cdot n)\f$ - one \f$O(n)\f$ modular
 * reduction per prime
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(1)\f$
 *
 * @param[in,out] n Number to trial-divide (reduced in place).
 * @param[in]     g Upper bound on the trial primes.
 *
 * @return true  If no prime \f$p < g\f$ divides \f$n\f$.
 * @return false If a prime \f$p < g\f$ divides \f$n\f$.
 */
bool trialdiv(bignum* n, u64 g);

/**
 * @brief Trial-divide \f$n\f$ against the prime table up to the bound \f$g\f$,
 * storing the first divisor found in \f$f\f$.
 *
 * Returns false and sets \f$f\f$ to the smallest prime \f$p < g\f$ (up to 70000
 * table entries) that divides \f$n\f$. Returns true if no such divisor was
 * found (\f$f\f$ is left unchanged).
 *
 * Complexity:
 *   - Time: \f$O(\min(70000, \pi(g)) \cdot n)\f$
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(1)\f$ limbs for f
 *
 * @param[out]    f Result storing the first divisor found.
 * @param[in,out] n Number to trial-divide (reduced in place).
 * @param[in]     g Upper bound on the trial primes.
 *
 * @return true  If no prime \f$p < g\f$ divides \f$n\f$ (\f$f\f$ left
 * unchanged).
 * @return false If a prime \f$p < g\f$ divides \f$n\f$ (\f$f\f$ set to it).
 */
bool trialdiv_factor(bignum* f, bignum* n, u64 g);

/**
 * @brief One attempt at Pollard's rho with Floyd's cycle detection.
 *
 * Let \f$n =\f$ n->size, measured in 64-bit limbs.
 *
 * Iterates \f$x \leftarrow x^2 + c\f$ (\f$\bmod n\f$) with two pointers (one
 * step and two steps per round) and computes \f$\gcd(|x - y|, n)\f$ each round.
 * A nontrivial \f$\gcd\f$ (\f$1 < d < n\f$) is a factor of \f$n\f$.
 *
 * Returns true and stores the factor in \f$f\f$ on success. Returns false if
 * the cycle closes without a factor (\f$d = n\f$) or after 1000000
 * iterations.
 *
 * Complexity:
 *   - Time: \f$O(n^{1/4} \cdot n^2)\f$ expected - about \f$n^{1/4}\f$
 * iterations, each costing a few n-limb multiplications and a \f$\gcd\f$
 *   - Auxiliary memory: \f$O(n)\f$ limbs for temporaries
 *   - Output memory: \f$O(n)\f$ limbs for f
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
 * @brief Pollard's rho factorization of \f$n\f$.
 *
 * Let \f$n =\f$ n->size, measured in 64-bit limbs.
 *
 * Tries bn_pollard_rho_inner() up to 64 times with fresh random
 * constants. Returns true and stores a nontrivial factor in \f$f\f$ on
 * success, false if all attempts failed.
 *
 * Complexity:
 *   - Time: \f$O(n^{1/4} \cdot n^2)\f$ expected per attempt, up to 64 attempts
 *   - Auxiliary memory: \f$O(n)\f$ limbs
 *   - Output memory: \f$O(n)\f$ limbs for f
 *
 * @param[out] f Result storing the factor found.
 * @param[in]  n Number to factor.
 *
 * @return true  On success (a nontrivial factor is stored in f).
 * @return false If all 64 attempts failed.
 */
bool bn_pollard_rho(bignum* f, const bignum* n);

/**
 * @brief Completely factorize \f$n\f$, appending the factors to the vector
 * \f$v\f$.
 *
 * Let \f$n =\f$ n->size, measured in 64-bit limbs.
 *
 * Strategy, applied recursively:
 *
 *   1. stop at \f$0\f$ or \f$1\f$,
 *   2. append \f$2\f$ or \f$3\f$ directly,
 *   3. append \f$n\f$ if BPSW says it is prime,
 *   4. peel off all factors of \f$2\f$,
 *   5. trial division up to 50000 (peel off all powers of the factor),
 *   6. Pollard's rho (up to 16 attempts), then Pollard \f$p - 1\f$,
 *   7. if the found factor \f$f\f$ is prime (trial division up to 10000),
 *      peel off all powers of \f$f\f$ and recurse on the cofactor; otherwise
 *      recurse on both \f$f\f$ and the cofactor.
 *
 * \f$n\f$ is modified (it is a working copy); \f$v\f$ receives the factors in
 * no particular order, with consecutive duplicates suppressed.
 *
 * Complexity:
 *   - Time: subexponential in practice; dominated by the Pollard rho /
 *           \f$p - 1\f$ attempts on the hardest composite branch
 *   - Auxiliary memory: \f$O(n)\f$ limbs for temporaries
 *   - Output memory: \f$O(n)\f$ limbs for the factor list
 *
 * @param[out]    v Vector receiving the factors.
 * @param[in,out] n Number to factorize (modified in place).
 */
void bn_factorize(bigvector* v, bignum* n);

/**
 * @brief One stage of Pollard's \f$p - 1\f$ method.
 *
 * Let \f$n =\f$ n->size, measured in 64-bit limbs.
 *
 * Repeatedly (up to `iterations` times) picks a random base \f$a\f$ in
 * \f$[2, n - 1]\f$ and computes \f$a^M \bmod n\f$ where \f$M\f$ is the product
 * of the highest prime powers \f$q^e \le B\f$ for all primes \f$q \le B\f$. If
 * \f$p - 1\f$ is
 * \f$B\f$-smooth for some prime factor \f$p\f$ of \f$n\f$, then \f$a^M = 1
 * \bmod p\f$, so
 * \f$\gcd(a^M - 1, n)\f$ reveals \f$p\f$.
 *
 * Returns true and stores the factor in \f$f\f$ on success, false otherwise.
 *
 * Complexity:
 *   - Time: \f$O(iterations \cdot \pi(B) \cdot n^2)\f$ - one modular
 * exponentiation per
 *           \f$prime \le B\f$ per iteration
 *   - Auxiliary memory: \f$O(n)\f$ limbs for temporaries
 *   - Output memory: \f$O(n)\f$ limbs for f
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
 * @brief Pollard's \f$p - 1\f$ factorization of \f$n\f$, with three escalating
 * stages.
 *
 * Let \f$n =\f$ n->size, measured in 64-bit limbs.
 *
 * Runs bn_pollard_p_minus_one_stage_1() with increasing smoothness
 * bounds and iteration counts:
 *
 *   1. \f$B = 10000\f$, 1000 rounds
 *   2. \f$B = 10000\f$, 100 rounds
 *   3. \f$B = 70000\f$, 100 rounds
 *
 * Returns true and stores a nontrivial factor in \f$f\f$ as soon as any
 * stage succeeds.
 *
 * Complexity:
 *   - Time: sum of the three stages, \f$O(\pi(B) \cdot n^2)\f$ per iteration
 *   - Auxiliary memory: \f$O(n)\f$ limbs
 *   - Output memory: \f$O(n)\f$ limbs for f
 *
 * @param[out] f Result storing the factor found.
 * @param[in]  n Number to factor.
 *
 * @return true  On success (a nontrivial factor is stored in f).
 * @return false If all stages failed.
 */
bool bn_pollard_p_minus_one(bignum* f, bignum* n);

#endif
