#ifndef BIGPRIME_H
#define BIGPRIME_H

/*===========================================================================
 *  bigprime.h
 *
 *  Prime generation and proving: BPSW / strong Lucas tests, proven
 *  primes with Pocklington certificates, Proth primes, and RNS
 *  prime-table construction.
 *
 *  Layout:
 *    - tests                   (bn_is_perfect_square, bn_stronglucas,
 *                               bn_bpsw, checkLemma1)
 *    - generation              (bn_gen_prime, bn_gen_safe_prime,
 *                               gen_rel_size, bn_gen_proth_primes,
 *                               gen_rns_primes)
 *    - proving                 (bn_provable_prime, bn_provable_prime_inner,
 *                               optimal_table_len,
 *                               gen_provable_primes_arithmetic)
 *===========================================================================*/

#include "bigcert.h"
#include "bigcore.h"

/**
 * @brief Generate a random prime of the given bit length into p.
 *
 * Let k = bits.
 *
 * Repeatedly draws random k-bit candidates from /dev/urandom with the
 * top and bottom bits set, rejects candidates divisible by 3, 5, 7,
 * 11, 13, or 17, trial-divides against the first 1500 primes, and
 * accepts the first candidate that passes BPSW.
 *
 * Returns true on success, false if /dev/urandom cannot be opened or
 * reading from it fails.
 *
 * Complexity:
 *   Time: O(k log k) expected - about k / ln(k) candidates, each
 *         costing trial division plus one BPSW test
 *   Auxiliary memory: O(k/64) limbs
 *   Output memory: O(k/64) limbs
 *
 * @param[out] p    Result storing the generated prime.
 * @param[in]  bits Desired bit length of the prime.
 *
 * @return true  On success.
 * @return false If /dev/urandom cannot be opened or reading fails.
 */
bool bn_gen_prime(bignum* p, int bits);

/**
 * @brief Generates a safe prime number p of a specified bit length.
 *
 * A prime p is safe if p = 2q + 1, where q is also a prime (a Sophie Germain
 * prime). This function repeatedly generates a random prime q of length (bits -
 * 1) until the calculated p passes the BPSW primality test.
 *
 * @param[out] p    Pointer to the bignum structure where the generated safe
 * prime will be stored.
 * @param[in]  bits The desired total bit length of the safe prime p.
 *
 * @return True  If the safe prime was successfully generated.
 * @return False If prime generation failed (e.g., maximum iterations reached or
 * invalid bit size).
 */
bool bn_gen_safe_prime(bignum* p, int bits);

/**
 * @brief Test whether n is a perfect square.
 *
 * Let n_l = n->size, measured in 64-bit limbs.
 *
 * Computes the integer square root with Newton-Raphson iteration
 * (x <- (x + n/x) / 2, seeded at 2^(bits/2)) and checks whether
 * x^2 == n. Zero is considered a square; negative numbers are not.
 *
 * Complexity:
 *   Time: O(n_l^2 log n_l) - O(log n_l) iterations, each dominated by
 *         a division
 *   Auxiliary memory: O(n_l) limbs for temporaries
 *   Output memory: O(1)
 *
 * @param[in] n Number to test.
 *
 * @return true  If n is a perfect square.
 * @return false If n is not a perfect square (or is negative).
 */
bool bn_is_perfect_square(const bignum* n);

/**
 * @brief Strong Lucas test of n with Lucas parameters (P, Q).
 *
 * Let n_l = n->size, measured in 64-bit limbs.
 *
 * Writes n + 1 = d * 2^s with d odd, computes (U_d, V_d, Q^d) mod n,
 * and n passes if:
 *
 *   - U_d == 0 mod n, or
 *   - V_d == 0 mod n, or
 *   - V_{d * 2^r} == 0 mod n for some 1 <= r < s
 *
 * (V is iterated with the doubling identity V_2k = V_k^2 - 2 Q^k.)
 *
 * Returns true if n passes the test (probably prime), false otherwise.
 *
 * Complexity:
 *   Time: O(n_l^3) for the Lucas sequence computation (Montgomery
 *         context init), then O(s * n_l^2) for the V iteration
 *   Auxiliary memory: O(n_l) limbs for temporaries
 *   Output memory: O(1)
 *
 * @param[in] n Number to test.
 * @param[in] P Lucas parameter P.
 * @param[in] Q Lucas parameter Q.
 *
 * @return true  If n passes the test (probably prime).
 * @return false If n fails the test.
 */
bool bn_stronglucas(const bignum* n, bignum* P, bignum* Q);

/**
 * @brief Baillie-PSW primality test.
 *
 * Let n_l = n->size, measured in 64-bit limbs.
 *
 * n passes if all of the following hold:
 *
 *   1. n is not divisible by any prime below 10000 (n itself is
 *      accepted if it is one of them),
 *   2. n passes Miller-Rabin to base 2,
 *   3. n passes the strong Lucas test with parameters (P, Q) =
 *      (1, (1 - D)/4) or (1, (1 + D)/4), where D is the first value
 *      in the sequence 5, -7, 9, -11, ... (Selfridge's method A*)
 *      with Jacobi symbol (D / n) = -1.
 *
 * If no such D is found within 15 rounds, n is rejected (this also
 * catches perfect squares). No known Baillie-PSW pseudoprime exists.
 *
 * Returns true if n is probably prime, false if n is even or a
 * witness for compositeness was found.
 *
 * Complexity:
 *   Time: O(n_l^3) - dominated by the Miller-Rabin and strong Lucas
 *         tests (each paying a Montgomery context initialization)
 *   Auxiliary memory: O(n_l) limbs for temporaries
 *   Output memory: O(1)
 *
 * @param[in] n Number to test.
 *
 * @return true  If n is probably prime.
 * @return false If n is even or a witness for compositeness was found.
 */
bool bn_bpsw(const bignum* n);

/**
 * @brief Lemma 1 check of Maurer's algorithm for the case r = 1.
 *
 * Let n_l = n->size, measured in 64-bit limbs.
 *
 * Given n = 2 R q + 1 with q prime, this checks the primality
 * criterion for a random base a:
 *
 *   X = a^((n-1)/q) mod n
 *
 * n is certified if X != 1, gcd(X - 1, n) == 1, and X^q == 1 mod n.
 *
 * Returns true if the check passes (n is prime), false otherwise.
 *
 * Complexity:
 *   Time: O(n_l^2) - two modular exponentiations
 *   Auxiliary memory: O(n_l) limbs for temporaries
 *   Output memory: O(1)
 *
 * @param[in]      n      Number to certify (n = 2 R q + 1).
 * @param[in]      n_min1 Value n - 1.
 * @param[in]      a      Random base.
 * @param[in]      q      Prime factor (q in n = 2 R q + 1).
 *
 * @return true  If the check passes (n is prime).
 * @return false If the check fails.
 */
bool checkLemma1(bignum* n, bignum* n_min1, bignum* a, bignum* q);

/**
 * @brief Draw a relative size for Maurer's algorithm: 2^u with u uniform in
 * [0, 1), i.e. a value in [1/2, 1] biased toward 1.
 *
 * Complexity:
 *   Time: O(1)
 *   Auxiliary memory: O(1)
 *   Output memory: O(1)
 *
 * @return The relative size, a value in [1/2, 1].
 */
double gen_rel_size();

/**
 * @brief Generate a provable prime of k bits into p (Maurer's algorithm).
 *
 * Repeatedly invokes bn_provable_prime_inner() until it succeeds; each
 * failed attempt unwinds its recursion and frees all memory, so the
 * loop is safe to retry indefinitely.
 *
 * Complexity:
 *   Time: O(k^3) expected per successful attempt (recursive provable
 *         prime generation plus trial division and Lemma 1 checks)
 *   Auxiliary memory: O(k^2 / 64) limbs on the recursion stack
 *   Output memory: O(k/64) limbs
 *
 * @param[out] p Result storing the provable prime.
 * @param[in]  k Desired bit length of the prime.
 */
void bn_provable_prime(bignum* p, u64 k);

/**
 * @brief One attempt at Maurer's simpler algorithm for a provable k-bit prime.
 *
 * Let k = bit length of the target prime.
 *
 * Strategy:
 *
 *   - base case k <= 20: accept a random k-bit candidate that passes
 *     BPSW (treated as proven at this size),
 *   - otherwise: pick a relative size rel_size in [1/2, 1] with
 *     rel_size * k < k - k/6, recursively generate a provable prime q
 *     of that size, then search for n = 2 * R q + 1 with R in
 *     [2^(k-1)/(2q), 2^k/(2q)] that survives trial division up to
 *     0.1 k^2 + 1 and passes the Lemma 1 check for some random base a
 *     (up to 200 bases).
 *
 * Returns true and stores the prime in p on success. Returns false
 * (having freed all temporaries) if the recursion fails or the search
 * exceeds 1000 candidates.
 *
 * Complexity:
 *   Time: O(k^3) expected - dominated by the recursive call and the
 *         Lemma 1 modular exponentiations
 *   Auxiliary memory: O(k^2 / 64) limbs on the recursion stack
 *   Output memory: O(k/64) limbs
 *
 * @param[out] p Result storing the provable prime.
 * @param[in]  k Desired bit length of the prime.
 *
 * @return true  On success (the prime is stored in p).
 * @return false If the recursion fails or the search exceeds 1000
 *               candidates.
 */
bool bn_provable_prime_inner(bignum* p, u64 k);

/**
 * @brief Generate `count` Proth primes of the form p = c * 2^k + 1 and print
 * them as a C array initializer.
 *
 * Let k = exponent of the power of two.
 *
 * Scans odd c starting from the given c (rounded up to odd), testing
 * each candidate c * 2^k + 1 with BPSW, until `count` primes are
 * found or the 64-bit search space is exhausted. The output is a
 * `static const uint64_t RNS_PRIMES[]` table suitable for pasting
 * into primes.h.
 *
 * Complexity:
 *   Time: O(count * k^3) expected - one BPSW test per candidate
 *   Auxiliary memory: O(k/64) limbs
 *   Output memory: O(count) printed entries
 *
 * @param[in] count Number of Proth primes to generate.
 * @param[in]     k Exponent of the power of two.
 * @param[in]     c Starting value for the odd multiplier c.
 */
void bn_gen_proth_primes(u64 count, u64 k, u64 c);

/**
 * @brief Generate `count` primes just below 2^64 and print them as a C array
 * initializer.
 *
 * Scans downward from 2^64 - 47 in steps of 2 (so all candidates are
 * of the form 2^64 - 1 - 2i, i.e. close to the top of the 64-bit
 * range), testing each candidate with BPSW until `count` primes are
 * found. The output is a `static const u64 RNS_PRIMES[]` table
 * suitable for pasting into primes.h.
 *
 * Complexity:
 *   Time: O(count * 1024^3) expected - one BPSW test per candidate
 *   Auxiliary memory: O(1) limbs (single-limb candidates)
 *   Output memory: O(count) printed entries
 *
 * @param[in] count Number of primes to generate.
 */
void gen_rns_primes(u64 count);

/**
 * @brief Heuristic table length for the arithmetic-progression sieve used
 * by Pocklington prime generation.
 *
 * Returns the sieve length s (number of candidates N_0 + i*a, i in
 * [0, s)) for a target prime of n bits where each candidate is processed
 * in B-bit words. s comes from the heuristic
 *
 *   s = (0.4 * n * m) / log(n^2 * m),   m = n / B
 *
 * which balances the cost of the small-prime trial table against the
 * expected number of progression steps needed to find a candidate.
 *
 * Complexity:
 *   Time: O(1)
 *   Auxiliary memory: O(1)
 *   Output memory: O(0)
 *
 * @param[in] n Target prime size in bits.
 * @param[in] B Word size in bits used per candidate (e.g. 64).
 *
 * @return The heuristic sieve length s.
 */
u64 optimal_table_len(u64 n, u64 B);

/**
 * @brief  Generates a provable prime of approximately @p n bits.
 *
 * @details
 * This function recursively constructs a prime number @p p with roughly @p n
 * bits using a Pocklington-style primality proof.
 *
 * For small values of @p n, the function delegates directly to bn_gen_prime().
 *
 * For larger values of @p n, the function performs the following steps:
 *
 *  1. Recursively generates a smaller provable prime @c F with approximately
 *     @p n / 2 bits.
 *  2. Chooses a random integer @c t such that:
 *
 *         2^(n-2) / F < t < 2^(n-1) / F - s*n
 *
 *     where @c s is the sieve/table length returned by optimal_table_len().
 *
 *  3. Constructs an arithmetic progression of candidate integers:
 *
 *         N_i = N_0 + i*a
 *
 *     where:
 *
 *         a  = 2F
 *         N_0 = t*a + 1
 *         0 <= i <= s
 *
 *  4. Sieves candidates using small primes up to a bound proportional to @p n.
 *
 *  5. Applies a Rabin-Miller test, currently with base 2, as a fast
 *     compositeness filter.
 *
 *  6. For candidates that pass the probable-prime test, attempts to prove
 *     primality using Pocklington's lemma with the known large factor @c F
 *     of @c N-1.
 *
 * The function returns by storing the first verified prime found in @p p.
 *
 * @param[out] p  Destination bignum receiving the generated provable prime.
 *                The caller is responsible for managing its lifetime according
 *                to the conventions of the bignum library.
 *
 * @param[out] cert_out Destination pocklington_cert to store the certificate.
 *                      If set to NULL, function will create/store a
 * certificate.
 *
 * @param[in]  n  Desired approximate bit length of the output prime.
 *
 * @pre The bignum library must be initialized.
 *
 * @pre The global small-prime table used by this function must be valid and
 *      must contain enough primes to support the trial-division bound used
 *      internally.
 *
 * @pre The random-number subsystem must be initialized if the internal
 *      bn_gen_random_range() function depends on it.
 *
 * @post On successful completion, @p p contains a prime number of approximately
 *       @p n bits.
 *
 * @note This function is recursive. Its stack usage and runtime grow with @p n.
 *
 * @note The primality proof relies on @c F satisfying the Pocklington condition
 *       @c F > sqrt(N - 1), or an equivalent sufficient condition. If this
 *       condition is not guaranteed by the caller or by the size bounds,
 *       the generated number may be probable-prime but not formally proven
 *       prime by this routine.
 *
 * @warning The current implementation may use variable-length array allocations
 *          and repeated temporary bignum allocations. For large @p n, this may
 *          lead to high stack usage or degraded performance.
 *
 * @warning If no suitable prime is found in the generated arithmetic
 *          progression, the function retries with a new random @c t. It has
 *          no explicit iteration limit and may therefore run for an unbounded
 *          amount of time.
 *
 * @see optimal_table_len
 * @see bn_gen_prime
 * @see bn_rabin
 * @see bn_mod_exp
 */
void gen_provable_primes_arithmetic(bignum* p, u64 n,
                                    pocklington_cert** cert_out);

#endif
