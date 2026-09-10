#ifndef BIGPRIME_H
#define BIGPRIME_H

#include "bigcert.h"
#include "bigcore.h"

/**
 * @brief Generate a random prime of the given bit length into \f$p\f$.
 *
 * Let \f$k =\f$ bits.
 *
 * Repeatedly draws random \f$k\f$-bit candidates from /dev/urandom with the
 * top and bottom bits set, rejects candidates divisible by 3, 5, 7,
 * 11, 13, or 17, trial-divides against the first 1500 primes, and
 * accepts the first candidate that passes BPSW.
 *
 * Returns true on success, false if /dev/urandom cannot be opened or
 * reading from it fails.
 *
 * Complexity:
 *   - Time: \f$O(k \log k)\f$ expected - about \f$k / \ln(k)\f$ candidates,
 * each costing trial division plus one BPSW test
 *   - Auxiliary memory: \f$O(k/64)\f$ limbs
 *   - Output memory: \f$O(k/64)\f$ limbs
 *
 * @param[out] p    Result storing the generated prime.
 * @param[in]  bits Desired bit length of the prime.
 *
 * @return true  On success.
 * @return false If /dev/urandom cannot be opened or reading fails.
 */
bool bn_gen_prime(bignum* p, int bits);

/**
 * @brief Generates a safe prime number \f$p\f$ of a specified bit length.
 *
 * A prime \f$p\f$ is safe if \f$p = 2q + 1\f$, where \f$q\f$ is also a prime (a
 * Sophie Germain prime). This function repeatedly generates a random prime
 * \f$q\f$ of length \f$bits - 1\f$ until the calculated \f$p\f$ passes the BPSW
 * primality test.
 *
 * @param[out] p    Pointer to the bignum structure where the generated safe
 * prime will be stored.
 * @param[in]  bits The desired total bit length of the safe prime \f$p\f$.
 *
 * @return True  If the safe prime was successfully generated.
 * @return False If prime generation failed (e.g., maximum iterations reached or
 * invalid bit size).
 */
bool bn_gen_safe_prime(bignum* p, int bits);

/**
 * @brief Strong Lucas test of \f$n\f$ with Lucas parameters \f$(P, Q)\f$.
 *
 * Let \f$n_l =\f$ n->size, measured in 64-bit limbs.
 *
 * Writes \f$n + 1 = d \cdot 2^s\f$ with \f$d\f$ odd, computes \f$(U_d, V_d,
 * Q^d) \bmod n\f$, and \f$n\f$ passes if:
 *
 *   - \f$U_d = 0 \bmod n\f$, or
 *   - \f$V_d = 0 \bmod n\f$, or
 *   - \f$V_{d \cdot 2^r} = 0 \bmod n\f$ for some \f$1 \le r < s\f$
 *
 * (\f$V\f$ is iterated with the doubling identity \f$V_{2k} = V_k^2 - 2
 * Q^k\f$.)
 *
 * Returns true if \f$n\f$ passes the test (probably prime), false otherwise.
 *
 * Complexity:
 *   - Time: \f$O(n_l^3)\f$ for the Lucas sequence computation (Montgomery
 *           context init), then \f$O(s \cdot n_{l}^2)\f$ for the \f$V\f$
 * iteration
 *   - Auxiliary memory: \f$O(n_l)\f$ limbs for temporaries
 *   - Output memory: \f$O(1)\f$
 *
 * @param[in] n Number to test.
 * @param[in] P Lucas parameter P.
 * @param[in] Q Lucas parameter Q.
 *
 * @return true  If \f$n\f$ passes the test (probably prime).
 * @return false If \f$n\f$ fails the test.
 */
bool bn_stronglucas(const bignum* n, bignum* P, bignum* Q);

/**
 * @brief Baillie-PSW primality test.
 *
 * Let \f$n_l =\f$ n->size, measured in 64-bit limbs.
 *
 * \f$n\f$ passes if all of the following hold:
 *
 *   1. \f$n\f$ is not divisible by any prime below 10000 (\f$n\f$ itself is
 *      accepted if it is one of them),
 *   2. \f$n\f$ passes Miller-Rabin to base \f$2\f$,
 *   3. \f$n\f$ passes the strong Lucas test with parameters \f$(P, Q)\f$ =
 *      \f$(1, (1 - D)/4)\f$ or \f$(1, (1 + D)/4)\f$, where \f$D\f$ is the first
 * value in the sequence \f$5, -7, 9, -11, \ldots\f$ (Selfridge's method A*)
 *      with Jacobi symbol \f$\left(\frac{D}{n}\right) = -1\f$.
 *
 * If no such \f$D\f$ is found within 15 rounds, \f$n\f$ is rejected (this also
 * catches perfect squares). No known Baillie-PSW pseudoprime exists.
 *
 * Returns true if \f$n\f$ is probably prime, false if \f$n\f$ is even or a
 * witness for compositeness was found.
 *
 * Complexity:
 *   - Time: \f$O(n_l^3)\f$ - dominated by the Miller-Rabin and strong Lucas
 *           tests (each paying a Montgomery context initialization)
 *   - Auxiliary memory: \f$O(n_l)\f$ limbs for temporaries
 *   - Output memory: \f$O(1)\f$
 *
 * @param[in] n Number to test.
 *
 * @return true  If \f$n\f$ is probably prime.
 * @return false If \f$n\f$ is even or a witness for compositeness was found.
 */
bool bn_bpsw(const bignum* n);

/**
 * @brief Lemma 1 check of Maurer's algorithm for the case \f$r = 1\f$.
 *
 * Let \f$n_l =\f$ n->size, measured in 64-bit limbs.
 *
 * Given \f$n = 2 R q + 1\f$ with \f$q\f$ prime, this checks the primality
 * criterion for a random base \f$a\f$:
 *
 *   \f$X = a^{(n - 1)/q} \bmod n\f$
 *
 * \f$n\f$ is certified if \f$X \neq 1\f$, \f$\gcd(X - 1, n) = 1\f$, and \f$X^q
 * = 1 \bmod n\f$.
 *
 * Returns true if the check passes (\f$n\f$ is prime), false otherwise.
 *
 * Complexity:
 *   - Time: \f$O(n_l^2)\f$ - two modular exponentiations
 *   - Auxiliary memory: \f$O(n_l)\f$ limbs for temporaries
 *   - Output memory: \f$O(1)\f$
 *
 * @param[in]      n      Number to certify (\f$n = 2 R q + 1\f$).
 * @param[in]      n_min1 Value \f$n - 1\f$.
 * @param[in]      a      Random base.
 * @param[in]      q      Prime factor (\f$q\f$ in \f$n = 2 R q + 1\f$).
 *
 * @return true  If the check passes (\f$n\f$ is prime).
 * @return false If the check fails.
 */
bool checkLemma1(bignum* n, bignum* n_min1, bignum* a, bignum* q);

/**
 * @brief Draw a relative size for Maurer's algorithm: \f$2^u\f$ with \f$u\f$
 * uniform in
 * \f$[0, 1)\f$, i.e. a value in \f$[1/2, 1]\f$ biased toward 1.
 *
 * Complexity:
 *   - Time: \f$O(1)\f$
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(1)\f$
 *
 * @return The relative size, a value in \f$[1/2, 1]\f$.
 */
double gen_rel_size();

/**
 * @brief Generate a provable prime of \f$k\f$ bits into \f$p\f$ (Maurer's
 * algorithm).
 *
 * Repeatedly invokes bn_provable_prime_inner() until it succeeds; each
 * failed attempt unwinds its recursion and frees all memory, so the
 * loop is safe to retry indefinitely.
 *
 * Complexity:
 *   - Time: \f$O(k^3)\f$ expected per successful attempt (recursive provable
 *           prime generation plus trial division and Lemma 1 checks)
 *   - Auxiliary memory: \f$O(k^2 / 64)\f$ limbs on the recursion stack
 *   - Output memory: \f$O(k/64)\f$ limbs
 *
 * @param[out] p Result storing the provable prime.
 * @param[in]  k Desired bit length of the prime.
 */
void bn_provable_prime(bignum* p, u64 k);

/**
 * @brief One attempt at Maurer's simpler algorithm for a provable k-bit prime.
 *
 * Let \f$k =\f$ bit length of the target prime.
 *
 * Strategy:
 *
 *   - base case \f$k \le 20\f$: accept a random \f$k\f$-bit candidate that
 * passes BPSW (treated as proven at this size),
 *   - otherwise: pick a relative size \f$rel_{size}\f$ in \f$[1/2, 1]\f$ with
 *     \f$rel_{size} \cdot k < k - k/6\f$, recursively generate a provable prime
 * \f$q\f$ of that size, then search for \f$n = 2 \cdot R q + 1\f$ with \f$R\f$
 * in
 *     \f$[2^{(k - 1)}/(2q), 2^k/(2q)]\f$ that survives trial division up to
 *     \f$0.1 k^2 + 1\f$ and passes the Lemma 1 check for some random base
 * \f$a\f$ (up to 200 bases).
 *
 * Returns true and stores the prime in \f$p\f$ on success. Returns false
 * (having freed all temporaries) if the recursion fails or the search
 * exceeds 1000 candidates.
 *
 * Complexity:
 *   - Time: \f$O(k^3)\f$ expected - dominated by the recursive call and the
 *           Lemma 1 modular exponentiations
 *   - Auxiliary memory: \f$O(k^2 / 64)\f$ limbs on the recursion stack
 *   - Output memory: \f$O(k/64)\f$ limbs
 *
 * @param[out] p Result storing the provable prime.
 * @param[in]  k Desired bit length of the prime.
 *
 * @return true  On success (the prime is stored in \f$p\f$).
 * @return false If the recursion fails or the search exceeds 1000
 *               candidates.
 */
bool bn_provable_prime_inner(bignum* p, u64 k);

/**
 * @brief Generate `count` Proth primes of the form \f$p = c \cdot 2^k + 1\f$
 * and print them as a C array initializer.
 *
 * Let \f$k =\f$ exponent of the power of two.
 *
 * Scans odd \f$c\f$ starting from the given \f$c\f$ (rounded up to odd),
 * testing each candidate \f$c \cdot 2^k + 1\f$ with BPSW, until `count` primes
 * are found or the 64-bit search space is exhausted. The output is a `static
 * const uint64_t RNS_PRIMES[]` table suitable for pasting into primes.h.
 *
 * Complexity:
 *   - Time: \f$O(count \cdot k^3)\f$ expected - one BPSW test per candidate
 *   - Auxiliary memory: \f$O(k/64)\f$ limbs
 *   - Output memory: \f$O(count)\f$ printed entries
 *
 * @param[in] count Number of Proth primes to generate.
 * @param[in]     k Exponent of the power of two.
 * @param[in]     c Starting value for the odd multiplier \f$c\f$.
 */
void bn_gen_proth_primes(u64 count, u64 k, u64 c);

/**
 * @brief Generate `count` primes just below 2^64 and print them as a C array
 * initializer.
 *
 * Scans downward from \f$2^{64} - 47\f$ in steps of \f$2\f$ (so all candidates
 * are of the form \f$2^{64} - 1 - 2i\f$, i.e. close to the top of the 64-bit
 * range), testing each candidate with BPSW until `count` primes are
 * found. The output is a `static const u64 RNS_PRIMES[]` table
 * suitable for pasting into primes.h.
 *
 * Complexity:
 *   - Time: \f$O(count \cdot 1024^3)\f$ expected - one BPSW test per candidate
 *   - Auxiliary memory: \f$O(1)\f$ limbs (single-limb candidates)
 *   - Output memory: \f$O(count)\f$ printed entries
 *
 * @param[in] count Number of primes to generate.
 */
void gen_rns_primes(u64 count);

/**
 * @brief Heuristic table length for the arithmetic-progression sieve used
 * by Pocklington prime generation.
 *
 * Returns the sieve length \f$s\f$ (number of candidates \f$N_0 + i\cdota\f$,
 * \f$i\f$ in
 * \f$[0, s)\f$) for a target prime of \f$n\f$ bits where each candidate is
 * processed in \f$B\f$-bit words. \f$s\f$ comes from the heuristic
 *
 *   \f$s = (0.4 \cdot n \cdot m) / \log(n^2 \cdot m)\f$,   \f$m = n / B\f$
 *
 * which balances the cost of the small-prime trial table against the
 * expected number of progression steps needed to find a candidate.
 *
 * Complexity:
 *   - Time: \f$O(1)\f$
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(0)\f$
 *
 * @param[in] n Target prime size in bits.
 * @param[in] B Word size in bits used per candidate (e.g. 64).
 *
 * @return The heuristic sieve length \f$s\f$.
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
 *     \f$n / 2\f$ bits.
 *  2. Chooses a random integer @c t such that:
 *
 *         \f$2^{(n - 2)}/F < t < 2^{(n - 1)}/F - s \cdot n\f$
 *
 *     where @c s is the sieve/table length returned by optimal_table_len().
 *
 *  3. Constructs an arithmetic progression of candidate integers:
 *
 *         \f$N_i = N_0 + i\cdota\f$
 *
 *     where:
 *
 *         \f$a = 2F\f$
 *         \f$N_0 = t\cdota + 1\f$
 *         \f$0 \le i \le s\f$
 *
 *  4. Sieves candidates using small primes up to a bound proportional to @p n.
 *
 *  5. Applies a Rabin-Miller test, currently with base \f$2\f$, as a fast
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
 *       \f$F > \sqrt{N - 1}\f$, or an equivalent sufficient condition. If this
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
