#ifndef RZPRIME_H
#define RZPRIME_H

#include "rabin_errors.h"
#include "rz.h"
#include "rzcert.h"

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
 * @return RABIN_SUCCESS on success, RABIN_ERR_OUT_OF_MEMORY on allocation
 * failure, or RABIN_ERR_INVALID_ARG if /dev/urandom cannot be opened,
 * reading fails, or \f$bits < 2\f$.
 */
rabin_err_t rz_gen_prime(rz_t* p, int bits);

/**
 * @brief Generates a safe prime number \f$p\f$ of a specified bit length.
 *
 * A prime \f$p\f$ is safe if \f$p = 2q + 1\f$, where \f$q\f$ is also a prime (a
 * Sophie Germain prime). This function repeatedly generates a random prime
 * \f$q\f$ of length \f$bits - 1\f$ until the calculated \f$p\f$ passes the BPSW
 * primality test.
 *
 * @param[out] p    Pointer to the rz_t structure where the generated safe
 * prime will be stored.
 * @param[in]  bits The desired total bit length of the safe prime \f$p\f$.
 *
 * @return RABIN_SUCCESS on success, RABIN_ERR_OUT_OF_MEMORY on allocation
 * failure, or RABIN_ERR_INVALID_ARG if /dev/urandom cannot be opened or
 * \f$bits < 3\f$.
 */
rabin_err_t rz_gen_safe_prime(rz_t* p, int bits);

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
bool rz_stronglucas(const rz_t* n, const rz_t* P, const rz_t* Q);

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
bool rz_bpsw(const rz_t* n);

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
bool rz_check_lemma1(const rz_t* n, const rz_t* n_min1, const rz_t* a,
                     const rz_t* q);

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
double rz_gen_rel_size();

/**
 * @brief Generate a provable prime of \f$k\f$ bits into \f$p\f$ (Maurer's
 * algorithm).
 *
 * Repeatedly invokes rz_provable_prime_inner() until it succeeds; each
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
 * @param[in]  k Desired bit length of the prime (>= 2).
 *
 * @return RABIN_SUCCESS on success, RABIN_ERR_OUT_OF_MEMORY on allocation
 * failure, or RABIN_ERR_INVALID_ARG if /dev/urandom cannot be opened or
 * \f$k < 2\f$.
 */
rabin_err_t rz_provable_prime(rz_t* p, u64 k);

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
 * @return RABIN_SUCCESS on success (the prime is stored in \f$p\f$),
 * RABIN_ERR_INVALID_ARG if the recursion fails or the search exceeds 1000
 * candidates (caller should retry), or RABIN_ERR_OUT_OF_MEMORY.
 */
rabin_err_t rz_provable_prime_inner(rz_t* p, u64 k);

/**
 * @brief Generate `count` Proth primes of the form \f$p = c \cdot 2^k + 1\f$
 * and print them as a C array initializer.
 *
 * Let \f$k =\f$ exponent of the power of two.
 *
 * Scans odd \f$c\f$ starting from the given \f$c\f$ (rounded up to odd),
 * testing each candidate \f$c \cdot 2^k + 1\f$ with BPSW, until `count` primes
 * are found or the 64-bit search space is exhausted. The output is a `static
 * const uint64_t RNS_PRIMES[]` table suitable for pasting into rns_primes.h.
 *
 * Complexity:
 *   - Time: \f$O(count \cdot k^3)\f$ expected - one BPSW test per candidate
 *   - Auxiliary memory: \f$O(k/64)\f$ limbs
 *   - Output memory: \f$O(count)\f$ printed entries
 *
 * @param[in] count Number of Proth primes to generate.
 * @param[in]     k Exponent of the power of two.
 * @param[in]     c Starting value for the odd multiplier \f$c\f$.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_INVALID_ARG,
 * RABIN_ERR_OVERFLOW, or RABIN_ERR_OUT_OF_MEMORY.
 */
rabin_err_t rz_gen_proth_primes(u64 count, u64 k, u64 c);

/**
 * @brief Generate `count` primes just below 2^64 and print them as a C array
 * initializer.
 *
 * Scans downward from \f$2^{64} - 47\f$ in steps of \f$2\f$ (so all candidates
 * are of the form \f$2^{64} - 1 - 2i\f$, i.e. close to the top of the 64-bit
 * range), testing each candidate with BPSW until `count` primes are
 * found. The output is a `static const u64 RNS_PRIMES[]` table
 * suitable for pasting into rns_primes.h.
 *
 * Complexity:
 *   - Time: \f$O(count \cdot 1024^3)\f$ expected - one BPSW test per candidate
 *   - Auxiliary memory: \f$O(1)\f$ limbs (single-limb candidates)
 *   - Output memory: \f$O(count)\f$ printed entries
 *
 * @param[in] count Number of primes to generate.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_INVALID_ARG,
 * RABIN_ERR_OVERFLOW, or RABIN_ERR_OUT_OF_MEMORY.
 */
rabin_err_t rns_gen_primes(u64 count);

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
u64 rz_optimal_table_len(u64 n, u64 B);

/**
 * @brief Generates a provable prime of approximately @p n bits.
 *
 * @details
 * Recursively constructs a prime number @p p using a Pocklington-style
 * primality proof. For small @p n, this delegates to rz_gen_prime().
 *
 * For larger @p n, the function performs the following:
 *  1. Recursively generates a provable prime F of \f$ \approx n/2 \f$ bits.
 *  2. Chooses a random integer t such that:
 *     \f$ 2^{(n - 2)}/F < t < 2^{(n - 1)}/F - s \cdot n \f$
 *  3. Constructs an arithmetic progression \f$ N_i = N_0 + i \cdot a \f$
 *     where \f$ a = 2F \f$ and \f$ N_0 = t \cdot a + 1 \f$.
 *  4. Sieves candidates and applies a base-2 Rabin-Miller filter.
 *  5. Proves primality using Pocklington's lemma with the known factor F.
 *
 * @param[out] p        Initialized rz_t receiving the provable prime.
 * @param[out] cert_out Initialized certificate struct. If NULL, no cert is
 * stored.
 * @param[in]  n        Desired approximate bit length.
 *
 * @return RABIN_SUCCESS on success, or an appropriate error code (e.g.,
 *         RABIN_ERR_OUT_OF_MEMORY).
 *
 * @pre The rz_t library and random-number subsystem must be initialized.
 *      The global small-prime table must contain enough primes for trial
 * division.
 *
 * @warning This function is recursive. High @p n values require significant
 * stack depth and dynamic allocations. It loops until a prime is found, meaning
 *          execution time is theoretically unbounded.
 *
 * @par Algorithm Reference:
 * U. Maurer, "Fast Generation of Prime Numbers and Secure Public-Key
 * Cryptographic Parameters," Journal of Cryptology, vol. 8, pp. 123–155, 1995.
 *
 * @see rz_optimal_table_len()
 * @see rz_gen_prime()
 */
rabin_err_t rz_gen_provable_arithmetic(rz_t* p, u64 n, rzcert_t** cert_out);

#endif
