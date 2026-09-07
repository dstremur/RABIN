#ifndef BIGMATH_H
#define BIGMATH_H

#include "bigcore.h"

/**
 * @brief Compute the Jacobi symbol (a / m).
 *
 * Let n = m->size, measured in 64-bit limbs.
 *
 * Returns 1 if a is a quadratic residue mod m, -1 if it is a
 * nonresidue, and 0 if gcd(a, m) > 1. m must be positive and odd;
 * otherwise 0 is returned and a message is printed.
 *
 * Uses the binary algorithm: repeatedly strip factors of 2 from a
 * (flipping the sign when both the stripped power and m mod 8 are
 * odd), apply the mutual-reciprocity rule when both a and m are 3
 * mod 4, and reduce m mod a.
 *
 * Complexity:
 *   Time: O(n^2) - O(n) iterations of O(n) modular reductions
 *   Auxiliary memory: O(n) limbs for temporaries
 *   Output memory: O(1)
 *
 * @param[in] a Numerator of the Jacobi symbol.
 * @param[in] m Denominator (must be positive and odd).
 *
 * @return 1 If a is a quadratic residue mod m, -1 if a nonresidue, 0
 *           if gcd(a, m) > 1 or m is not positive and odd.
 */
i64 bn_jacobi(const bignum* a, const bignum* m);

/**
 * @brief Tonelli-Shanks: square root of n modulo the prime p.
 *
 * Let n_l = p->size, measured in 64-bit limbs.
 *
 * Inputs:
 * p, a prime
 * n, an element of Z / p Z such that solutions to the congruence
 * r^2 = n exist; when this is so we say that n is a quadratic
 * residue mod p.
 *
 * Computes r with r^2 = n (mod p) and stores it in r. The algorithm
 * factors p - 1 = Q * 2^S with Q odd, finds a quadratic nonresidue z,
 * and iteratively refines the candidate root R until t = n^Q * c^2
 * reaches 1.
 *
 * If n is zero, r is set to 0. If n is not a quadratic residue mod p
 * (Jacobi symbol != 1), a message is printed and r is left unchanged.
 *
 * Complexity:
 *   Time: O(n_l^2 * S^2) worst case - O(S) iterations, each with O(S)
 *         modular squarings, plus O(n_l^2) for the initial
 *         exponentiations and the nonresidue search
 *   Auxiliary memory: O(n_l) limbs for temporaries
 *   Output memory: O(n_l) limbs
 *
 * @param[out] r Result storing a square root of n mod p (if it exists).
 * @param[in]  n Value whose square root is computed (a quadratic
 *               residue mod p).
 * @param[in]  p Prime modulus.
 */
void tonelli_shanks(bignum* r, const bignum* n, const bignum* p);

/**
 * @brief  Computes the greatest common divisor (GCD) of two bignums.
 *
 * @details Uses the classic Euclidean algorithm,
 *          @c gcd(a, b) == gcd(b, a % b), iterating until the remainder
 *          is zero; the last non-zero value is the GCD. The loop rotates
 *          three internal buffers by pointer swap, so no temporary
 *          bignums are allocated per iteration.
 *
 *          Zero handling follows the standard conventions:
 *          - @p a is zero  -> result is @p b
 *          - @p b is zero  -> result is @p a
 *          - both zero     -> result is zero, i.e. `gcd(0, 0) == 0`
 *
 *          The result is the principal (non-negative) GCD; the sign of the
 *          inputs is ignored, so `bn_gcd(d, &a, &b) == bn_gcd(d, &a, &nb)`.
 *
 * @param[out] d  Receives `gcd(a, b)`. Must be initialized before the call.
 *                May alias @p a or @p b — both operands are copied into
 *                temporaries before any work is done — so
 *                `bn_gcd(&a, &a, &b)` is valid.
 * @param[in]  a  First operand. Not modified. Must not be @c NULL.
 * @param[in]  b  Second operand. Not modified. Must not be @c NULL.
 *
 * @pre @p d, @p a and @p b are initialized (e.g. via bn_init() /
 *      bn_init_multi()) and have valid size fields.
 *
 * @note The three temporaries created here are freed before returning;
 *       the caller only has to manage the lifetime of @p d.
 * @note bn_mod() is only called with a non-zero divisor, which is
 *       guaranteed by the `while (!bn_is_zero(v))` loop condition.
 *
 * @par Memory
 * Allocates 3 temporary bignums (each up to `max(size(a), size(b))`),
 * peak extra memory ≈ 3 operands. Fails silently on allocation error
 * if the underlying allocator does, so check @p d if that matters to you.
 *
 * @par Complexity
 * O(log min(a, b)) modulo operations; each modulo is O(n·m) word
 * divisions for n- and m-word operands.
 *
 * @warning <b>Not constant-time.</b> The number and shape of divisions
 *          depend on the operand values, which leaks information through
 *          execution time. Do not use on secret inputs (e.g. private keys
 *          in `invmod`/key-derivation paths) without a constant-time
 *          variant such as binary GCD.
 *
 * @par Example
 * @code
 * bignum a, b, g;
 * bn_init_multi(&a, &b, &g, NULL);
 *
 * bn_set_str(&a, "1071", 10);   // 1071 = 3 * 3 * 7 * 17
 * bn_set_str(&b, "462",  10);   //  462 = 2 * 3 * 7 * 11
 * bn_gcd(&g, &a, &b);           // g == 21
 *
 * bn_gcd(&a, &a, &b);           // also fine: a == 21, b unchanged
 *
 * bn_free_multi(&a, &b, &g, NULL);
 * @endcode
 *
 * @see bn_mod(), bn_copy(), bn_is_zero()
 */
void bn_gcd(bignum* d, const bignum* a, const bignum* b);

#endif
