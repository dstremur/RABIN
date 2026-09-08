#ifndef BIGMATH_H
#define BIGMATH_H

#include "bigcore.h"

/**
 * @brief Compute the Jacobi symbol \f$\left(\frac{a}{m}\right)\f$.
 *
 * Let \f$n =\f$ m->size, measured in 64-bit limbs.
 *
 * Returns \f$1\f$ if \f$a\f$ is a quadratic residue \f$\bmod m\f$, \f$-1\f$ if
 * it is a nonresidue, and \f$0\f$ if \f$\gcd(a, m) > 1\f$. \f$m\f$ must be
 * positive and odd; otherwise \f$0\f$ is returned and a message is printed.
 *
 * Uses the binary algorithm: repeatedly strip factors of \f$2\f$ from \f$a\f$
 * (flipping the sign when both the stripped power and \f$m \bmod 8\f$ are
 * odd), apply the mutual-reciprocity rule when both \f$a\f$ and \f$m\f$ are
 * \f$3
 * \bmod 4\f$, and reduce \f$m \bmod a\f$.
 *
 * Complexity:
 *   - Time: \f$O(n^2)\f$ - \f$O(n)\f$ iterations of \f$O(n)\f$ modular
 * reductions
 *   - Auxiliary memory: \f$O(n)\f$ limbs for temporaries
 *   - Output memory: \f$O(1)\f$
 *
 * @param[in] a Numerator of the Jacobi symbol.
 * @param[in] m Denominator (must be positive and odd).
 *
 * @return 1 If \f$a\f$ is a quadratic residue \f$\bmod m\f$, \f$-1\f$ if a
 * nonresidue, \f$0\f$ if \f$\gcd(a, m) > 1\f$ or \f$m\f$ is not positive and
 * odd.
 */
i64 bn_jacobi(const bignum* a, const bignum* m);

/**
 * @brief Tonelli-Shanks: square root of \f$n\f$ modulo the prime \f$p\f$.
 *
 * Let \f$n_l =\f$ p->size, measured in 64-bit limbs.
 *
 * Inputs:
 * \f$p\f$, a prime
 * \f$n\f$, an element of \f$Z / p Z\f$ such that solutions to the congruence
 * \f$r^2 = n\f$ exist; when this is so we say that \f$n\f$ is a quadratic
 * residue \f$\bmod p\f$.
 *
 * Computes \f$r\f$ with \f$r^2 = n\f$ (\f$\bmod p\f$) and stores it in \f$r\f$.
 * The algorithm factors \f$p - 1 = Q \cdot 2^S\f$ with \f$Q\f$ odd, finds a
 * quadratic nonresidue \f$z\f$, and iteratively refines the candidate root
 * \f$R\f$ until \f$t = n^Q \cdot c^2\f$ reaches \f$1\f$.
 *
 * If \f$n\f$ is zero, \f$r\f$ is set to \f$0\f$. If \f$n\f$ is not a quadratic
 * residue \f$\bmod p\f$ (Jacobi symbol \f$\neq 1\f$), a message is printed and
 * \f$r\f$ is left unchanged.
 *
 * Complexity:
 *   - Time: \f$O(n_l^2 \cdot S^2)\f$ worst case - \f$O(S)\f$ iterations, each
 * with \f$O(S)\f$ modular squarings, plus \f$O(n_l^2)\f$ for the initial
 *           exponentiations and the nonresidue search
 *   - Auxiliary memory: \f$O(n_l)\f$ limbs for temporaries
 *   - Output memory: \f$O(n_l)\f$ limbs
 *
 * @param[out] r Result storing a square root of \f$n \bmod p\f$ (if it exists).
 * @param[in]  n Value whose square root is computed (a quadratic
 *               residue \f$\bmod p\f$).
 * @param[in]  p Prime modulus.
 */
void tonelli_shanks(bignum* r, const bignum* n, const bignum* p);

/**
 * @brief  Computes the greatest common divisor (GCD) of two bignums.
 *
 * @details Uses the classic Euclidean algorithm,
 *          \f$\gcd(a, b) = \gcd(b, a \bmod b)\f$, iterating until the remainder
 *          is zero; the last non-zero value is the GCD. The loop rotates
 *          three internal buffers by pointer swap, so no temporary
 *          bignums are allocated per iteration.
 *
 *          Zero handling follows the standard conventions:
 *          - @p a is zero  -> result is @p b
 *          - @p b is zero  -> result is @p a
 *          - both zero     -> result is zero, i.e. \f$\gcd(0, 0) = 0\f$
 *
 *          The result is the principal (non-negative) GCD; the sign of the
 *          inputs is ignored, so `bn_gcd(d, &a, &b) == bn_gcd(d, &a, &nb)`.
 *
 * @param[out] d  Receives \f$\gcd(a, b)\f$. Must be initialized before the
 * call. May alias @p a or @p b — both operands are copied into temporaries
 * before any work is done — so `bn_gcd(&a, &a, &b)` is valid.
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
 * Allocates 3 temporary bignums (each up to \f$\max(size(a), size(b))\f$),
 * peak extra memory \f$\approx 3\f$ operands. Fails silently on allocation
 * error if the underlying allocator does, so check @p d if that matters to you.
 *
 * @par Complexity
 * \f$O(\log \min(a, b))\f$ modulo operations; each modulo is \f$O(n \cdot m)\f$
 * word divisions for \f$n\f$- and \f$m\f$-word operands.
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
