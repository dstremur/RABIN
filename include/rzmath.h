#ifndef RZMATH_H
#define RZMATH_H

#include "rabin_errors.h"
#include "rz.h"

/**
 * @brief Compute the Jacobi symbol \f$\left(\frac{a}{m}\right)\f$.
 *
 * Wrapper around rz_kronecker()
 *
 * Complexity:
 *   - Time: \f$O(\ln(n)^2)\f$
 *   - Auxiliary memory: \f$O(n)\f$ limbs for temporaries
 *   - Output memory: \f$O(1)\f$
 *
 * @param[in] a Numerator of the Jacobi symbol.
 * @param[in] m Denominator (must be positive and odd).
 *
 * @return 1 If \f$a\f$ is a quadratic residue \f$\bmod m\f$, \f$-1\f$ if a
 * nonresidue, \f$0\f$ if \f$\gcd(a, m) > 1\f$ or \f$m\f$ is not positive and
 * odd.
 * @see rz_kronecker()
 */
i64 rz_jacobi(const rz_t* a, const rz_t* m);

/**
 * @brief Compute the Kronecker symbol \f$\left(\frac{a}{b}\right)\f$.
 *
 * Uses the standard algorithm based around quadratic reciprocity.
 *
 * Complexity:
 *   - Time: \f$O(\ln(n)^2)\f$
 *   - Auxiliary memory: \f$O(n)\f$ limbs for temporaries
 *   - Output memory: \f$O(1)\f$
 *
 * @param[in] a Numerator of the Kronecker symbol.
 * @param[in] b Denominator.
 *
 * @return 1 If \f$a\f$ is a quadratic residue \f$\bmod b\f$, \f$-1\f$ if a
 * nonresidue, \f$0\f$ if \f$\gcd(a, b) > 1\f$ or \f$b\f$ is not positive and
 * odd.
 *
 * @par Algorithm Reference:
 * H. Cohen, "A Course in Computational Algebraic Number Theory,"
 * Springer-Verlag, Berlin, 1993, p. 29, Algorithm 1.4.10.
 * @see rz_jacobi()
 */
i64 rz_kronecker(const rz_t* a, const rz_t* b);

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
 * residue \f$\bmod p\f$ (Jacobi symbol \f$\neq 1\f$), the call fails and
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
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR,
 * RABIN_ERR_INVALID_ARG (if \f$n\f$ is not a quadratic residue \f$\bmod
 * p\f$), or RABIN_ERR_OUT_OF_MEMORY.
 *
 * @par Algorithm Reference:
 * A. Tonelli, "Solution generale de l'equation \f$i^{i} = n \pmod{p}\f$,"
 * Rendiconti del Circolo Matematico di Palermo, vol. 11, 1891.
 * D. E. Shanks, "Solving the Congruence \f$g^x \equiv a \pmod{p}\f$," in
 * Proceedings of the First Symposium on Symbolic and Algebraic
 * Computation, ACM, 1971.
 * @see rz_mod_exp(), rz_jacobi()
 */
rabin_err_t rz_tonelli_shanks(rz_t* r, const rz_t* n, const rz_t* p);

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
 *          inputs is ignored, so `rz_gcd(d, &a, &b) == rz_gcd(d, &a, &nb)`.
 *
 * @param[out] d  Receives \f$\gcd(a, b)\f$. Must be initialized before the
 * call. May alias @p a or @p b — both operands are copied into temporaries
 * before any work is done — so `rz_gcd(&a, &a, &b)` is valid.
 * @param[in]  a  First operand. Not modified. Must not be @c NULL.
 * @param[in]  b  Second operand. Not modified. Must not be @c NULL.
 *
 * @pre @p d, @p a and @p b are initialized (e.g. via rz_init() /
 *      rz_init_multi()) and have valid size fields.
 *
 * @note rz_mod() is only called with a non-zero divisor, which is
 *       guaranteed by the `while (!rz_is_zero(v))` loop condition.
 *
 * Complexity:
 *   - Time: \f$O(\log \min(a, b))\f$ modulo operations; each modulo is
 *     \f$O(n \cdot m)\f$ word divisions for \f$n\f$- and \f$m\f$-word operands
 *   - Auxiliary memory: \f$O(\max(a, b))\f$ limbs in three rotating buffers
 *   - Output memory: \f$O(\max(a, b))\f$ limbs
 *
 * @warning <b>Not constant-time.</b> The number and shape of divisions
 *          depend on the operand values, which leaks information through
 *          execution time. Do not use on secret inputs (e.g. private keys
 *          in `invmod`/key-derivation paths) without a constant-time
 *          variant such as binary GCD.
 *
 * @par Example
 * @code
 * rz_t a, b, g;
 * rz_init_multi(&a, &b, &g, NULL);
 *
 * rz_set_str(&a, "1071", 10);   // 1071 = 3 * 3 * 7 * 17
 * rz_set_str(&b, "462",  10);   //  462 = 2 * 3 * 7 * 11
 * rz_gcd(&g, &a, &b);           // g == 21
 *
 * rz_gcd(&a, &a, &b);           // also fine: a == 21, b unchanged
 *
 * rz_clear_multi(&a, &b, &g, NULL);
 * @endcode
 *
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR /
 * RABIN_ERR_OUT_OF_MEMORY.
 * @see rz_mod(), rz_copy(), rz_is_zero()
 */
rabin_err_t rz_gcd(rz_t* d, const rz_t* a, const rz_t* b);

rabin_err_t rz_gcd_binary(rz_t* d, const rz_t* a, const rz_t* b);

rabin_err_t rz_gcd_lehmer(rz_t* d, const rz_t* a, const rz_t* b);

rabin_err_t rz_lcm(rz_t* l, const rz_t* a, const rz_t* b);
/**
 * @brief Extended Euclidean algorithm: Bézout coefficients of a and b.
 *
 * Computes integers \f$u\f$ and \f$v\f$ with \f$u \cdot a + v \cdot b =
 * d = \gcd(a, b)\f$, with \f$d \ge 0\f$ (GMP convention). \p d is
 * initialized from \f$a\f$ before the iteration, \f$v\f$ is derived as
 * \f$v = (d - u \cdot a) / b\f$, and neither operand is modified.
 *
 * Complexity:
 *   - Time: \f$O(\log \min(a, b))\f$ modulo operations, see rz_gcd()
 *   - Auxiliary memory: \f$O(\max(a, b))\f$ limbs
 *   - Output memory: \f$O(\max(a, b))\f$ limbs per result
 *
 * @param[out] u Receives the Bézout coefficient of \f$a\f$.
 * @param[out] v Receives the Bézout coefficient of \f$b\f$.
 * @param[out] d Receives \f$\gcd(a, b)\f$ (\f$\ge 0\f$).
 * @param[in]  a First operand.
 * @param[in]  b Second operand.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR /
 * RABIN_ERR_OUT_OF_MEMORY.
 * @see rz_gcd()
 */
rabin_err_t rz_gcd_extended(rz_t* u, rz_t* v, rz_t* d, const rz_t* a,
                            const rz_t* b);

rabin_err_t rz_gcd_extended_lehmer(rz_t* u, rz_t* v, rz_t* d, const rz_t* a,
                                   const rz_t* b);

/**
 * @brief Compute a solution to the Diophantine equation \f$x^2 + d \cdot y^2
 * = p\f$.
 *
 * Uses the algorithm of Cornacchia: finds \f$x_0 = \sqrt{-d} \bmod p\f$
 * (via rz_tonelli_shanks()), runs the Euclidean algorithm on
 * \f$(p, x_0)\f$ down to \f$\lfloor\sqrt{p}\rfloor\f$, and checks that
 * \f$d \mid (p - b^2)\f$ and that \f$(p - b^2)/d\f$ is a perfect square.
 *
 * Complexity:
 *   - Time: \f$O(\log(n)^2)\f$
 *   - Auxiliary memory: \f$O(n)\f$ limbs
 *   - Output memory: \f$O(n)\f$ limbs
 *
 * @param[out] x Receives the solution \f$x\f$.
 * @param[out] y Receives the solution \f$y\f$.
 * @param[in]  p Prime modulus.
 * @param[in]  d Coefficient \f$d > 0\f$.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR,
 * RABIN_ERR_INVALID_ARG (if \f$-d\f$ is a quadratic nonresidue \f$\bmod
 * p\f$ or no solution exists), or RABIN_ERR_OUT_OF_MEMORY.
 *
 * @par Algorithm Reference:
 * H. Cohen, "A Course in Computational Algebraic Number Theory,"
 * Springer-Verlag, Berlin, 1993, p. 34, Algorithm 1.3.5.
 * @see rz_tonelli_shanks()
 */
rabin_err_t rz_cornacchia(rz_t* x, rz_t* y, const rz_t* p, const rz_t* d);

#endif
