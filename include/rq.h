#ifndef RQ_H
#define RQ_H

#include "rabin.h"
#include "rabin_errors.h"

/**
 * @brief A rational number \f$r = num / den\f$.
 *
 * Invariants of a normalized rational:
 *
 *   - \f$den > 0\f$
 *   - \f$\gcd(num, den) = 1\f$
 *   - \f$den = 1\f$ when \f$num = 0\f$
 *
 * All arithmetic operations below normalize their result; use
 * rq_normalize() to canonicalize a manually constructed value.
 */
typedef struct rq_t {
  rz_t num;
  rz_t den;
} rq_t;

/**
 * @brief Initialize a rational to the empty (zero) state.
 *
 * Both \f$num\f$ and \f$den\f$ are set to zero; the value is not
 * normalized (the \f$den > 0\f$ invariant holds only after rq_normalize()
 * or any arithmetic operation).
 *
 * Complexity:
 *   - Time: \f$O(1)\f$
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(1)\f$
 *
 * @param[out] r Rational to initialize.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR.
 */
rabin_err_t rq_init(rq_t* r);

/**
 * @brief Initialize a NULL-terminated list of rationals.
 *
 * @param[out] r First rational to initialize.
 * @param      ... Further rationals, terminated by a NULL pointer.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR.
 */
rabin_err_t rq_init_multi(rq_t* r, ...);

/**
 * @brief Free the storage of a rational.
 *
 * Complexity:
 *   - Time: \f$O(1)\f$
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(1)\f$
 *
 * @param[in,out] r Rational to free.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR.
 */
rabin_err_t rq_clear(rq_t* r);

/**
 * @brief Free a NULL-terminated list of rationals.
 *
 * @param[in,out] r First rational to free.
 * @param      ... Further rationals, terminated by a NULL pointer.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR.
 */
rabin_err_t rq_clear_multi(rq_t* r, ...);

/**
 * @brief Sum of two rationals: \f$r = a + b\f$.
 *
 * Let \f$n =\f$ the size of the numerators/denominators in 64-bit limbs.
 *
 * Computes \f$(a_{num} b_{den} + b_{num} a_{den}) / (a_{den} b_{den})\f$
 * and normalizes the result.
 *
 * Complexity:
 *   - Time: \f$O(n^2 \log n)\f$ (two multiplications plus normalization)
 *   - Auxiliary memory: \f$O(n)\f$ limbs
 *   - Output memory: \f$O(n)\f$ limbs
 *
 * @param[out] r Result. May alias \p a or \p b.
 * @param[in] a First operand.
 * @param[in] b Second operand.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR /
 * RABIN_ERR_OUT_OF_MEMORY.
 */
rabin_err_t rq_add(rq_t* r, const rq_t* a, const rq_t* b);

/**
 * @brief Difference of two rationals: \f$r = a - b\f$.
 *
 * Let \f$n =\f$ the size of the numerators/denominators in 64-bit limbs.
 *
 * Computes \f$(a_{num} b_{den} - b_{num} a_{den}) / (a_{den} b_{den})\f$
 * and normalizes the result.
 *
 * Complexity:
 *   - Time: \f$O(n^2 \log n)\f$ (two multiplications plus normalization)
 *   - Auxiliary memory: \f$O(n)\f$ limbs
 *   - Output memory: \f$O(n)\f$ limbs
 *
 * @param[out] r Result. May alias \p a or \p b.
 * @param[in] a First operand.
 * @param[in] b Second operand.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR /
 * RABIN_ERR_OUT_OF_MEMORY.
 */
rabin_err_t rq_sub(rq_t* r, const rq_t* a, const rq_t* b);

/**
 * @brief Product of two rationals: \f$r = a \cdot b\f$.
 *
 * Let \f$n =\f$ the size of the numerators/denominators in 64-bit limbs.
 *
 * Cross-cancels \f$\gcd(a_{num}, b_{den})\f$ and \f$\gcd(a_{den},
 * b_{num})\f$ before multiplying to keep the intermediates small, then
 * normalizes.
 *
 * Complexity:
 *   - Time: \f$O(n^2 \log n)\f$
 *   - Auxiliary memory: \f$O(n)\f$ limbs
 *   - Output memory: \f$O(n)\f$ limbs
 *
 * @param[out] r Result. May alias \p a or \p b.
 * @param[in] a First operand.
 * @param[in] b Second operand.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR /
 * RABIN_ERR_OUT_OF_MEMORY.
 */
rabin_err_t rq_mul(rq_t* r, const rq_t* a, const rq_t* b);

/**
 * @brief Quotient of two rationals: \f$r = a / b\f$.
 *
 * Let \f$n =\f$ the size of the numerators/denominators in 64-bit limbs.
 *
 * Computes \f$a \cdot b^{-1} = (a_{num} b_{den}) / (a_{den} b_{num})\f$,
 * cross-cancelling \f$\gcd(a_{num}, b_{num})\f$ and \f$\gcd(a_{den},
 * b_{den})\f$, then normalizes.
 *
 * Complexity:
 *   - Time: \f$O(n^2 \log n)\f$
 *   - Auxiliary memory: \f$O(n)\f$ limbs
 *   - Output memory: \f$O(n)\f$ limbs
 *
 * @param[out] r Result. May alias \p a or \p b.
 * @param[in] a Dividend.
 * @param[in] b Divisor; must be nonzero.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR,
 * RABIN_ERR_DIV_BY_ZERO, or RABIN_ERR_OUT_OF_MEMORY.
 */
rabin_err_t rq_div(rq_t* r, const rq_t* a, const rq_t* b);

/**
 * @brief Negate a rational in place: \f$r = -r\f$.
 *
 * Flips the sign of the numerator only.
 *
 * Complexity:
 *   - Time: \f$O(n)\f$ where \f$n =\f$ the numerator size in limbs
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(n)\f$ limbs
 *
 * @param[in,out] r Rational to negate.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR /
 * RABIN_ERR_OUT_OF_MEMORY.
 */
rabin_err_t rq_neg(rq_t* r);

/**
 * @brief Normalize a rational in place.
 *
 * Sets \f$den = 1\f$ when \f$num = 0\f$, moves a negative denominator to
 * the numerator, and divides both by their gcd, restoring the invariants
 * \f$den > 0\f$ and \f$\gcd(num, den) = 1\f$.
 *
 * Complexity:
 *   - Time: \f$O(n^2 \log n)\f$ where \f$n =\f$ the size in 64-bit limbs
 *   - Auxiliary memory: \f$O(n)\f$ limbs
 *   - Output memory: \f$O(n)\f$ limbs
 *
 * @param[in,out] r Rational to normalize.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR,
 * RABIN_ERR_INVALID_ARG (if \f$den = 0\f$), or RABIN_ERR_OUT_OF_MEMORY.
 */
rabin_err_t rq_normalize(rq_t* r);

/**
 * @brief Print a rational to stdout.
 *
 * Prints \f$0\f$ for the zero value, \f$num\f$ alone when \f$den = 1\f$,
 * and \f$num/den\f$ otherwise.
 *
 * Complexity:
 *   - Time: \f$O(n)\f$ where \f$n =\f$ the size in 64-bit limbs
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(n)\f$ characters written
 *
 * @param[in] r Rational to print.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR.
 */
rabin_err_t rq_print(const rq_t* r);

#endif
