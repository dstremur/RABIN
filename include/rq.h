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

/**
 * @brief Set a rational to an integer value: \f$r = a\f$.
 *
 * The denominator is set to \f$1\f$, so the result is normalized.
 *
 * Complexity:
 *   - Time: \f$O(n)\f$ where \f$n =\f$ the size of \p a in limbs
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(n)\f$ limbs
 *
 * @param[out] r Rational to set.
 * @param[in]  a Integer value.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR /
 * RABIN_ERR_OUT_OF_MEMORY.
 */
rabin_err_t rq_set(rq_t* r, const rz_t* a);

/**
 * @brief Test whether a rational is zero.
 *
 * Complexity:
 *   - Time: \f$O(1)\f$
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(1)\f$
 *
 * @param[in] r Rational to test.
 * @return true if \f$r = 0\f$, false otherwise.
 */
bool rq_is_zero(const rq_t* r);

/**
 * @brief Compare two rationals.
 *
 * Both values are compared via cross-multiplication
 * \f$a_{\text{num}} \cdot b_{\text{den}}\f$ vs \f$b_{\text{num}} \cdot
 * a_{\text{den}}\f$; each denominator must be nonzero (any value
 * produced by the arithmetic operations above satisfies the \f$den >
 * 0\f$ invariant).
 *
 * Complexity:
 *   - Time: \f$O(n \cdot m)\f$ where \f$n, m\f$ are the limb sizes
 *   - Auxiliary memory: \f$O(n + m)\f$ limbs
 *   - Output memory: \f$O(1)\f$
 *
 * @param[in] a First rational.
 * @param[in] b Second rational.
 * @return -1, 0, or 1 if \f$a < b\f$, \f$a = b\f$, or \f$a > b\f$.
 */
int rq_cmp(const rq_t* a, const rq_t* b);

/**
 * @brief Copy a rational: \f$r = a\f$.
 *
 * The result is normalized, so copying a value in the empty (zero)
 * state yields the normalized zero \f$0/1\f$.
 *
 * Complexity:
 *   - Time: \f$O(n^2 \log n)\f$ where \f$n =\f$ the size of \p a in limbs
 *   - Auxiliary memory: \f$O(n)\f$ limbs
 *   - Output memory: \f$O(n)\f$ limbs
 *
 * @param[out] r Destination. May alias \p a.
 * @param[in]  a Source rational.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR,
 * RABIN_ERR_INVALID_ARG, or RABIN_ERR_OUT_OF_MEMORY.
 */
rabin_err_t rq_copy(rq_t* r, const rq_t* a);

/**
 * @brief Square a rational: \f$r = a^2\f$.
 *
 * Complexity:
 *   - Time: \f$O(n^2 \log n)\f$ where \f$n =\f$ the size of \p a in limbs
 *   - Auxiliary memory: \f$O(n)\f$ limbs
 *   - Output memory: \f$O(n)\f$ limbs
 *
 * @param[out] r Result. May alias \p a.
 * @param[in]  a Rational to square.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR or
 * RABIN_ERR_OUT_OF_MEMORY.
 */
rabin_err_t rq_sqr(rq_t* r, const rq_t* a);

/**
 * @brief Set a rational to an unsigned integer value: \f$r = a\f$.
 *
 * The denominator is set to \f$1\f$, so the result is normalized.
 *
 * Complexity:
 *   - Time: \f$O(1)\f$
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(1)\f$ limbs
 *
 * @param[out] r Rational to set.
 * @param[in]  a Unsigned integer value.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR or
 * RABIN_ERR_OUT_OF_MEMORY.
 */
rabin_err_t rq_set_u64(rq_t* r, u64 a);

/**
 * @brief Parse a decimal rational from a string: \f$p\f$ or \f$p/q\f$.
 *
 * Accepted grammar (whitespace only at the ends of the string):
 * \f$[\pm]\f$ digits, optionally followed by \f$/q\f$ with \f$q > 0\f$.
 * Examples: \c "3/4", \c "-7/2", \c "+5", \c "0/9". The result is
 * normalized, so \c "3/6" yields \f$1/2\f$.
 *
 * Complexity:
 *   - Time: \f$O(n^2)\f$ where \f$n =\f$ the number of digits
 *   - Auxiliary memory: \f$O(n)\f$ limbs
 *   - Output memory: \f$O(n)\f$ limbs
 *
 * @param[out] r Rational to set. Must be initialized with rq_init();
 *              left at \f$0/0\f$ on parse failure.
 * @param[in]  str String to parse.
 * @return RABIN_SUCCESS on success; RABIN_ERR_NULL_PTR if \p r or
 * \p str is NULL; RABIN_ERR_INVALID_ARG if \p str is malformed;
 * RABIN_ERR_DIV_BY_ZERO if the denominator is zero;
 * RABIN_ERR_OUT_OF_MEMORY on allocation failure.
 * @see rq_to_string(), rq_init()
 */
rabin_err_t rq_init_val(rq_t* r, const char* str);

/**
 * @brief Write a rational as a decimal string \f$p/q\f$ (omitting
 * \f$/q\f$ when the denominator is \f$1\f$) into a caller buffer.
 *
 * The output is exactly NUL-terminated and needs rq_to_string_len()
 * bytes including the terminator; a smaller \p cap leaves \p out
 * untouched and fails.
 *
 * Complexity:
 *   - Time: \f$O(n^2)\f$ decimal conversion, see rz_to_string()
 *   - Auxiliary memory: \f$O(n)\f$ limbs
 *   - Output memory: \f$O(n)\f$ characters
 *
 * @param[out] out Caller buffer receiving the string.
 * @param[in]  r   Rational to print. Expected normalized (see the
 *                 struct invariants); unnormalized values are printed
 *                 as stored.
 * @param[in]  cap Size of \p out in bytes, NUL included.
 * @return RABIN_SUCCESS on success; RABIN_ERR_NULL_PTR if \p out or
 * \p r is NULL; RABIN_ERR_INVALID_ARG if \p cap is too small;
 * RABIN_ERR_OUT_OF_MEMORY on allocation failure.
 * @see rq_to_string_len(), rq_init_val()
 */
rabin_err_t rq_to_string(const rq_t* r, char* out, size_t cap);

/**
 * @brief Number of bytes rq_to_string() needs for \p r, NUL included.
 *
 * Complexity:
 *   - Time: \f$O(n^2)\f$ decimal conversion, see rz_to_string()
 *   - Auxiliary memory: \f$O(n)\f$ limbs
 *
 * @param[in] r Rational to measure.
 * @return Required buffer size including the NUL terminator, or 0 if
 * \p r is NULL or conversion fails.
 * @see rq_to_string()
 */
size_t rq_to_string_len(const rq_t* r);

/**
 * @brief Set a rational to a signed integer value: \f$r = a\f$.
 *
 * The denominator is set to \f$1\f$, so the result is normalized.
 *
 * Complexity:
 *   - Time: \f$O(1)\f$
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(1)\f$ limbs
 *
 * @param[out] r Rational to set.
 * @param[in]  a Signed integer value.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR or
 * RABIN_ERR_OUT_OF_MEMORY.
 */
rabin_err_t rq_set_i64(rq_t* r, i64 a);

/**
 * @brief Reciprocal: \f$r = 1/a\f$.
 *
 * The result is normalized: a negative \f$a_\text{num}\f$ lands in the
 * numerator, keeping \f$den > 0\f$.
 *
 * Complexity:
 *   - Time: \f$O(n^2 \log n)\f$ where \f$n =\f$ the size of \p a in limbs
 *   - Auxiliary memory: \f$O(n)\f$ limbs
 *   - Output memory: \f$O(n)\f$ limbs
 *
 * @param[out] r Result. May alias \p a.
 * @param[in]  a Rational to invert.
 * @return RABIN_SUCCESS on success; RABIN_ERR_NULL_PTR if \p r or
 * \p a is NULL; RABIN_ERR_DIV_BY_ZERO if \p a is zero;
 * RABIN_ERR_INVALID_ARG if \p a has a zero denominator;
 * RABIN_ERR_OUT_OF_MEMORY on allocation failure.
 * @see rq_div()
 */
rabin_err_t rq_inv(rq_t* r, const rq_t* a);

/**
 * @brief Copy the numerator of \p r into \p out.
 *
 * Complexity:
 *   - Time: \f$O(n)\f$
 *   - Output memory: \f$O(n)\f$ limbs
 *
 * @param[out] out Destination; must be initialized.
 * @param[in]  r   Source rational.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR or
 * RABIN_ERR_OUT_OF_MEMORY.
 * @see rq_get_den()
 */
rabin_err_t rq_get_num(rz_t* out, const rq_t* r);

/**
 * @brief Copy the denominator of \p r into \p out.
 *
 * Complexity:
 *   - Time: \f$O(n)\f$
 *   - Output memory: \f$O(n)\f$ limbs
 *
 * @param[out] out Destination; must be initialized.
 * @param[in]  r   Source rational.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR or
 * RABIN_ERR_OUT_OF_MEMORY.
 * @see rq_get_num()
 */
rabin_err_t rq_get_den(rz_t* out, const rq_t* r);

/**
 * @brief Compare a normalized rational \p r with a machine integer \p a.
 *
 * Compares \f$r\f$ and \f$a\f$ as \f$\mathrm{num}(r)\f$ against
 * \f$a \cdot \mathrm{den}(r)\f$, avoiding any rounding.
 *
 * Complexity:
 *   - Time: \f$O(n)\f$ where \f$n =\f$ the size of \p r in limbs
 *   - Auxiliary memory: \f$O(n)\f$ limbs
 *
 * @param[in] r Rational operand. Expected normalized.
 * @param[in] a Integer operand.
 * @return Negative, zero, or positive if \f$r < a\f$, \f$r = a\f$, or
 * \f$r > a\f$; 0 if \p r is NULL or allocation fails.
 * @see rq_cmp()
 */
int rq_cmp_i64(const rq_t* r, i64 a);

#endif
