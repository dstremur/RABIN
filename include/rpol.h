#ifndef RPOL_H
#define RPOL_H

#include "rabin_errors.h"
#include "rz.h"

typedef struct {
  rz_t* coeff;
  u64 size;
  u64 deg;
} rpol_t;

/**
 * @brief Initialize a polynomial to the empty state (no coefficients).
 *
 * Complexity:
 *   - Time: \f$O(1)\f$
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(1)\f$
 *
 * @param[out] p Polynomial to initialize.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR /
 * RABIN_ERR_OUT_OF_MEMORY.
 */
rabin_err_t rpol_init(rpol_t* p);

/**
 * @brief Set a polynomial from an array of rz_t coefficients.
 *
 * Let \f$d =\f$ deg.
 *
 * Allocates \f$d + 1\f$ coefficient slots and deep-copies \f$coeff[0 \ldots
 * d]\f$ into the polynomial, setting the degree to \f$d\f$.
 *
 * Complexity:
 *   - Time: \f$O(d \cdot n)\f$ where \f$n =\f$ the size of the coefficients in
 * limbs
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(d)\f$ bignums
 *
 * @param[out] p Polynomial to set.
 * @param[in]  coeff Array of coefficients.
 * @param[in]  deg Degree of the polynomial.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR /
 * RABIN_ERR_OUT_OF_MEMORY.
 */
rabin_err_t rpol_set(rpol_t* p, const rz_t* coeff, u64 deg);

/**
 * @brief Set a polynomial from an array of 64-bit signed coefficients.
 *
 * Let \f$d =\f$ deg.
 *
 * Allocates coefficient slots and stores \f$coeff[0 \ldots d]\f$ as single-limb
 * bignums, then trims leading zero coefficients.
 *
 * Complexity:
 *   - Time: \f$O(d)\f$
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(d)\f$ bignums
 *
 * @param[out] p Polynomial to set.
 * @param[in]  coeff Array of 64-bit signed coefficients.
 * @param[in]  deg Degree of the polynomial.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR /
 * RABIN_ERR_OUT_OF_MEMORY.
 */
rabin_err_t rpol_set_i64(rpol_t* p, const i64* coeff, u64 deg);

/**
 * @brief Copy the coefficients of \f$q\f$ into \f$p\f$ (up to q->deg).
 *
 * Let \f$d =\f$ q->deg.
 *
 * \f$p\f$ must have capacity for at least \f$d + 1\f$ coefficients. The degree
 * of
 * \f$p\f$ is left unchanged.
 *
 * Complexity:
 *   - Time: \f$O(d \cdot n)\f$ where \f$n =\f$ the size of the coefficients in
 * limbs
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(d)\f$ bignums
 *
 * @param[out] p Destination polynomial.
 * @param[in]  q Source polynomial.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR, RABIN_ERR_OVERFLOW,
 * or RABIN_ERR_OUT_OF_MEMORY.
 */
rabin_err_t rpol_copy(rpol_t* p, const rpol_t* q);

/**
 * @brief Test whether two polynomials are equal.
 *
 * Let \f$d =\f$ max(a->deg, b->deg).
 *
 * Returns true iff the degrees match and all coefficients up to the
 * degree are equal (deep rz_t comparison).
 *
 * Complexity:
 *   - Time: \f$O(d \cdot n)\f$ where \f$n =\f$ the size of the coefficients in
 * limbs
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(1)\f$
 *
 * @param[in] a First polynomial.
 * @param[in] b Second polynomial.
 *
 * @return true  If the polynomials are equal.
 * @return false If the polynomials differ.
 */
bool rpol_equal(const rpol_t* a, const rpol_t* b);

/**
 * @brief Free all coefficient storage of a polynomial.
 *
 * Complexity:
 *   - Time: \f$O(size)\f$
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(1)\f$
 *
 * @param[in,out] p Polynomial to free.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR.
 */
rabin_err_t rpol_clear(rpol_t* p);

/**
 * @brief Grow the coefficient array of a polynomial to at least deg + 1
 * slots.
 *
 * If the capacity is already sufficient, nothing happens. Otherwise
 * the array is reallocated with exponential growth and the new slots
 * are initialized to zero bignums.
 *
 * Complexity:
 *   - Time: \f$O(new_cap)\f$ for the initialization of the new slots
 *   - Auxiliary memory: \f$O(new_cap)\f$ during the realloc
 *   - Output memory: \f$O(new_cap)\f$ bignums
 *
 * @param[in,out] p   Polynomial to grow.
 * @param[in]  deg Minimum number of coefficient slots minus one.
 *
 * @return RABIN_SUCCESS on success, or RABIN_ERR_OUT_OF_MEMORY.
 */
rabin_err_t rpol_alloc(rpol_t* p, u64 deg);

/**
 * @brief Remove leading zero coefficients from a polynomial.
 *
 * Lowers p->deg until the top coefficient is nonzero (or the degree
 * is 0). The capacity is left unchanged.
 *
 * Complexity:
 *   - Time: \f$O(number of trimmed coefficients)\f$
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(1)\f$
 *
 * @param[in,out] p Polynomial to trim.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR.
 */
rabin_err_t rpol_trim(rpol_t* p);

/**
 * @brief Print a polynomial to stdout in descending powers of \f$x\f$.
 *
 * Complexity:
 *   - Time: \f$O(d \cdot n)\f$ where \f$d =\f$ the degree and \f$n =\f$ the
 * coefficient size
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(d \cdot n)\f$ characters written
 *
 * @param[in] p Polynomial to print.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_OUT_OF_MEMORY.
 */
rabin_err_t rpol_print(const rpol_t* p);

/**
 * @brief Add two polynomials: \f$r = p + q\f$.
 *
 * Let \f$d =\f$ max(p->deg, q->deg).
 *
 * Coefficients up to the smaller degree are added with rz_add(); the
 * remaining coefficients of the longer polynomial are copied. The
 * result is trimmed.
 *
 * Complexity:
 *   - Time: \f$O(d \cdot n)\f$ where \f$n =\f$ the size of the coefficients in
 * limbs
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(d)\f$ bignums
 *
 * @param[out] r Result polynomial.
 * @param[in]  p First polynomial.
 * @param[in]  q Second polynomial.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR /
 * RABIN_ERR_OUT_OF_MEMORY.
 */
rabin_err_t rpol_add(rpol_t* r, const rpol_t* p, const rpol_t* q);

/**
 * @brief Subtract two polynomials: \f$r = p - q\f$.
 *
 * Let \f$d =\f$ max(p->deg, q->deg).
 *
 * Coefficients up to the smaller degree are subtracted with rz_sub();
 * the remaining coefficients of the longer polynomial are copied with
 * their sign flipped. The result is trimmed.
 *
 * Complexity:
 *   - Time: \f$O(d \cdot n)\f$ where \f$n =\f$ the size of the coefficients in
 * limbs
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(d)\f$ bignums
 *
 * @param[out] r Result polynomial.
 * @param[in]  p First polynomial.
 * @param[in]  q Second polynomial.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR /
 * RABIN_ERR_OUT_OF_MEMORY.
 */
rabin_err_t rpol_sub(rpol_t* r, const rpol_t* p, const rpol_t* q);

/**
 * @brief Multiply two polynomials with the rz_t NTT: \f$r = p \cdot q\f$.
 *
 * Let \f$d =\f$ p->deg + q->deg + 1 (length of the true product) and
 * \f$n = 2^k\f$ the next power of two \f$\ge d\f$.
 *
 * Pads both polynomials to length \f$n\f$, applies the forward rz_t NTT
 * (Goldilocks field) to each, multiplies pointwise in the Montgomery
 * domain, and applies the inverse NTT. The first \f$d\f$ coefficients of the
 * cyclic convolution equal the true product coefficients (no
 * wrap-around, since the Goldilocks prime is large enough for the
 * coefficient sizes used).
 *
 * A fresh NTT context is initialized per call.
 *
 * Complexity:
 *   - Time: \f$O(n \log n \cdot k_{l}^2)\f$ where \f$k_{l} =\f$ the size of the
 * modulus in limbs, plus \f$O(n \cdot k_{l}^2)\f$ for the context tables
 *   - Auxiliary memory: \f$O(n)\f$ bignums for the transformed copies
 *   - Output memory: \f$O(d)\f$ bignums
 *
 * @param[out] r Result polynomial.
 * @param[in]  p First polynomial.
 * @param[in]  q Second polynomial.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR, RABIN_ERR_OVERFLOW,
 * or RABIN_ERR_OUT_OF_MEMORY.
 *
 * @par Algorithm Reference:
 * C. M. Cooley and J. W. Tukey, "An Algorithm for the Machine
 * Calculation of Complex Fourier Series," Mathematics of Computation,
 * vol. 19, no. 90, 1965.
 * @see rpol_mul_ntt_u64(), rpol_mul_school()
 */
rabin_err_t rpol_mul_ntt(rpol_t* r, const rpol_t* p, const rpol_t* q);

/**
 * @brief Multiply two polynomials with the u64 NTT: \f$r = p \cdot q\f$.
 *
 * Let \f$d =\f$ p->deg + q->deg + 1 (length of the true product) and
 * \f$n = 2^k\f$ the next power of two \f$\ge d\f$.
 *
 * Only works for coefficients that fit in a single u64 limb: each
 * coefficient is reduced to its low limb, the convolution is computed
 * with the fast u64 Goldilocks NTT (flat arrays, Montgomery
 * butterflies), and the result limbs are lifted back to bignums.
 *
 * The Goldilocks NTT context for the transform size is taken from a
 * shared per-size cache (built once per k, see
 * u64_ntt_ctx_golden_cached()), and the flat u64 arrays come from the
 * thread-local scratch arena.
 *
 * Complexity:
 *   - Time: \f$O(n \log n)\f$ for the NTTs
 *   - Auxiliary memory: \f$O(n)\f$ u64s for the flat arrays (scratch arena)
 *   - Output memory: \f$O(d)\f$ bignums
 *
 * @param[out] r Result polynomial.
 * @param[in]  p First polynomial.
 * @param[in]  q Second polynomial.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR, RABIN_ERR_OVERFLOW,
 * or RABIN_ERR_OUT_OF_MEMORY.
 *
 * @par Algorithm Reference:
 * C. M. Cooley and J. W. Tukey, "An Algorithm for the Machine
 * Calculation of Complex Fourier Series," Mathematics of Computation,
 * vol. 19, no. 90, 1965.
 * @see rpol_mul_ntt(), rpol_mul_school()
 */
rabin_err_t rpol_mul_ntt_u64(rpol_t* r, const rpol_t* p, const rpol_t* q);

/**
 * @brief Multiply two polynomials with the schoolbook algorithm: \f$r = p \cdot
 * q\f$.
 *
 * Let \f$d_p =\f$ p->deg, \f$d_q =\f$ q->deg, and \f$d = d_p + d_q\f$.
 *
 * Each pair of coefficients is multiplied and accumulated into the
 * corresponding slot of a temporary polynomial, which is then copied
 * into r and trimmed.
 *
 * Complexity:
 *   - Time: \f$O((d_p + 1) \cdot (d_q + 1) \cdot n^2)\f$ where \f$n =\f$ the
 * size of the coefficients in limbs
 *   - Auxiliary memory: \f$O(d)\f$ bignums for the temporary
 *   - Output memory: \f$O(d)\f$ bignums
 *
 * @param[out] r Result polynomial.
 * @param[in]  p First polynomial.
 * @param[in]  q Second polynomial.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR, RABIN_ERR_OVERFLOW,
 * or RABIN_ERR_OUT_OF_MEMORY.
 */
rabin_err_t rpol_mul_school(rpol_t* r, const rpol_t* p, const rpol_t* q);

/**
 * @brief Multiply two polynomials: \f$r = p \cdot q\f$.
 *
 * Currently a thin wrapper around rpol_mul_school().
 *
 * Complexity:
 *   - Time: see rpol_mul_school()
 *   - Auxiliary memory: \f$O(d)\f$ bignums
 *   - Output memory: \f$O(d)\f$ bignums
 *
 * @param[out] r Result polynomial.
 * @param[in]  p First polynomial.
 * @param[in]  q Second polynomial.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR, RABIN_ERR_OVERFLOW,
 * or RABIN_ERR_OUT_OF_MEMORY.
 */
rabin_err_t rpol_mul(rpol_t* r, const rpol_t* p, const rpol_t* q);

/**
 * @brief Pseudo-division over \f$\Z[x]\f$: \f$\delta \cdot a = q \cdot b +
 * r\f$.
 *
 * Let \f$d_a = \f$ a->deg, \f$d_b = \f$ b->deg, and let the scaling factor be
 * \f$\delta = lc(b)^{d_a - d_b + 1}\f$ (with \f$\delta = 1\f$ when
 * \f$d_a < d_b\f$).
 *
 * Multiplying \f$a\f$ by \f$\delta\f$ makes every leading-coefficient
 * division of ordinary long division exact, so \f$q\f$ and \f$r\f$ stay in
 * \f$\Z[x]\f$ without ever introducing rationals. The remainder is zero or
 * has degree \f$< d_b\f$. When \f$a = 0\f$ or \f$d_a < d_b\f$ the result is
 * \f$q = 0\f$, \f$r = a\f$.
 *
 * Complexity:
 *   - Time: \f$O((d_a - d_b + 1) \cdot d_b \cdot n^2)\f$ where \f$n =\f$ the
 * size of the coefficients in limbs (coefficients grow by up to
 * \f$(d_a - d_b + 1) \cdot \log_2 |lc(b)|\f$ bits)
 *   - Auxiliary memory: \f$O(d_a)\f$ bignums for the local copies
 *   - Output memory: \f$O(d_a)\f$ bignums
 *
 * @param[out] q Pseudo-quotient (degree \f$\le d_a - d_b\f$).
 * @param[out] r Pseudo-remainder (zero, or degree \f$< d_b\f$).
 * @param[in]  a Dividend.
 * @param[in]  b Divisor (nonzero).
 * @pre q and r must be initialized (rpol_init()) and distinct.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR,
 * RABIN_ERR_DIV_BY_ZERO, or RABIN_ERR_OUT_OF_MEMORY. On error q and r are
 * left unchanged.
 *
 * @par Algorithm Reference:
 * D. E. Knuth, The Art of Computer Programming, Vol. 2: Seminumerical
 * Algorithms, 3rd ed., Addison-Wesley, 1997, section 4.6.1,
 * Algorithm P (pseudodivision).
 * @see rpol_div_exact(), rpol_gcd()
 */
rabin_err_t rpol_pseudo_div(rpol_t* q, rpol_t* r, const rpol_t* a,
                            const rpol_t* b);

/**
 * @brief Exact division over \f$\Z[x]\f$: \f$q = a / b\f$ (requires
 * \f$b \mid a\f$).
 *
 * Runs ordinary long division over \f$\Z\f$ and aborts at the first
 * leading coefficient not divisible by \f$lc(b)\f$. Fails with
 * RABIN_ERR_INVALID_ARG iff no \f$q \in \Z[x]\f$ satisfies \f$a = q \cdot
 * b\f$ (in particular when the remainder is nonzero).
 *
 * Complexity:
 *   - Time: \f$O((d_a - d_b + 1) \cdot d_b \cdot n^2)\f$ where \f$n =\f$ the
 * size of the coefficients in limbs
 *   - Auxiliary memory: \f$O(d_a)\f$ bignums for the local copies
 *   - Output memory: \f$O(d_a - d_b + 1)\f$ bignums
 *
 * @pre q must be initialized (rpol_init()).
 * @param[out] q Quotient \f$a / b\f$.
 * @param[in]  a Dividend.
 * @param[in]  b Divisor (nonzero, must divide \f$a\f$ exactly).
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR,
 * RABIN_ERR_DIV_BY_ZERO, RABIN_ERR_INVALID_ARG, or
 * RABIN_ERR_OUT_OF_MEMORY. On error q is left unchanged.
 * @see rpol_pseudo_div(), rpol_mul()
 */
rabin_err_t rpol_div_exact(rpol_t* q, const rpol_t* a, const rpol_t* b);

/**
 * @brief Content of \f$p\f$: the nonnegative gcd of its coefficients.
 *
 * The zero polynomial has content 0.
 *
 * Complexity:
 *   - Time: \f$O(d \cdot \log M)\f$ gcd steps where \f$d = \f$ p->deg and
 * \f$M =\f$ the largest coefficient
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: one bignum
 *
 * @param[out] g Content \f$\mathrm{cont}(p)\f$.
 * @param[in]  p Polynomial.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR /
 * RABIN_ERR_OUT_OF_MEMORY.
 * @see rpol_primitive_part(), rpol_gcd()
 */
rabin_err_t rpol_content(rz_t* g, const rpol_t* p);

/**
 * @brief Primitive part of \f$p\f$: \f$p / \mathrm{cont}(p)\f$ normalized to
 * a positive leading coefficient.
 *
 * The zero polynomial maps to itself. For nonzero \f$p\f$ the result has
 * content 1 and \f$lc > 0\f$, and satisfies
 * \f$p = \pm \mathrm{cont}(p) \cdot \mathrm{pp}(p)\f$.
 *
 * Complexity:
 *   - Time: see rpol_content(), plus \f$O(d \cdot n^2)\f$ divisions
 *   - Auxiliary memory: \f$O(d)\f$ bignums for the local copy
 *   - Output memory: \f$O(d)\f$ bignums
 *
 * @pre pp must be initialized (rpol_init()).
 * @param[out] pp Primitive part \f$\mathrm{pp}(p)\f$.
 * @param[in]  p Polynomial.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR /
 * RABIN_ERR_OUT_OF_MEMORY.
 * @see rpol_content()
 */
rabin_err_t rpol_primitive_part(rpol_t* pp, const rpol_t* p);

/**
 * @brief Greatest common divisor over \f$\Z[x]\f$, normalized to a positive
 * leading coefficient.
 *
 * Let \f$d = \f$ max(a->deg, b->deg).
 *
 * By Gauss's lemma \f$\gcd(a, b) = \gcd(\mathrm{cont}\, a,
 * \mathrm{cont}\, b) \cdot \gcd(\mathrm{pp}\, a, \mathrm{pp}\, b)\f$;
 * the primitive part is computed with a reduced pseudo-remainder sequence
 * (pseudo-division via rpol_pseudo_div(), with each remainder's content
 * stripped to keep coefficients small). Conventions: \f$\gcd(0, 0) =
 * 0\f$ and \f$\gcd(0, b) = \pm b\f$ normalized to \f$lc > 0\f$.
 *
 * Complexity:
 *   - Time: \f$O(d^2)\f$ pseudo-division steps, each polynomial in the
 * coefficient sizes
 *   - Auxiliary memory: \f$O(d)\f$ bignums for the local copies
 *   - Output memory: \f$O(d)\f$ bignums
 *
 * @pre g must be initialized (rpol_init()).
 * @param[out] g Greatest common divisor (zero, or leading coefficient > 0).
 * @param[in]  a First polynomial.
 * @param[in]  b Second polynomial.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR /
 * RABIN_ERR_OUT_OF_MEMORY.
 *
 * @par Algorithm Reference:
 * D. E. Knuth, The Art of Computer Programming, Vol. 2: Seminumerical
 * Algorithms, 3rd ed., Addison-Wesley, 1997, section 4.6.1 (gcd in
 * \f$\Z[x]\f$); H. Cohen, A Course in Computational Algebraic Number
 * Theory, Springer, 1993, Algorithm 3.5.4 (reduced PRS).
 * @warning The reduced PRS has slightly worse coefficient growth than the
 * subresultant PRS (Cohen, Algorithm 3.5.9); swap it in if this becomes a
 * hot path.
 * @see rpol_pseudo_div(), rpol_content(), rpol_primitive_part()
 */
rabin_err_t rpol_gcd(rpol_t* g, const rpol_t* a, const rpol_t* b);

/**
 * @brief Evaluate \f$p\f$ at \f$x\f$ with Horner's scheme: \f$r = p(x)\f$.
 *
 * Let \f$d = \f$ p->deg. Evaluates
 * \f$r = (\ldots((c_d x + c_{d-1}) x + \ldots) x + c_0\f$).
 *
 * Complexity:
 *   - Time: \f$O(d \cdot n^2)\f$ where \f$n =\f$ the size of the
 * accumulated result in limbs
 *   - Auxiliary memory: \f$O(1)\f$ bignums for the accumulator
 *   - Output memory: one bignum
 *
 * @param[out] r Result \f$p(x)\f$.
 * @param[in]  p Polynomial.
 * @param[in]  x Evaluation point.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR /
 * RABIN_ERR_OUT_OF_MEMORY.
 * @see rpol_derive()
 */
rabin_err_t rpol_eval(rz_t* r, const rpol_t* p, const rz_t* x);

/**
 * @brief Formal derivative: \f$\frac{d}{dx} \sum_i c_i x^i = \sum_i i
 * \cdot c_i x^{i-1}\f$.
 *
 * Let \f$d = \f$ p->deg. The derivative of a constant (and of the zero
 * polynomial) is the zero polynomial.
 *
 * Complexity:
 *   - Time: \f$O(d \cdot n^2)\f$ where \f$n =\f$ the size of the
 * coefficients in limbs
 *   - Auxiliary memory: \f$O(d)\f$ bignums for the local copy
 *   - Output memory: \f$O(d)\f$ bignums
 *
 * @pre dp must be initialized (rpol_init()).
 * @param[out] dp Derivative polynomial.
 * @param[in]  p Polynomial.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR /
 * RABIN_ERR_OUT_OF_MEMORY.
 * @see rpol_eval()
 */
rabin_err_t rpol_derive(rpol_t* dp, const rpol_t* p);

/**
 * @brief Decompose a rz_t into a polynomial of fixed-width chunks.
 *
 * Let \f$n =\f$ n->size, measured in 64-bit limbs, and \f$W =\f$ width.
 *
 * Slices \f$n\f$ into chunks of \f$W\f$ bits (base \f$2^W\f$) and stores them
 * as the coefficients of \f$r\f$ in ascending order:
 *
 *   \f$n = \sum_i r->coeff[i] \cdot (2^W)^i\f$
 *
 * The degree is \f$\lceil n / W \rceil - 1\f$ (\f$0\f$ for \f$n = 0\f$).
 *
 * Complexity:
 *   - Time: \f$O(n \cdot \lceil n/W \rceil)\f$ - one \f$O(n)\f$ shift per chunk
 *   - Auxiliary memory: \f$O(n)\f$ limbs for temporaries
 *   - Output memory: \f$O(n/W)\f$ bignums
 *
 * @param[out]    r    Result polynomial storing the chunks.
 * @param[in]     n    Bignum to decompose.
 * @param[in] width Chunk width in bits.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR,
 * RABIN_ERR_INVALID_ARG, or RABIN_ERR_OUT_OF_MEMORY.
 */
rabin_err_t rz_decompose(rpol_t* r, const rz_t* n, u64 width);

/**
 * @brief Propagate carries between the fixed-width coefficient slots of a
 * polynomial, in place.
 *
 * Let \f$d =\f$ r->deg and \f$W =\f$ bit_width.
 *
 * After an NTT-based convolution the coefficient slots may exceed
 * \f$2^W - 1\f$. This normalizes them: for each slot, \f$total = coeff +
 * carry\f$,
 * \f$coeff = total \bmod 2^W\f$, \f$carry = total / 2^W\f$, rippling upward and
 * growing the polynomial if the carry outlives the current degree.
 *
 * Complexity:
 *   - Time: \f$O(d \cdot n)\f$ where \f$n =\f$ the size of the coefficients in
 * limbs
 *   - Auxiliary memory: \f$O(n)\f$ limbs for temporaries
 *   - Output memory: \f$O(d)\f$ bignums (possibly grown by the final carry)
 *
 * @param[in,out]      r       Polynomial to normalize (modified in place).
 * @param[in] bit_width Width of each coefficient slot in bits.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR,
 * RABIN_ERR_INVALID_ARG, or RABIN_ERR_OUT_OF_MEMORY.
 */
rabin_err_t rpol_carry_propagation(rpol_t* r, u64 bit_width);

/**
 * @brief Recompose a polynomial of fixed-width chunks into a rz_t.
 *
 * Let \f$d =\f$ p->deg and \f$W =\f$ bit_width.
 *
 * The inverse of rz_decompose():
 *
 *   \f$n = \sum_i p->coeff[i] \cdot 2^{i \cdot W}\f$
 *
 * Complexity:
 *   - Time: \f$O(d \cdot n)\f$ where n is the size of the result in limbs
 *   - Auxiliary memory: \f$O(n)\f$ limbs for the running term
 *   - Output memory: \f$O(n)\f$ limbs
 *
 * @param[out]        n         Result rz_t.
 * @param[in]         p         Polynomial of chunks.
 * @param[in] bit_width Width of each chunk in bits.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR /
 * RABIN_ERR_OUT_OF_MEMORY.
 */
rabin_err_t rz_recompose(rz_t* n, const rpol_t* p, u64 bit_width);

/**
 * @brief Number of bytes rpol_to_string() needs for \p p, NUL included.
 *
 * Complexity:
 *   - Time: \f$O((d+1) \cdot n^2)\f$ decimal conversion, see rz_to_string()
 *   - Auxiliary memory: \f$O(n)\f$ limbs
 *
 * @param[in] p Polynomial to measure.
 * @return Required buffer size including the NUL terminator, or 0 if
 * \p p is NULL or conversion fails.
 * @see rpol_to_string()
 */
size_t rpol_to_string_len(const rpol_t* p);

/**
 * @brief Write a polynomial as a bracketed ascending coefficient list
 * \f$[c_0,\dots,c_d]\f$ into a caller buffer.
 *
 * Coefficients are printed from \f$c_0\f$ up to \f$c_d\f$ of the
 * stored (trimmed) degree; the empty state prints \c []. The output
 * is exactly NUL-terminated and needs rpol_to_string_len() bytes
 * including the terminator; a smaller \p cap leaves \p out untouched
 * and fails.
 *
 * Complexity:
 *   - Time: \f$O((d+1) \cdot n^2)\f$ decimal conversion, see rz_to_string()
 *   - Auxiliary memory: \f$O(n)\f$ limbs
 *   - Output memory: \f$O(d)\f$ characters
 *
 * @param[out] out Caller buffer receiving the string.
 * @param[in]  p   Polynomial to print.
 * @param[in]  cap Size of \p out in bytes, NUL included.
 * @return RABIN_SUCCESS on success; RABIN_ERR_NULL_PTR if \p out or
 * \p p is NULL; RABIN_ERR_INVALID_ARG if \p cap is too small;
 * RABIN_ERR_OUT_OF_MEMORY on allocation failure.
 * @see rpol_init_val()
 */
rabin_err_t rpol_to_string(const rpol_t* p, char* out, size_t cap);

/**
 * @brief Parse a bracketed ascending coefficient list \f$[c_0,\dots,c_d]\f$
 * into a polynomial.
 *
 * Coefficients are signed decimal integers; whitespace is allowed
 * between tokens. Trailing high-degree zeros are trimmed, e.g. \c
 * "[5,0,0]" stores degree \f$0\f$; \c [] yields the empty zero
 * polynomial. \p p must be initialized; a failed parse leaves \p p
 * unchanged. The trimmed result round-trips exactly with
 * rpol_to_string().
 *
 * Complexity:
 *   - Time: \f$O((d+1) \cdot n^2)\f$ where \f$n =\f$ the digit count
 *   - Auxiliary memory: \f$O((d+1) \cdot n)\f$ limbs
 *   - Output memory: \f$O((d+1) \cdot n)\f$ limbs
 *
 * @param[out] p   Polynomial to set; must be initialized.
 * @param[in]  str String to parse.
 * @return RABIN_SUCCESS on success; RABIN_ERR_NULL_PTR if \p p or
 * \p str is NULL; RABIN_ERR_INVALID_ARG if \p str is malformed;
 * RABIN_ERR_OUT_OF_MEMORY on allocation failure.
 * @see rpol_to_string()
 */
rabin_err_t rpol_init_val(rpol_t* p, const char* str);

#endif
