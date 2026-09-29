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
 * @brief Self-test of the polynomial operations.
 *
 * Builds \f$p(x) = 2x^2 + 3x + 1\f$ and \f$q(x) = 4x + 5\f$, prints them, and
 * prints \f$p + q\f$ and \f$p \cdot q\f$ (schoolbook and u64-NTT). Intended for
 * manual verification, not part of the library API.
 *
 * Complexity:
 *   - Time: \f$O(1)\f$ (fixed-size inputs)
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(1)\f$
 *
 * @return RABIN_SUCCESS on success, or RABIN_ERR_OVERFLOW /
 * RABIN_ERR_OUT_OF_MEMORY.
 */
rabin_err_t rpol_test();

#endif
