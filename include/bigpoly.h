#ifndef BIGPOLY_H
#define BIGPOLY_H

#include "bigcore.h"

typedef struct {
  bignum* coeff;
  u64 size;
  u64 deg;
} bigpoly;

/**
 * @brief Initialize a polynomial to the empty state (no coefficients).
 *
 * Complexity:
 *   Time: O(1)
 *   Auxiliary memory: O(1)
 *   Output memory: O(1)
 *
 * @param[out] p Polynomial to initialize.
 */
void bigpoly_init(bigpoly* p);

/**
 * @brief Bitwise AND of two bignums: r = a & mask.
 *
 * Let n = a->size, measured in 64-bit limbs.
 *
 * Copies a into r and ANDs each limb with the corresponding limb of
 * mask (limbs beyond mask->size are zeroed), then trims. Used to
 * extract fixed-width chunks of a bignum.
 *
 * Complexity:
 *   Time: O(n)
 *   Auxiliary memory: O(1)
 *   Output memory: O(n) limbs
 *
 * @param[out]   r    Result storing a & mask.
 * @param[in]    a    First operand.
 * @param[in] mask Mask.
 */
void bn_and(bignum* r, const bignum* a, const bignum* mask);

/**
 * @brief Set a polynomial from an array of bignum coefficients.
 *
 * Let d = deg.
 *
 * Allocates d + 1 coefficient slots and deep-copies coeff[0..d] into
 * the polynomial, setting the degree to d.
 *
 * Complexity:
 *   Time: O(d * n) where n is the size of the coefficients in limbs
 *   Auxiliary memory: O(1)
 *   Output memory: O(d) bignums
 *
 * @param[out]    p     Polynomial to set.
 * @param[in] coeff Array of coefficients.
 * @param[in]   deg   Degree of the polynomial.
 */
void bigpoly_set(bigpoly* p, bignum* coeff, u64 deg);

/**
 * @brief Set a polynomial from an array of 64-bit signed coefficients.
 *
 * Let d = deg.
 *
 * Allocates coefficient slots and stores coeff[0..d] as single-limb
 * bignums, then trims leading zero coefficients.
 *
 * Complexity:
 *   Time: O(d)
 *   Auxiliary memory: O(1)
 *   Output memory: O(d) bignums
 *
 * @param[out]    p     Polynomial to set.
 * @param[in] coeff Array of 64-bit signed coefficients.
 * @param[in]   deg   Degree of the polynomial.
 */
void bigpoly_set_i64(bigpoly* p, i64* coeff, u64 deg);

/**
 * @brief Copy the coefficients of q into p (up to q->deg).
 *
 * Let d = q->deg.
 *
 * p must have capacity for at least d + 1 coefficients. The degree of
 * p is left unchanged.
 *
 * Complexity:
 *   Time: O(d * n) where n is the size of the coefficients in limbs
 *   Auxiliary memory: O(1)
 *   Output memory: O(d) bignums
 *
 * @param[out] p Destination polynomial.
 * @param[in]  q Source polynomial.
 */
void bigpoly_copy(bigpoly* p, bigpoly* q);

/**
 * @brief Test whether two polynomials are equal.
 *
 * Let d = max(a->deg, b->deg).
 *
 * Returns true iff the degrees match and all coefficients up to the
 * degree are equal (deep bignum comparison).
 *
 * Complexity:
 *   Time: O(d * n) where n is the size of the coefficients in limbs
 *   Auxiliary memory: O(1)
 *   Output memory: O(1)
 *
 * @param[in] a First polynomial.
 * @param[in] b Second polynomial.
 *
 * @return true  If the polynomials are equal.
 * @return false If the polynomials differ.
 */
bool bigpoly_equal(bigpoly* a, bigpoly* b);

/**
 * @brief Free all coefficient storage of a polynomial.
 *
 * Complexity:
 *   Time: O(size)
 *   Auxiliary memory: O(1)
 *   Output memory: O(1)
 *
 * @param[in,out] p Polynomial to free.
 */
void bigpoly_free(bigpoly* p);

/**
 * @brief Grow the coefficient array of a polynomial to at least deg + 1
 * slots.
 *
 * If the capacity is already sufficient, nothing happens. Otherwise
 * the array is reallocated with exponential growth and the new slots
 * are initialized to zero bignums.
 *
 * Returns true on success, false if the realloc fails.
 *
 * Complexity:
 *   Time: O(new_cap) for the initialization of the new slots
 *   Auxiliary memory: O(new_cap) during the realloc
 *   Output memory: O(new_cap) bignums
 *
 * @param[in,out] p   Polynomial to grow.
 * @param[in]  deg Minimum number of coefficient slots minus one.
 *
 * @return true  On success.
 * @return false If the realloc fails.
 */
bool bigpoly_alloc(bigpoly* p, u64 deg);

/**
 * @brief Remove leading zero coefficients from a polynomial.
 *
 * Lowers p->deg until the top coefficient is nonzero (or the degree
 * is 0). The capacity is left unchanged.
 *
 * Complexity:
 *   Time: O(number of trimmed coefficients)
 *   Auxiliary memory: O(1)
 *   Output memory: O(1)
 *
 * @param[in,out] p Polynomial to trim.
 */
void bigpoly_trim(bigpoly* p);

/**
 * @brief Print a polynomial to stdout in descending powers of x.
 *
 * Complexity:
 *   Time: O(d * n) where d is the degree and n the coefficient size
 *   Auxiliary memory: O(1)
 *   Output memory: O(d * n) characters written
 *
 * @param[in] p Polynomial to print.
 */
void bigpoly_print(const bigpoly* p);

/**
 * @brief Add two polynomials: r = p + q.
 *
 * Let d = max(p->deg, q->deg).
 *
 * Coefficients up to the smaller degree are added with bn_add(); the
 * remaining coefficients of the longer polynomial are copied. The
 * result is trimmed.
 *
 * Complexity:
 *   Time: O(d * n) where n is the size of the coefficients in limbs
 *   Auxiliary memory: O(1)
 *   Output memory: O(d) bignums
 *
 * @param[out] r Result polynomial.
 * @param[in]  p First polynomial.
 * @param[in]  q Second polynomial.
 */
void bigpoly_add(bigpoly* r, const bigpoly* p, const bigpoly* q);

/**
 * @brief Subtract two polynomials: r = p - q.
 *
 * Let d = max(p->deg, q->deg).
 *
 * Coefficients up to the smaller degree are subtracted with bn_sub();
 * the remaining coefficients of the longer polynomial are copied with
 * their sign flipped. The result is trimmed.
 *
 * Complexity:
 *   Time: O(d * n) where n is the size of the coefficients in limbs
 *   Auxiliary memory: O(1)
 *   Output memory: O(d) bignums
 *
 * @param[out] r Result polynomial.
 * @param[in]  p First polynomial.
 * @param[in]  q Second polynomial.
 */
void bigpoly_sub(bigpoly* r, const bigpoly* p, const bigpoly* q);

/**
 * @brief Multiply two polynomials with the bignum NTT: r = p * q.
 *
 * Let d = p->deg + q->deg + 1 (length of the true product) and
 * n = 2^k the next power of two >= d.
 *
 * Pads both polynomials to length n, applies the forward bignum NTT
 * (Goldilocks field) to each, multiplies pointwise in the Montgomery
 * domain, and applies the inverse NTT. The first d coefficients of the
 * cyclic convolution equal the true product coefficients (no
 * wrap-around, since the Goldilocks prime is large enough for the
 * coefficient sizes used).
 *
 * A fresh NTT context is initialized per call.
 *
 * Complexity:
 *   Time: O(n log n * k_l^2) where k_l is the size of the modulus in
 *         limbs, plus O(n * k_l^2) for the context tables
 *   Auxiliary memory: O(n) bignums for the transformed copies
 *   Output memory: O(d) bignums
 *
 * @param[out] r Result polynomial.
 * @param[in]  p First polynomial.
 * @param[in]  q Second polynomial.
 */
void bigpoly_mul_ntt(bigpoly* r, const bigpoly* p, const bigpoly* q);

/**
 * @brief Multiply two polynomials with the u64 NTT: r = p * q.
 *
 * Let d = p->deg + q->deg + 1 (length of the true product) and
 * n = 2^k the next power of two >= d.
 *
 * Only works for coefficients that fit in a single u64 limb: each
 * coefficient is reduced to its low limb, the convolution is computed
 * with the fast u64 Goldilocks NTT (flat arrays, Montgomery
 * butterflies), and the result limbs are lifted back to bignums.
 *
 * The Goldilocks NTT context for the transform size is taken from a
 * shared per-size cache (built once per k, see
 * ntt_ctx_u64_golden_cached), and the flat u64 arrays come from the
 * thread-local scratch arena.
 *
 * Complexity:
 *   Time: O(n log n) for the NTTs
 *   Auxiliary memory: O(n) u64s for the flat arrays (scratch arena)
 *   Output memory: O(d) bignums
 *
 * @param[out] r Result polynomial.
 * @param[in]  p First polynomial.
 * @param[in]  q Second polynomial.
 */
void bigpoly_mul_ntt_u64(bigpoly* r, const bigpoly* p, const bigpoly* q);

/**
 * @brief Multiply two polynomials with the schoolbook algorithm: r = p * q.
 *
 * Let d = p->deg + q->deg.
 *
 * Each pair of coefficients is multiplied and accumulated into the
 * corresponding slot of a temporary polynomial, which is then copied
 * into r and trimmed.
 *
 * Complexity:
 *   Time: O((p->deg + 1) * (q->deg + 1) * n^2) where n is the size of
 *         the coefficients in limbs
 *   Auxiliary memory: O(d) bignums for the temporary
 *   Output memory: O(d) bignums
 *
 * @param[out] r Result polynomial.
 * @param[in]  p First polynomial.
 * @param[in]  q Second polynomial.
 */
void bigpoly_mul_school(bigpoly* r, const bigpoly* p, const bigpoly* q);

/**
 * @brief Multiply two polynomials: r = p * q.
 *
 * Currently a thin wrapper around bigpoly_mul_school().
 *
 * Complexity:
 *   Time: see bigpoly_mul_school()
 *   Auxiliary memory: O(d) bignums
 *   Output memory: O(d) bignums
 *
 * @param[out] r Result polynomial.
 * @param[in]  p First polynomial.
 * @param[in]  q Second polynomial.
 */
void bigpoly_mul(bigpoly* r, const bigpoly* p, const bigpoly* q);

/**
 * @brief Decompose a bignum into a polynomial of fixed-width chunks.
 *
 * Let n = n->size, measured in 64-bit limbs, and W = width.
 *
 * Slices n into chunks of W bits (base 2^W) and stores them as the
 * coefficients of r in ascending order:
 *
 *   n = sum_i r->coeff[i] * (2^W)^i
 *
 * The degree is ceil(n / W) - 1 (0 for n == 0).
 *
 * Complexity:
 *   Time: O(n * ceil(n/W)) - one O(n) shift per chunk
 *   Auxiliary memory: O(n) limbs for temporaries
 *   Output memory: O(n/W) bignums
 *
 * @param[out]    r    Result polynomial storing the chunks.
 * @param[in]     n    Bignum to decompose.
 * @param[in] width Chunk width in bits.
 */
void bn_decompose(bigpoly* r, const bignum* n, u64 width);

/**
 * @brief Propagate carries between the fixed-width coefficient slots of a
 * polynomial, in place.
 *
 * Let d = r->deg and W = bit_width.
 *
 * After an NTT-based convolution the coefficient slots may exceed
 * 2^W - 1. This normalizes them: for each slot, total = coeff + carry,
 * coeff = total mod 2^W, carry = total / 2^W, rippling upward and
 * growing the polynomial if the carry outlives the current degree.
 *
 * Complexity:
 *   Time: O(d * n) where n is the size of the coefficients in limbs
 *   Auxiliary memory: O(n) limbs for temporaries
 *   Output memory: O(d) bignums (possibly grown by the final carry)
 *
 * @param[in,out]      r       Polynomial to normalize (modified in place).
 * @param[in] bit_width Width of each coefficient slot in bits.
 */
void poly_carry_propagation(bigpoly* r, u64 bit_width);

/**
 * @brief Recompose a polynomial of fixed-width chunks into a bignum.
 *
 * Let d = p->deg and W = bit_width.
 *
 * The inverse of bn_decompose():
 *
 *   n = sum_i p->coeff[i] << (i * W)
 *
 * Complexity:
 *   Time: O(d * n) where n is the size of the result in limbs
 *   Auxiliary memory: O(n) limbs for the running term
 *   Output memory: O(n) limbs
 *
 * @param[out]        n         Result bignum.
 * @param[in]         p         Polynomial of chunks.
 * @param[in] bit_width Width of each chunk in bits.
 */
void bn_recompose(bignum* n, const bigpoly* p, u64 bit_width);

/**
 * @brief Self-test of the polynomial operations.
 *
 * Builds p(x) = 2x^2 + 3x + 1 and q(x) = 4x + 5, prints them, and
 * prints p + q and p * q (schoolbook and u64-NTT). Intended for
 * manual verification, not part of the library API.
 *
 * Complexity:
 *   Time: O(1) (fixed-size inputs)
 *   Auxiliary memory: O(1)
 *   Output memory: O(1)
 */
void bigpoly_test();

#endif
