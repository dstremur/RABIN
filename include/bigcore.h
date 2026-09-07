#ifndef BIGCORE_H
#define BIGCORE_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>

typedef uint64_t u64;
typedef int64_t i64;

typedef struct bignum {
  uint64_t* limbs;  // 64-bit limbs (Base 2^64)
  size_t size;
  size_t capacity;
  bool is_neg;
} bignum;

typedef struct {
  bignum n;        /* modulus (odd, > 1)               */
  uint64_t n_inv;  /* -n^{-1} mod 2^64                  */
  bignum r_square; /* R^2 mod N   (R = 2^{64 * n.limbs})*/
  bignum one_mont; /* R   mod N   (representation of 1) */
  bignum tmp;      // scratchpad
} bn_mont_ctx;

extern bignum BN_ONE;
extern bignum BN_TWO;
extern bignum BN_ZERO;

#define MAX(a, b) ((a) > (b) ? (a) : (b));
#define MIN(a, b) ((a) < (b) ? (a) : (b));

// bigconst.c
/**
 * @brief Initialize the global constants BN_ZERO, BN_ONE, and BN_TWO.
 *
 * Must be called once before any of the constants is used. Calling it
 * more than once without an intervening bn_free_constants() leaks the
 * previously allocated limb storage.
 *
 * Complexity:
 *   Time: O(1)
 *   Auxiliary memory: O(1)
 *   Output memory: O(1) limbs (one limb per constant)
 */
void bn_init_constants();

/**
 * @brief Free the limb storage of the global constants.
 *
 * After this call the constants are uninitialized and must not be used
 * until bn_init_constants() is called again.
 *
 * Complexity:
 *   Time: O(1)
 *   Auxiliary memory: O(1)
 *   Output memory: O(1)
 */
void bn_free_constants();

// bnscratch.c (thread-local bump arena; one get + one release per call)
/**
 * @brief Reserve n u64s of thread-local scratch space.
 *
 * Returns a pointer to at least n u64s. The memory is NOT zeroed.
 * The block stays valid until bn_scratch_release() is called.
 *
 * Complexity:
 *   Time: O(1) amortized, O(1) realloc on arena growth
 *   Auxiliary memory: O(n) u64s (grow-only per thread)
 *   Output memory: O(n) u64s
 *
 * @param[in] n Number of u64s to reserve.
 *
 * @return A pointer to at least n u64s of scratch (not zeroed).
 */
u64* bn_scratch_get(u64 n);

/**
 * @brief Release the current thread's scratch block(s) by rewinding the bump
 * pointer. The arena memory itself is kept for reuse.
 *
 * Complexity:
 *   Time: O(1)
 *   Auxiliary memory: O(1)
 *   Output memory: O(1)
 */
void bn_scratch_release(void);

// bignum.c
/**
 * @brief Initialize a bignum to zero, releasing no memory (r must not have
 * been allocated yet).
 *
 * Complexity:
 *   Time: O(1)
 *   Auxiliary memory: O(1)
 *   Output memory: O(1)
 *
 * @param[out] r Bignum to initialize to zero.
 */
void bn_init(bignum* r);

/**
 * @brief Initialize multiple bignums to zero.
 *
 * Takes a NULL-terminated list of bignum pointers; the first argument
 * is the first bignum and the variadic arguments continue the list.
 * IMPORTANT: the list must be terminated with NULL.
 *
 * Complexity:
 *   Time: O(k) for k bignums
 *   Auxiliary memory: O(1)
 *   Output memory: O(1)
 *
 * @param[out] r First bignum to initialize; further bignums follow as
 *               variadic arguments, terminated by NULL.
 */
void bn_init_multi(bignum* r, ...);

/**
 * @brief Swap the contents of two bignums by exchanging their structs.
 *
 * Only the struct fields (pointers, sizes, sign) are swapped, so this
 * is O(1) regardless of size.
 *
 * Complexity:
 *   Time: O(1)
 *   Auxiliary memory: O(1)
 *   Output memory: O(1)
 *
 * @param[in,out] a First bignum.
 * @param[in,out] b Second bignum.
 */
void bn_swap(bignum* a, bignum* b);

/**
 * @brief Test whether a bignum is even.
 *
 * Returns 1 if n is even (including zero), 0 otherwise. Only the least
 * significant bit of the lowest limb is inspected.
 *
 * Complexity:
 *   Time: O(1)
 *   Auxiliary memory: O(1)
 *   Output memory: O(1)
 *
 * @param[in] n Bignum to test.
 *
 * @return 1 If n is even (including zero), 0 otherwise.
 */
int bn_is_even(const bignum* n);

/**
 * @brief Test whether a bignum is zero.
 *
 * Complexity:
 *   Time: O(1)
 *   Auxiliary memory: O(1)
 *   Output memory: O(1)
 *
 * @param[in] a Bignum to test.
 *
 * @return true If a is zero, false otherwise.
 */
bool bn_is_zero(const bignum* a);

/**
 * @brief Test whether a bignum equals a 64-bit signed integer.
 *
 * Returns true iff n has exactly one limb and that limb equals a.
 *
 * Complexity:
 *   Time: O(1)
 *   Auxiliary memory: O(1)
 *   Output memory: O(1)
 *
 * @param[in] n Bignum to test.
 * @param[in] a 64-bit signed integer to compare against.
 *
 * @return true If n equals a, false otherwise.
 */
bool bn_is_eq_i64(const bignum* n, i64 a);

/**
 * @brief Grow the limb storage of r to at least capacity limbs.
 *
 * If r->capacity is already sufficient, nothing happens. Otherwise the
 * storage is reallocated with exponential growth (at least double the
 * old capacity, or the requested capacity, whichever is larger) and
 * the newly added limbs are zeroed. Existing limb contents are
 * preserved.
 *
 * Returns true on success, false if the realloc fails (r is left
 * unchanged).
 *
 * Complexity:
 *   Time: O(capacity) for the zero-fill and the realloc copy
 *   Auxiliary memory: O(capacity) during the realloc
 *   Output memory: O(capacity) limbs
 *
 * @param[in,out] r        Bignum whose storage is grown.
 * @param[in]     capacity Minimum number of limbs required.
 *
 * @return true  On success.
 * @return false If the realloc fails (r is left unchanged).
 */
bool bn_alloc(bignum* r, u64 capacity);

/**
 * @brief Parse a base-10 decimal string into a bignum.
 *
 * Let d = number of decimal digits in str.
 *
 * Frees any previous contents of n and sets it to the value of str,
 * which may start with a '-' sign. Non-digit characters are skipped.
 * Digits are processed one at a time with a multiply-by-10-and-add
 * carry loop over the limbs.
 *
 * Complexity:
 *   Time: O(d * n) where n is the number of limbs, i.e. O(d^2 / 64)
 *   Auxiliary memory: O(1)
 *   Output memory: O(n) limbs
 *
 * @param[out] n   Bignum to set (previous contents are freed).
 * @param[in]  str Base-10 decimal string to parse (may start with '-').
 */
void bn_init_val(bignum* n, const char* str);

/**
 * @brief Convert a bignum to a base-10 decimal string.
 *
 * Let n = a->size, measured in 64-bit limbs.
 *
 * Returns a newly allocated string (caller must free it) containing
 * the decimal representation of n, with a leading '-' for negative
 * values. The magnitude is reduced repeatedly by 10^19 (which fits in
 * a u64), collecting 19-digit chunks; the most significant chunk is
 * printed without padding and the rest with zero padding.
 *
 * Returns NULL on allocation failure.
 *
 * Complexity:
 *   Time: O(n^2) - O(n) divisions by a 64-bit divisor, each O(n)
 *   Auxiliary memory: O(n) limbs for temporaries, O(n) for the chunk
 *                    array and the output string
 *   Output memory: O(n) characters
 *
 * @param[in] n Bignum to convert.
 *
 * @return A newly allocated decimal string (caller must free it), or
 *         NULL on allocation failure.
 */
char* bn_to_string(const bignum* n);

/**
 * @brief Print a bignum to stdout in base 10.
 *
 * Let n = a->size, measured in 64-bit limbs.
 *
 * Same chunking strategy as bn_to_string(): the magnitude is reduced
 * repeatedly by 10^19 and the chunks are printed most-significant
 * first, with a leading '-' for negative values.
 *
 * Complexity:
 *   Time: O(n^2) - O(n) divisions by a 64-bit divisor, each O(n)
 *   Auxiliary memory: O(n) limbs for temporaries, O(n) for the chunk
 *                    array
 *   Output memory: O(n) characters written
 *
 * @param[in] n Bignum to print.
 */
void bn_print(const bignum* n);

/**
 * @brief Print a bignum to stdout in base 10 followed by a newline.
 *
 * Complexity:
 *   Time: O(n^2), see bn_print()
 *   Auxiliary memory: O(n)
 *   Output memory: O(n) characters written
 *
 * @param[in] n Bignum to print.
 */
void bn_println(const bignum* n);

/**
 * @brief Free the limb storage of a bignum and reset it to the zero state.
 *
 * Safe to call on a NULL pointer or on an already-freed bignum.
 *
 * Complexity:
 *   Time: O(1)
 *   Auxiliary memory: O(1)
 *   Output memory: O(1)
 *
 * @param[in,out] r Bignum to free (may be NULL).
 */
void bn_free(bignum* r);

/**
 * @brief Free multiple bignums.
 *
 * Takes a NULL-terminated list of bignum pointers; the first argument
 * is the first bignum and the variadic arguments continue the list.
 * IMPORTANT: the list must be terminated with NULL.
 *
 * Complexity:
 *   Time: O(k) for k bignums
 *   Auxiliary memory: O(1)
 *   Output memory: O(1)
 *
 * @param[in,out] r First bignum to free; further bignums follow as
 *                  variadic arguments, terminated by NULL.
 */
void bn_free_multi(bignum* r, ...);

/**
 * @brief Remove trailing zero limbs from a bignum.
 *
 * Reduces r->size so that the most significant limb is nonzero (or
 * size is 1). The capacity is left unchanged.
 *
 * Complexity:
 *   Time: O(number of trimmed limbs)
 *   Auxiliary memory: O(1)
 *   Output memory: O(1)
 *
 * @param[in,out] r Bignum to trim.
 */
void bn_trim(bignum* r);

/**
 * @brief Return the value of the i-th bit of a (0 or 1).
 *
 * Bit 0 is the least significant bit. Bits beyond the current size
 * read as 0.
 *
 * Complexity:
 *   Time: O(1)
 *   Auxiliary memory: O(1)
 *   Output memory: O(1)
 *
 * @param[in] a Bignum to read the bit from.
 * @param[in] i Bit index (0 = least significant).
 *
 * @return The value of bit i (0 or 1).
 */
int bn_get_bit(const bignum* a, int i);

/**
 * @brief Compare two signed bignums.
 *
 * Returns 1 if a > b, -1 if a < b, 0 if a == b.
 *
 * Signs are handled first (a positive number is greater than a
 * negative one); for equal signs the magnitudes are compared, most
 * significant limb first, and the result is negated when both
 * operands are negative.
 *
 * Complexity:
 *   Time: O(n) worst case, where n = max(a->size, b->size)
 *   Auxiliary memory: O(1)
 *   Output memory: O(1)
 *
 * @param[in] a First bignum.
 * @param[in] b Second bignum.
 *
 * @return 1 If a > b, -1 if a < b, 0 if a == b.
 */
int bn_cmp(const bignum* a, const bignum* b);

/**
 * @brief Compare the absolute values (magnitudes) of two bignums.
 *
 * Returns 1 if |a| > |b|, -1 if |a| < |b|, 0 if |a| == |b|. Signs are
 * ignored. The comparison is by size first, then most significant
 * limb first.
 *
 * Complexity:
 *   Time: O(n) worst case, where n = max(a->size, b->size)
 *   Auxiliary memory: O(1)
 *   Output memory: O(1)
 *
 * @param[in] a First bignum.
 * @param[in] b Second bignum.
 *
 * @return 1 If |a| > |b|, -1 if |a| < |b|, 0 if |a| == |b|.
 */
int bn_cmp_abs(const bignum* a, const bignum* b);

/**
 * @brief Deep copy a bignum: r = a.
 *
 * Copies the limb storage and the sign. r may alias a (a no-op in
 * that case). If a is zero-sized, r is reset to size 0.
 *
 * Complexity:
 *   Time: O(n) where n = a->size
 *   Auxiliary memory: O(1)
 *   Output memory: O(n) limbs
 *
 * @param[out] r Destination of the copy.
 * @param[in]  a Source bignum.
 */
void bn_copy(bignum* r, const bignum* a);

/**
 * @brief Set a bignum to an unsigned 64-bit value.
 *
 * Frees any previous contents and stores val in a single limb with a
 * positive sign.
 *
 * Complexity:
 *   Time: O(1)
 *   Auxiliary memory: O(1)
 *   Output memory: O(1) limbs
 *
 * @param[out] n   Bignum to set.
 * @param[in]  val Unsigned 64-bit value to store.
 */
void bn_set_u64(bignum* n, uint64_t val);

/**
 * @brief Set a bignum to a signed 64-bit value.
 *
 * For negative val the magnitude is computed as -val (with overflow
 * protection for INT64_MIN) and the sign flag is set.
 *
 * Complexity:
 *   Time: O(1)
 *   Auxiliary memory: O(1)
 *   Output memory: O(1) limbs
 *
 * @param[out] n   Bignum to set.
 * @param[in]  val Signed 64-bit value to store.
 */
void bn_set_i64(bignum* n, int64_t val);

/**
 * @brief Set the i-th bit of a to 1.
 *
 * Bit 0 is the least significant bit. If the bit lies beyond the
 * current size, the bignum is grown (with zeroed limbs) to reach it.
 *
 * Complexity:
 *   Time: O(1) amortized (O(limb) for the zero-fill on growth)
 *   Auxiliary memory: O(1)
 *   Output memory: O(1) limbs, possibly grown
 *
 * @param[in,out] a Bignum to modify.
 * @param[in]     i Bit index (0 = least significant).
 */
void bn_set_bit(bignum* a, int i);

/**
 * @brief Clear the i-th bit of a to 0.
 *
 * Bit 0 is the least significant bit. Bits beyond the current size
 * are already 0 and are left untouched.
 *
 * Complexity:
 *   Time: O(1)
 *   Auxiliary memory: O(1)
 *   Output memory: O(1)
 *
 * @param[in,out] a Bignum to modify.
 * @param[in]     i Bit index (0 = least significant).
 */
void bn_clear_bit(bignum* a, int i);

/**
 * @brief Return the bit length of a: the index of the highest set bit plus
 * one (0 for zero).
 *
 * Only the most significant limb is inspected.
 *
 * Complexity:
 *   Time: O(1)
 *   Auxiliary memory: O(1)
 *   Output memory: O(1)
 *
 * @param[in] a Bignum to measure.
 *
 * @return The bit length of a (index of highest set bit + 1, or 0).
 */
int bn_bit_length(const bignum* a);

/**
 * @brief Return the number of trailing zero bits of a.
 *
 * Counts whole zero limbs (64 bits each) plus the trailing zeros of
 * the first nonzero limb. Returns 0 for a == 0.
 *
 * Complexity:
 *   Time: O(number of zero limbs)
 *   Auxiliary memory: O(1)
 *   Output memory: O(1)
 *
 * @param[in] a Bignum to count trailing zero bits of.
 *
 * @return The number of trailing zero bits (0 if a is zero).
 */
u64 bn_cnt_trailing_zeros(const bignum* a);

// bigadd.c

/**
 * @brief Add two signed bignums.
 *
 * Let n = max(a->size, b->size), measured in 64-bit limbs.
 *
 * If a and b have the same sign, their magnitudes are added.
 *
 * If a and b have different signs, the smaller magnitude is subtracted
 * from the larger magnitude, and the result takes the sign of the operand
 * with the larger magnitude.
 *
 * If the magnitudes are equal, the result is normalized to positive zero.
 *
 * Complexity:
 *   Time: O(n)
 *   Auxiliary memory: O(1), except for any temporary storage used by
 *                     bn_add_abs(), bn_sub_abs(), bn_copy(), or bn_alloc()
 *   Output memory: O(n) limbs
 *
 * @param[out] r Result of the signed addition a + b.
 * @param[in]  a First operand.
 * @param[in]  b Second operand.
 */
void bn_add(bignum* r, const bignum* a, const bignum* b);

/**
 * @brief Add the absolute values of a and b into r.
 *
 * Let n = max(a->size, b->size), measured in 64-bit limbs.
 *
 * This function only adds magnitudes. It does not interpret or modify
 * the sign of the operands or the result. The caller is responsible for
 * setting r->is_neg.
 *
 * r may alias a or b.
 *
 * Complexity:
 *   Time: O(n)
 *   Auxiliary memory: O(1) in the normal case,
 *                     O(n) if r aliases a or b and a temporary is used
 *   Output memory: O(n) limbs, at most n + 1 limbs
 *
 * @param[out] r Result storing the magnitude sum |a| + |b|.
 * @param[in]  a First operand (magnitude only).
 * @param[in]  b Second operand (magnitude only).
 */
void bn_add_abs(bignum* r, const bignum* a, const bignum* b);

/**
 * @brief Add an unsigned 64-bit value to a signed bignum.
 *
 * Let n = a->size, measured in 64-bit limbs.
 *
 * This computes:
 *
 *   r = a + b
 *
 * where b is a nonnegative u64.
 *
 * If a is positive, this is ordinary magnitude addition.
 *
 * If a is negative, this is magnitude subtraction:
 *
 *   r = -|a| + b
 *
 * which is equivalent to subtracting b from |a| and preserving the
 * correct sign.
 *
 * r may alias a, assuming bn_copy() supports aliasing.
 *
 * Complexity:
 *   Time: O(n) worst case.
 *         The carry/borrow loop may stop early if propagation ends
 *         before reaching the most significant limb.
 *   Auxiliary memory: O(1), except for any temporary storage used by
 *                     bn_copy() or bn_alloc()
 *   Output memory: O(n) limbs, or O(n + 1) if a new carry limb is created
 *
 * @param[out] r Result of the addition a + b.
 * @param[in]  a Signed bignum operand.
 * @param[in]  b Nonnegative 64-bit value to add.
 */
void bn_add_u64(bignum* r, const bignum* a, u64 b);

/**
 * @brief Add a to r at the given limb offset.
 *
 * Let n = a->size, measured in 64-bit limbs.
 * Let o = offset.
 * Let m = o + n + 1 be the maximum resulting size in limbs.
 *
 * This effectively performs:
 *
 *   r = r + (a << (offset * 64))
 *
 * for unsigned magnitudes.
 *
 * The addition loop itself only touches about n + 1 limbs, but allocation
 * and trimming may need to consider the full output range up to m limbs.
 *
 * Complexity:
 *   Time: O(m) worst case, where m = offset + a->size + 1.
 *         The inner addition loop is O(n), but bn_alloc() and bn_trim()
 *         may make the total worst case proportional to the output size.
 *   Auxiliary memory: O(1)
 *   Output memory: O(m) limbs
 *
 * Note:
 *   The current implementation does not explicitly handle r == a safely
 *   for all nonzero offsets. If aliasing is not supported, document that
 *   here and enforce it in the API.
 *
 * @param[out] r Accumulator; receives r + (a << (offset * 64)).
 * @param[in]  a Value to add.
 * @param[in]  offset Limb offset at which a is added.
 */
void bn_add_at_offset(bignum* r, const bignum* a, u64 offset);

// bigsub.c

/**
 * @brief Subtract two signed bignums: r = a - b.
 *
 * Let n = max(a->size, b->size), measured in 64-bit limbs.
 *
 * Subtraction is reduced to magnitude addition or subtraction:
 *
 *   a - (-b) = a + b          (opposite signs, b negative)
 *   (-a) - b = -(a + b)       (opposite signs, a negative)
 *   same signs: subtract the smaller magnitude from the larger one and
 *               take the sign of the operand with the larger magnitude
 *
 * If the magnitudes are equal, the result is zero.
 *
 * Complexity:
 *   Time: O(n)
 *   Auxiliary memory: O(1), except for any temporary storage used by
 *                     bn_add_abs(), bn_sub_abs(), or bn_alloc()
 *   Output memory: O(n) limbs
 *
 * @param[out] r Result of the signed subtraction a - b.
 * @param[in]  a Minuend.
 * @param[in]  b Subtrahend.
 */
void bn_sub(bignum* r, const bignum* a, const bignum* b);

/**
 * @brief Subtract the absolute value of b from the absolute value of a into r.
 *
 * Let n = a->size, measured in 64-bit limbs.
 *
 * This computes:
 *
 *   r = |a| - |b|
 *
 * and requires |a| >= |b|; if |b| > |a| the result wraps around
 * (two's-complement style) and is meaningless. It does not interpret or
 * modify the sign of the operands or the result. The caller is
 * responsible for setting r->is_neg.
 *
 * r may alias a or b.
 *
 * Complexity:
 *   Time: O(n)
 *   Auxiliary memory: O(1) in the normal case,
 *                     O(n) if r aliases a or b and a temporary is used
 *   Output memory: O(n) limbs
 *
 * @param[out] r Result storing |a| - |b| (requires |a| >= |b|).
 * @param[in]  a Minuend (magnitude only).
 * @param[in]  b Subtrahend (magnitude only).
 */
void bn_sub_abs(bignum* r, const bignum* a, const bignum* b);

// bigmul.c

/**
 * @brief Multiply two signed bignums: r = a * b.
 *
 * Let n = max(a->size, b->size), measured in 64-bit limbs.
 *
 * This computes the signed product of a and b. The sign of the result
 * is the xor of the operand signs; the magnitude is computed with:
 *
 *   - the squaring path (bn_sqr) when a == b
 *   - Karatsuba (O(n^1.585)) when both operands have at least
 *     KARATSUBA_LIMIT limbs
 *   - schoolbook (O(n^2)) otherwise
 *
 * r may alias a or b.
 *
 * Complexity:
 *   Time: O(n^1.585) for large operands, O(n^2) for small ones
 *   Auxiliary memory: O(n) limbs of scratch, plus O(n) for the
 *                     zero-padding buffers on the Karatsuba path and
 *                     O(n) if r aliases an operand
 *   Output memory: O(n) limbs (at most a->size + b->size)
 *
 * @param[out] r Result of the signed product a * b.
 * @param[in]  a First operand.
 * @param[in]  b Second operand.
 */
void bn_mul(bignum* r, const bignum* a, const bignum* b);
void bn_mul_raw(bignum* r, const bignum* a, const bignum* b);
void bn_mul_karatsuba(bignum* r, const bignum* a, const bignum* b);

/**
 * @brief Multiply two signed bignums with the schoolbook (grade-school)
 * algorithm: r = a * b.
 *
 * Let n = a->size and m = b->size, measured in 64-bit limbs.
 *
 * Each limb of a is multiplied by the whole of b and accumulated into
 * r at the corresponding offset (bn_mul_add_inner), skipping zero
 * limbs. The sign of the result is the xor of the operand signs.
 *
 * Complexity:
 *   Time: O(n * m)
 *   Auxiliary memory: O(1)
 *   Output memory: O(n + m) limbs
 *
 * @param[out] r Result of the schoolbook product a * b.
 * @param[in]  a First operand.
 * @param[in]  b Second operand.
 */
void bn_mul_school(bignum* r, const bignum* a, const bignum* b);

/**
 * @brief Multiply two bignums with the NTT-based fast path: res = a * b.
 *
 * Let n = max(a->size, b->size), measured in 64-bit limbs.
 *
 * This computes the product of the magnitudes of a and b by:
 *
 *   1. decomposing each operand into a polynomial whose coefficients
 *      are 16-bit chunks (base 2^16),
 *   2. multiplying the polynomials with a cyclic NTT
 *      (bigpoly_mul_ntt),
 *   3. propagating carries between the 16-bit coefficient slots,
 *   4. recomposing the coefficients back into a bignum.
 *
 * Unlike bn_mul() this path is not wired into the general dispatch;
 * it is a standalone fast path for very large operands.
 *
 * Complexity:
 *   Time: O(n log n) for the NTT, plus O(n) for decompose/carry/
 *         recompose
 *   Auxiliary memory: O(n) limbs for the polynomial arrays
 *   Output memory: O(n) limbs (at most a->size + b->size)
 *
 * @param[out] res Result of the NTT-based product a * b.
 * @param[in]  a   First operand.
 * @param[in]  b   Second operand.
 */
void bn_mul_fast(bignum* res, const bignum* a, const bignum* b);

/**
 * @brief Square a bignum: r = a * a.
 *
 * Let n = a->size, measured in 64-bit limbs.
 *
 * This computes the square of the magnitude of a using the Karatsuba
 * squaring kernel (limbs_sqr_karatsuba), which exploits the symmetry
 * of squaring to save about a quarter of the multiplications. The
 * result is always nonnegative.
 *
 * If a is zero, r is set to zero.
 *
 * Complexity:
 *   Time: O(n^1.585) for n >= KARATSUBA_LIMIT, O(n^2) below it
 *   Auxiliary memory: O(n) limbs of scratch
 *   Output memory: O(n) limbs (at most 2n)
 *
 * @param[out] r Result storing a * a.
 * @param[in]  a Value to square.
 */
void bn_sqr(bignum* r, const bignum* a);

// bigdiv.c

/**
 * @brief q = a / b (truncated division, C semantics).
 *
 * Let n = a->size and m = b->size, measured in 64-bit limbs.
 *
 * Thin wrapper around bn_divmod() that discards the remainder.
 *
 * Complexity:
 *   Time: O(n) for a single-limb divisor, O((n - m + 1) * m) otherwise
 *   Auxiliary memory: O(n + m) limbs
 *   Output memory: O(n - m + 1) limbs
 *
 * @param[out] q Quotient a / b.
 * @param[in]  a Dividend.
 * @param[in]  b Divisor (must be nonzero).
 */
void bn_div(bignum* q, const bignum* a, const bignum* b);

/**
 * @brief q = a / b, r = a % b  (either result may be NULL).
 *
 * Let n = a->size and m = b->size, measured in 64-bit limbs.
 *
 * Truncated division (C semantics): q truncates toward 0 and
 * sign(r) == sign(a).
 *
 * If |a| < |b| the result is q = 0, r = a. Otherwise the magnitude
 * division dispatches to the single-limb path (limbs_divrem_1) or to
 * Knuth's Algorithm D (bn_divmod_limbs). q and r may each alias a or
 * b (handled via temporaries); q and r must not alias each other.
 *
 * Scratch for the multi-limb path is taken from the stack when it fits
 * in BN_DIV_STACK_LIMBS limbs, otherwise from the heap.
 *
 * Complexity:
 *   Time: O(n) for a single-limb divisor, O((n - m + 1) * m) for a
 *         multi-limb divisor
 *   Auxiliary memory: O(n + m) limbs of scratch (stack or heap), plus
 *                     O(n) if q or r aliases an operand
 *   Output memory: O(n - m + 1) limbs for q, O(m) for r
 *
 * @param[out] q Quotient a / b (may be NULL).
 * @param[out] r Remainder a % b (may be NULL).
 * @param[in]  a Dividend.
 * @param[in]  b Divisor (must be nonzero).
 */
void bn_divmod(bignum* q, bignum* r, const bignum* a, const bignum* b);

/**
 * @brief Computes the division a / d using a Newton-iterated reciprocal: q = a
 * / d.
 *
 * Let n = a->size and m = d->size, measured in 64-bit limbs.
 *
 * Works with P = bit_length(a) + 32 bits of precision. The reciprocal
 * x = 2^P / d is seeded with one real division, then refined with the
 * Newton iteration
 *
 *   x <- x + (x * (2^P - d*x)) >> P
 *
 * which roughly doubles the correct bits each step. The quotient is
 * then q = (a * x) >> P, followed by a rare off-by-one fix-up using
 * the remainder r = a - q*d.
 *
 * If d is zero, q is left unchanged. If |a| < |d|, q is set to 0.
 *
 * Complexity:
 *   Time: O(n^2) - O(log m) Newton iterations of n-limb multiplications
 *         plus one seeding division
 *   Auxiliary memory: O(n) limbs for temporaries
 *   Output memory: O(n) limbs
 *
 * @param[out] q Quotient a / d.
 * @param[in]  a Dividend.
 * @param[in]  d Divisor.
 */
void bn_newton_div(bignum* q, const bignum* a, const bignum* d);

/**
 * @brief q = a / b, assuming b divides a exactly (Jebelean's exact division).
 *
 * Let n = a->size and m = b->size, measured in 64-bit limbs.
 *
 * Strips the common power of two from a and b (so the low limb of the
 * divisor becomes odd), then computes the quotient limbs from the
 * bottom up: with D[0] odd, each quotient limb qi = A[i] * D[0]^{-1}
 * (mod 2^64) is forced, because it must cancel limb i of the running
 * remainder. No estimation or correction steps are needed, which makes
 * this about 2x faster than a real division.
 *
 * The behaviour is undefined if b does not divide a. If b is zero, q
 * is left unchanged (asserts in debug builds). If a is zero, q is set
 * to 0. q may alias a or b.
 *
 * Complexity:
 *   Time: O(n * m) - one O(m) submul per quotient limb
 *   Auxiliary memory: O(n + m) limbs for the shifted copies
 *   Output memory: O(n - m + 1) limbs
 *
 * @param[out] q Quotient a / b (requires b | a).
 * @param[in]  a Dividend.
 * @param[in]  b Divisor (must divide a exactly).
 */
void bn_div_exact(bignum* q, const bignum* a, const bignum* b);

// bigmod.c

/**
 * @brief r = a % b (truncated remainder, C semantics: sign(r) == sign(a)).
 *
 * Let n = a->size and m = b->size, measured in 64-bit limbs.
 *
 * Thin wrapper around bn_divmod() that discards the quotient.
 *
 * Complexity:
 *   Time: O(n) for a single-limb divisor, O((n - m + 1) * m) otherwise
 *   Auxiliary memory: O(n + m) limbs
 *   Output memory: O(m) limbs
 *
 * @param[out] r Remainder a % b.
 * @param[in]  a Dividend.
 * @param[in]  b Divisor (must be nonzero).
 */
void bn_mod(bignum* r, const bignum* a, const bignum* b);

/**
 * @brief Divide a by a 64-bit divisor d, storing the quotient in q and
 * returning the remainder.
 *
 * Let n = a->size, measured in 64-bit limbs.
 *
 * This computes:
 *
 *   q = a / d,  return value = a mod d
 *
 * by processing the limbs of a most-significant-first with a running
 * 128-bit dividend: each step yields one quotient limb and the new
 * remainder. The quotient is truncated toward zero (C semantics) for
 * negative a.
 *
 * q may alias a.
 *
 * Complexity:
 *   Time: O(n)
 *   Auxiliary memory: O(1) in the normal case,
 *                     O(n) if q aliases a and a temporary is used
 *   Output memory: O(n) limbs for q
 *
 * @param[out] q Quotient a / d (truncated toward zero).
 * @param[in]  a Dividend.
 * @param[in]  d 64-bit divisor (must be nonzero).
 *
 * @return The remainder a mod d.
 */
uint64_t bn_divmod_u64(bignum* q, const bignum* a, uint64_t d);

/**
 * @brief Reduce a modulo a 64-bit divisor d.
 *
 * Let n = a->size, measured in 64-bit limbs.
 *
 * This computes:
 *
 *   a mod d
 *
 * by processing the limbs of a most-significant-first with a running
 * 128-bit remainder: rem = (rem * 2^64 + limb) mod d. The result is
 * returned in [0, d).
 *
 * For negative a the result is the nonnegative residue: if a < 0 and
 * the magnitude remainder is nonzero, d - rem is returned (C-style
 * truncating remainder mapped into [0, d)).
 *
 * Complexity:
 *   Time: O(n)
 *   Auxiliary memory: O(1)
 *   Output memory: O(1)
 *
 * @param[in] a Bignum to reduce.
 * @param[in] d 64-bit divisor (must be nonzero).
 *
 * @return The residue a mod d, in the range [0, d).
 */
uint64_t bn_mod_u64(const bignum* a, uint64_t d);

/**
 * @brief Compute the modular multiplicative inverse of a modulo m.
 *
 * Let n = m->size, measured in 64-bit limbs.
 *
 * This computes:
 *
 *   res = a^(-1) mod m
 *
 * i.e. the value in [0, m) such that (a * res) mod m == 1, using the
 * binary extended GCD (Stein's algorithm) with coefficient tracking.
 * The coefficients are kept reduced modulo m at every step so they
 * never grow beyond n limbs.
 *
 * Returns true and stores the inverse in res if it exists (i.e. gcd(a,
 * m) == 1). Returns false and leaves res unchanged if a is zero, m is
 * zero or one, or a and m are not coprime.
 *
 * Complexity:
 *   Time: O(n^2) - O(n) iterations of shifts and subtractions of
 *         n-limb values (binary GCD), each O(n)
 *   Auxiliary memory: O(n) limbs for temporaries
 *   Output memory: O(n) limbs
 *
 * @param[out] res Receives a^(-1) mod m if the inverse exists.
 * @param[in]  a   Value to invert.
 * @param[in]  m   Modulus.
 *
 * @return true  If the inverse exists (gcd(a, m) == 1); res is set.
 * @return false If a is zero, m is zero or one, or a and m are not
 *               coprime; res is left unchanged.
 */
bool bn_mod_inverse(bignum* res, const bignum* a, const bignum* m);

// bigexp.c
/**
 * @brief Calculate a^b into r using binary exponentiation.
 *
 * Let n = a->size and e = b->size, measured in 64-bit limbs.
 *
 * This computes:
 *
 *   r = a^b
 *
 * using left-to-right square-and-multiply over the bits of b. The
 * result has roughly e * n limbs, so this is only practical for small
 * exponents; use bn_mod_exp() for large ones.
 *
 * If b is zero, r is set to 1 (including 0^0).
 *
 * Complexity:
 *   Time: O(e * n^2) - one squaring per exponent bit (64*e of them)
 *         plus one multiplication per set bit, each an O(n^2) bignum
 *         multiply of growing operands
 *   Auxiliary memory: O(n) limbs for temporaries
 *   Output memory: O(e * n) limbs
 *
 * @param[out] r Result of a^b.
 * @param[in]  a Base.
 * @param[in]  b Exponent.
 */
void bn_pow(bignum* r, const bignum* a, const bignum* b);

/**
 * @brief Calculate a_bar^d in the Montgomery domain.
 *
 * Let n = ctx->n.size and e = d->size, measured in 64-bit limbs.
 *
 * This computes:
 *
 *   r_bar = a_bar^d (mod n)
 *
 * where a_bar is already in Montgomery form (a_bar = a * R mod n) and
 * r_bar is returned in Montgomery form.
 *
 * For short exponents (fewer than 64 bits) it uses plain left-to-right
 * square-and-multiply with bn_mont_mul(). For longer exponents it uses
 * a fixed window of width w = 4: the odd powers a^1, a^3, ..., a^15 are
 * precomputed in Montgomery form, then the exponent is scanned from the
 * most significant bit; each run is consumed as a window of up to w bits
 * ending in a 1-bit, costing w squarings plus one multiplication by the
 * precomputed window value. This uses roughly 15-20% fewer Montgomery
 * multiplications than plain binary.
 *
 * Complexity:
 *   Time: O(e * n^2) - one Montgomery squaring per exponent bit plus
 *         one Montgomery multiplication per window (plus the O(1)
 *         precompute for the window table)
 *   Auxiliary memory: O(n) limbs for the window table
 *   Output memory: O(n) limbs
 *
 * @param[out] r_bar Result in Montgomery form: a_bar^d (mod n).
 * @param[in]  a_bar Base in Montgomery form (a * R mod n).
 * @param[in]  d     Exponent.
 * @param[in]  ctx   Initialized Montgomery context for the modulus n.
 */
void bn_mont_exp(bignum* r_bar, const bignum* a_bar, const bignum* d,
                 bn_mont_ctx* ctx);

/**
 * @brief Calculate a^b mod m into r using plain (non-Montgomery) arithmetic.
 *
 * Let n = m->size and e = b->size, measured in 64-bit limbs.
 *
 * This computes:
 *
 *   r = a^b mod m
 *
 * using right-to-left binary exponentiation: the base is squared and
 * reduced modulo m for every bit of b, and the accumulator is
 * multiplied by the base and reduced for every set bit.
 *
 * This is the slow path: every step costs a full bignum multiplication
 * plus a full Knuth division. bn_mod_exp() uses Montgomery arithmetic
 * for odd moduli and is much faster.
 *
 * Complexity:
 *   Time: O(e * n^2) multiplications plus O(e * n^2) divisions, i.e.
 *         O(e * n^2) with a large constant
 *   Auxiliary memory: O(n) limbs for temporaries
 *   Output memory: O(n) limbs
 *
 * @param[out] r Result of a^b mod m.
 * @param[in]  a Base.
 * @param[in]  b Exponent.
 * @param[in]  m Modulus.
 */
void bn_mod_exp_slow(bignum* r, const bignum* a, const bignum* b,
                     const bignum* m);

/**
 * @brief Calculate a^b mod m into r.
 *
 * Let n = m->size and e = b->size, measured in 64-bit limbs.
 *
 * This computes:
 *
 *   r = a^b mod m
 *
 * For odd m it uses the fast Montgomery path (bn_mod_exp_mont) with a
 * freshly initialized context. For even m it falls back to
 * bn_mod_exp_slow(), since Montgomery reduction requires an odd
 * modulus.
 *
 * Complexity:
 *   Time: O(e * n^2) - O(e) Montgomery multiplications (or plain
 *         multiply+divide pairs for even m)
 *   Auxiliary memory: O(n) limbs for the context and temporaries
 *   Output memory: O(n) limbs
 *
 * @param[out] r Result of a^b mod m.
 * @param[in]  a Base.
 * @param[in]  b Exponent.
 * @param[in]  m Modulus.
 */
void bn_mod_exp(bignum* r, const bignum* a, const bignum* b, const bignum* m);

/**
 * @brief Calculate a^b mod m into r using a caller-provided Montgomery context.
 *
 * Let n = m->size and e = b->size, measured in 64-bit limbs.
 *
 * This computes:
 *
 *   r = a^b mod m
 *
 * by converting a into the Montgomery domain (bn_mont_in), exponentiating
 * with bn_mont_exp(), and converting the result back (bn_mont_out).
 *
 * Precondition: m is nonzero and odd, and ctx was initialized with
 * bn_mont_ctx_init() for this m.
 *
 * Complexity:
 *   Time: O(e * n^2) - O(e) Montgomery multiplications plus two
 *         conversions, each one Montgomery multiplication
 *   Auxiliary memory: O(n) limbs for temporaries
 *   Output memory: O(n) limbs
 *
 * @param[out] r     Result of a^b mod m.
 * @param[in]  a     Base.
 * @param[in]  b     Exponent.
 * @param[in]  m     Modulus (nonzero and odd).
 * @param[in]  ctx   Montgomery context initialized for m.
 */
void bn_mod_exp_mont(bignum* r, const bignum* a, const bignum* b,
                     const bignum* m, bn_mont_ctx* ctx);

// bigshift.c
/**
 * @brief Shift r left by one bit, in place: r = r * 2.
 *
 * Let n = r->size, measured in 64-bit limbs.
 *
 * The bits are propagated from the least significant limb to the most
 * significant one; if the top bit overflows, a new limb is appended.
 * The sign is preserved.
 *
 * Complexity:
 *   Time: O(n)
 *   Auxiliary memory: O(1)
 *   Output memory: O(n) limbs, or O(n + 1) if a carry limb is created
 *
 * @param[in,out] r Value to shift left by one bit (modified in place).
 */
void bn_lshift1(bignum* r);

/**
 * @brief Shift r right by one bit, in place: r = r / 2 (truncated).
 *
 * Let n = r->size, measured in 64-bit limbs.
 *
 * The bits are propagated from the most significant limb to the least
 * significant one; leading zero limbs are trimmed afterwards. The sign
 * is preserved.
 *
 * Complexity:
 *   Time: O(n)
 *   Auxiliary memory: O(1)
 *   Output memory: O(n) limbs
 *
 * @param[in,out] r Value to shift right by one bit (modified in place).
 */
void bn_rshift1(bignum* r);

/**
 * @brief Shift a left by shift bits into r: r = a << shift.
 *
 * Let n = a->size, measured in 64-bit limbs.
 * Let w = shift / 64 (whole limbs) and b = shift % 64 (remaining bits).
 *
 * Each limb of a is shifted left by b bits and placed w limbs higher
 * in r, with the overflow bits carried into the next limb. The sign is
 * preserved.
 *
 * r may alias a.
 *
 * Complexity:
 *   Time: O(n)
 *   Auxiliary memory: O(1) in the normal case,
 *                     O(n) if r aliases a and a temporary is used
 *   Output memory: O(n) limbs, at most n + w + 1 limbs
 *
 * @param[out] r     Result of a << shift.
 * @param[in]  a     Value to shift.
 * @param[in]  shift Number of bits to shift left.
 */
void bn_lshift(bignum* r, const bignum* a, int shift);

/**
 * @brief Shift a right by shift bits into r: r = a >> shift (truncated).
 *
 * Let n = a->size, measured in 64-bit limbs.
 * Let w = shift / 64 (whole limbs) and b = shift % 64 (remaining bits).
 *
 * The top w limbs are dropped and the remaining limbs are shifted
 * right by b bits, pulling in the low bits of the next limb. If the
 * shift is at least the size of a, the result is zero. The sign is
 * preserved.
 *
 * r may alias a.
 *
 * Complexity:
 *   Time: O(n)
 *   Auxiliary memory: O(1) in the normal case,
 *                     O(n) if r aliases a and a temporary is used
 *   Output memory: O(n) limbs, at most n - w limbs
 *
 * @param[out] r     Result of a >> shift.
 * @param[in]  a     Value to shift.
 * @param[in]  shift Number of bits to shift right.
 */
void bn_rshift(bignum* r, const bignum* a, int shift);

/**
 * @brief Shift r left by one bit and add bit (0 or 1), in place:
 * r = (r << 1) + bit.
 *
 * Let n = r->size, measured in 64-bit limbs.
 *
 * This is the primitive used when building numbers bit by bit (e.g.
 * during parsing or in exponentiation loops). The added bit enters at
 * the least significant position and the carry propagates upward; if
 * the top bit overflows, a new limb is appended. The sign is
 * preserved.
 *
 * Complexity:
 *   Time: O(n)
 *   Auxiliary memory: O(1)
 *   Output memory: O(n) limbs, or O(n + 1) if a carry limb is created
 *
 * @param[in,out] r   Value to shift and add into (modified in place).
 * @param[in]     bit Bit (0 or 1) to add at the least significant position.
 */
void bn_lshift1_add(bignum* r, int bit);

// bigsqrt.c
/**
 * @brief Calculate the integer square root of a using Heron's method.
 *
 * Let n = a->size, measured in 64-bit limbs.
 *
 * This computes:
 *
 *   r = floor(sqrt(a))
 *
 * The iteration starts at x0 = 2^(ceil(k/2) + 1) where k is the bit
 * length of a, and repeatedly applies:
 *
 *   x_{i+1} = (x_i + a / x_i) / 2
 *
 * until the sequence stops decreasing. Each iteration roughly doubles
 * the number of correct bits, so O(log n) iterations suffice.
 *
 * If a is negative, r is left unchanged (no real square root exists).
 *
 * Complexity:
 *   Time: O(n^2 log n) - O(log n) iterations, each dominated by a
 *         division of an n-limb value by an n-limb value
 *   Auxiliary memory: O(n) limbs for temporaries
 *   Output memory: O(n) limbs (at most n/2 + 1 limbs)
 *
 * @param[out] r Result storing floor(sqrt(a)).
 * @param[in]  a Value to take the square root of.
 */
void bn_isqrt_heron(bignum* r, bignum* a);

/**
 * @brief Calculate the integer square root of a.
 *
 * Let n = a->size, measured in 64-bit limbs.
 *
 * This computes:
 *
 *   r = floor(sqrt(a))
 *
 * Currently a thin wrapper around bn_isqrt_heron().
 *
 * Complexity:
 *   Time: O(n^2 log n), see bn_isqrt_heron()
 *   Auxiliary memory: O(n) limbs
 *   Output memory: O(n) limbs
 *
 * @param[out] r Result storing floor(sqrt(a)).
 * @param[in]  a Value to take the square root of.
 */
void bn_isqrt(bignum* r, bignum* a);

// biglog.c
/**
 * @brief Calculate the integer natural logarithm of a, truncated.
 *
 * Let n = a->size, measured in 64-bit limbs.
 *
 * This computes:
 *
 *   r = floor(ln(a))
 *
 * using the identity ln(a) = ln(2) * log_2(a). The base-2 logarithm is
 * exact (see bn_log_2), and ln(2) is approximated by the fixed-point
 * rational 69314718 / 100000000 (8 decimal digits), so the result is
 * accurate to within roughly 1 for large a.
 *
 * If a is zero, r is left unchanged (ln(0) is undefined).
 *
 * Complexity:
 *   Time: O(n) for the bit length, plus O(n^2) for the bn_mul() and
 *         O(n^2) for the bn_div() of two n-limb values
 *   Auxiliary memory: O(n) limbs for temporaries
 *   Output memory: O(1) limbs (the result fits in a few limbs)
 *
 * @param[out] r Result storing floor(ln(a)).
 * @param[in]  a Value to take the natural logarithm of.
 */
void bn_ln(bignum* r, bignum* a);

/**
 * @brief Calculate the integer base-2 logarithm of a.
 *
 * Let n = a->size, measured in 64-bit limbs.
 *
 * This computes:
 *
 *   r = floor(log_2(a))
 *
 * which is exact and equals the bit length of a minus one.
 *
 * If a is zero, r is left unchanged (log_2(0) is undefined).
 *
 * Complexity:
 *   Time: O(1) (only the most significant limb is inspected)
 *   Auxiliary memory: O(1)
 *   Output memory: O(1) limbs
 *
 * @param[out] r Result storing floor(log_2(a)).
 * @param[in]  a Value to take the base-2 logarithm of.
 */
void bn_log_2(bignum* r, bignum* a);

// bigmont.c
/**
 * @brief Initialize a Montgomery context for the modulus n.
 *
 * Let n_l = n->size, measured in 64-bit limbs.
 *
 * Computes the context fields:
 *
 *   ctx->n        = n
 *   ctx->n_inv    = -n^{-1} mod 2^64   (from the low limb of n)
 *   ctx->one_mont = R mod n            (R = 2^(64 * n_l))
 *   ctx->r_square = R^2 mod n
 *
 * and preallocates ctx->tmp with 2 * n_l + 1 limbs of scratch space.
 *
 * Precondition: n is odd and greater than 1.
 *
 * Complexity:
 *   Time: O(n_l^2) - one_mont = R mod n costs one division, r_square
 *         costs one multiply plus one division
 *   Auxiliary memory: O(n_l) limbs
 *   Output memory: O(n_l) limbs per context field
 *
 * @param[out] ctx Context to initialize.
 * @param[in]  n   Modulus (odd and greater than 1).
 */
void bn_mont_ctx_init(bn_mont_ctx* ctx, const bignum* n);

/**
 * @brief Free all bignum storage held by a Montgomery context.
 *
 * After this call the context must not be used until re-initialized.
 *
 * Complexity:
 *   Time: O(1)
 *   Auxiliary memory: O(1)
 *   Output memory: O(1)
 *
 * @param[in,out] ctx Context to free.
 */
void bn_mont_ctx_free(bn_mont_ctx* ctx);

/**
 * @brief Montgomery reduction (REDC).
 *
 * Let n_l = ctx->n.size, measured in 64-bit limbs.
 *
 * Given t with 0 <= t < n * 2^(64 * n_l), this computes:
 *
 *   r = t * R^{-1} mod n
 *
 * in O(n_l^2) time. For each limb i it forms the multiple
 * m = t[i] * (-n^{-1} mod 2^64) and adds m * n shifted by i limbs to
 * t, which zeroes out limb i (mod 2^64). After n_l steps the lower
 * half of t is zero, so the result is the upper half, followed by a
 * final conditional subtraction of n to bring r into [0, n).
 *
 * t is destroyed (overwritten) and must have at least 2 * n_l + 1
 * limbs of capacity.
 *
 * Complexity:
 *   Time: O(n_l^2) for the reduction loop, plus O(n_l) for the final
 *         correction
 *   Auxiliary memory: O(1)
 *   Output memory: O(n_l) limbs
 *
 * @param[out]    r   Result t * R^{-1} mod n, in [0, n).
 * @param[in,out] t   Value to reduce (destroyed; needs 2 * n_l + 1 limbs).
 * @param[in]     ctx Initialized Montgomery context.
 */
void bn_mont_redc(bignum* r, bignum* t, bn_mont_ctx* ctx);

/**
 * @brief Convert a value from the normal domain into the Montgomery domain.
 *
 * Let n_l = ctx->n.size, measured in 64-bit limbs.
 *
 * This computes:
 *
 *   A_bar = A * R mod n
 *
 * as one Montgomery multiplication of A by R^2 mod n:
 * REDC(A * R^2) = A * R mod n.
 *
 * Complexity:
 *   Time: O(n_l^2)
 *   Auxiliary memory: O(n_l) limbs (ctx->tmp)
 *   Output memory: O(n_l) limbs
 *
 * @param[out] A_bar Result in Montgomery form: A * R mod n.
 * @param[in]  A     Value in the normal domain.
 * @param[in]  ctx   Initialized Montgomery context.
 */
void bn_mont_in(bignum* A_bar, const bignum* A, bn_mont_ctx* ctx);

/**
 * @brief Convert a value from the Montgomery domain back to the normal domain.
 *
 * Let n_l = ctx->n.size, measured in 64-bit limbs.
 *
 * This computes:
 *
 *   A = A_bar * R^{-1} mod n
 *
 * as one Montgomery multiplication of A_bar by 1.
 *
 * Complexity:
 *   Time: O(n_l^2)
 *   Auxiliary memory: O(n_l) limbs (ctx->tmp and a temporary for 1)
 *   Output memory: O(n_l) limbs
 *
 * @param[out] A     Result in the normal domain: A_bar * R^{-1} mod n.
 * @param[in]  A_bar Value in the Montgomery domain.
 * @param[in]  ctx   Initialized Montgomery context.
 */
void bn_mont_out(bignum* A, const bignum* A_bar, bn_mont_ctx* ctx);

/**
 * @brief Montgomery multiplication of two values in the Montgomery domain.
 *
 * Let n_l = ctx->n.size, measured in 64-bit limbs.
 *
 * This computes:
 *
 *   r = a_bar * b_bar * R^{-1} mod n
 *
 * so that if a_bar = A * R mod n and b_bar = B * R mod n, then
 * r = A * B * R mod n (the Montgomery form of A * B).
 *
 * Thin wrapper around bn_mont_mul_raw().
 *
 * Complexity:
 *   Time: O(n_l^2)
 *   Auxiliary memory: O(n_l) limbs (ctx->tmp)
 *   Output memory: O(n_l) limbs
 *
 * @param[out] r     Result a_bar * b_bar * R^{-1} mod n.
 * @param[in]  a_bar First operand in Montgomery form.
 * @param[in]  b_bar Second operand in Montgomery form.
 * @param[in]  ctx   Initialized Montgomery context.
 */
void bn_mont_mul(bignum* r, const bignum* a_bar, const bignum* b_bar,
                 bn_mont_ctx* ctx);

/**
 * @brief Montgomery multiplication using the context's scratch buffer.
 *
 * Let n_l = ctx->n.size, measured in 64-bit limbs.
 *
 * This computes:
 *
 *   result = A_bar * B_bar * R^{-1} mod n
 *
 * by first forming the full product T = A_bar * B_bar (at most
 * 2 * n_l limbs) in ctx->tmp, zero-padding it to 2 * n_l + 1 limbs,
 * and applying bn_mont_redc().
 *
 * Note: uses ctx->tmp as scratch, so it is not reentrant and must not
 * be called with A_bar or B_bar aliasing ctx->tmp.
 *
 * Complexity:
 *   Time: O(n_l^2) - one bignum multiply plus one REDC
 *   Auxiliary memory: O(n_l) limbs (ctx->tmp)
 *   Output memory: O(n_l) limbs
 *
 * @param[out] result  Result A_bar * B_bar * R^{-1} mod n.
 * @param[in]  A_bar   First operand in Montgomery form.
 * @param[in]  B_bar   Second operand in Montgomery form.
 * @param[in]  ctx     Initialized Montgomery context (provides scratch).
 */
void bn_mont_mul_raw(bignum* result, const bignum* A_bar, const bignum* B_bar,
                     bn_mont_ctx* ctx);

#endif
