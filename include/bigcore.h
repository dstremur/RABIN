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
 *   - Time: \f$O(1)\f$
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(1)\f$ limbs
 */
void bn_init_constants();

/**
 * @brief Free the limb storage of the global constants.
 *
 * After this call the constants are uninitialized and must not be used
 * until bn_init_constants() is called again.
 *
 * Complexity:
 *   - Time: \f$O(1)\f$
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(1)\f$
 */
void bn_free_constants();

// bnscratch.c
/**
 * @brief Reserve n u64s of thread-local scratch space.
 *
 * Returns a pointer to at least n u64s. The memory is NOT zeroed.
 * The block stays valid until bn_scratch_release() is called.
 *
 * Complexity:
 *   - Time: \f$O(1)\f$ amortized, \f$O(1)\f$ realloc on arena growth
 *   - Auxiliary memory: \f$O(n)\f$ u64s (grow-only per thread)
 *   - Output memory: \f$O(n)\f$ u64s
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
 *   - Time: \f$O(1)\f$
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(1)\f$
 */
void bn_scratch_release(void);

// bignum.c
/**
 * @brief Initialize a bignum to zero, releasing no memory (r must not have
 * been allocated yet).
 *
 * Complexity:
 *   - Time: \f$O(1)\f$
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(1)\f$
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
 *   - Time: \f$O(k)\f$ for k bignums
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(1)\f$
 *
 * @param[out] r First bignum to initialize; further bignums follow as
 *               variadic arguments, terminated by NULL.
 */
void bn_init_multi(bignum* r, ...);

/**
 * @brief Swap the contents of two bignums by exchanging their structs.
 *
 * Only the struct fields (pointers, sizes, sign) are swapped, so this
 * is \f$O(1)\f$ regardless of size.
 *
 * Complexity:
 *   - Time: \f$O(1)\f$
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(1)\f$
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
 *   - Time: \f$O(1)\f$
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(1)\f$
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
 *   - Time: \f$O(1)\f$
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(1)\f$
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
 *   - Time: \f$O(1)\f$
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(1)\f$
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
 * storage is reallocated with exponential growth and
 * the newly added limbs are zeroed. Existing limb contents are
 * preserved.
 *
 * Returns true on success, false if the realloc fails (r is left
 * unchanged).
 *
 * Complexity:
 *   - Time: \f$O(capacity)\f$ for the zero-fill and the realloc copy
 *   - Auxiliary memory: \f$O(capacity)\f$ during the realloc
 *   - Output memory: \f$O(capacity)\f$ limbs
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
 * Let \f$d =\f$ number of decimal digits in str.
 *
 * Frees any previous contents of n and sets it to the value of str,
 * which may start with a '-' sign. Non-digit characters are skipped.
 * Digits are processed one at a time with a multiply-by-10-and-add
 * carry loop over the limbs.
 *
 * Complexity:
 *   - Time: \f$O(d \cdot n)\f$ where \f$n =\f$ the number of limbs, i.e.
 * \f$O(d^2 / 64)\f$
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(n)\f$ limbs
 *
 * @param[out] n   Bignum to set (previous contents are freed).
 * @param[in]  str Base-10 decimal string to parse (may start with '-').
 */
void bn_init_val(bignum* n, const char* str);

/**
 * @brief Convert a bignum to a base-10 decimal string.
 *
 * Let \f$n =\f$ a->size, measured in 64-bit limbs.
 *
 * Returns a newly allocated string (caller must free it) containing
 * the decimal representation of n, with a leading '-' for negative
 * values. The magnitude is reduced repeatedly by \f$10^{19}\f$ (which fits in
 * a u64), collecting 19-digit chunks; the most significant chunk is
 * printed without padding and the rest with zero padding.
 *
 * Returns NULL on allocation failure.
 *
 * Complexity:
 *   - Time: \f$O(n^2)\f$ - \f$O(n)\f$ divisions by a 64-bit divisor, each
 * \f$O(n)\f$
 *   - Auxiliary memory: \f$O(n)\f$ limbs for temporaries, \f$O(n)\f$ for the
 * chunk array and the output string
 *   - Output memory: \f$O(n)\f$ characters
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
 * Let \f$n =\f$ a->size, measured in 64-bit limbs.
 *
 * Same chunking strategy as bn_to_string(): the magnitude is reduced
 * repeatedly by \f$10^{19}\f$ and the chunks are printed most-significant
 * first, with a leading '-' for negative values.
 *
 * Complexity:
 *   - Time: \f$O(n^2)\f$ - \f$O(n)\f$ divisions by a 64-bit divisor, each
 * \f$O(n)\f$
 *   - Auxiliary memory: \f$O(n)\f$ limbs for temporaries, \f$O(n)\f$ for the
 * chunk array
 *   - Output memory: \f$O(n)\f$ characters written
 *
 * @param[in] n Bignum to print.
 */
void bn_print(const bignum* n);

/**
 * @brief Print a bignum to stdout in base 10 followed by a newline.
 *
 * Complexity:
 *   - Time: \f$O(n^2)\f$, see bn_print()
 *   - Auxiliary memory: \f$O(n)\f$
 *   - Output memory: \f$O(n)\f$ characters written
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
 *   - Time: \f$O(1)\f$
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(1)\f$
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
 *   - Time: \f$O(k)\f$ for k bignums
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(1)\f$
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
 *   - Time: \f$O(number of trimmed limbs)\f$
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(1)\f$
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
 *   - Time: \f$O(1)\f$
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(1)\f$
 *
 * @param[in] a Bignum to read the bit from.
 * @param[in] i Bit index (\f$0 =\f$ least significant).
 *
 * @return The value of bit i (0 or 1).
 */
int bn_get_bit(const bignum* a, int i);

/**
 * @brief Compare two signed bignums.
 *
 * Returns 1 if \f$a > b\f$, -1 if \f$a < b\f$, 0 if \f$a = b\f$.
 *
 * Signs are handled first; for equal signs the magnitudes are compared, most
 * significant limb first, and the result is negated when both
 * operands are negative.
 *
 * Complexity:
 *   - Time: \f$O(n)\f$ worst case, where \f$n =\f$ max(a->size, b->size)
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(1)\f$
 *
 * @param[in] a First bignum.
 * @param[in] b Second bignum.
 *
 * @return 1 If \f$a > b\f$, -1 if \f$a < b\f$, 0 if \f$a = b\f$.
 */
int bn_cmp(const bignum* a, const bignum* b);

/**
 * @brief Compare two signed bignums in constant time
 *
 * Returns 1 if \f$a > b\f$, -1 if \f$a < b\f$, 0 if \f$a = b\f$.
 *
 * Same semantics as bn_cmp, but this function always compares
 * all limbs to ensure constant time.
 *
 * Complexity:
 *   - Time: \f$O(n)\f$ always, where \f$n =\f$ max(a->size, b->size)
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(1)\f$
 *
 * @param[in] a First bignum.
 * @param[in] b Second bignum.
 *
 * @return 1 If \f$a > b\f$, -1 if \f$a < b\f$, 0 if \f$a = b\f$.
 */
int bn_cmp_const_time(const bignum* a, const bignum* b);

/**
 * @brief Compare the absolute values (magnitudes) of two bignums.
 *
 * Returns 1 if \f$|a| > |b|\f$, -1 if \f$|a| < |b|\f$, 0 if \f$|a| = |b|\f$.
 * Signs are ignored. The comparison is by size first, then most significant
 * limb first.
 *
 * Complexity:
 *   - Time: \f$O(n)\f$ worst case, where \f$n =\f$ max(a->size, b->size)
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(1)\f$
 *
 * @param[in] a First bignum.
 * @param[in] b Second bignum.
 *
 * @return 1 If \f$|a| > |b|\f$, -1 if \f$|a| < |b|\f$, 0 if \f$|a| = |b|\f$.
 */
int bn_cmp_abs(const bignum* a, const bignum* b);

/**
 * @brief Deep copy a bignum: \f$r = a\f$.
 *
 * Copies the limb storage and the sign. r may alias a (a no-op in
 * that case). If a is zero-sized, r is reset to size 0.
 *
 * Complexity:
 *   - Time: \f$O(n)\f$ where \f$n =\f$ a->size
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(n)\f$ limbs
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
 *   - Time: \f$O(1)\f$
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(1)\f$ limbs
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
 *   - Time: \f$O(1)\f$
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(1)\f$ limbs
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
 *   - Time: \f$O(1)\f$ amortized (\f$O(limb)\f$ for the zero-fill on growth)
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(1)\f$ limbs, possibly grown
 *
 * @param[in,out] a Bignum to modify.
 * @param[in]     i Bit index (\f$0 =\f$ least significant).
 */
void bn_set_bit(bignum* a, int i);

/**
 * @brief Clear the i-th bit of a to 0.
 *
 * Bit 0 is the least significant bit. Bits beyond the current size
 * are already 0 and are left untouched.
 *
 * Complexity:
 *   - Time: \f$O(1)\f$
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(1)\f$
 *
 * @param[in,out] a Bignum to modify.
 * @param[in]     i Bit index (\f$0 =\f$ least significant).
 */
void bn_clear_bit(bignum* a, int i);

/**
 * @brief Return the bit length of a: the index of the highest set bit plus
 * one (0 for zero).
 *
 * Only the most significant limb is inspected.
 *
 * Complexity:
 *   - Time: \f$O(1)\f$
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(1)\f$
 *
 * @param[in] a Bignum to measure.
 *
 * @return The bit length of \f$a\f$ (index of highest set bit + \f$1\f$, or
 * \f$0\f$).
 */
int bn_bit_length(const bignum* a);

/**
 * @brief Return the number of trailing zero bits of a.
 *
 * Counts whole zero limbs (64 bits each) plus the trailing zeros of
 * the first nonzero limb. Returns 0 for \f$a = 0\f$.
 *
 * Complexity:
 *   - Time: \f$O(number of zero limbs)\f$
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(1)\f$
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
 * Let \f$n =\f$ max(a->size, b->size), measured in 64-bit limbs.
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
 *   - Time: \f$O(n)\f$
 *   - Auxiliary memory: \f$O(1)\f$, except for any temporary storage used by
 *                       bn_add_abs(), bn_sub_abs(), bn_copy(), or bn_alloc()
 *   - Output memory: \f$O(n)\f$ limbs
 *
 * @param[out] r Result of the signed addition \f$a + b\f$.
 * @param[in]  a First operand.
 * @param[in]  b Second operand.
 */
void bn_add(bignum* r, const bignum* a, const bignum* b);

/**
 * @brief Add the absolute values of a and b into r.
 *
 * Let \f$n =\f$ max(a->size, b->size), measured in 64-bit limbs.
 *
 * This function only adds magnitudes. It does not interpret or modify
 * the sign of the operands or the result. The caller is responsible for
 * setting r->is_neg.
 *
 * r may alias a or b.
 *
 * Complexity:
 *   - Time: \f$O(n)\f$
 *   - Auxiliary memory: \f$O(1)\f$ in the normal case,
 *                       \f$O(n)\f$ if r aliases a or b and a temporary is used
 *   - Output memory: \f$O(n)\f$ limbs, at most \f$n + 1\f$ limbs
 *
 * @param[out] r Result storing the magnitude sum \f$|a| + |b|\f$.
 * @param[in]  a First operand (magnitude only).
 * @param[in]  b Second operand (magnitude only).
 */
void bn_add_abs(bignum* r, const bignum* a, const bignum* b);

/**
 * @brief Add an unsigned 64-bit value to a signed bignum.
 *
 * Let \f$n =\f$ a->size, measured in 64-bit limbs.
 *
 * This computes:
 *
 *   \f$r = a + b\f$
 *
 * where \f$b\f$ is a nonnegative u64.
 *
 * If a is positive, this is ordinary magnitude addition.
 *
 * If a is negative, this is magnitude subtraction:
 *
 *   \f$r = -|a| + b\f$
 *
 * which is equivalent to subtracting \f$b\f$ from \f$|a|\f$ and preserving the
 * correct sign.
 *
 * r may alias a, assuming bn_copy() supports aliasing.
 *
 * Complexity:
 *   - Time: \f$O(n)\f$ worst case.
 *           The carry/borrow loop may stop early if propagation ends
 *           before reaching the most significant limb.
 *   - Auxiliary memory: \f$O(1)\f$, except for any temporary storage used by
 *                       bn_copy() or bn_alloc()
 *   - Output memory: \f$O(n)\f$ limbs, or \f$O(n + 1)\f$ if a new carry limb is
 * created
 *
 * @param[out] r Result of the addition \f$a + b\f$.
 * @param[in]  a Signed bignum operand.
 * @param[in]  b Nonnegative 64-bit value to add.
 */
void bn_add_u64(bignum* r, const bignum* a, u64 b);

/**
 * @brief Add a to r at the given limb offset.
 *
 * Let \f$n =\f$ a->size, measured in 64-bit limbs.
 * Let \f$o =\f$ offset.
 * Let \f$m =\f$ o + n + 1 be the maximum resulting size in limbs.
 *
 * This effectively performs:
 *
 *   \f$r = r + (a \ll (offset \cdot 64))\f$
 *
 * for unsigned magnitudes.
 *
 * The addition loop itself only touches about \f$n + 1\f$ limbs, but allocation
 * and trimming may need to consider the full output range up to m limbs.
 *
 * Complexity:
 *   - Time: \f$O(m)\f$ worst case, where \f$m =\f$ offset + a->size + 1.
 *           The inner addition loop is \f$O(n)\f$, but bn_alloc() and bn_trim()
 *           may make the total worst case proportional to the output size.
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(m)\f$ limbs
 *
 * Note:
 *   The current implementation does not explicitly handle \f$r = a\f$ safely
 *   for all nonzero offsets.
 *
 * @param[out] r Accumulator; receives \f$r + (a \ll (offset \cdot 64))\f$.
 * @param[in]  a Value to add.
 * @param[in]  offset Limb offset at which a is added.
 */
void bn_add_at_offset(bignum* r, const bignum* a, u64 offset);

// bigsub.c

/**
 * @brief Subtract two signed bignums: \f$r = a - b\f$.
 *
 * Let \f$n =\f$ max(a->size, b->size), measured in 64-bit limbs.
 *
 * Subtraction is reduced to magnitude addition or subtraction:
 *
 *   \f$a - (-b) = a + b\f$          (opposite signs, \f$b\f$ negative)
 *   \f$(-a) - b = -(a + b)\f$       (opposite signs, \f$a\f$ negative)
 *   same signs: subtract the smaller magnitude from the larger one and
 *               take the sign of the operand with the larger magnitude
 *
 * If the magnitudes are equal, the result is zero.
 *
 * Complexity:
 *   - Time: \f$O(n)\f$
 *   - Auxiliary memory: \f$O(1)\f$, except for any temporary storage used by
 *                       bn_add_abs(), bn_sub_abs(), or bn_alloc()
 *   - Output memory: \f$O(n)\f$ limbs
 *
 * @param[out] r Result of the signed subtraction \f$a - b\f$.
 * @param[in]  a Minuend.
 * @param[in]  b Subtrahend.
 */
void bn_sub(bignum* r, const bignum* a, const bignum* b);

/**
 * @brief Subtract the absolute value of b from the absolute value of a into r.
 *
 * Let \f$n =\f$ a->size, measured in 64-bit limbs.
 *
 * This computes:
 *
 *   \f$r = |a| - |b|\f$
 *
 * and requires \f$|a| \ge |b|\f$; if \f$|b| > |a|\f$ the result wraps around
 * (two's-complement style) and is meaningless. It does not interpret or
 * modify the sign of the operands or the result. The caller is
 * responsible for setting r->is_neg.
 *
 * r may alias a or b.
 *
 * Complexity:
 *   - Time: \f$O(n)\f$
 *   - Auxiliary memory: \f$O(1)\f$ in the normal case,
 *                       \f$O(n)\f$ if r aliases a or b and a temporary is used
 *   - Output memory: \f$O(n)\f$ limbs
 *
 * @param[out] r Result storing \f$|a| - |b|\f$ (requires \f$|a| \ge |b|\f$).
 * @param[in]  a Minuend (magnitude only).
 * @param[in]  b Subtrahend (magnitude only).
 */
void bn_sub_abs(bignum* r, const bignum* a, const bignum* b);

// bigmul.c

/**
 * @brief Multiply two signed bignums: \f$r = a \cdot b\f$.
 *
 * Let \f$n =\f$ max(a->size, b->size), measured in 64-bit limbs.
 *
 * This computes the signed product of a and b. The sign of the result
 * is the xor of the operand signs; the magnitude is computed with:
 *
 *   - the squaring path (bn_sqr) when \f$a = b\f$
 *   - Karatsuba (\f$O(n^{1.585})\f$) when both operands have at least
 *     KARATSUBA_LIMIT limbs
 *   - schoolbook (\f$O(n^2)\f$) otherwise
 *
 * r may alias a or b.
 *
 * Complexity:
 *   - Time: \f$O(n^1.585)\f$ for large operands, \f$O(n^2)\f$ for small ones
 *   - Auxiliary memory: \f$O(n)\f$ limbs of scratch, plus \f$O(n)\f$ for the
 *                       zero-padding buffers on the Karatsuba path and
 *                       \f$O(n)\f$ if r aliases an operand
 *   - Output memory: \f$O(n)\f$ limbs (at most a->size + b->size)
 *
 * @param[out] r Result of the signed product \f$a \cdot b\f$.
 * @param[in]  a First operand.
 * @param[in]  b Second operand.
 */
void bn_mul(bignum* r, const bignum* a, const bignum* b);

/**
 * @brief Multiply a bignum a with a 64-bit integer c
 *
 * @param r
 * @param a
 * @param c
 */
void bn_mul_i64(bignum* r, const bignum* a, const i64 c);

/**
 * @brief Multiply two signed bignums with the schoolbook (grade-school)
 * algorithm: \f$r = a \cdot b\f$.
 *
 * Let \f$n =\f$ a->size and \f$m =\f$ b->size, measured in 64-bit limbs.
 *
 * Each limb of a is multiplied by the whole of b and accumulated into
 * r at the corresponding offset (bn_mul_add_inner), skipping zero
 * limbs. The sign of the result is the xor of the operand signs.
 *
 * Complexity:
 *   - Time: \f$O(n \cdot m)\f$
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(n + m)\f$ limbs
 *
 * @param[out] r Result of the schoolbook product \f$a \cdot b\f$.
 * @param[in]  a First operand.
 * @param[in]  b Second operand.
 */
void bn_mul_school(bignum* r, const bignum* a, const bignum* b);

/**
 * @brief Multiply two bignums with the NTT-based fast path: \f$res = a \cdot
 * b\f$.
 *
 * Let \f$n =\f$ max(a->size, b->size), measured in 64-bit limbs.
 *
 * This computes the product of the magnitudes of a and b by:
 *
 *   1. decomposing each operand into a polynomial whose coefficients
 *      are 16-bit chunks (base \f$2^{16}\f$),
 *   2. multiplying the polynomials with a cyclic NTT
 *      (bigpoly_mul_ntt),
 *   3. propagating carries between the 16-bit coefficient slots,
 *   4. recomposing the coefficients back into a bignum.
 *
 * Unlike bn_mul() this path is not wired into the general dispatch;
 * it is a standalone fast path for very large operands.
 *
 * Complexity:
 *   - Time: \f$O(n \log n)\f$ for the NTT, plus \f$O(n)\f$ for decompose/carry/
 *           recompose
 *   - Auxiliary memory: \f$O(n)\f$ limbs for the polynomial arrays
 *   - Output memory: \f$O(n)\f$ limbs (at most a->size + b->size)
 *
 * @param[out] res Result of the NTT-based product \f$a \cdot b\f$.
 * @param[in]  a   First operand.
 * @param[in]  b   Second operand.
 */
void bn_mul_fast(bignum* res, const bignum* a, const bignum* b);

/**
 * @brief Square a bignum: \f$r = a \cdot a\f$.
 *
 * Let \f$n =\f$ a->size, measured in 64-bit limbs.
 *
 * This computes the square of the magnitude of a using the Karatsuba
 * squaring kernel (limbs_sqr_karatsuba), which exploits the symmetry
 * of squaring to save about a quarter of the multiplications. The
 * result is always nonnegative.
 *
 * If a is zero, r is set to zero.
 *
 * Complexity:
 *   - Time: \f$O(n^1.585)\f$ for \f$n \ge\f$ KARATSUBA_LIMIT, \f$O(n^2)\f$
 * below it
 *   - Auxiliary memory: \f$O(n)\f$ limbs of scratch
 *   - Output memory: \f$O(n)\f$ limbs (at most 2n)
 *
 * @param[out] r Result storing \f$a \cdot a\f$.
 * @param[in]  a Value to square.
 */
void bn_sqr(bignum* r, const bignum* a);

// bigdiv.c

/**
 * @brief \f$q = a / b\f$ (truncated division, C semantics).
 *
 * Let \f$n =\f$ a->size and \f$m =\f$ b->size, measured in 64-bit limbs.
 *
 * Thin wrapper around bn_divmod() that discards the remainder.
 *
 * Complexity:
 *   - Time: \f$O(n)\f$ for a single-limb divisor, \f$O((n - m + 1) \cdot m)\f$
 * otherwise
 *   - Auxiliary memory: \f$O(n + m)\f$ limbs
 *   - Output memory: \f$O(n - m + 1)\f$ limbs
 *
 * @param[out] q Quotient \f$a / b\f$.
 * @param[in]  a Dividend.
 * @param[in]  b Divisor (must be nonzero).
 */
void bn_div(bignum* q, const bignum* a, const bignum* b);

/**
 * @brief \f$q = a / b\f$, \f$r = a \% b\f$  (either result may be NULL).
 *
 * Let \f$n =\f$ a->size and \f$m =\f$ b->size, measured in 64-bit limbs.
 *
 * Truncated division (C semantics): \f$q\f$ truncates toward 0 and
 * \f$sign(r) = sign(a)\f$.
 *
 * If \f$|a| < |b|\f$ the result is \f$q = 0\f$, \f$r = a\f$. Otherwise the
 * magnitude division dispatches to the single-limb path (limbs_divrem_1) or to
 * Knuth's Algorithm D (bn_divmod_limbs). \f$q\f$ and \f$r\f$ may each alias
 * \f$a\f$ or
 * \f$b\f$ (handled via temporaries); \f$q\f$ and \f$r\f$ must not alias each
 * other.
 *
 * Scratch for the multi-limb path is taken from the stack when it fits
 * in BN_DIV_STACK_LIMBS limbs, otherwise from the heap.
 *
 * Complexity:
 *   - Time: \f$O(n)\f$ for a single-limb divisor, \f$O((n - m + 1) \cdot m)\f$
 * for a multi-limb divisor
 *   - Auxiliary memory: \f$O(n + m)\f$ limbs of scratch (stack or heap), plus
 *                       \f$O(n)\f$ if \f$q\f$ or \f$r\f$ aliases an operand
 *   - Output memory: \f$O(n - m + 1)\f$ limbs for \f$q\f$, \f$O(m)\f$ for
 * \f$r\f$
 *
 * @param[out] q Quotient \f$a / b\f$ (may be NULL).
 * @param[out] r Remainder \f$a \% b\f$ (may be NULL).
 * @param[in]  a Dividend.
 * @param[in]  b Divisor (must be nonzero).
 */
void bn_divmod(bignum* q, bignum* r, const bignum* a, const bignum* b);

/**
 * @brief Computes the division \f$a / d\f$ using a Newton-iterated
 * reciprocal: \f$q = a / d\f$.
 *
 * Let \f$n =\f$ a->size and \f$m =\f$ d->size, measured in 64-bit limbs.
 *
 * Works with \f$P =\f$ bit_length(a) + 32 bits of precision. The reciprocal
 * \f$x = 2^P / d\f$ is seeded with one real division, then refined with the
 * Newton iteration
 *
 *   \f$x \leftarrow x + (x \cdot (2^P - d\cdotx)) \gg P\f$
 *
 * which roughly doubles the correct bits each step. The quotient is
 * then \f$q = (a \cdot x) \gg P\f$, followed by a rare off-by-one fix-up using
 * the remainder \f$r = a - q\cdotd\f$.
 *
 * If \f$d\f$ is zero, \f$q\f$ is left unchanged. If \f$|a| < |d|\f$, \f$q\f$ is
 * set to 0.
 *
 * Complexity:
 *   - Time: \f$O(n^2)\f$ - \f$O(\log m)\f$ Newton iterations of n-limb
 * multiplications plus one seeding division
 *   - Auxiliary memory: \f$O(n)\f$ limbs for temporaries
 *   - Output memory: \f$O(n)\f$ limbs
 *
 * @param[out] q Quotient \f$a / d\f$.
 * @param[in]  a Dividend.
 * @param[in]  d Divisor.
 */
void bn_newton_div(bignum* q, const bignum* a, const bignum* d);

/**
 * @brief \f$q = a / b\f$, assuming \f$b\f$ divides \f$a\f$ exactly (Jebelean's
 * exact division).
 *
 * Let \f$n =\f$ a->size and \f$m =\f$ b->size, measured in 64-bit limbs.
 *
 * Strips the common power of two from a and b (so the low limb of the
 * divisor becomes odd), then computes the quotient limbs from the
 * bottom up: with \f$D[0]\f$ odd, each quotient limb \f$q_i = A[i] \cdot
 * D[0]^{-1}\f$
 * (\f$\bmod 2^{64}\f$) is forced, because it must cancel limb \f$i\f$ of the
 * running remainder. No estimation or correction steps are needed, which makes
 * this about 2x faster than a real division.
 *
 * The behaviour is undefined if \f$b\f$ does not divide \f$a\f$. If \f$b\f$ is
 * zero,
 * \f$q\f$ is left unchanged (asserts in debug builds). If \f$a\f$ is zero,
 * \f$q\f$ is set to 0. \f$q\f$ may alias \f$a\f$ or \f$b\f$.
 *
 * Complexity:
 *   - Time: \f$O(n \cdot m)\f$ - one \f$O(m)\f$ submul per quotient limb
 *   - Auxiliary memory: \f$O(n + m)\f$ limbs for the shifted copies
 *   - Output memory: \f$O(n - m + 1)\f$ limbs
 *
 * @param[out] q Quotient \f$a / b\f$ (requires \f$b \mid a\f$).
 * @param[in]  a Dividend.
 * @param[in]  b Divisor (must divide \f$a\f$ exactly).
 */
void bn_div_exact(bignum* q, const bignum* a, const bignum* b);

// bigmod.c

/**
 * @brief \f$r = a \% b\f$ (truncated remainder, C semantics: \f$sign(r) =
 * sign(a)\f$).
 *
 * Let \f$n =\f$ a->size and \f$m =\f$ b->size, measured in 64-bit limbs.
 *
 * Thin wrapper around bn_divmod() that discards the quotient.
 *
 * Complexity:
 *   - Time: \f$O(n)\f$ for a single-limb divisor, \f$O((n - m + 1) \cdot m)\f$
 * otherwise
 *   - Auxiliary memory: \f$O(n + m)\f$ limbs
 *   - Output memory: \f$O(m)\f$ limbs
 *
 * @param[out] r Remainder \f$a \% b\f$.
 * @param[in]  a Dividend.
 * @param[in]  b Divisor (must be nonzero).
 */
void bn_mod(bignum* r, const bignum* a, const bignum* b);

/**
 * @brief Divide a by a 64-bit divisor d, storing the quotient in q and
 * returning the remainder.
 *
 * Let \f$n =\f$ a->size, measured in 64-bit limbs.
 *
 * This computes:
 *
 *   \f$q = a / d\f$,  return value = \f$a \bmod d\f$
 *
 * by processing the limbs of a most-significant-first with a running
 * 128-bit dividend: each step yields one quotient limb and the new
 * remainder. The quotient is truncated toward zero (C semantics) for
 * negative \f$a\f$.
 *
 * \f$q\f$ may alias \f$a\f$.
 *
 * Complexity:
 *   - Time: \f$O(n)\f$
 *   - Auxiliary memory: \f$O(1)\f$ in the normal case,
 *                       \f$O(n)\f$ if \f$q\f$ aliases \f$a\f$ and a temporary
 * is used
 *   - Output memory: \f$O(n)\f$ limbs for \f$q\f$
 *
 * @param[out] q Quotient \f$a / d\f$ (truncated toward zero).
 * @param[in]  a Dividend.
 * @param[in]  d 64-bit divisor (must be nonzero).
 *
 * @return The remainder \f$a \bmod d\f$.
 */
uint64_t bn_divmod_u64(bignum* q, const bignum* a, uint64_t d);

/**
 * @brief Reduce a modulo a 64-bit divisor d.
 *
 * Let \f$n =\f$ a->size, measured in 64-bit limbs.
 *
 * This computes:
 *
 *   \f$a \bmod d\f$
 *
 * by processing the limbs of \f$a\f$ most-significant-first with a running
 * 128-bit remainder: \f$rem = (rem \cdot 2^{64} + limb) \bmod d\f$. The result
 * is returned in \f$[0, d)\f$.
 *
 * For negative \f$a\f$ the result is the nonnegative residue: if \f$a < 0\f$
 * and the magnitude remainder is nonzero, \f$d - rem\f$ is returned (C-style
 * truncating remainder mapped into \f$[0, d)\f$).
 *
 * Complexity:
 *   - Time: \f$O(n)\f$
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(1)\f$
 *
 * @param[in] a Bignum to reduce.
 * @param[in] d 64-bit divisor (must be nonzero).
 *
 * @return The residue \f$a \bmod d\f$, in the range \f$[0, d)\f$.
 */
uint64_t bn_mod_u64(const bignum* a, uint64_t d);

/**
 * @brief Compute the modular multiplicative inverse of a modulo m.
 *
 * Let \f$n =\f$ m->size, measured in 64-bit limbs.
 *
 * This computes:
 *
 *   \f$res = a^{-1} \bmod m\f$
 *
 * i.e. the value in \f$[0, m)\f$ such that \f$(a \cdot res) \bmod m = 1\f$,
 * using the binary extended GCD (Stein's algorithm) with coefficient tracking.
 * The coefficients are kept reduced modulo \f$m\f$ at every step so they
 * never grow beyond \f$n\f$ limbs.
 *
 * Returns true and stores the inverse in \f$res\f$ if it exists (i.e.
 * \f$\gcd(a, m) = 1\f$). Returns false and leaves \f$res\f$ unchanged if
 * \f$a\f$ is zero, \f$m\f$ is zero or one, or \f$a\f$ and \f$m\f$ are not
 * coprime.
 *
 * Complexity:
 *   - Time: \f$O(n^2)\f$ - \f$O(n)\f$ iterations of shifts and subtractions of
 *           n-limb values (binary GCD), each \f$O(n)\f$
 *   - Auxiliary memory: \f$O(n)\f$ limbs for temporaries
 *   - Output memory: \f$O(n)\f$ limbs
 *
 * @param[out] res Receives \f$a^{-1} \bmod m\f$ if the inverse exists.
 * @param[in]  a   Value to invert.
 * @param[in]  m   Modulus.
 *
 * @return true  If the inverse exists (\f$\gcd(a, m) = 1\f$); \f$res\f$ is set.
 * @return false If \f$a\f$ is zero, \f$m\f$ is zero or one, or \f$a\f$ and
 * \f$m\f$ are not coprime; \f$res\f$ is left unchanged.
 */
bool bn_mod_inverse(bignum* res, const bignum* a, const bignum* m);

// bigexp.c
/**
 * @brief Calculate \f$a^b\f$ into \f$r\f$ using binary exponentiation.
 *
 * Let \f$n =\f$ a->size and \f$e =\f$ b->size, measured in 64-bit limbs.
 *
 * This computes:
 *
 *   \f$r = a^b\f$
 *
 * using left-to-right square-and-multiply over the bits of \f$b\f$. The
 * result has roughly \f$e \cdot n\f$ limbs, so this is only practical for small
 * exponents; use bn_mod_exp() for large ones.
 *
 * If \f$b\f$ is zero, \f$r\f$ is set to 1 (including \f$0^0\f$).
 *
 * Complexity:
 *   - Time: \f$O(e \cdot n^2)\f$ - one squaring per exponent bit (\f$64 \cdot
 * e\f$ of them) plus one multiplication per set bit, each an \f$O(n^2)\f$
 * bignum multiply of growing operands
 *   - Auxiliary memory: \f$O(n)\f$ limbs for temporaries
 *   - Output memory: \f$O(e \cdot n)\f$ limbs
 *
 * @param[out] r Result of \f$a^b\f$.
 * @param[in]  a Base.
 * @param[in]  b Exponent.
 */
void bn_pow(bignum* r, const bignum* a, const bignum* b);

/**
 * @brief Calculate \f$a_{bar}^d\f$ in the Montgomery domain.
 *
 * Let \f$n =\f$ ctx->n.size and \f$e =\f$ d->size, measured in 64-bit limbs.
 *
 * This computes:
 *
 *   \f$r_{bar} = a_{bar}^d \bmod n\f$
 *
 * where \f$a_{bar}\f$ is already in Montgomery form (\f$a_{bar} = a \cdot R
 * \bmod n\f$) and \f$r_{bar}\f$ is returned in Montgomery form.
 *
 * For short exponents (fewer than 64 bits) it uses plain left-to-right
 * square-and-multiply with bn_mont_mul(). For longer exponents it uses
 * a fixed window of width \f$w = 4\f$: the odd powers \f$a^1, a^3, \ldots,
 * a^{15}\f$ are precomputed in Montgomery form, then the exponent is scanned
 * from the most significant bit; each run is consumed as a window of up to
 * \f$w\f$ bits ending in a 1-bit, costing \f$w\f$ squarings plus one
 * multiplication by the
 * precomputed window value. This uses roughly 15-20% fewer Montgomery
 * multiplications than plain binary.
 *
 * Complexity:
 *   - Time: \f$O(e \cdot n^2)\f$ - one Montgomery squaring per exponent bit
 * plus one Montgomery multiplication per window (plus the \f$O(1)\f$ precompute
 * for the window table)
 *   - Auxiliary memory: \f$O(n)\f$ limbs for the window table
 *   - Output memory: \f$O(n)\f$ limbs
 *
 * @param[out] r_bar Result in Montgomery form: \f$a_{bar}^d \bmod n\f$.
 * @param[in]  a_bar Base in Montgomery form (\f$a \cdot R \bmod n\f$).
 * @param[in]  d     Exponent.
 * @param[in]  ctx   Initialized Montgomery context for the modulus \f$n\f$.
 */
void bn_mont_exp(bignum* r_bar, const bignum* a_bar, const bignum* d,
                 const bn_mont_ctx* ctx);

/**
 * @brief Calculate \f$a^b \bmod m\f$ into \f$r\f$ using plain (non-Montgomery)
 * arithmetic.
 *
 * Let \f$n =\f$ m->size and \f$e =\f$ b->size, measured in 64-bit limbs.
 *
 * This computes:
 *
 *   \f$r = a^b \bmod m\f$
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
 *   - Time: \f$O(e \cdot n^2)\f$ multiplications plus \f$O(e \cdot n^2)\f$
 * divisions, i.e.
 *           \f$O(e \cdot n^2)\f$ with a large constant
 *   - Auxiliary memory: \f$O(n)\f$ limbs for temporaries
 *   - Output memory: \f$O(n)\f$ limbs
 *
 * @param[out] r Result of \f$a^b \bmod m\f$.
 * @param[in]  a Base.
 * @param[in]  b Exponent.
 * @param[in]  m Modulus.
 */
void bn_mod_exp_slow(bignum* r, const bignum* a, const bignum* b,
                     const bignum* m);

/**
 * @brief Calculate \f$a^b \bmod m\f$ into \f$r\f$.
 *
 * Let \f$n =\f$ m->size and \f$e =\f$ b->size, measured in 64-bit limbs.
 *
 * This computes:
 *
 *   \f$r = a^b \bmod m\f$
 *
 * For odd \f$m\f$ it uses the fast Montgomery path (bn_mod_exp_mont) with a
 * freshly initialized context. For even \f$m\f$ it falls back to
 * bn_mod_exp_slow(), since Montgomery reduction requires an odd
 * modulus.
 *
 * Complexity:
 *   - Time: \f$O(e \cdot n^2)\f$ - \f$O(e)\f$ Montgomery multiplications (or
 * plain multiply+divide pairs for even \f$m\f$)
 *   - Auxiliary memory: \f$O(n)\f$ limbs for the context and temporaries
 *   - Output memory: \f$O(n)\f$ limbs
 *
 * @param[out] r Result of \f$a^b \bmod m\f$.
 * @param[in]  a Base.
 * @param[in]  b Exponent.
 * @param[in]  m Modulus.
 */
void bn_mod_exp(bignum* r, const bignum* a, const bignum* b, const bignum* m);

/**
 * @brief Calculate \f$a^b \bmod m\f$ into \f$r\f$ using a caller-provided
 * Montgomery context.
 *
 * Let \f$n =\f$ m->size and \f$e =\f$ b->size, measured in 64-bit limbs.
 *
 * This computes:
 *
 *   \f$r = a^b \bmod m\f$
 *
 * by converting a into the Montgomery domain (bn_mont_in), exponentiating
 * with bn_mont_exp(), and converting the result back (bn_mont_out).
 *
 * Precondition: \f$m\f$ is nonzero and odd, and ctx was initialized with
 * bn_mont_ctx_init() for this \f$m\f$.
 *
 * Complexity:
 *   - Time: \f$O(e \cdot n^2)\f$ - \f$O(e)\f$ Montgomery multiplications plus
 * two conversions, each one Montgomery multiplication
 *   - Auxiliary memory: \f$O(n)\f$ limbs for temporaries
 *   - Output memory: \f$O(n)\f$ limbs
 *
 * @param[out] r     Result of \f$a^b \bmod m\f$.
 * @param[in]  a     Base.
 * @param[in]  b     Exponent.
 * @param[in]  m     Modulus (nonzero and odd).
 * @param[in]  ctx   Montgomery context initialized for \f$m\f$.
 */
void bn_mod_exp_mont(bignum* r, const bignum* a, const bignum* b,
                     const bignum* m, const bn_mont_ctx* ctx);

// bigshift.c
/**
 * @brief Shift \f$r\f$ left by one bit, in place: \f$r = r \cdot 2\f$.
 *
 * Let \f$n =\f$ r->size, measured in 64-bit limbs.
 *
 * The bits are propagated from the least significant limb to the most
 * significant one; if the top bit overflows, a new limb is appended.
 * The sign is preserved.
 *
 * Complexity:
 *   - Time: \f$O(n)\f$
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(n)\f$ limbs, or \f$O(n + 1)\f$ if a carry limb is
 * created
 *
 * @param[in,out] r Value to shift left by one bit (modified in place).
 */
void bn_lshift1(bignum* r);

/**
 * @brief Shift \f$r\f$ right by one bit, in place: \f$r = r / 2\f$ (truncated).
 *
 * Let \f$n =\f$ r->size, measured in 64-bit limbs.
 *
 * The bits are propagated from the most significant limb to the least
 * significant one; leading zero limbs are trimmed afterwards. The sign
 * is preserved.
 *
 * Complexity:
 *   - Time: \f$O(n)\f$
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(n)\f$ limbs
 *
 * @param[in,out] r Value to shift right by one bit (modified in place).
 */
void bn_rshift1(bignum* r);

/**
 * @brief Shift \f$a\f$ left by \f$shift\f$ bits into \f$r\f$: \f$r = a \ll
 * shift\f$.
 *
 * Let \f$n =\f$ a->size, measured in 64-bit limbs.
 * Let \f$w =\f$ shift / 64 (whole limbs) and \f$b =\f$ shift % 64 (remaining
 * bits).
 *
 * Each limb of \f$a\f$ is shifted left by \f$b\f$ bits and placed \f$w\f$ limbs
 * higher in \f$r\f$, with the overflow bits carried into the next limb. The
 * sign is preserved.
 *
 * r may alias a.
 *
 * Complexity:
 *   - Time: \f$O(n)\f$
 *   - Auxiliary memory: \f$O(1)\f$ in the normal case,
 *                       \f$O(n)\f$ if \f$r\f$ aliases \f$a\f$ and a temporary
 * is used
 *   - Output memory: \f$O(n)\f$ limbs, at most \f$n + w + 1\f$ limbs
 *
 * @param[out] r     Result of \f$a \ll shift\f$.
 * @param[in]  a     Value to shift.
 * @param[in]  shift Number of bits to shift left.
 */
void bn_lshift(bignum* r, const bignum* a, int shift);

/**
 * @brief Shift \f$a\f$ right by \f$shift\f$ bits into \f$r\f$: \f$r = a \gg
 * shift\f$ (truncated).
 *
 * Let \f$n =\f$ a->size, measured in 64-bit limbs.
 * Let \f$w =\f$ shift / 64 (whole limbs) and \f$b =\f$ shift % 64 (remaining
 * bits).
 *
 * The top \f$w\f$ limbs are dropped and the remaining limbs are shifted
 * right by \f$b\f$ bits, pulling in the low bits of the next limb. If the
 * shift is at least the size of \f$a\f$, the result is zero. The sign is
 * preserved.
 *
 * r may alias a.
 *
 * Complexity:
 *   - Time: \f$O(n)\f$
 *   - Auxiliary memory: \f$O(1)\f$ in the normal case,
 *                       \f$O(n)\f$ if \f$r\f$ aliases \f$a\f$ and a temporary
 * is used
 *   - Output memory: \f$O(n)\f$ limbs, at most \f$n - w\f$ limbs
 *
 * @param[out] r     Result of \f$a \gg shift\f$.
 * @param[in]  a     Value to shift.
 * @param[in]  shift Number of bits to shift right.
 */
void bn_rshift(bignum* r, const bignum* a, int shift);

/**
 * @brief Shift \f$r\f$ left by one bit and add \f$bit\f$ (0 or 1), in place:
 * \f$r = (r \ll 1) + bit\f$.
 *
 * Let \f$n =\f$ r->size, measured in 64-bit limbs.
 *
 * This is the primitive used when building numbers bit by bit (e.g.
 * during parsing or in exponentiation loops). The added bit enters at
 * the least significant position and the carry propagates upward; if
 * the top bit overflows, a new limb is appended. The sign is
 * preserved.
 *
 * Complexity:
 *   - Time: \f$O(n)\f$
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(n)\f$ limbs, or \f$O(n + 1)\f$ if a carry limb is
 * created
 *
 * @param[in,out] r   Value to shift and add into (modified in place).
 * @param[in]     bit Bit (0 or 1) to add at the least significant position.
 */
void bn_lshift1_add(bignum* r, int bit);

// bigsqrt.c
/**
 * @brief Calculate the integer square root of a using Heron's method.
 *
 * Let \f$n =\f$ a->size, measured in 64-bit limbs.
 *
 * This computes:
 *
 *   \f$r = \lfloor\sqrt{a}\rfloor\f$
 *
 * The iteration starts at \f$x_0 = 2^{\lceil k/2 \rceil + 1}\f$ where \f$k\f$
 * is the bit length of \f$a\f$, and repeatedly applies:
 *
 *   \f$x_{i+1} = (x_i + a / x_i) / 2\f$
 *
 * until the sequence stops decreasing. Each iteration roughly doubles
 * the number of correct bits, so \f$O(\log n)\f$ iterations suffice.
 *
 * If \f$a\f$ is negative, \f$r\f$ is left unchanged (no real square root
 * exists).
 *
 * Complexity:
 *   - Time: \f$O(n^2 \log n)\f$ - \f$O(\log n)\f$ iterations, each dominated by
 * a division of an n-limb value by an n-limb value
 *   - Auxiliary memory: \f$O(n)\f$ limbs for temporaries
 *   - Output memory: \f$O(n)\f$ limbs (at most \f$n/2 + 1\f$ limbs)
 *
 * @param[out] r Result storing \f$\lfloor\sqrt{a}\rfloor\f$.
 * @param[in]  a Value to take the square root of.
 */
void bn_isqrt_heron(bignum* r, const bignum* a);

/**
 * @brief Calculate the integer square root of a.
 *
 * Let \f$n =\f$ a->size, measured in 64-bit limbs.
 *
 * This computes:
 *
 *   \f$r = \lfloor\sqrt{a}\rfloor\f$
 *
 * Currently a thin wrapper around bn_isqrt_heron().
 *
 * Complexity:
 *   - Time: \f$O(n^2 \log n)\f$, see bn_isqrt_heron()
 *   - Auxiliary memory: \f$O(n)\f$ limbs
 *   - Output memory: \f$O(n)\f$ limbs
 *
 * @param[out] r Result storing \f$\lfloor\sqrt{a}\rfloor\f$.
 * @param[in]  a Value to take the square root of.
 */
void bn_isqrt(bignum* r, const bignum* a);

// biglog.c
/**
 * @brief Calculate the integer natural logarithm of a, truncated.
 *
 * Let \f$n =\f$ a->size, measured in 64-bit limbs.
 *
 * This computes:
 *
 *   \f$r = \lfloor\ln(a)\rfloor\f$
 *
 * using the identity \f$\ln(a) = \ln(2) \cdot \log_2(a)\f$. The base-2
 * logarithm is exact (see bn_log_2), and \f$\ln(2)\f$ is approximated by the
 * fixed-point rational 69314718 / 100000000 (8 decimal digits), so the result
 * is accurate to within roughly 1 for large \f$a\f$.
 *
 * If \f$a\f$ is zero, \f$r\f$ is left unchanged (\f$\ln(0)\f$ is undefined).
 *
 * Complexity:
 *   - Time: \f$O(n)\f$ for the bit length, plus \f$O(n^2)\f$ for the bn_mul()
 * and
 *           \f$O(n^2)\f$ for the bn_div() of two n-limb values
 *   - Auxiliary memory: \f$O(n)\f$ limbs for temporaries
 *   - Output memory: \f$O(1)\f$ limbs (the result fits in a few limbs)
 *
 * @param[out] r Result storing \f$\lfloor\ln(a)\rfloor\f$.
 * @param[in]  a Value to take the natural logarithm of.
 */
void bn_ln(bignum* r, bignum* a);

/**
 * @brief Calculate the integer base-2 logarithm of a.
 *
 * Let \f$n =\f$ a->size, measured in 64-bit limbs.
 *
 * This computes:
 *
 *   \f$r = \lfloor\log_2(a)\rfloor\f$
 *
 * which is exact and equals the bit length of \f$a\f$ minus one.
 *
 * If \f$a\f$ is zero, \f$r\f$ is left unchanged (\f$\log_2(0)\f$ is undefined).
 *
 * Complexity:
 *   - Time: \f$O(1)\f$ (only the most significant limb is inspected)
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(1)\f$ limbs
 *
 * @param[out] r Result storing \f$\lfloor\log_2(a)\rfloor\f$.
 * @param[in]  a Value to take the base-2 logarithm of.
 */
void bn_log_2(bignum* r, bignum* a);

// bigmont.c
/**
 * @brief Initialize a Montgomery context for the modulus \f$n\f$.
 *
 * Let \f$n_l =\f$ n->size, measured in 64-bit limbs.
 *
 * Computes the context fields:
 *
 *   ctx->n        = \f$n\f$
 *   ctx->n_inv    = \f$-n^{-1} \bmod 2^{64}\f$   (from the low limb of \f$n\f$)
 *   ctx->one_mont = \f$R \bmod n\f$            (\f$R = 2^{64 \cdot n_l}\f$)
 *   ctx->r_square = \f$R^2 \bmod n\f$
 *
 * and preallocates ctx->tmp with \f$2 \cdot n_l + 1\f$ limbs of scratch space.
 *
 * Precondition: \f$n\f$ is odd and greater than \f$1\f$.
 *
 * Complexity:
 *   - Time: \f$O(n_l^2)\f$ - one_mont = \f$R \bmod n\f$ costs one division,
 * r_square = \f$R^2 \bmod n\f$ costs one multiply plus one division
 *   - Auxiliary memory: \f$O(n_l)\f$ limbs
 *   - Output memory: \f$O(n_l)\f$ limbs per context field
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
 *   - Time: \f$O(1)\f$
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(1)\f$
 *
 * @param[in,out] ctx Context to free.
 */
void bn_mont_ctx_free(bn_mont_ctx* ctx);

/**
 * @brief Montgomery reduction (REDC).
 *
 * Let \f$n_l =\f$ ctx->n.size, measured in 64-bit limbs.
 *
 * Given \f$t\f$ with \f$0 \le t < n \cdot 2^{64 \cdot n_l}\f$, this computes:
 *
 *   \f$r = t \cdot R^{-1} \bmod n\f$
 *
 * in \f$O(n_l^2)\f$ time. For each limb \f$i\f$ it forms the multiple
 * \f$m = t[i] \cdot (-n^{-1} \bmod 2^{64})\f$ and adds \f$m \cdot n\f$ shifted
 * by \f$i\f$ limbs to \f$t\f$, which zeroes out limb \f$i\f$ (\f$\bmod
 * 2^{64}\f$). After \f$n_l\f$ steps the lower half of \f$t\f$ is zero, so the
 * result is the upper half, followed by a final conditional subtraction of
 * \f$n\f$ to bring \f$r\f$ into \f$[0, n)\f$.
 *
 * \f$t\f$ is destroyed (overwritten) and must have at least \f$2 \cdot n_l +
 * 1\f$ limbs of capacity.
 *
 * Complexity:
 *   - Time: \f$O(n_l^2)\f$ for the reduction loop, plus \f$O(n_l)\f$ for the
 * final correction
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(n_l)\f$ limbs
 *
 * @param[out]    r   Result \f$t \cdot R^{-1} \bmod n\f$, in \f$[0, n)\f$.
 * @param[in,out] t   Value to reduce (destroyed; needs \f$2 \cdot n_l + 1\f$
 * limbs).
 * @param[in]     ctx Initialized Montgomery context.
 */
void bn_mont_redc(bignum* r, bignum* t, const bn_mont_ctx* ctx);

/**
 * @brief Convert a value from the normal domain into the Montgomery domain.
 *
 * Let \f$n_l =\f$ ctx->n.size, measured in 64-bit limbs.
 *
 * This computes:
 *
 *   \f$A_{bar} = A \cdot R \bmod n\f$
 *
 * as one Montgomery multiplication of \f$A\f$ by \f$R^2 \bmod n\f$:
 * REDC(\f$A \cdot R^2\f$) = \f$A \cdot R \bmod n\f$.
 *
 * Complexity:
 *   - Time: \f$O(n_l^2)\f$
 *   - Auxiliary memory: \f$O(n_l)\f$ limbs (ctx->tmp)
 *   - Output memory: \f$O(n_l)\f$ limbs
 *
 * @param[out] A_bar Result in Montgomery form: \f$A \cdot R \bmod n\f$.
 * @param[in]  A     Value in the normal domain.
 * @param[in]  ctx   Initialized Montgomery context.
 */
void bn_mont_in(bignum* A_bar, const bignum* A, const bn_mont_ctx* ctx);

/**
 * @brief Convert a value from the Montgomery domain back to the normal domain.
 *
 * Let \f$n_l =\f$ ctx->n.size, measured in 64-bit limbs.
 *
 * This computes:
 *
 *   \f$A = A_{bar} \cdot R^{-1} \bmod n\f$
 *
 * as one Montgomery multiplication of \f$A_{bar}\f$ by 1.
 *
 * Complexity:
 *   - Time: \f$O(n_l^2)\f$
 *   - Auxiliary memory: \f$O(n_l)\f$ limbs (ctx->tmp and a temporary for 1)
 *   - Output memory: \f$O(n_l)\f$ limbs
 *
 * @param[out] A     Result in the normal domain: \f$A_{bar} \cdot R^{-1} \bmod
 * n\f$.
 * @param[in]  A_bar Value in the Montgomery domain.
 * @param[in]  ctx   Initialized Montgomery context.
 */
void bn_mont_out(bignum* A, const bignum* A_bar, const bn_mont_ctx* ctx);

/**
 * @brief Montgomery multiplication of two values in the Montgomery domain.
 *
 * Let \f$n_l =\f$ ctx->n.size, measured in 64-bit limbs.
 *
 * This computes:
 *
 *   \f$r = a_{bar} \cdot b_{bar} \cdot R^{-1} \bmod n\f$
 *
 * so that if \f$a_{bar} = A \cdot R \bmod n\f$ and \f$b_{bar} = B \cdot R \bmod
 * n\f$, then
 * \f$r = A \cdot B \cdot R \bmod n\f$ (the Montgomery form of \f$A \cdot B\f$).
 *
 * Thin wrapper around bn_mont_mul_raw().
 *
 * Complexity:
 *   - Time: \f$O(n_l^2)\f$
 *   - Auxiliary memory: \f$O(n_l)\f$ limbs (ctx->tmp)
 *   - Output memory: \f$O(n_l)\f$ limbs
 *
 * @param[out] r     Result \f$a_{bar} \cdot b_{bar} \cdot R^{-1} \bmod n\f$.
 * @param[in]  a_bar First operand in Montgomery form.
 * @param[in]  b_bar Second operand in Montgomery form.
 * @param[in]  ctx   Initialized Montgomery context.
 */
void bn_mont_mul(bignum* r, const bignum* a_bar, const bignum* b_bar,
                 const bn_mont_ctx* ctx);

/**
 * @brief Montgomery multiplication using the context's scratch buffer.
 *
 * Let \f$n_l =\f$ ctx->n.size, measured in 64-bit limbs.
 *
 * This computes:
 *
 *   \f$result = A_{bar} \cdot B_{bar} \cdot R^{-1} \bmod n\f$
 *
 * by first forming the full product \f$T = A_{bar} \cdot B_{bar}\f$ (at most
 * \f$2 \cdot n_l\f$ limbs) in ctx->tmp, zero-padding it to \f$2 \cdot n_l +
 * 1\f$ limbs, and applying bn_mont_redc().
 *
 * Note: uses ctx->tmp as scratch, so it is not reentrant and must not
 * be called with \f$A_{bar}\f$ or \f$B_{bar}\f$ aliasing ctx->tmp.
 *
 * Complexity:
 *   - Time: \f$O(n_l^2)\f$ - one bignum multiply plus one REDC
 *   - Auxiliary memory: \f$O(n_l)\f$ limbs (ctx->tmp)
 *   - Output memory: \f$O(n_l)\f$ limbs
 *
 * @param[out] result  Result \f$A_{bar} \cdot B_{bar} \cdot R^{-1} \bmod n\f$.
 * @param[in]  A_bar   First operand in Montgomery form.
 * @param[in]  B_bar   Second operand in Montgomery form.
 * @param[in]  ctx     Initialized Montgomery context (provides scratch).
 */
void bn_mont_mul_raw(bignum* result, const bignum* A_bar, const bignum* B_bar,
                     const bn_mont_ctx* ctx);

#endif
