#ifndef BIGRAND_H
#define BIGRAND_H

#include <stdio.h>

#include "bigcore.h"

/**
 * @brief Generate a random bits-long number.
 *
 * Opens /dev/urandom, delegates to bn_gen_random_with_fd(), and
 * closes the descriptor.
 *
 * Returns true on success, false if /dev/urandom cannot be opened or
 * the read fails.
 *
 * Complexity:
 *   - Time: \f$O(bits / 64)\f$
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(bits / 64)\f$ limbs
 *
 * @param[out] r    Result storing the random number.
 * @param[in]  bits Desired bit length of the number.
 *
 * @return true  On success.
 * @return false If /dev/urandom cannot be opened or the read fails.
 */
bool bn_gen_random(bignum* r, u64 bits);

/**
 * @brief Generate a random bits-long odd number using a given urandom fd.
 *
 * Reads \f$\lceil bits / 64 \rceil\f$ random limbs from \f$fd\f$, masks the top
 * limb to the exact bit length, sets the MSB so the number is exactly
 * \f$bits\f$ long, and sets the LSB so it is odd.
 *
 * Returns true on success, false on allocation or read failure.
 *
 * Complexity:
 *   - Time: \f$O(bits / 64)\f$
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(bits / 64)\f$ limbs
 *
 * @param[out] r    Result storing the random odd number.
 * @param[in]  bits Desired bit length of the number.
 * @param[in]  fd   Open file descriptor for /dev/urandom.
 *
 * @return true  On success.
 * @return false On allocation or read failure.
 */
bool bn_gen_random_odd_with_fd(bignum* r, u64 bits, int fd);

/**
 * @brief Generate a random bits-long number using a given urandom fd.
 *
 * Returns true on success, false on allocation or read failure.
 *
 * Complexity:
 *   - Time: \f$O(bits / 64)\f$
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(bits / 64)\f$ limbs
 *
 * @param[out] r    Result storing the random number.
 * @param[in]  bits Desired bit length of the number.
 * @param[in]  fd   Open file descriptor for /dev/urandom.
 *
 * @return true  On success.
 * @return false On allocation or read failure.
 */
bool bn_gen_random_with_fd(bignum* r, u64 bits, int fd);

/**
 * @brief Generate a random number in the closed range \f$[low, high]\f$.
 *
 * Let \f$b =\f$ bit length of \f$(high - low)\f$.
 *
 * Rejection sampling: draws random \f$b\f$-bit numbers until one is \f$\le
 * (high - low)\f$, then adds \f$low\f$. If \f$low = high\f$, copies \f$low\f$
 * directly.
 *
 * Complexity:
 *   - Time: \f$O(b / 64)\f$ per draw, expected \f$O(b / 64)\f$ overall
 * (geometric number of draws), plus \f$O(b)\f$ for the range computation
 *   - Auxiliary memory: \f$O(b)\f$ limbs for the range
 *   - Output memory: \f$O(b)\f$ limbs
 *
 * @param[out] r    Result storing the random number in \f$[low, high]\f$.
 * @param[in]  low  Lower bound (inclusive).
 * @param[in]  high Upper bound (inclusive).
 */
void bn_gen_random_range(bignum* r, const bignum* low, const bignum* high);

#endif
