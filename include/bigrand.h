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
 *   Time: O(bits / 64)
 *   Auxiliary memory: O(1)
 *   Output memory: O(bits / 64) limbs
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
 * Reads ceil(bits / 64) random limbs from fd, masks the top limb to
 * the exact bit length, sets the MSB so the number is exactly
 * bits long, and sets the LSB so it is odd.
 *
 * Returns true on success, false on allocation or read failure.
 *
 * Complexity:
 *   Time: O(bits / 64)
 *   Auxiliary memory: O(1)
 *   Output memory: O(bits / 64) limbs
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
 *   Time: O(bits / 64)
 *   Auxiliary memory: O(1)
 *   Output memory: O(bits / 64) limbs
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
 * @brief Generate a random number in the closed range [low, high].
 *
 * Let b = bit length of (high - low).
 *
 * Rejection sampling: draws random b-bit numbers until one is <=
 * (high - low), then adds low. If low == high, copies low directly.
 *
 * Complexity:
 *   Time: O(b / 64) per draw, expected O(b / 64) overall (geometric
 *         number of draws), plus O(b) for the range computation
 *   Auxiliary memory: O(b) limbs for the range
 *   Output memory: O(b) limbs
 *
 * @param[out] r    Result storing the random number in [low, high].
 * @param[in]  low  Lower bound (inclusive).
 * @param[in]  high Upper bound (inclusive).
 */
void bn_gen_random_range(bignum* r, const bignum* low, const bignum* high);

#endif
