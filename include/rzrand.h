#ifndef RZRAND_H
#define RZRAND_H

#include <stdio.h>

#include "rabin_errors.h"
#include "rz.h"

/**
 * @brief Generate a random bits-long number.
 *
 * Opens /dev/urandom, delegates to rz_gen_random_with_fd(), and
 * closes the descriptor.
 *
 * Complexity:
 *   - Time: \f$O(bits / 64)\f$
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(bits / 64)\f$ limbs
 *
 * @param[out] r    Result storing the random number.
 * @param[in]  bits Desired bit length of the number.
 *
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR,
 * RABIN_ERR_INVALID_ARG (if /dev/urandom cannot be opened or the read
 * fails), or RABIN_ERR_OUT_OF_MEMORY.
 * @see rz_gen_random_with_fd()
 */
rabin_err_t rz_gen_random(rz_t* r, u64 bits);

/**
 * @brief Generate a random bits-long odd number using a given urandom fd.
 *
 * Reads \f$\lceil bits / 64 \rceil\f$ random limbs from \f$fd\f$, masks the top
 * limb to the exact bit length, sets the MSB so the number is exactly
 * \f$bits\f$ long, and sets the LSB so it is odd.
 *
 * Complexity:
 *   - Time: \f$O(bits / 64)\f$
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(bits / 64)\f$ limbs
 *
 * @param[out] r    Result storing the random odd number.
 * @param[in]  bits Desired bit length of the number (>= 1).
 * @param[in]  fd   Open file descriptor for /dev/urandom.
 *
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR,
 * RABIN_ERR_INVALID_ARG (if \f$bits < 1\f$ or the read fails), or
 * RABIN_ERR_OUT_OF_MEMORY.
 * @see rz_gen_random_with_fd()
 */
rabin_err_t rz_gen_random_odd_with_fd(rz_t* r, u64 bits, int fd);

/**
 * @brief Generate a random bits-long number using a given urandom fd.
 *
 * Reads \f$\lceil bits / 64 \rceil\f$ random limbs from \f$fd\f$ and masks the
 * top limb to the exact bit length; \f$bits = 0\f$ stores zero.
 *
 * Complexity:
 *   - Time: \f$O(bits / 64)\f$
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(bits / 64)\f$ limbs
 *
 * @param[out] r    Result storing the random number.
 * @param[in]  bits Desired bit length of the number (0 yields zero).
 * @param[in]  fd   Open file descriptor for /dev/urandom.
 *
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR,
 * RABIN_ERR_INVALID_ARG (if the read fails), or RABIN_ERR_OUT_OF_MEMORY.
 * @see rz_gen_random(), rz_gen_random_odd_with_fd()
 */
rabin_err_t rz_gen_random_with_fd(rz_t* r, u64 bits, int fd);

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
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR,
 * RABIN_ERR_INVALID_ARG, or RABIN_ERR_OUT_OF_MEMORY.
 */
rabin_err_t rz_gen_random_range(rz_t* r, const rz_t* low, const rz_t* high);

#endif
