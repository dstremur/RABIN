#ifndef BIGPSEUDO_H
#define BIGPSEUDO_H

/*===========================================================================
 *  bigpseudo.h
 *
 *  Generation of strong pseudoprimes to base 2.
 *
 *  Layout:
 *    - strong pseudoprime generator (bn_gen_strps)
 *===========================================================================*/

#include "bigcore.h"

/**
 * @brief Generate a k-bit strong pseudoprime to base 2 into p.
 *
 * Repeatedly draws random k-bit candidates from /dev/urandom until one
 * is found that is composite (fails BPSW) yet passes the base-2
 * Miller-Rabin test (bn_rabin_mont). Such numbers are strong
 * pseudoprimes to base 2.
 *
 * Returns true on success, false if /dev/urandom cannot be opened.
 *
 * Complexity:
 *   Time: unbounded expected - each attempt costs one BPSW test plus
 *         one base-2 Miller-Rabin test; the hit rate depends on k
 *   Auxiliary memory: O(k/64) limbs
 *   Output memory: O(k/64) limbs
 *
 * @param[out] p Result storing the strong pseudoprime.
 * @param[in]  k Desired bit length of the pseudoprime.
 *
 * @return true  On success.
 * @return false If /dev/urandom cannot be opened.
 */
bool bn_gen_strps(bignum* p, u64 k);

#endif
