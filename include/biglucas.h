#ifndef BIGLUCAS_H
#define BIGLUCAS_H

/*===========================================================================
 *  biglucas.h
 *
 *  Lucas sequences and Lucas primality tests.
 *
 *  Layout:
 *    - sequence computation   (bn_lucas, bn_lucas_mod)
 *    - exponentiation of the   (bn_lucas_solve, bn_lucas_solve_mod)
 *      sequence pair
 *===========================================================================*/

#include "bigcore.h"

/**
 * @brief Compute (U_n, V_n, Q^n) exactly, by recursion on the binary
 * expansion of n.
 *
 * Let n_l = n->size, measured in 64-bit limbs.
 *
 * Base cases: n = 0 gives (U, V, Q^n) = (0, 2, 1); n = 1 gives
 * (1, P, Q). The recursion halves n and applies the doubling
 * identities above, so the sequence of values is built from the
 * halved result.
 *
 * Complexity:
 *   Time: O(n_l^2 log n) - O(log n) levels of recursion, each with a
 *         few multiplications of up to n_l-limb values
 *   Auxiliary memory: O(n_l) limbs per recursion level, O(n_l log n)
 *                     total on the stack of temporaries
 *   Output memory: O(n_l) limbs per result
 *
 * @param[out] u  Result storing U_n.
 * @param[out] v  Result storing V_n.
 * @param[in]  p  Lucas parameter P.
 * @param[in]  q  Lucas parameter Q.
 * @param[out] qn Result storing Q^n.
 * @param[in]  n  Index of the Lucas sequence term.
 */
void bn_lucas_solve(bignum* u, bignum* v, const bignum* p, const bignum* q,
                    bignum* qn, const bignum* n);

/**
 * @brief Compute (U_n mod m, V_n mod m, Q^n mod m) iteratively over the bits
 * of n, in the Montgomery domain.
 *
 * Let n_l = m->size, measured in 64-bit limbs.
 *
 * Same doubling identities as bn_lucas_solve(), but applied
 * left-to-right over the bits of n (starting from (U_1, V_1, Q^1) =
 * (1, P, Q)) with all multiplications as Montgomery multiplications
 * modulo m. Subtractions that would go negative are brought back into
 * [0, m) by adding m. Division by 2 is done by adding m when the
 * value is odd, then shifting.
 *
 * The results are converted out of the Montgomery domain before
 * return.
 *
 * Complexity:
 *   Time: O(n_l^3) for the context initialization, then O(log n *
 *         n_l^2) for the bit loop
 *   Auxiliary memory: O(n_l) limbs for the context and temporaries
 *   Output memory: O(n_l) limbs per result
 *
 * @param[out] u  Result storing U_n mod m.
 * @param[out] v  Result storing V_n mod m.
 * @param[in]  p  Lucas parameter P.
 * @param[in]  q  Lucas parameter Q.
 * @param[out] qn Result storing Q^n mod m.
 * @param[in]  n  Index of the Lucas sequence term.
 * @param[in]  m  Modulus.
 */
void bn_lucas_solve_mod(bignum* u, bignum* v, const bignum* p, const bignum* q,
                        bignum* qn, const bignum* n, const bignum* m);

/**
 * @brief Compute the n-th terms of the Lucas sequences U_n(P, Q) and
 * V_n(P, Q) exactly.
 *
 * Let n_l = n->size, measured in 64-bit limbs.
 *
 * Thin wrapper around bn_lucas_solve() that discards Q^n.
 *
 * Complexity:
 *   Time: O(n_l^2 log n), see bn_lucas_solve()
 *   Auxiliary memory: O(n_l log n) limbs
 *   Output memory: O(n_l) limbs per result
 *
 * @param[out] u Result storing U_n.
 * @param[out] v Result storing V_n.
 * @param[in]  p Lucas parameter P.
 * @param[in]  q Lucas parameter Q.
 * @param[in]  n Index of the Lucas sequence term.
 */
void bn_lucas(bignum* u, bignum* v, const bignum* p, const bignum* q,
              const bignum* n);

/**
 * @brief Compute U_n(P, Q) mod m and V_n(P, Q) mod m.
 *
 * Let n_l = m->size, measured in 64-bit limbs.
 *
 * Thin wrapper around bn_lucas_solve_mod() that discards Q^n mod m.
 *
 * Complexity:
 *   Time: O(n_l^3) for the context initialization, then O(log n *
 *         n_l^2), see bn_lucas_solve_mod()
 *   Auxiliary memory: O(n_l) limbs
 *   Output memory: O(n_l) limbs per result
 *
 * @param[out] u Result storing U_n mod m.
 * @param[out] v Result storing V_n mod m.
 * @param[in]  p Lucas parameter P.
 * @param[in]  q Lucas parameter Q.
 * @param[in]  n Index of the Lucas sequence term.
 * @param[in]  m Modulus.
 */
void bn_lucas_mod(bignum* u, bignum* v, const bignum* p, const bignum* q,
                  const bignum* n, const bignum* m);

#endif
