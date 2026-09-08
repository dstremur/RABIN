#ifndef BIGLUCAS_H
#define BIGLUCAS_H

#include "bigcore.h"

/**
 * @brief Compute \f$(U_n, V_n, Q^n)\f$ exactly, by recursion on the binary
 * expansion of \f$n\f$.
 *
 * Let \f$n_l =\f$ n->size, measured in 64-bit limbs.
 *
 * Base cases: \f$n = 0\f$ gives \f$(U, V, Q^n) = (0, 2, 1)\f$; \f$n = 1\f$
 * gives
 * \f$(1, P, Q)\f$. The recursion halves \f$n\f$ and applies the doubling
 * identities above, so the sequence of values is built from the
 * halved result.
 *
 * Complexity:
 *   - Time: \f$O(n_l^2 \log n)\f$ - \f$O(\log n)\f$ levels of recursion, each
 * with a few multiplications of up to \f$n_l\f$-limb values
 *   - Auxiliary memory: \f$O(n_l)\f$ limbs per recursion level, \f$O(n_l \log
 * n)\f$ total on the stack of temporaries
 *   - Output memory: \f$O(n_l)\f$ limbs per result
 *
 * @param[out] u  Result storing \f$U_n\f$.
 * @param[out] v  Result storing \f$V_n\f$.
 * @param[in]  p  Lucas parameter P.
 * @param[in]  q  Lucas parameter Q.
 * @param[out] qn Result storing \f$Q^n\f$.
 * @param[in]  n  Index of the Lucas sequence term.
 */
void bn_lucas_solve(bignum* u, bignum* v, const bignum* p, const bignum* q,
                    bignum* qn, const bignum* n);

/**
 * @brief Compute \f$(U_n \bmod m, V_n \bmod m, Q^n \bmod m)\f$ iteratively over
 * the bits of \f$n\f$, in the Montgomery domain.
 *
 * Let \f$n_l =\f$ m->size, measured in 64-bit limbs.
 *
 * Same doubling identities as bn_lucas_solve(), but applied
 * left-to-right over the bits of \f$n\f$ (starting from \f$(U_1, V_1, Q^1)\f$ =
 * \f$(1, P, Q)\f$) with all multiplications as Montgomery multiplications
 * modulo \f$m\f$. Subtractions that would go negative are brought back into
 * \f$[0, m)\f$ by adding \f$m\f$. Division by 2 is done by adding \f$m\f$ when
 * the value is odd, then shifting.
 *
 * The results are converted out of the Montgomery domain before
 * return.
 *
 * Complexity:
 *   - Time: \f$O(n_l^3)\f$ for the context initialization, then
 *           \f$O(\log n \cdot n_l^2)\f$ for the bit loop
 *   - Auxiliary memory: \f$O(n_l)\f$ limbs for the context and temporaries
 *   - Output memory: \f$O(n_l)\f$ limbs per result
 *
 * @param[out] u  Result storing \f$U_n \bmod m\f$.
 * @param[out] v  Result storing \f$V_n \bmod m\f$.
 * @param[in]  p  Lucas parameter P.
 * @param[in]  q  Lucas parameter Q.
 * @param[out] qn Result storing \f$Q^n \bmod m\f$.
 * @param[in]  n  Index of the Lucas sequence term.
 * @param[in]  m  Modulus.
 */
void bn_lucas_solve_mod(bignum* u, bignum* v, const bignum* p, const bignum* q,
                        bignum* qn, const bignum* n, const bignum* m);

/**
 * @brief Compute the \f$n\f$-th terms of the Lucas sequences \f$U_n(P, Q)\f$
 * and
 * \f$V_n(P, Q)\f$ exactly.
 *
 * Let \f$n_l =\f$ n->size, measured in 64-bit limbs.
 *
 * Thin wrapper around bn_lucas_solve() that discards \f$Q^n\f$.
 *
 * Complexity:
 *   - Time: \f$O(n_l^2 \log n)\f$, see bn_lucas_solve()
 *   - Auxiliary memory: \f$O(n_l \log n)\f$ limbs
 *   - Output memory: \f$O(n_l)\f$ limbs per result
 *
 * @param[out] u Result storing \f$U_n\f$.
 * @param[out] v Result storing \f$V_n\f$.
 * @param[in]  p Lucas parameter P.
 * @param[in]  q Lucas parameter Q.
 * @param[in]  n Index of the Lucas sequence term.
 */
void bn_lucas(bignum* u, bignum* v, const bignum* p, const bignum* q,
              const bignum* n);

/**
 * @brief Compute \f$U_n(P, Q) \bmod m\f$ and \f$V_n(P, Q) \bmod m\f$.
 *
 * Let \f$n_l =\f$ m->size, measured in 64-bit limbs.
 *
 * Thin wrapper around bn_lucas_solve_mod() that discards \f$Q^n \bmod m\f$.
 *
 * Complexity:
 *   - Time: \f$O(n_l^3)\f$ for the context initialization, then
 *           \f$O(\log n \cdot n_l^2)\f$, see bn_lucas_solve_mod()
 *   - Auxiliary memory: \f$O(n_l)\f$ limbs
 *   - Output memory: \f$O(n_l)\f$ limbs per result
 *
 * @param[out] u Result storing \f$U_n \bmod m\f$.
 * @param[out] v Result storing \f$V_n \bmod m\f$.
 * @param[in]  p Lucas parameter P.
 * @param[in]  q Lucas parameter Q.
 * @param[in]  n Index of the Lucas sequence term.
 * @param[in]  m Modulus.
 */
void bn_lucas_mod(bignum* u, bignum* v, const bignum* p, const bignum* q,
                  const bignum* n, const bignum* m);

#endif
