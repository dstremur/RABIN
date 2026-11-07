#ifndef RZLUCAS_H
#define RZLUCAS_H

#include "rabin_errors.h"
#include "rz.h"

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
 * @param[in]  n  Index of the Lucas sequence term (non-negative).
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR /
 * RABIN_ERR_INVALID_ARG (negative \f$n\f$) / RABIN_ERR_OUT_OF_MEMORY.
 *
 * @par Algorithm Reference:
 * E. Lucas, "Theorie des fonctions numeriques periodiques, des
 * fonctions entieres et des nombres qui se rapportent aux nombres
 * premiers," Journal de Mathematiques Pures et Appliees, vol. 3, 1878.
 * @see rz_lucas_solve_mod()
 */
rabin_err_t rz_lucas_solve(rz_t* u, rz_t* v, const rz_t* p, const rz_t* q,
                           rz_t* qn, const rz_t* n);

/**
 * @brief Compute \f$(U_n \bmod m, V_n \bmod m, Q^n \bmod m)\f$ iteratively over
 * the bits of \f$n\f$, in the Montgomery domain.
 *
 * Let \f$n_l =\f$ m->size, measured in 64-bit limbs.
 *
 * Same doubling identities as rz_lucas_solve(), but applied
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
 * @param[in]  p  Lucas parameter P (may be negative).
 * @param[in]  q  Lucas parameter Q (may be negative).
 * @param[out] qn Result storing \f$Q^n \bmod m\f$.
 * @param[in]  n  Index of the Lucas sequence term (non-negative).
 * @param[in]  m  Modulus.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR /
 * RABIN_ERR_INVALID_ARG (negative \f$n\f$, or \f$m \le 1\f$ or \f$m\f$ even),
 * or RABIN_ERR_OUT_OF_MEMORY.
 *
 * @pre \f$n \ge 0\f$ and \f$m\f$ is odd with \f$m > 1\f$.
 */
rabin_err_t rz_lucas_solve_mod(rz_t* u, rz_t* v, const rz_t* p, const rz_t* q,
                               rz_t* qn, const rz_t* n, const rz_t* m);

/**
 * @brief Run the strong Lucas test on \f$n\f$ without ever leaving the
 * Montgomery domain.
 *
 * Let \f$n_l\f$ = n->size, measured in 64-bit limbs. Requires \f$n + 1
 * = d \cdot 2^s\f$ with \f$d \ge 1\f$ odd. One Montgomery context is
 * initialized for \f$n\f$; the double-and-add pass over the bits of
 * \f$d\f$ and the subsequent \f$2^s\f$ doubling chain run back to back
 * in the Montgomery domain, and the zero tests are evaluated directly
 * on the Montgomery residues (\f$0 \cdot R = 0\f$ and all residues are
 * canonical in \f$[0, n)\f$), so no rz_mont_out() conversion is needed.
 * In the doubling chain \f$U\f$ is not tracked, costing two Montgomery
 * multiplications per step.
 *
 * Complexity:
 *   - Time: \f$O(n_l^3)\f$ for the context initialization, then
 *           \f$O(\log n \cdot n_l^2)\f$ for both phases
 *   - Auxiliary memory: \f$O(n_l)\f$ limbs for the context and temporaries
 *
 * @param[out] res Set to true if n passes the strong Lucas test for
 *                 (P, Q), false otherwise; set to false on error.
 * @param[in]  n   Odd modulus \f$> 1\f$ under test.
 * @param[in]  P   Lucas parameter P (may be negative).
 * @param[in]  Q   Lucas parameter Q (may be negative).
 * @param[in]  d   Odd cofactor of \f$n + 1\f$ (\f$d \ge 1\f$).
 * @param[in]  s   Exponent with \f$n + 1 = d \cdot 2^s\f$.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR /
 * RABIN_ERR_INVALID_ARG / RABIN_ERR_OUT_OF_MEMORY.
 *
 * @pre \f$n\f$ is odd with \f$n > 1\f$ and \f$n + 1 = d \cdot 2^s\f$
 * with \f$d \ge 1\f$ odd.
 *
 * @see rz_lucas_solve_mod(), rz_stronglucas()
 */
rabin_err_t rz_strong_lucas_solve_mod(bool* res, const rz_t* n, const rz_t* P,
                                      const rz_t* Q, const rz_t* d, u64 s);

/**
 * @brief Compute the \f$n\f$-th terms of the Lucas sequences \f$U_n(P, Q)\f$
 * and
 * \f$V_n(P, Q)\f$ exactly.
 *
 * Let \f$n_l =\f$ n->size, measured in 64-bit limbs.
 *
 * Thin wrapper around rz_lucas_solve() that discards \f$Q^n\f$.
 *
 * Complexity:
 *   - Time: \f$O(n_l^2 \log n)\f$, see rz_lucas_solve()
 *   - Auxiliary memory: \f$O(n_l \log n)\f$ limbs
 *   - Output memory: \f$O(n_l)\f$ limbs per result
 *
 * @param[out] u Result storing \f$U_n\f$.
 * @param[out] v Result storing \f$V_n\f$.
 * @param[in]  p Lucas parameter P.
 * @param[in]  q Lucas parameter Q.
 * @param[in]  n Index of the Lucas sequence term.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR /
 * RABIN_ERR_OUT_OF_MEMORY.
 */
rabin_err_t rz_lucas(rz_t* u, rz_t* v, const rz_t* p, const rz_t* q,
                     const rz_t* n);

/**
 * @brief Compute \f$U_n(P, Q) \bmod m\f$ and \f$V_n(P, Q) \bmod m\f$.
 *
 * Let \f$n_l =\f$ m->size, measured in 64-bit limbs.
 *
 * Thin wrapper around rz_lucas_solve_mod() that discards \f$Q^n \bmod m\f$.
 *
 * Complexity:
 *   - Time: \f$O(n_l^3)\f$ for the context initialization, then
 *           \f$O(\log n \cdot n_l^2)\f$, see rz_lucas_solve_mod()
 *   - Auxiliary memory: \f$O(n_l)\f$ limbs
 *   - Output memory: \f$O(n_l)\f$ limbs per result
 *
 * @param[out] u Result storing \f$U_n \bmod m\f$.
 * @param[out] v Result storing \f$V_n \bmod m\f$.
 * @param[in]  p Lucas parameter P.
 * @param[in]  q Lucas parameter Q.
 * @param[in]  n Index of the Lucas sequence term.
 * @param[in]  m Modulus.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR /
 * RABIN_ERR_OUT_OF_MEMORY.
 */
rabin_err_t rz_lucas_mod(rz_t* u, rz_t* v, const rz_t* p, const rz_t* q,
                         const rz_t* n, const rz_t* m);

#endif
