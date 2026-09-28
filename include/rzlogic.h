#ifndef RZLOGIC_H
#define RZLOGIC_H

#include <stdbool.h>

#include "rabin_errors.h"
#include "rz.h"

/**
 * @brief Report whether the current CPU supports AVX-512F.
 *
 * The bitwise routines in this module dispatch to AVX-512F intrinsics when
 * both the compiler was built with `__AVX512F__` and the runtime CPU
 * exposes the feature.
 *
 * @return True if AVX-512F is available at runtime, false otherwise.
 */
bool rz_supports_avx512(void);

/**
 * @brief Bitwise AND of two bignums: \f$r = a \& m\f$.
 *
 * Let \f$n =\f$ a->size, measured in 64-bit limbs.
 *
 * Copies \f$a\f$ into \f$r\f$ and ANDs each limb with the corresponding limb
 * of \f$m\f$ (limbs beyond m->size are zeroed), then trims. Used to extract
 * fixed-width chunks of a rz_t.
 *
 * Complexity:
 *   - Time: \f$O(n)\f$
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(n)\f$ limbs
 *
 * @param[out] r Result storing \f$a \& m\f$. May alias \p a or \p m.
 * @param[in] a First operand.
 * @param[in] m Second operand (mask).
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR.
 * @pre r, a, and m must be non-NULL.
 */
rabin_err_t rz_and(rz_t* r, const rz_t* a, const rz_t* m);

/**
 * @brief Bitwise OR of two bignums: \f$r = a \mid m\f$.
 *
 * Let \f$n =\f$ max(a->size, m->size), measured in 64-bit limbs.
 *
 * Computes the limb-wise OR of \f$a\f$ and \f$m\f$, zero-extending the
 * shorter operand, then trims the result.
 *
 * Complexity:
 *   - Time: \f$O(n)\f$
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(n)\f$ limbs
 *
 * @param[out] r Result storing \f$a \mid m\f$. May alias \p a or \p m.
 * @param[in] a First operand.
 * @param[in] m Second operand.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR.
 * @pre r, a, and m must be non-NULL.
 */
rabin_err_t rz_or(rz_t* r, const rz_t* a, const rz_t* m);

/**
 * @brief Bitwise XOR of two bignums: \f$r = a \oplus m\f$.
 *
 * Let \f$n =\f$ max(a->size, m->size), measured in 64-bit limbs.
 *
 * Computes the limb-wise XOR of \f$a\f$ and \f$m\f$, zero-extending the
 * shorter operand, then trims the result.
 *
 * Complexity:
 *   - Time: \f$O(n)\f$
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(n)\f$ limbs
 *
 * @param[out] r Result storing \f$a \oplus m\f$. May alias \p a or \p m.
 * @param[in] a First operand.
 * @param[in] m Second operand.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR.
 * @pre r, a, and m must be non-NULL.
 */
rabin_err_t rz_xor(rz_t* r, const rz_t* a, const rz_t* m);

/**
 * @brief Bitwise NOT (complement) of a bignum: \f$r = \sim a\f$.
 *
 * Let \f$n =\f$ a->size, measured in 64-bit limbs.
 *
 * Flips every bit of the \f$n\f$ limbs of \f$a\f$ into \f$r\f$, then trims
 * (the result is the least-significant-limb complement, i.e. \f$r \equiv
 * 2^{64n} - 1 - a\f$). A zero input yields a zero result.
 *
 * Complexity:
 *   - Time: \f$O(n)\f$
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(n)\f$ limbs
 *
 * @param[out] r Result storing \f$\sim a\f$. May alias \p a.
 * @param[in] a Operand.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR.
 * @pre r and a must be non-NULL.
 */
rabin_err_t rz_not(rz_t* r, const rz_t* a);

#endif
