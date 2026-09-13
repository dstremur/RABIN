#ifndef BIGVECTOR_H
#define BIGVECTOR_H

#include <assert.h>

#include "bigcore.h"

typedef struct {
  bignum* data;
  u64 size;
  u64 capacity;
  bool dynamic;
} bigvector;

/**
 * @brief Initialize a static vector of d zero bignums.
 *
 * Allocates exactly \f$d\f$ slots (\f$capacity = size = d\f$) and initializes
 * each to a zero bignum. The vector cannot be grown afterwards.
 *
 * Complexity:
 *   - Time: \f$O(d)\f$
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(d)\f$ bignums
 *
 * @param[out] a Vector to initialize.
 * @param[in]  d Number of elements.
 */
void bigvector_init(bigvector* a, u64 d);

/**
 * @brief Initialize an empty dynamic vector.
 *
 * The vector starts with \f$size = capacity = 0\f$ and may be grown with
 * bigvector_append().
 *
 * Complexity:
 *   - Time: \f$O(1)\f$
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(1)\f$
 *
 * @param[out] a Vector to initialize.
 */
void bigvector_init_dynamic(bigvector* a);

/**
 * @brief Free all storage of a vector.
 *
 * Complexity:
 *   - Time: \f$O(size)\f$
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(1)\f$
 *
 * @param[in,out] a Vector to free.
 */
void bigvector_free(bigvector* a);

/**
 * @brief Append a copy of a bignum to a dynamic vector.
 *
 * Only allowed on vectors created with bigvector_init_dynamic();
 * otherwise an error is printed and nothing happens. If the capacity
 * is exhausted it is doubled (starting at 4) via realloc and the new
 * slots are initialized to zero bignums.
 *
 * Complexity:
 *   - Time: \f$O(n)\f$ where \f$n =\f$ the size of \f$a\f$ in limbs, amortized
 * \f$O(1)\f$ for the growth
 *   - Auxiliary memory: \f$O(capacity)\f$ during the realloc
 *   - Output memory: \f$O(size + 1)\f$ bignums
 *
 * @param[in,out] v Dynamic vector to append to.
 * @param[in]     a Bignum to append.
 */
void bigvector_append(bigvector* v, bignum* a);

/**
 * @brief Copy vector b into a
 *
 * @param [out] a
 * @param [in] b
 */
void bigvector_copy(bigvector* a, bigvector* b);

/**
 * @brief Set element i of a vector to a copy of v.
 *
 * No-op if i is out of range.
 *
 * Complexity:
 *   - Time: \f$O(n)\f$ where \f$n =\f$ the size of \f$v\f$ in limbs
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(n)\f$ limbs
 *
 * @param[in,out] a Vector to modify.
 * @param[in]     v Value to store.
 * @param[in]     i Element index.
 */
void bigvector_set(bigvector* a, bignum* v, u64 i);

/**
 * @brief Component-wise addition of two vectors: \f$r = a + b\f$.
 *
 * Let \f$d =\f$ a->size.
 *
 * No-op if the sizes differ. \f$r\f$ must have capacity for at least \f$d\f$
 * elements.
 *
 * Complexity:
 *   - Time: \f$O(d \cdot n)\f$ where \f$n =\f$ the size of the elements in
 * limbs
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(d)\f$ bignums
 *
 * @param[out] r Result vector.
 * @param[in]  a First vector.
 * @param[in]  b Second vector.
 */
void bigvector_add(bigvector* r, const bigvector* a, const bigvector* b);

/**
 * @brief Component-wise subtraction of two vectors: \f$r = a - b\f$.
 *
 * Let \f$d =\f$ a->size.
 *
 * No-op if the sizes differ. \f$r\f$ must have capacity for at least \f$d\f$
 * elements.
 *
 * Complexity:
 *   - Time: \f$O(d \cdot n)\f$ where \f$n =\f$ the size of the elements in
 * limbs
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(d)\f$ bignums
 *
 * @param[out] r Result vector.
 * @param[in]  a First vector.
 * @param[in]  b Second vector.
 */
void bigvector_sub(bigvector* r, const bigvector* a, const bigvector* b);

/**
 * @brief Euclidean norm of a vector: \f$r = \sqrt{\sum_i a_i^2}\f$.
 *
 * Let \f$d =\f$ a->size.
 *
 * Only defined for static vectors (an error is printed for dynamic
 * ones). Computes the dot product of \f$a\f$ with itself and takes the
 * integer square root.
 *
 * Complexity:
 *   - Time: \f$O(d \cdot n^2)\f$ where n is the size of the elements in limbs,
 *           plus \f$O(n^2 \log n)\f$ for the square root
 *   - Auxiliary memory: \f$O(n)\f$ limbs
 *   - Output memory: \f$O(n)\f$ limbs
 *
 * @param[out] r Result storing the norm.
 * @param[in]  a Vector.
 */
void bigvector_norm(bignum* r, const bigvector* a);

/**
 * @brief Dot product of two vectors: \f$r = \sum_i a_i \cdot b_i\f$.
 *
 * Let \f$d =\f$ a->size.
 *
 * Only defined for static vectors (an error is printed for dynamic
 * ones); no-op if the sizes differ.
 *
 * Complexity:
 *   - Time: \f$O(d \cdot n^2)\f$ where n is the size of the elements in limbs
 *   - Auxiliary memory: \f$O(n)\f$ limbs
 *   - Output memory: \f$O(n)\f$ limbs
 *
 * @param[out] r Result storing the dot product.
 * @param[in]  a First vector.
 * @param[in]  b Second vector.
 */
void bigvector_dot(bignum* r, const bigvector* a, const bigvector* b);

/**
 * @brief Print a vector to stdout as \f$[a_0, a_1, \ldots]\f$.
 *
 * Complexity:
 *   - Time: \f$O(d \cdot n)\f$ where \f$d =\f$ the size and \f$n =\f$ the
 * element size
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(d \cdot n)\f$ characters written
 *
 * @param[in] a Vector to print.
 */
void bigvector_print(bigvector* a);

/**
 * @brief Print a vector to stdout followed by a newline.
 *
 * Complexity:
 *   - Time: \f$O(d \cdot n)\f$ where \f$d =\f$ the size and \f$n =\f$ the
 * element size
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(d \cdot n)\f$ characters written
 *
 * @param[in] a Vector to print.
 */
void bigvector_println(bigvector* a);

void bigvector_neg(bigvector* r, const bigvector* a);

#endif
