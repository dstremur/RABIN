#ifndef BIGVECTOR_H
#define BIGVECTOR_H

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
 * Allocates exactly d slots (capacity = size = d) and initializes
 * each to a zero bignum. The vector cannot be grown afterwards.
 *
 * Complexity:
 *   Time: O(d)
 *   Auxiliary memory: O(1)
 *   Output memory: O(d) bignums
 *
 * @param[out] a Vector to initialize.
 * @param[in]  d Number of elements.
 */
void bigvector_init(bigvector* a, u64 d);

/**
 * @brief Initialize an empty dynamic vector.
 *
 * The vector starts with size = capacity = 0 and may be grown with
 * bigvector_append().
 *
 * Complexity:
 *   Time: O(1)
 *   Auxiliary memory: O(1)
 *   Output memory: O(1)
 *
 * @param[out] a Vector to initialize.
 */
void bigvector_init_dynamic(bigvector* a);

/**
 * @brief Free all storage of a vector.
 *
 * Complexity:
 *   Time: O(size)
 *   Auxiliary memory: O(1)
 *   Output memory: O(1)
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
 *   Time: O(n) where n is the size of a in limbs, amortized O(1) for
 *         the growth
 *   Auxiliary memory: O(capacity) during the realloc
 *   Output memory: O(size + 1) bignums
 *
 * @param[in,out] v Dynamic vector to append to.
 * @param[in]     a Bignum to append.
 */
void bigvector_append(bigvector* v, bignum* a);

/**
 * @brief Set element i of a vector to a copy of v.
 *
 * No-op if i is out of range.
 *
 * Complexity:
 *   Time: O(n) where n is the size of v in limbs
 *   Auxiliary memory: O(1)
 *   Output memory: O(n) limbs
 *
 * @param[in,out] a Vector to modify.
 * @param[in]     v Value to store.
 * @param[in]     i Element index.
 */
void bigvector_set(bigvector* a, bignum* v, u64 i);

/**
 * @brief Component-wise addition of two vectors: r = a + b.
 *
 * Let d = a->size.
 *
 * No-op if the sizes differ. r must have capacity for at least d
 * elements.
 *
 * Complexity:
 *   Time: O(d * n) where n is the size of the elements in limbs
 *   Auxiliary memory: O(1)
 *   Output memory: O(d) bignums
 *
 * @param[out] r Result vector.
 * @param[in]  a First vector.
 * @param[in]  b Second vector.
 */
void bigvector_add(bigvector* r, const bigvector* a, const bigvector* b);

/**
 * @brief Component-wise subtraction of two vectors: r = a - b.
 *
 * Let d = a->size.
 *
 * No-op if the sizes differ. r must have capacity for at least d
 * elements.
 *
 * Complexity:
 *   Time: O(d * n) where n is the size of the elements in limbs
 *   Auxiliary memory: O(1)
 *   Output memory: O(d) bignums
 *
 * @param[out] r Result vector.
 * @param[in]  a First vector.
 * @param[in]  b Second vector.
 */
void bigvector_sub(bigvector* r, const bigvector* a, const bigvector* b);

/**
 * @brief Euclidean norm of a vector: r = sqrt(sum_i a_i^2).
 *
 * Let d = a->size.
 *
 * Only defined for static vectors (an error is printed for dynamic
 * ones). Computes the dot product of a with itself and takes the
 * integer square root.
 *
 * Complexity:
 *   Time: O(d * n^2) where n is the size of the elements in limbs,
 *         plus O(n^2 log n) for the square root
 *   Auxiliary memory: O(n) limbs
 *   Output memory: O(n) limbs
 *
 * @param[out] r Result storing the norm.
 * @param[in]  a Vector.
 */
void bigvector_norm(bignum* r, const bigvector* a);

/**
 * @brief Dot product of two vectors: r = sum_i a_i * b_i.
 *
 * Let d = a->size.
 *
 * Only defined for static vectors (an error is printed for dynamic
 * ones); no-op if the sizes differ.
 *
 * Complexity:
 *   Time: O(d * n^2) where n is the size of the elements in limbs
 *   Auxiliary memory: O(n) limbs
 *   Output memory: O(n) limbs
 *
 * @param[out] r Result storing the dot product.
 * @param[in]  a First vector.
 * @param[in]  b Second vector.
 */
void bigvector_dot(bignum* r, const bigvector* a, const bigvector* b);

/**
 * @brief Print a vector to stdout as [a_0, a_1, ...].
 *
 * Complexity:
 *   Time: O(d * n) where d is the size and n the element size
 *   Auxiliary memory: O(1)
 *   Output memory: O(d * n) characters written
 *
 * @param[in] a Vector to print.
 */
void bigvector_print(bigvector* a);

/**
 * @brief Print a vector to stdout followed by a newline.
 *
 * Complexity:
 *   Time: O(d * n) where d is the size and n the element size
 *   Auxiliary memory: O(1)
 *   Output memory: O(d * n) characters written
 *
 * @param[in] a Vector to print.
 */
void bigvector_println(bigvector* a);

#endif
