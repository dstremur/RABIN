#ifndef BIGMATRIX_H
#define BIGMATRIX_H

#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>

#include "bigcore.h"
#include "bigvector.h"

// use a flat array structure
// index = (row * nr_cols) + col
typedef struct bigmatrix {
  bignum* data;
  u64 r_size;
  u64 c_size;
} bigmatrix;

#define GET(M, r, c) (&(M)->data[(r) * (M)->c_size + (c)])

/**
 * @brief Initialize a matrix with r rows and c columns of zero bignums.
 *
 * Complexity:
 *   Time: O(r * c)
 *   Auxiliary memory: O(1)
 *   Output memory: O(r * c) bignums
 *
 * @param[out] M Matrix to initialize.
 * @param[in]  r Number of rows.
 * @param[in]  c Number of columns.
 */
void bigmatrix_init(bigmatrix* M, u64 r, u64 c);

/**
 * @brief Free all storage of a matrix.
 *
 * Complexity:
 *   Time: O(r * c)
 *   Auxiliary memory: O(1)
 *   Output memory: O(1)
 *
 * @param[in,out] A Matrix to free.
 */
void bigmatrix_free(bigmatrix* A);

/**
 * @brief Print a matrix to stdout in Python list-of-lists syntax.
 *
 * Complexity:
 *   Time: O(r * c * n) where n is the size of the entries in limbs
 *   Auxiliary memory: O(1)
 *   Output memory: O(r * c * n) characters written
 *
 * @param[in] A Matrix to print.
 */
void bigmatrix_print_python(const bigmatrix* A);

/**
 * @brief Copy a matrix: R = A.
 *
 * Let r = A->r_size, c = A->c_size.
 *
 * No-op if the dimensions do not match.
 *
 * Complexity:
 *   Time: O(r * c * n) where n is the size of the entries in limbs
 *   Auxiliary memory: O(1)
 *   Output memory: O(r * c) bignums
 *
 * @param[out] R Destination matrix.
 * @param[in]  A Source matrix.
 */
void bigmatrix_copy(bigmatrix* R, bigmatrix* A);

/**
 * @brief Get a single element: R = A[r][c].
 *
 * No-op if (r, c) is out of range.
 *
 * Complexity:
 *   Time: O(n) where n is the size of the entry in limbs
 *   Auxiliary memory: O(1)
 *   Output memory: O(n) limbs
 *
 * @param[out] R Result storing the element.
 * @param[in]  A Matrix.
 * @param[in]  r Row index.
 * @param[in]  c Column index.
 */
void bigmatrix_get(bignum* R, const bigmatrix* A, u64 r, u64 c);

/**
 * @brief Set a single element: A[r][c] = a.
 *
 * No-op if (r, c) is out of range.
 *
 * Complexity:
 *   Time: O(n) where n is the size of a in limbs
 *   Auxiliary memory: O(1)
 *   Output memory: O(n) limbs
 *
 * @param[in,out] A Matrix.
 * @param[in]    a Value to store.
 * @param[in]    r Row index.
 * @param[in]    c Column index.
 */
void bigmatrix_set(bigmatrix* A, const bignum* a, u64 r, u64 c);

/**
 * @brief Extract a column into a vector: c = A[:, col].
 *
 * The vector c must already be allocated with r_size elements; a
 * size mismatch is reported but not fatal.
 *
 * Complexity:
 *   Time: O(r * n) where n is the size of the entries in limbs
 *   Auxiliary memory: O(1)
 *   Output memory: O(r) bignums
 *
 * @param[out]  c   Result vector storing the column.
 * @param[in]   A   Matrix.
 * @param[in] col Column index.
 */
void bigmatrix_get_col(bigvector* c, const bigmatrix* A, u64 col);

/**
 * @brief Extract a row into a vector: r = A[row, :].
 *
 * The vector r must already be allocated with c_size elements; a
 * size mismatch is reported but not fatal.
 *
 * Complexity:
 *   Time: O(c * n) where n is the size of the entries in limbs
 *   Auxiliary memory: O(1)
 *   Output memory: O(c) bignums
 *
 * @param[out]  r   Result vector storing the row.
 * @param[in]   A   Matrix.
 * @param[in] row Row index.
 */
void bigmatrix_get_row(bigvector* r, const bigmatrix* A, u64 row);

/**
 * @brief Matrix-vector product: r = A * v.
 *
 * Let r = A->r_size, c = A->c_size.
 *
 * Each output element is the dot product of the corresponding row of
 * A with v. Sizes must match (r has r_size elements, v has c_size).
 *
 * Complexity:
 *   Time: O(r * c * n^2) where n is the size of the entries in limbs
 *   Auxiliary memory: O(c) bignums for the row buffer
 *   Output memory: O(r) bignums
 *
 * @param[out] r Result vector storing the product.
 * @param[in]  A Matrix.
 * @param[in]  v Input vector.
 */
void bigmatrix_mv(bigvector* r, const bigmatrix* A, bigvector* v);

/**
 * @brief Vector-matrix product: r = v * A.
 *
 * Let r = A->r_size, c = A->c_size.
 *
 * Each output element is the dot product of v with the corresponding
 * column of A. Sizes must match (v has r_size elements, r has
 * c_size).
 *
 * Complexity:
 *   Time: O(r * c * n^2) where n is the size of the entries in limbs
 *   Auxiliary memory: O(r) bignums for the column buffer
 *   Output memory: O(c) bignums
 *
 * @param[out] r Result vector storing the product.
 * @param[in]  A Matrix.
 * @param[in]  v Input vector.
 */
void bigmatrix_vm(bigvector* r, const bigmatrix* A, bigvector* v);

/**
 * @brief Hadamard bound on the determinant: r = prod_i ||col_i||.
 *
 * Let n = A->c_size.
 *
 * Sheldon Axler: let c be the max entry, then
 * |det A| <= c^n * n^{n / 2}. The tighter bound used here is the
 * product of the Euclidean norms of the columns. Since the norm is
 * computed with an integer square root, 1 is added to each norm to
 * keep the bound valid.
 *
 * Complexity:
 *   Time: O(n^2 * k^2) where k is the size of the entries in limbs
 *   Auxiliary memory: O(n) bignums for the column buffer
 *   Output memory: O(n * k) limbs
 *
 * @param[out] r Result storing the Hadamard bound.
 * @param[in]  A Matrix.
 */
void bigmatrix_hadamard(bignum* r, const bigmatrix* A);

/**
 * @brief Component-wise addition of two matrices: R = A + B.
 *
 * Let r = A->r_size, c = A->c_size.
 *
 * No-op if the dimensions do not match.
 *
 * Complexity:
 *   Time: O(r * c * n) where n is the size of the entries in limbs
 *   Auxiliary memory: O(1)
 *   Output memory: O(r * c) bignums
 *
 * @param[out] R Result matrix.
 * @param[in]  A First matrix.
 * @param[in]  B Second matrix.
 */
void bigmatrix_add(bigmatrix* R, const bigmatrix* A, const bigmatrix* B);

/**
 * @brief Print a matrix to stdout, one row per line.
 *
 * Complexity:
 *   Time: O(r * c * n) where n is the size of the entries in limbs
 *   Auxiliary memory: O(n) limbs for the temporary
 *   Output memory: O(r * c * n) characters written
 *
 * @param[in] A Matrix to print.
 */
void bigmatrix_print(const bigmatrix* A);

/**
 * @brief Schoolbook matrix multiplication: R = A * B.
 *
 * Let r = A->r_size, k = A->c_size, c = B->c_size.
 *
 * No-op if A->c_size != B->r_size.
 *
 * Complexity:
 *   Time: O(r * k * c * n^2) where n is the size of the entries in
 *         limbs
 *   Auxiliary memory: O(n) limbs for temporaries
 *   Output memory: O(r * c) bignums
 *
 * @param[out] R Result matrix.
 * @param[in]  A First matrix.
 * @param[in]  B Second matrix.
 */
void bigmatrix_mul(bigmatrix* R, const bigmatrix* A, const bigmatrix* B);

/**
 * @brief Determinant of a square bignum matrix.
 *
 * Let n = A->c_size.
 *
 * Currently delegates to the RNS path: estimates the number of primes
 * from the Hadamard bound, builds an RNS context over the first k
 * primes, and reconstructs the determinant with the CRT.
 *
 * The code after the early return is a direct Bareiss (fraction-free
 * Gaussian elimination) implementation, currently disabled ("broken").
 *
 * Optimization ideas for the Bareiss path: use exact division, better
 * cache locality, preallocate using the Hadamard bound, OpenMP,
 * compute mod primes larger than the Hadamard bound and reconstruct
 * with the CRT, or use Jebelean's algorithm.
 *
 * Complexity:
 *   Time: O(k * n^3) for the RNS determinants (parallel over k),
 *         plus O(k * n_b^2) for the CRT where n_b is the size of the
 *         product in limbs
 *   Auxiliary memory: O(n^2) bignums for the copy, O(n^2) u64s per
 *         thread in the RNS path
 *   Output memory: O(n_b) limbs
 *
 * @param[out] d Result storing the determinant.
 * @param[in]  A Square matrix.
 */
void bigmatrix_det(bignum* d, const bigmatrix* A);

#endif
