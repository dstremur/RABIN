#ifndef RMAT_H
#define RMAT_H

#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>

#include "rabin_errors.h"
#include "rpol.h"
#include "rvec.h"
#include "rz.h"

// use a flat array structure
// index = (row * nr_cols) + col
typedef struct rmat_t {
  rz_t* data;
  u64 rows;
  u64 cols;
} rmat_t;

#define RMAT_GET(M, r, c) (&(M)->data[(r) * (M)->cols + (c)])

/**
 * @brief Initialize a matrix with \f$r\f$ rows and \f$c\f$ columns of zero
 * bignums.
 *
 * Complexity:
 *   - Time: \f$O(r \cdot c)\f$
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(r \cdot c)\f$ bignums
 *
 * @param[out] M Matrix to initialize.
 * @param[in]  r Number of rows.
 * @param[in]  c Number of columns.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR, RABIN_ERR_OVERFLOW,
 * or RABIN_ERR_OUT_OF_MEMORY.
 */
rabin_err_t rmat_init(rmat_t* M, u64 r, u64 c);

/**
 * @brief Free all storage of a matrix.
 *
 * Complexity:
 *   - Time: \f$O(r \cdot c)\f$
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(1)\f$
 *
 * @param[in,out] A Matrix to free.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR.
 */
rabin_err_t rmat_clear(rmat_t* A);

/**
 * @brief Print a matrix to stdout in Python list-of-lists syntax.
 *
 * Complexity:
 *   - Time: \f$O(r \cdot c \cdot n)\f$ where \f$n =\f$ the size of the entries
 * in limbs
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(r \cdot c \cdot n)\f$ characters written
 *
 * @param[in] A Matrix to print.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR.
 */
rabin_err_t rmat_print_python(const rmat_t* A);

/**
 * @brief Copy a matrix: \f$R = A\f$.
 *
 * Let \f$r =\f$ A->rows, \f$c =\f$ A->cols.
 *
 * No-op if the dimensions do not match.
 *
 * Complexity:
 *   - Time: \f$O(r \cdot c \cdot n)\f$ where \f$n =\f$ the size of the entries
 * in limbs
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(r \cdot c)\f$ bignums
 *
 * @param[out] R Destination matrix.
 * @param[in]  A Source matrix.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR,
 * RABIN_ERR_OUT_OF_MEMORY, or RABIN_ERR_MATRIX_DIM.
 */
rabin_err_t rmat_copy(rmat_t* R, const rmat_t* A);

/**
 * @brief Get a single element: \f$R = A[r][c]\f$.
 *
 * No-op if (r, c) is out of range.
 *
 * Complexity:
 *   - Time: \f$O(n)\f$ where \f$n =\f$ the size of the entry in limbs
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(n)\f$ limbs
 *
 * @param[out] R Result storing the element.
 * @param[in]  A Matrix.
 * @param[in]  r Row index.
 * @param[in]  c Column index.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR,
 * RABIN_ERR_INVALID_ARG, or RABIN_ERR_OUT_OF_MEMORY.
 */
rabin_err_t rmat_get(rz_t* R, const rmat_t* A, u64 r, u64 c);

/**
 * @brief Set a single element: \f$A[r][c] = a\f$.
 *
 * No-op if (r, c) is out of range.
 *
 * Complexity:
 *   - Time: \f$O(n)\f$ where \f$n =\f$ the size of \f$a\f$ in limbs
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(n)\f$ limbs
 *
 * @param[in,out] A Matrix.
 * @param[in]    a Value to store.
 * @param[in]    r Row index.
 * @param[in]    c Column index.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR,
 * RABIN_ERR_INVALID_ARG, or RABIN_ERR_OUT_OF_MEMORY.
 */
rabin_err_t rmat_set(rmat_t* A, const rz_t* a, u64 r, u64 c);

/**
 * @brief Extract a column into a vector: \f$c = A[:, col]\f$.
 *
 * The vector \f$c\f$ must already be allocated with rows elements; a
 * size mismatch is reported but not fatal.
 *
 * Complexity:
 *   - Time: \f$O(r \cdot n)\f$ where \f$n =\f$ the size of the entries in limbs
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(r)\f$ bignums
 *
 * @param[out]  c   Result vector storing the column.
 * @param[in]   A   Matrix.
 * @param[in] col Column index.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR,
 * RABIN_ERR_INVALID_ARG, or RABIN_ERR_MATRIX_DIM.
 */
rabin_err_t rmat_get_col(rvec_t* c, const rmat_t* A, u64 col);

/**
 * @brief Extract a row into a vector: \f$r = A[row, :]\f$.
 *
 * The vector \f$r\f$ must already be allocated with cols elements; a
 * size mismatch is reported but not fatal.
 *
 * Complexity:
 *   - Time: \f$O(c \cdot n)\f$ where \f$n =\f$ the size of the entries in limbs
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(c)\f$ bignums
 *
 * @param[out]  r   Result vector storing the row.
 * @param[in]   A   Matrix.
 * @param[in] row Row index.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR,
 * RABIN_ERR_INVALID_ARG, or RABIN_ERR_MATRIX_DIM.
 */
rabin_err_t rmat_get_row(rvec_t* r, const rmat_t* A, u64 row);

/**
 * @brief Test whether a matrix is square (\f$r = c\f$).
 *
 * Complexity:
 *   - Time: \f$O(1)\f$
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(1)\f$
 *
 * @param[in] A Matrix.
 *
 * @return true If the number of rows equals the number of columns.
 */
bool rmat_is_square(const rmat_t* A);

/**
 * @brief Test whether all entries of a matrix are zero.
 *
 * Complexity:
 *   - Time: \f$O(r \cdot c)\f$ where \f$r =\f$ A->rows, \f$c =\f$ A->cols,
 *           stopping at the first nonzero entry
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(1)\f$
 *
 * @param[in] A Matrix.
 *
 * @return true If every entry is zero.
 */
bool rmat_is_zero(const rmat_t* A);

/**
 * @brief Test whether a matrix is the identity matrix.
 *
 * The matrix must be square with 1s on the main diagonal and 0s
 * elsewhere. Non-square matrices are rejected before any entry is
 * inspected.
 *
 * Complexity:
 *   - Time: \f$O(n^2)\f$ where \f$n =\f$ A->rows, stopping at the first
 *   violation
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(1)\f$
 *
 * @param[in] A Matrix.
 *
 * @return true If A is the identity matrix.
 */
bool rmat_is_identity(const rmat_t* A);

/**
 * @brief Test whether all off-diagonal entries of a matrix are zero.
 *
 * Works for rectangular matrices: the diagonal is the set of entries
 * \f$A[i][i]\f$ with \f$i < \min(r, c)\f$; the diagonal entries
 * themselves are not inspected.
 *
 * Complexity:
 *   - Time: \f$O(r \cdot c)\f$ in the worst case, stopping at the first
 *   nonzero off-diagonal entry
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(1)\f$
 *
 * @param[in] A Matrix.
 *
 * @return true If all off-diagonal entries are zero.
 */
bool rmat_is_diagonal(const rmat_t* A);

/**
 * @brief Test whether all entries below the main diagonal are zero.
 *
 * Only the sub-regions \f$A[i][j]\f$ with \f$i > j\f$ are inspected.
 *
 * Complexity:
 *   - Time: \f$O(\min(r, c)^2 / 2)\f$ in the worst case, stopping at the
 *   first nonzero entry below the diagonal
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(1)\f$
 *
 * @param[in] A Matrix.
 *
 * @return true If all entries with row > col are zero.
 */
bool rmat_is_upper_triangular(const rmat_t* A);

/**
 * @brief Test whether all entries above the main diagonal are zero.
 *
 * Only the sub-regions \f$A[i][j]\f$ with \f$i < j\f$ are inspected.
 *
 * Complexity:
 *   - Time: \f$O(\min(r, c)^2 / 2)\f$ in the worst case, stopping at the
 *   first nonzero entry above the diagonal
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(1)\f$
 *
 * @param[in] A Matrix.
 *
 * @return true If all entries with row < col are zero.
 */
bool rmat_is_lower_triangular(const rmat_t* A);

/**
 * @brief Test whether a matrix is symmetric (\f$A[i][j] = A[j][i]\f$).
 *
 * The matrix must be square; only the strict upper triangle is walked,
 * comparing each mirror pair once. Non-square matrices are rejected
 * before any entry is inspected.
 *
 * Complexity:
 *   - Time: \f$O(n^2 \cdot n_{b})\f$ where \f$n =\f$ A->rows and
 *   \f$n_{b} =\f$ the size of the entries in limbs, stopping at the
 *   first mismatching pair
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(1)\f$
 *
 * @param[in] A Matrix.
 *
 * @return true If A is square and A[i][j] == A[j][i] for all (i, j).
 */
bool rmat_is_symmetric(const rmat_t* A);

/**
 * @brief Scalar multiplication of a matrix: \f$R = A \cdot a\f$.
 *
 * Let \f$r =\f$ A->rows, \f$c =\f$ A->cols.
 *
 * Every entry is multiplied by \f$a\f$. R must already be initialized
 * with the same dimensions as A.
 *
 * Complexity:
 *   - Time: \f$O(r \cdot c \cdot n)\f$ where \f$n =\f$ the size of the
 * entries in limbs
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(r \cdot c)\f$ bignums
 *
 * @param[out] R Result matrix.
 * @param[in]  A Matrix to scale.
 * @param[in]  a Scalar factor.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR,
 * RABIN_ERR_MATRIX_DIM, or RABIN_ERR_OUT_OF_MEMORY.
 */
rabin_err_t rmat_scalar(rmat_t* R, const rmat_t* A, const rz_t* a);

/**
 * @brief Exact scalar division of a matrix: \f$R = A / a\f$.
 *
 * Let \f$r =\f$ A->rows, \f$c =\f$ A->cols.
 *
 * Every entry is divided by \f$a\f$; the division must be exact. R must
 * already be initialized with the same dimensions as A.
 *
 * Complexity:
 *   - Time: \f$O(r \cdot c \cdot n)\f$ where \f$n =\f$ the size of the
 * entries in limbs
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(r \cdot c)\f$ bignums
 *
 * @param[out] R Result matrix.
 * @param[in]  A Matrix to divide.
 * @param[in]  a Divisor scalar; must be nonzero and must divide every
 * entry of A.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR,
 * RABIN_ERR_DIV_BY_ZERO, RABIN_ERR_MATRIX_DIM, or
 * RABIN_ERR_OUT_OF_MEMORY.
 */
rabin_err_t rmat_div_exact_scalar(rmat_t* R, const rmat_t* A, const rz_t* a);

/**
 * @brief Matrix-vector product: \f$r = A \cdot v\f$.
 *
 * Let \f$r =\f$ A->rows, \f$c =\f$ A->cols.
 *
 * Each output element is the dot product of the corresponding row of
 * \f$A\f$ with \f$v\f$. Sizes must match (\f$r\f$ has rows elements, \f$v\f$
 * has cols).
 *
 * Complexity:
 *   - Time: \f$O(r \cdot c \cdot n^2)\f$ where \f$n =\f$ the size of the
 * entries in limbs
 *   - Auxiliary memory: \f$O(c)\f$ bignums for the row buffer
 *   - Output memory: \f$O(r)\f$ bignums
 *
 * @param[out] r Result vector storing the product.
 * @param[in]  A Matrix.
 * @param[in]  v Input vector.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR,
 * RABIN_ERR_OUT_OF_MEMORY, or RABIN_ERR_MATRIX_DIM.
 */
rabin_err_t rmat_mv(rvec_t* r, const rmat_t* A, rvec_t* v);

/**
 * @brief Vector-matrix product: \f$r = v \cdot A\f$.
 *
 * Let \f$r =\f$ A->rows, \f$c =\f$ A->cols.
 *
 * Each output element is the dot product of \f$v\f$ with the corresponding
 * column of \f$A\f$. Sizes must match (\f$v\f$ has rows elements, \f$r\f$ has
 * cols).
 *
 * Complexity:
 *   - Time: \f$O(r \cdot c \cdot n^2)\f$ where \f$n =\f$ the size of the
 * entries in limbs
 *   - Auxiliary memory: \f$O(r)\f$ bignums for the column buffer
 *   - Output memory: \f$O(c)\f$ bignums
 *
 * @param[out] r Result vector storing the product.
 * @param[in]  A Matrix.
 * @param[in]  v Input vector.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR,
 * RABIN_ERR_OUT_OF_MEMORY, or RABIN_ERR_MATRIX_DIM.
 */
rabin_err_t rmat_vm(rvec_t* r, const rmat_t* A, rvec_t* v);

/**
 * @brief Hadamard bound on the determinant: \f$r = \prod_i \lVert col_i
 * \rVert\f$.
 *
 * Let \f$n =\f$ A->cols.
 *
 * Sheldon Axler: let \f$c\f$ be the max entry, then
 * \f$|\det A| \le c^n \cdot n^{n / 2}\f$. The tighter bound used here is the
 * product of the Euclidean norms of the columns. Since the norm is
 * computed with an integer square root, \f$1\f$ is added to each norm to
 * keep the bound valid.
 *
 * Complexity:
 *   - Time: \f$O(n^2 \cdot k^2)\f$ where \f$k =\f$ the size of the entries in
 * limbs
 *   - Auxiliary memory: \f$O(n)\f$ bignums for the column buffer
 *   - Output memory: \f$O(n \cdot k)\f$ limbs
 *
 * @param[out] r Result storing the Hadamard bound.
 * @param[in]  A Matrix.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR /
 * RABIN_ERR_OUT_OF_MEMORY.
 */
rabin_err_t rmat_hadamard(rz_t* r, const rmat_t* A);

/**
 * @brief Component-wise addition of two matrices: \f$R = A + B\f$.
 *
 * Let \f$r =\f$ A->rows, \f$c =\f$ A->cols.
 *
 * No-op if the dimensions do not match.
 *
 * Complexity:
 *   - Time: \f$O(r \cdot c \cdot n)\f$ where \f$n =\f$ the size of the entries
 * in limbs
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(r \cdot c)\f$ bignums
 *
 * @param[out] R Result matrix.
 * @param[in]  A First matrix.
 * @param[in]  B Second matrix.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR,
 * RABIN_ERR_OUT_OF_MEMORY, or RABIN_ERR_MATRIX_DIM.
 */
rabin_err_t rmat_add(rmat_t* R, const rmat_t* A, const rmat_t* B);

/**
 * @brief Component-wise subtraction of two matrices: \f$R = A - B\f$.
 *
 * Let \f$r =\f$ A->rows, \f$c =\f$ A->cols.
 *
 * R must already be initialized with the same dimensions as A and B.
 *
 * Complexity:
 *   - Time: \f$O(r \cdot c \cdot n)\f$ where \f$n =\f$ the size of the entries
 * in limbs
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(r \cdot c)\f$ bignums
 *
 * @param[out] R Result matrix.
 * @param[in]  A First matrix.
 * @param[in]  B Second matrix.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR,
 * RABIN_ERR_MATRIX_DIM, or RABIN_ERR_OUT_OF_MEMORY.
 */
rabin_err_t rmat_sub(rmat_t* R, const rmat_t* A, const rmat_t* B);

/**
 * @brief Print a matrix to stdout, one row per line, with right-aligned
 * columns.
 *
 * The width of each column is derived from its widest entry. Entries whose
 * decimal form is wider than 16 characters are truncated to their leading
 * digits followed by an ellipsis, so large matrices stay compact. A single
 * oversized entry (e.g. the last invariant factor of a Smith form sitting
 * in a column of zeros) does not inflate its column: such outliers overflow
 * the column to the right instead. Use rmat_print_full() for
 * untruncated output and rmat_print_tail() for truncation to the
 * trailing digits instead.
 *
 * Complexity:
 *   - Time: \f$O(r \cdot c \cdot n)\f$ where \f$n =\f$ the size of the entries
 * in limbs
 *   - Auxiliary memory: \f$O(c)\f$ size_t for the column widths, plus one
 *   temporary string per entry
 *   - Output memory: \f$O(r \cdot c \cdot \min(n, 16))\f$ characters written
 *
 * @param[in] A Matrix to print.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR /
 * RABIN_ERR_OUT_OF_MEMORY.
 */
rabin_err_t rmat_print(const rmat_t* A);

/**
 * @brief Print a matrix to stdout, one row per line, with right-aligned
 * columns and untruncated (full) entries.
 *
 * Complexity:
 *   - Time: \f$O(r \cdot c \cdot n)\f$ where \f$n =\f$ the size of the entries
 * in limbs
 *   - Auxiliary memory: \f$O(c)\f$ size_t for the column widths, plus one
 *   temporary string per entry
 *   - Output memory: \f$O(r \cdot c \cdot n)\f$ characters written
 *
 * @param[in] A Matrix to print.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR /
 * RABIN_ERR_OUT_OF_MEMORY.
 */
rabin_err_t rmat_print_full(const rmat_t* A);

/**
 * @brief Print a matrix to stdout, one row per line, with right-aligned
 * columns.
 *
 * Like rmat_print(), entries wider than 16 characters are truncated,
 * but to an ellipsis followed by their trailing digits.
 *
 * Complexity:
 *   - Time: \f$O(r \cdot c \cdot n)\f$ where \f$n =\f$ the size of the entries
 * in limbs
 *   - Auxiliary memory: \f$O(c)\f$ size_t for the column widths, plus one
 *   temporary string per entry
 *   - Output memory: \f$O(r \cdot c \cdot \min(n, 16))\f$ characters written
 *
 * @param[in] A Matrix to print.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR /
 * RABIN_ERR_OUT_OF_MEMORY.
 */
rabin_err_t rmat_print_tail(const rmat_t* A);

/**
 * @brief Schoolbook matrix multiplication: \f$R = A \cdot B\f$.
 *
 * Let \f$r =\f$ A->rows, \f$k =\f$ A->cols, \f$c =\f$ B->cols.
 *
 * No-op if A->cols != B->rows.
 *
 * Complexity:
 *   - Time: \f$O(r \cdot k \cdot c \cdot n^2)\f$ where \f$n =\f$ the size of
 * the entries in limbs
 *   - Auxiliary memory: \f$O(n)\f$ limbs for temporaries
 *   - Output memory: \f$O(r \cdot c)\f$ bignums
 *
 * @param[out] R Result matrix.
 * @param[in]  A First matrix.
 * @param[in]  B Second matrix.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR,
 * RABIN_ERR_OUT_OF_MEMORY, or RABIN_ERR_MATRIX_DIM.
 */
rabin_err_t rmat_mul(rmat_t* R, const rmat_t* A, const rmat_t* B);

/**
 * @brief Determinant of a square rz_t matrix.
 *
 * Let \f$n =\f$ A->cols.
 *
 * Currently delegates to the RNS path: estimates the number of primes
 * from the Hadamard bound, builds an RNS context over the first \f$k\f$
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
 *   - Time: \f$O(k \cdot n^3)\f$ for the RNS determinants (parallel over k),
 *           plus \f$O(k \cdot n_{b}^2)\f$ for the CRT where \f$n_{b} =\f$ the
 * size of the product in limbs
 *   - Auxiliary memory: \f$O(n^2)\f$ bignums for the copy, \f$O(n^2)\f$ u64s
 * per thread in the RNS path
 *   - Output memory: \f$O(n_{b})\f$ limbs
 *
 * @param[out] d Result storing the determinant.
 * @param[in]  A Square matrix.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR,
 * RABIN_ERR_OUT_OF_MEMORY, or RABIN_ERR_MATRIX_DIM.
 */
rabin_err_t rmat_det(rz_t* d, const rmat_t* A);

/**
 * @brief Determinant of a square matrix via signed fraction-free (Bareiss)
 * elimination.
 *
 * Let \f$n =\f$ A->rows.
 *
 * Runs the exact division variant of Bareiss elimination on a copy of A,
 * so the input is not modified; all intermediate values stay integral and
 * grow only to the size of the leading principal minors. The elimination
 * is parallelized across threads with OpenMP; pivot search and row swaps
 * remain sequential. \f$det = 1\f$ for the \f$0 \times 0\f$ matrix.
 *
 * Complexity:
 *   - Time: \f$O(n^3 \cdot d^2)\f$ where \f$d =\f$ the entry size in limbs
 *   - Auxiliary memory: \f$O(n^2)\f$ bignums (working copy and temporaries)
 *   - Output memory: \f$O(d)\f$ limbs
 *
 * @param[out] det Determinant.
 * @param[in]  A Square matrix.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR,
 * RABIN_ERR_MATRIX_DIM, or RABIN_ERR_OUT_OF_MEMORY.
 * @see rmat_det()
 */
rabin_err_t rmat_det_bareiss(rz_t* det, const rmat_t* A);

/**
 * @brief Returns the identity matrix
 *
 * @param I
 * @param n
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR,
 * RABIN_ERR_OUT_OF_MEMORY, or RABIN_ERR_MATRIX_DIM.
 */
rabin_err_t rmat_id(rmat_t* I, const u64 n);

/**
 * @brief Greatest common divisor of all entries of a matrix.
 *
 * Let \f$d =\f$ A->rows \f$\cdot\f$ A->cols.
 *
 * Computes \f$g = \gcd(A[0][0], A[0][1], \ldots)\f$, skipping zero
 * entries. \f$g = 1\f$ for an empty matrix.
 *
 * Complexity:
 *   - Time: \f$O(d \cdot n^2)\f$ where \f$n =\f$ the entry size in limbs
 *   - Auxiliary memory: \f$O(n)\f$ limbs
 *   - Output memory: \f$O(n)\f$ limbs
 *
 * @param[out] g GCD of all entries.
 * @param[in]  A Matrix.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR /
 * RABIN_ERR_OUT_OF_MEMORY.
 */
rabin_err_t rmat_gcd_all(rz_t* g, const rmat_t* A);

/**
 * @brief Swap the headers of two initialized matrices.
 *
 * Exchanges the dimension fields and the data pointer, so the matrices
 * trade their contents. Both must already be initialized; the caller
 * keeps owning both.
 *
 * Complexity:
 *   - Time: \f$O(1)\f$
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(1)\f$
 *
 * @param[in,out] a First matrix.
 * @param[in,out] b Second matrix.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR.
 */
rabin_err_t rmat_swap(rmat_t* a, rmat_t* b);

/**
 * @brief Trace of a square matrix: \f$t = \sum_i A[i][i]\f$.
 *
 * Let \f$n =\f$ A->rows.
 *
 * @param[out] t Trace.
 * @param[in]  A Square matrix.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR,
 * RABIN_ERR_MATRIX_DIM, or RABIN_ERR_OUT_OF_MEMORY.
 */
rabin_err_t rmat_trace(rz_t* t, const rmat_t* A);

/**
 * @brief Characteristic polynomial and adjugate of a square matrix
 * (Danilevsky's algorithm).
 *
 * Let \f$n =\f$ A->rows.
 *
 * Computes the monic characteristic polynomial \f$p(x) =
 * \det(xI - A)\f$ of degree \f$n\f$ into \p p, and, if \p J is not NULL,
 * the adjugate matrix \f$\operatorname{adj}(A)\f$ into \p J.
 *
 * Complexity:
 *   - Time: \f$O(n^4 \cdot d^2)\f$ where \f$d =\f$ the entry size in limbs
 *   - Auxiliary memory: \f$O(n^2)\f$ bignums
 *   - Output memory: \f$O(n)\f$ bignums (polynomial),
 * \f$O(n^2)\f$ bignums (adjugate)
 *
 * @param[out] p Characteristic polynomial (must be zero-initialized).
 * @param[out] J Adjugate matrix; ignored if NULL, otherwise must already
 * be initialized with the same dimensions as A.
 * @param[in]  A Square matrix.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR,
 * RABIN_ERR_MATRIX_DIM, or RABIN_ERR_OUT_OF_MEMORY.
 */
rabin_err_t rmat_charpoly_adj(rpol_t* p, rmat_t* J, const rmat_t* A);

/**
 * @brief Compute the LLL algorithm on the basis matrix B
 *
 * @note TODO: Not yet implemented (placeholder declaration).
 *
 * @param B
 * @param n
 * @param delta
 * @param H
 */
rabin_err_t rmat_LLL(rmat_t* B, u64 n, double delta, rmat_t* H);

/**
 * @brief Hermite normal form via the exact (coefficient-blowup prone)
 * column reduction (Cohen's Algorithm 2.4.4).
 *
 * Let \f$m =\f$ A->rows, \f$n =\f$ A->cols.
 *
 * Handles zero entries in the pivot positions \f$a_{i,k}\f$ (the extended
 * GCD is safe on a zero argument). This is the reference implementation:
 * use it to cross-check rmat_hermite_mod_d().
 *
 * @param[out] W Result storing the HNF of A.
 * @param[in]  A Input matrix.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR /
 * RABIN_ERR_OUT_OF_MEMORY.
 */
rabin_err_t rmat_hermite(rmat_t* W, const rmat_t* A);

/**
 * @brief Hermite normal form, Cohen's Algorithm 2.4.5 (exact integer
 * arithmetic).
 *
 * Let \f$m =\f$ A->rows, \f$n =\f$ A->cols.
 *
 * Works entirely with integers and Euclidean divisions, but is subject to
 * the coefficient explosion phenomenon; prefer
 * rmat_hermite_mod_d() when a multiple of the module determinant is
 * known.
 *
 * Handles zero entries in the pivot positions, and each Euclidean step uses
 * the Bezout pair normalized to the canonical range of Cohen's "Important
 * Remark" (v = 0 when a_{i,k} | a_{i,j}, |v| <= |a_{i,k} / d| / 2
 * otherwise), which keeps the intermediate coefficients small.
 *
 * @param[out] W Result storing the HNF of A.
 * @param[in]  A Input matrix.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR /
 * RABIN_ERR_OUT_OF_MEMORY.
 */
rabin_err_t rmat_hermite_gcd(rmat_t* W, const rmat_t* A);

/**
 * @brief Hermite normal form modulo D, Cohen's Algorithm 2.4.8 (essentially
 * due to Domich et al. [DKT]).
 *
 * Let \f$m =\f$ A->rows, \f$n =\f$ A->cols.
 *
 * Column operations are reduced modulo a running modulus \f$R\f$ (starting
 * at \f$D\f$, shrinking by a factor \f$\gcd(a_{i,k}, R)\f$ per row) with
 * residues taken in \f$(-R/2, R/2]\f$, which keeps the coefficients bounded
 * and avoids the coefficient explosion of the exact algorithms. The lift
 * \f$W_i = u A_k \bmod R\f$ is taken in \f$[0, R)\f$.
 *
 * Preconditions: \f$m \le n\f$, A has rank \f$m\f$, and \f$D > 0\f$ is a
 * multiple of the determinant \f$\Delta\f$ of the
 * \f$\mathbb{Z}\f$-module generated by the columns of A, i.e. the GCD of
 * all \f$m \times m\f$ minors of A. If \f$D\f$ is not such a multiple the
 * result is incorrect. For square A, \f$\Delta = |\det A|\f$, so
 * \f$D = |\det A|\f$ (e.g. from rmat_det()) always works. A is allowed
 * to have zero entries, including at the trailing pivot \f$a_{m - 1, n - 1}\f$.
 * If \f$m > n\f$, the computation falls back to
 * rmat_hermite_gcd().
 *
 * Each Euclidean step uses the Bezout pair normalized to the canonical range
 * of Cohen's "Important Remark", which keeps the lift coefficients bounded
 * (and is required for the \f$W_i = u A_k \bmod R\f$ lift to be correct).
 *
 * Complexity:
 *  - Time: \f$O(n^2 \cdot m \cdot T(R))\f$ where \f$T(R)\f$ is the cost of a
 *          modular reduction of a \f$\mathcal{O}(R)\f$-sized integer, and
 *          \f$R\f$ shrinks by at least a factor 2 on average per row
 *  - Auxiliary memory: \f$O(n \cdot m)\f$ bignums for the working copy
 *  - Output memory: \f$O(m^2)\f$ bignums
 *
 * @param[out] W Result storing the \f$m \times m\f$ HNF of A.
 * @param[in]  A Input matrix of rank \f$m\f$.
 * @param[in]  D Positive multiple of the module determinant.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR,
 * RABIN_ERR_DIV_BY_ZERO, or RABIN_ERR_OUT_OF_MEMORY.
 */
rabin_err_t rmat_hermite_mod_d(rmat_t* W, const rmat_t* A, const rz_t* D);

/**
 * @brief Smith normal form of a square matrix.
 *
 * Let \f$n =\f$ A->rows.
 *
 * Reduces a copy of A to diagonal Smith form \f$S = U A V\f$ with
 * diagonal entries \f$d_1 \mid d_2 \mid \ldots \mid d_n\f$ using integer
 * row and column operations; the input is not modified. S must already be
 * initialized with the same dimensions as A.
 *
 * Complexity:
 *   - Time: \f$O(n^3 \cdot d^2)\f$ where \f$d =\f$ the entry size in limbs
 * (worst case, dominated by the gcd reductions)
 *   - Auxiliary memory: \f$O(n^2)\f$ bignums (working copy and
 * temporaries)
 *   - Output memory: \f$O(n^2)\f$ bignums
 *
 * @param[out] S Smith normal form of A.
 * @param[in]  A Square matrix.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR,
 * RABIN_ERR_MATRIX_DIM, or RABIN_ERR_OUT_OF_MEMORY.
 */
rabin_err_t rmat_smith(rmat_t* S, const rmat_t* A);

/**
 * @brief Test two matrices for entrywise equality.
 *
 * Complexity:
 *   - Time: \f$O(r \cdot c \cdot n)\f$ where \f$n =\f$ the size of the
 *   entries in limbs, stopping at the first differing entry
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(1)\f$
 *
 * @param[in] A First matrix.
 * @param[in] B Second matrix.
 *
 * @return true If A and B have the same shape and equal entries.
 */
bool rmat_equal(const rmat_t* A, const rmat_t* B);

/**
 * @brief Validate that H strictly satisfies the canonical Hermite normal
 * form structure.
 *
 * Let \f$r =\f$ H->rows, \f$c =\f$ H->cols.
 *
 * Checks, with early exit on the first violation:
 *   - Upper triangular: \f$H[i][j] = 0\f$ for all \f$i > j\f$.
 *   - Positive pivots: the diagonal pivot of every nonzero row is
 *     strictly positive.
 *   - Reduced entries: for every nonzero row \f$i\f$, the entries right
 *     of the pivot, \f$H[i][j]\f$ for \f$j > i\f$, satisfy
 *     \f$0 \le H[i][j] < H[i][i]\f$. (This is the reduction invariant
 *     of this library's column-operation HNF algorithms: each row is
 *     reduced modulo its own pivot.)
 *   - Zero-row order: all zero rows are grouped strictly at the bottom.
 *
 * Complexity:
 *   - Time: \f$O(r \cdot c \cdot n)\f$ where \f$n =\f$ the size of the
 *   entries in limbs, stopping at the first violation
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(1)\f$
 *
 * @param[in] H Matrix to validate.
 *
 * @return true If H is in canonical Hermite normal form.
 */
bool rmat_hnf_check_structure(const rmat_t* H);

/**
 * @brief Verify the unimodular transformation \f$A \cdot U = H\f$.
 *
 * This library's HNF routines act by column operations, so the
 * transformation matrix \f$U\f$ right-multiplies: \f$H = A \cdot U\f$.
 *
 * The product is computed with rmat_mul() into a temporary matrix
 * and compared entrywise with H; the temporary is freed before
 * returning.
 *
 * Complexity:
 *   - Time: \f$O(r \cdot c \cdot n \cdot n_{b}^2)\f$ where \f$r =\f$ A->rows,
 *   \f$c =\f$ A->cols, \f$n =\f$ U->cols, \f$n_{b} =\f$ the size of
 *   the entries in limbs
 *   - Auxiliary memory: \f$O(r \cdot n)\f$ bignums for the product
 *   - Output memory: \f$O(1)\f$
 *
 * @param[in] A Original matrix.
 * @param[in] H Claimed normal form of A.
 * @param[in] U Transformation matrix (expected square, \f$c =\f$ A->cols).
 *
 * @return true If the dimensions are compatible and A * U == H entrywise.
 */
bool rmat_hnf_check_transformation(const rmat_t* A, const rmat_t* H,
                                   const rmat_t* U);

/**
 * @brief Verify that U is unimodular (integer matrix with
 * \f$|\det U| = 1\f$).
 *
 * Complexity:
 *   - Time: same as rmat_det() for an \f$n \times n\f$ matrix
 *   - Auxiliary memory: as in rmat_det()
 *   - Output memory: \f$O(1)\f$
 *
 * @param[in] U Matrix to validate.
 *
 * @return true If U is square and |det(U)| == 1.
 */
bool rmat_hnf_check_unimodular(const rmat_t* U);

/**
 * @brief Master HNF verification: structure, transformation, and
 * unimodularity.
 *
 * Runs rmat_hnf_check_structure(),
 * rmat_hnf_check_transformation(), and
 * rmat_hnf_check_unimodular() in that order, returning early on the
 * first failure.
 *
 * Note: the HNF routines of this library only return H; the
 * transformation matrix U must be known independently (e.g. computed
 * alongside the elimination). For non-square A the true U is rectangular
 * (a column subset of a unimodular matrix), so the unimodularity check
 * only applies to square inputs.
 *
 * @param[in] A Original matrix.
 * @param[in] H Claimed Hermite normal form of A.
 * @param[in] U Transformation matrix with H = A * U.
 *
 * @return true Only if all three invariants hold.
 */
bool rmat_hnf_verify(const rmat_t* A, const rmat_t* H, const rmat_t* U);
#endif
