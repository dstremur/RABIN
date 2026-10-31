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

/// Rational matrix with a common denominator; defined in rmatq.h.
typedef struct rmatq_t rmatq_t;

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
 * Returns RABIN_ERR_MATRIX_DIM if the dimensions do not match.
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
 * Returns RABIN_ERR_INVALID_ARG if (r, c) is out of range.
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
 * Returns RABIN_ERR_INVALID_ARG if (r, c) is out of range.
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
 * The vector \f$c\f$ must already be allocated with \f$A->rows\f$ elements;
 * a size mismatch returns RABIN_ERR_MATRIX_DIM.
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
 * The vector \f$r\f$ must already be allocated with \f$A->cols\f$ elements;
 * a size mismatch returns RABIN_ERR_MATRIX_DIM.
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
 *
 * @par Algorithm Reference:
 * J. Hadamard, "Note sur les Determinants," Comptes Rendus de l'Academie
 * des Sciences, vol. 117, 1893.
 */
rabin_err_t rmat_hadamard(rz_t* r, const rmat_t* A);

/**
 * @brief Component-wise addition of two matrices: \f$R = A + B\f$.
 *
 * Let \f$r =\f$ A->rows, \f$c =\f$ A->cols.
 *
 * Returns RABIN_ERR_MATRIX_DIM if the dimensions do not match.
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
 * Returns RABIN_ERR_MATRIX_DIM if \f$A->cols \neq B->rows\f$.
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
 * Delegates to the RNS path: estimates the number of primes from the
 * Hadamard bound, builds an RNS context over the first \f$k\f$ primes
 * (rns_ctx_init()), and reconstructs the determinant with the CRT
 * (rmat_det_rns()). For a direct exact-arithmetic determinant see
 * rmat_det_bareiss().
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
 *
 * @par Algorithm Reference:
 * D. H. Bareiss, "A Division-Free Analogue of the Gaussian Algorithm
 * for Solving Linear Systems of Linear Equations," Numerische
 * Mathematik, vol. 22, no. 2, 1974.
 * @see rmat_det()
 */
rabin_err_t rmat_det_bareiss(rz_t* det, const rmat_t* A);

/**
 * @brief Set an \f$n \times n\f$ matrix to the identity matrix.
 *
 * Stores 1 on the main diagonal and 0 elsewhere. The matrix must already
 * be initialized with dimensions \f$n \times n\f$.
 *
 * Complexity:
 *   - Time: \f$O(n^2)\f$
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(n^2)\f$ bignums
 *
 * @param[in,out] I Matrix to set (already initialized to \f$n \times n\f$).
 * @param[in]     n Matrix dimension.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR,
 * RABIN_ERR_MATRIX_DIM, or RABIN_ERR_OUT_OF_MEMORY.
 */
rabin_err_t rmat_id(rmat_t* I, const u64 n);

/**
 * @brief Set all entries of a matrix to zero, in place.
 *
 * Complexity:
 *   - Time: \f$O(r \cdot c)\f$
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(1)\f$
 *
 * @param[in,out] A Matrix to clear in place.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR.
 * @see rmat_ones(), rmat_id()
 */
rabin_err_t rmat_zero(rmat_t* A);

/**
 * @brief Set all entries of a matrix to one, in place.
 *
 * Complexity:
 *   - Time: \f$O(r \cdot c)\f$
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(1)\f$
 *
 * @param[in,out] A Matrix to set in place.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR,
 * or RABIN_ERR_OUT_OF_MEMORY.
 * @see rmat_zero()
 */
rabin_err_t rmat_ones(rmat_t* A);

/**
 * @brief Build the diagonal matrix \f$D = \mathrm{diag}(v)\f$.
 *
 * \f$D\f$ must already be initialized as \f$|v| \times |v|\f$; every
 * entry of \f$D\f$ is overwritten.
 *
 * Complexity:
 *   - Time: \f$O(n^2 \cdot d)\f$ where \f$n = |v|\f$ and \f$d =\f$ the
 * entry size in limbs
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(n^2)\f$ bignums
 *
 * @param[out] D Diagonal result (\f$|v| \times |v|\f$).
 * @param[in]  v Vector of diagonal entries.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR,
 * RABIN_ERR_MATRIX_DIM, or RABIN_ERR_OUT_OF_MEMORY.
 * @see rmat_id()
 */
rabin_err_t rmat_diag(rmat_t* D, const rvec_t* v);

/**
 * @brief Swap rows \f$i\f$ and \f$j\f$ of a matrix, in place.
 *
 * Complexity:
 *   - Time: \f$O(c)\f$ pointer swaps
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(1)\f$
 *
 * @param[in,out] A Matrix.
 * @param[in]     i First row index.
 * @param[in]     j Second row index.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR,
 * or RABIN_ERR_INVALID_ARG (index out of range).
 * @see rmat_swap_cols()
 */
rabin_err_t rmat_swap_rows(rmat_t* A, u64 i, u64 j);

/**
 * @brief Swap columns \f$i\f$ and \f$j\f$ of a matrix, in place.
 *
 * Complexity:
 *   - Time: \f$O(r)\f$ pointer swaps
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(1)\f$
 *
 * @param[in,out] A Matrix.
 * @param[in]     i First column index.
 * @param[in]     j Second column index.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR,
 * or RABIN_ERR_INVALID_ARG (index out of range).
 * @see rmat_swap_rows()
 */
rabin_err_t rmat_swap_cols(rmat_t* A, u64 i, u64 j);

/**
 * @brief Horizontal concatenation: \f$R = [A \,|\, B]\f$.
 *
 * The matrices must have the same number of rows; \f$R\f$ must already
 * be initialized as \f$A\text{.rows} \times (A\text{.cols} + B\text{.cols})\f$.
 *
 * Complexity:
 *   - Time: \f$O(r \cdot (c_1 + c_2) \cdot d)\f$
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(r \cdot (c_1 + c_2))\f$ bignums
 *
 * @param[out] R Result matrix.
 * @param[in]  A Left block.
 * @param[in]  B Right block.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR,
 * RABIN_ERR_MATRIX_DIM, RABIN_ERR_OVERFLOW, or RABIN_ERR_OUT_OF_MEMORY.
 * @see rmat_vstack()
 */
rabin_err_t rmat_hstack(rmat_t* R, const rmat_t* A, const rmat_t* B);

/**
 * @brief Vertical concatenation: \f$R = \begin{bmatrix} A \\ B
 * \end{bmatrix}\f$.
 *
 * The matrices must have the same number of columns; \f$R\f$ must
 * already be initialized as \f$(A\text{.rows} + B\text{.rows}) \times
 * A\text{.cols}\f$.
 *
 * Complexity:
 *   - Time: \f$O((r_1 + r_2) \cdot c \cdot d)\f$
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O((r_1 + r_2) \cdot c)\f$ bignums
 *
 * @param[out] R Result matrix.
 * @param[in]  A Top block.
 * @param[in]  B Bottom block.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR,
 * RABIN_ERR_MATRIX_DIM, RABIN_ERR_OVERFLOW, or RABIN_ERR_OUT_OF_MEMORY.
 * @see rmat_hstack()
 */
rabin_err_t rmat_vstack(rmat_t* R, const rmat_t* A, const rmat_t* B);

/**
 * @brief Integer matrix power: \f$R = A^e\f$ by binary exponentiation.
 *
 * \f$A\f$ must be square and \f$R\f$ must be initialized with the same
 * dimensions; \f$e = 0\f$ produces the identity. The result is computed
 * in a temporary and copied out, so \f$R\f$ may alias \f$A\f$.
 *
 * Complexity:
 *   - Time: \f$O(n^3 \cdot d^2 \cdot \log e)\f$
 *   - Auxiliary memory: \f$O(n^2)\f$ bignums
 *   - Output memory: \f$O(n^2)\f$ bignums
 *
 * @param[out] R Result matrix (\f$A^0 = I\f$).
 * @param[in]  A Square base matrix.
 * @param[in]  e Exponent.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR,
 * RABIN_ERR_MATRIX_DIM, or RABIN_ERR_OUT_OF_MEMORY.
 * @see rmat_mul(), rmat_id()
 */
rabin_err_t rmat_pow(rmat_t* R, const rmat_t* A, u64 e);

/**
 * @brief Transpose of a matrix: \f$T = A^t\f$.
 *
 * \f$T\f$ must already be initialized as \f$c \times r\f$ where
 * \f$r \times c\f$ are the dimensions of \f$A\f$; \f$T\f$ may not alias
 * \f$A\f$.
 *
 * Complexity:
 *   - Time: \f$O(r \cdot c \cdot d)\f$ where \f$d =\f$ the entry size in limbs
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(r \cdot c)\f$ bignums
 *
 * @param[out] T Result matrix (\f$c \times r\f$).
 * @param[in]  A Source matrix (\f$r \times c\f$).
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR,
 * RABIN_ERR_MATRIX_DIM, RABIN_ERR_OVERFLOW, or RABIN_ERR_OUT_OF_MEMORY.
 * @see rmat_copy()
 */
rabin_err_t rmat_transpose(rmat_t* T, const rmat_t* A);

/**
 * @brief Exact inverse of a unimodular integer matrix: \f$B = A^{-1}\f$.
 *
 * Verifies \f$\det(A) = \pm 1\f$ with rmat_det(), then runs fraction-free
 * Gauss-Jordan elimination on the augmented system \f$[A | I]\f$ (every
 * intermediate value stays integral) and normalizes each row by its pivot,
 * checking the product against the identity before storing it. \f$B\f$ must
 * already be initialized with the same dimensions as \f$A\f$ and is left
 * untouched when the inverse does not exist or is not integral.
 *
 * Complexity:
 *   - Time: \f$O(n^3 \cdot d^2)\f$ where \f$n \times n\f$ is the dimension
 * of \f$A\f$ and \f$d =\f$ the entry size in limbs
 *   - Auxiliary memory: \f$O(n^2)\f$ bignums (augmented system)
 *   - Output memory: \f$O(n^2)\f$ bignums
 *
 * @param[out] B Result inverse matrix.
 * @param[in]  A Square integer matrix with \f$\det(A) = \pm 1\f$.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR,
 * RABIN_ERR_MATRIX_DIM, RABIN_ERR_OVERFLOW, RABIN_ERR_OUT_OF_MEMORY, or
 * RABIN_ERR_INVALID_ARG (not unimodular or verification failed).
 * @see rmat_det(), rmat_mul()
 */
rabin_err_t rmat_inv_unimodular(rmat_t* B, const rmat_t* A);

/**
 * @brief Rank of an integer matrix over \f$\mathbb{Q}\f$.
 *
 * Runs exact division-free elimination (Euclid descent pivoting) on a
 * working copy of \f$A\f$; the rank is the number of nonzero pivots.
 *
 * Complexity:
 *   - Time: \f$O(m \cdot n \cdot \min(m,n) \cdot d^2)\f$
 *   - Auxiliary memory: \f$O(m \cdot n)\f$ bignums (working copy)
 *   - Output memory: \f$O(1)\f$
 *
 * @param[out] rank Result rank.
 * @param[in]  A Matrix (\f$m \times n\f$).
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR,
 * or RABIN_ERR_OUT_OF_MEMORY.
 * @par Algorithm Reference:
 * E. G. Bareiss, "Recursive determinants", Numer. Math. 13 (1968);
 * H. Cohen, "A Course in Computational Algebraic Number Theory",
 * Alg. 2.3.1 (exact elimination).
 * @see rmat_rref(), rmat_det_bareiss()
 */
rabin_err_t rmat_rank(u64* rank, const rmat_t* A);

/**
 * @brief Reduced row echelon form over \f$\mathbb{Q}\f$: exact RREF of
 * the integer matrix \f$A\f$ stored in \f$R\f$.
 *
 * Gauss-Jordan elimination on an exact rational copy of \f$A\f$ with
 * partial pivoting on nonzero entries; every pivot is scaled to
 * \f$1\f$ and its column cleared in all other rows. \f$R\f$ must
 * already be initialized with the dimensions of \f$A\f$ and is
 * normalized to a minimal common denominator on output. The rank
 * out-param is optional (\f$NULL\f$ allowed).
 *
 * Complexity:
 *   - Time: \f$O(m \cdot n \cdot \min(m,n) \cdot M(d))\f$ with \f$M(d)\f$
 * the rational multiply cost
 *   - Auxiliary memory: \f$O(m \cdot n)\f$ rationals
 *   - Output memory: \f$O(m \cdot n)\f$ bignums
 *
 * @param[out] R Result matrix in RREF (rational).
 * @param[out] rank Optional result rank (may be NULL).
 * @param[in]  A Integer matrix (\f$m \times n\f$).
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR,
 * RABIN_ERR_MATRIX_DIM, RABIN_ERR_OVERFLOW, or RABIN_ERR_OUT_OF_MEMORY.
 * @par Algorithm Reference:
 * H. Cohen, "A Course in Computational Algebraic Number Theory",
 * Alg. 2.7 (Gauss-Jordan elimination, exact).
 * @see rmat_rank(), rmat_solve()
 */
rabin_err_t rmat_rref(rmatq_t* R, u64* rank, const rmat_t* A);

/**
 * @brief Exact integer LU decomposition with row pivoting.
 *
 * Division-free elimination (Euclid descent pivoting) on \f$A\f$
 * produces a permutation \f$P\f$ (returned through the array
 * \f$P[0..n-1]\f$, where \f$P[i]\f$ is the original row index sitting
 * at position \f$i\f$ after the swaps of step \f$i\f$), a unit lower
 * triangular integer matrix \f$L\f$ and an upper triangular integer
 * matrix \f$U\f$ with
 *
 * \f$(P A) = L \cdot U\f$ exactly.
 *
 * A zero pivot on the diagonal of \f$U\f$ signals that \f$A\f$ is
 * singular; the factors are still produced and satisfy \f$P A = L U\f$
 * exactly, and the function returns RABIN_SUCCESS for singular inputs
 * as well.
 *
 * Complexity:
 *   - Time: \f$O(n^3 \cdot d^2)\f$
 *   - Auxiliary memory: \f$O(n^2)\f$ bignums (augmented system)
 *   - Output memory: \f$O(n^2)\f$ bignums
 *
 * @param[out] L Unit lower triangular factor (caller-allocated \f$n
 * \times n\f$).
 * @param[out] U Upper triangular factor (caller-allocated \f$n \times
 * n\f$); a zero diagonal entry flags singularity.
 * @param[out] P Permutation array of \f$n\f$ u64 slots allocated by the
 * caller.
 * @param[in]  A Square integer matrix.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR,
 * RABIN_ERR_MATRIX_DIM, RABIN_ERR_OVERFLOW, or RABIN_ERR_OUT_OF_MEMORY.
 * @par Algorithm Reference:
 * E. G. Bareiss, "Recursive determinants", Numer. Math. 13 (1968);
 * exact-preserving variant of Gauss elimination with row exchanges.
 * @see rmat_det(), rmat_mul(), rmat_solve()
 */
rabin_err_t rmat_lu(rmat_t* L, rmat_t* U, u64* P, const rmat_t* A);

/**
 * @brief Solve \f$A X = B\f$ over \f$\mathbb{Q}\f$.
 *
 * \f$A\f$ is \f$n \times n\f$, \f$B\f$ is \f$n \times k\f$, and \f$X\f$
 * must already be initialized as \f$n \times k\f$. The system is
 * solved from the RREF of the augmented matrix \f$[A | B]\f$: the
 * solution is unique when the rank of \f$A\f$ is \f$n\f$.
 *
 * Error semantics: RABIN_ERR_INCONSISTENT when the augmented rank
 * exceeds the rank of \f$A\f$ (a zero row of the left block with a
 * nonzero right-hand side), RABIN_ERR_SINGULAR when \f$A\f$ has a free
 * column (rank below \f$n\f$) but the system is consistent.
 *
 * Complexity:
 *   - Time: \f$O(n^2 (n + k))\f$ rational operations
 *   - Auxiliary memory: \f$O(n (n + k))\f$ rationals
 *   - Output memory: \f$O(n k)\f$ bignums
 *
 * @param[out] X Solution matrix (rational, \f$n \times k\f$).
 * @param[in]  A Square coefficient matrix.
 * @param[in]  B Right-hand side (\f$B\text{.rows} = A\text{.rows}\f$).
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR,
 * RABIN_ERR_MATRIX_DIM, RABIN_ERR_OVERFLOW, RABIN_ERR_OUT_OF_MEMORY,
 * RABIN_ERR_SINGULAR, or RABIN_ERR_INCONSISTENT.
 * @see rmat_rref()
 */
rabin_err_t rmat_solve(rmatq_t* X, const rmat_t* A, const rmat_t* B);

/**
 * @brief Cofactor matrix: \f$C[i][j] = (-1)^{i+j} M_{ij}\f$.
 *
 * Each entry is the determinant of the \f$(n-1) \times (n-1)\f$ minor
 * obtained by deleting row \f$i\f$ and column \f$j\f$, computed with
 * rmat_det_bareiss(). For \f$n = 1\f$ the cofactor is \f$1\f$.
 *
 * Complexity:
 *   - Time: \f$O(n^2)\f$ determinants of size \f$n-1\f$
 *   - Auxiliary memory: \f$O(n^2)\f$ bignums (minor)
 *   - Output memory: \f$O(n^2)\f$ bignums
 *
 * @param[out] C Cofactor matrix (caller-allocated like \f$A\f$).
 * @param[in]  A Square integer matrix.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR,
 * RABIN_ERR_MATRIX_DIM, or RABIN_ERR_OUT_OF_MEMORY.
 * @see rmat_adj(), rmat_det_bareiss()
 */
rabin_err_t rmat_cofactor(rmat_t* C, const rmat_t* A);

/**
 * @brief Adjugate matrix: \f$\mathrm{adj}(A) = C^t\f$, the transpose of
 * the cofactor matrix.
 *
 * Satisfies \f$A \cdot \mathrm{adj}(A) = \det(A) \cdot I\f$.
 *
 * Complexity:
 *   - Time: \f$O(n^2)\f$ determinants of size \f$n-1\f$ plus one
 * transpose
 *   - Auxiliary memory: \f$O(n^2)\f$ bignums
 *   - Output memory: \f$O(n^2)\f$ bignums
 *
 * @param[out] ADJ Adjugate matrix (caller-allocated \f$n \times n\f$).
 * @param[in]  A Square integer matrix.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR,
 * RABIN_ERR_MATRIX_DIM, or RABIN_ERR_OUT_OF_MEMORY.
 * @see rmat_cofactor(), rmat_det()
 */
rabin_err_t rmat_adj(rmat_t* ADJ, const rmat_t* A);

/**
 * @brief Singularity test: \f$\det(A) = 0\f$.
 *
 * Returns false for NULL and for non-square inputs.
 *
 * Complexity:
 *   - Time: as rmat_det()
 *   - Auxiliary memory: as rmat_det()
 *   - Output memory: \f$O(1)\f$
 *
 * @param[in] A Matrix to test.
 * @return true if \f$A\f$ is square and singular, false otherwise.
 * @see rmat_det()
 */
bool rmat_is_singular(const rmat_t* A);

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
 *
 * @par Algorithm Reference:
 * V. P. Danilevsky, "On a Method of Calculating the Eigenvectors of a
 * Matrix," Uch. Zap. LGU, ser. Mat. Fiz. Khim. Nauk, no. 17, 1966.
 */
rabin_err_t rmat_charpoly_adj(rpol_t* p, rmat_t* J, const rmat_t* A);

/**
 * @brief Compute the LLL-reduced basis of the lattice generated by the
 * columns of B (full rank, all-integer arithmetic).
 *
 * Let \f$m =\f$ B->rows, \f$n =\f$ B->cols. The columns of B are the
 * generating vectors of a lattice in \f$\mathbb{Z}^m\f$. On success the
 * columns of B are replaced by the LLL-reduced basis and H (re-
 * initialized to \f$n \times n\f$) holds the unimodular transformation:
 * \f$B_{\text{new}} = B_{\text{old}} \cdot H\f$.
 *
 * The Lovasz condition is tested with exact integer arithmetic for every
 * \f$\delta \in (0,1]\f$: any double \f$\delta\f$ is exactly a binary
 * rational \f$\delta = num / 2^e\f$, so the comparison becomes the
 * integer inequality \f$2^e d_{k+1} d_{k-1} < num \cdot d_k^2 -
 * 2^e \lambda^2\f$ and the reduced output is independent of any
 * floating-point rounding.
 *
 * @pre The columns of B must be linearly independent, otherwise
 * RABIN_ERR_INVALID_ARG is returned (see rmat_mlll() for the general
 * rank-deficient case); on failure B is left unchanged.
 *
 * Complexity:
 *   - Time: worst case \f$O(n^2)\f$ RED/SWAP steps, each
 *     \f$O(n \cdot d)\f$ bignum ops where \f$d =\f$ B->rows
 *   - Auxiliary memory: \f$O(n^2)\f$ bignums (H, \f$\lambda\f$) plus
 *     \f$O(n \cdot d)\f$ bignums (vector work copy)
 *   - Output memory: \f$O(n \cdot d)\f$ bignums
 *
 * @param[in,out] B Basis matrix, \f$m \times n\f$, columns = vectors.
 * @param[in]     n Number of basis vectors (must equal B->cols).
 * @param[in]     delta LLL parameter, \f$0 < \delta \le 1\f$ (typically
 *                     \f$0.75\f$).
 * @param[out]    H Transformation matrix; re-initialized to \f$n \times
 *                  n\f$ on success and reset to an empty (zero-initialized)
 *                  matrix on failure.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR,
 * RABIN_ERR_INVALID_ARG (invalid delta, \f$n \ne B \text{-}>\text{cols}\f$,
 * or dependent columns), or RABIN_ERR_OUT_OF_MEMORY.
 *
 * @par Algorithm Reference:
 * H. Cohen, "A Course in Computational Algebraic Number Theory,"
 * Springer-Verlag, Berlin, 1993, Algorithm 2.6.7;
 * H. W. Lenstra, H. W. Lenstra, and L. Lovasz, "Factoring Polynomials
 * in Polynomial Time," Mathematics with the IBM 3090 Computer, 1982.
 * @see rmat_mlll()
 */
rabin_err_t rmat_lll(rmat_t* B, u64 n, double delta, rmat_t* H);

/**
 * @brief MLLL: LLL reduction for possibly linearly dependent generating
 * vectors (exact rational arithmetic).
 *
 * Let \f$m =\f$ B->rows, \f$n =\f$ B->cols. The columns of B generate a
 * lattice of rank \f$p \le n\f$ in \f$\mathbb{Z}^m\f$; zero columns and
 * dependent columns are allowed. On success:
 *
 *   - \f$*rank = p\f$,
 *   - the first \f$n - p\f$ columns of B are zero and the last \f$p\f$
 *     columns form an LLL-reduced basis of the lattice,
 *   - H (re-initialized to \f$n \times n\f$) is unimodular with
 *     \f$B_{\text{new}} = B_{\text{old}} \cdot H\f$, and its first
 *     \f$n - p\f$ columns form a \f$\mathbb{Z}\f$-basis of the relation
 *     space \f$\{x \in \mathbb{Z}^n : B_{\text{old}} \cdot x = 0\}\f$.
 *
 * The Lovasz condition is tested with exact integer arithmetic for every
 * \f$\delta \in (0,1]\f$: any double \f$\delta\f$ is exactly a binary
 * rational \f$\delta = num / 2^e\f$, so the comparison becomes an
 * integer inequality and the reduced output is independent of any
 * floating-point rounding.
 *
 * Complexity:
 *   - Time: one HNF column reduction (Cohen 2.4.5) plus one
 *     \f$O(n^2)\f$-step LLL reduction on the \f$p\f$-dimensional
 *     lattice basis
 *   - Auxiliary memory: \f$O(n^2)\f$ bignums (U, H) plus
 *     \f$O(n \cdot m)\f$ bignums (work copies)
 *   - Output memory: \f$O(n^2 + n \cdot m)\f$ bignums
 *
 * @pre \f$B \text{-}>\text{rows} \le n\f$ (the column-reduction
 * machinery handles wide / square input only, as in Cohen 2.6.8);
 * B is left unchanged when an error is returned.
 *
 * @param[in,out] B Generator matrix, \f$m \times n\f$, columns =
 *                  vectors.
 * @param[in]     n Number of generating vectors (must equal B->cols).
 * @param[in]     delta LLL parameter, \f$0 < \delta \le 1\f$ (typically
 *                     \f$0.75\f$).
 * @param[out]    H Transformation matrix; re-initialized to \f$n \times
 *                  n\f$ on success.
 * @param[out]    rank Lattice rank \f$p\f$.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR,
 * RABIN_ERR_INVALID_ARG (invalid delta, \f$n \ne B \text{-}>\text{cols}\f$,
 * or \f$B \text{-}>\text{rows} > n\f$), or RABIN_ERR_OUT_OF_MEMORY.
 *
 * @par Algorithm Reference:
 * H. Cohen, "A Course in Computational Algebraic Number Theory,"
 * Springer-Verlag, Berlin, 1993, Algorithm 2.4.5 (rank-revealing HNF
 * with transformation) and Algorithm 2.6.7 (LLL), in the spirit of
 * Algorithm 2.6.8 (MLLL, after P. Pohst).
 * @see rmat_lll()
 */
rabin_err_t rmat_mlll(rmat_t* B, u64 n, double delta, rmat_t* H, u64* rank);

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
 *
 * @par Algorithm Reference:
 * H. Cohen, "A Course in Computational Algebraic Number Theory,"
 * Springer-Verlag, Berlin, 1993, Algorithm 2.4.4.
 * @see rmat_hermite_gcd(), rmat_hermite_mod_d()
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
 *
 * @par Algorithm Reference:
 * H. Cohen, "A Course in Computational Algebraic Number Theory,"
 * Springer-Verlag, Berlin, 1993, Algorithm 2.4.5.
 * @see rmat_hermite(), rmat_hermite_mod_d()
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
 *
 * @par Algorithm Reference:
 * H. Cohen, "A Course in Computational Algebraic Number Theory,"
 * Springer-Verlag, Berlin, 1993, Algorithm 2.4.8 (based on the method
 * of Domich, Klee, and Wets).
 * @see rmat_hermite(), rmat_hermite_gcd()
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

/**
 * @brief Integer kernel of a matrix: an LLL-reduced basis of
 * \f$\ker_{\mathbb{Z}}(A)\f$.
 *
 * Computes a Z-basis of the kernel \f$\{x \in \mathbb{Z}^n : A x = 0\}\f$
 * from the transformation matrix of the Hermite normal form of \f$A\f$,
 * then LLL-reduces it with \f$\delta = 0.75\f$.
 *
 * @param[in] A \f$m \times n\f$ integer matrix.
 * @param[out] K \f$n \times r\f$ matrix whose columns are an LLL-reduced
 * Z-basis of the kernel (\f$r = n - \mathrm{rank}(A)\f$).
 *
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR if A or K is
 * NULL, or a memory allocation error.
 *
 * @pre K must be uninitialized or previously cleared.
 *
 * @par Algorithm Reference:
 * Cohen, A Course in Computational Algebraic Number Theory, Alg. 2.7.3
 * (kernel via HNF transform) followed by rmat_lll().
 *
 * @see rmat_kernel_lll(), rmat_kernel_image(), rmat_lll()
 */
rabin_err_t rmat_kernel(const rmat_t* A, rmat_t* K);

/**
 * @brief Integer kernel of a matrix via standalone kernel-LLL (Cohen
 * 2.7.2).
 *
 * Same output contract as rmat_kernel(): an LLL-reduced Z-basis of
 * \f$\ker_{\mathbb{Z}}(A)\f$ of shape \f$n \times r\f$. Runs a modified
 * LLL in which Gram-Schmidt-zero vectors are pushed to the front by a
 * special swap, so the kernel basis is produced directly without first
 * computing an HNF.
 *
 * @param[in] A \f$m \times n\f$ integer matrix.
 * @param[out] K \f$n \times r\f$ matrix whose columns are an LLL-reduced
 * Z-basis of the kernel (\f$r = n - \mathrm{rank}(A)\f$).
 *
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR if A or K is
 * NULL, or a memory allocation error.
 *
 * @pre K must be uninitialized or previously cleared.
 *
 * @par Algorithm Reference:
 * Cohen, A Course in Computational Algebraic Number Theory, Alg. 2.7.2.
 *
 * @see rmat_kernel()
 */
rabin_err_t rmat_kernel_lll(const rmat_t* A, rmat_t* K);

/**
 * @brief Kernel and image of a matrix via an LLL-reduced change of basis.
 *
 * Produces a unimodular \f$n \times n\f$ matrix \f$H\f$ such that:
 *   - the first \f$r = n - p\f$ columns are an LLL-reduced Z-basis of
 *     \f$\ker_{\mathbb{Z}}(A)\f$,
 *   - the last \f$p = \mathrm{rank}(A)\f$ columns give a basis of a
 *     complement on which \f$A\f$ acts as an LLL-reduced image lattice,
 *     size-reduced against the kernel columns.
 *
 * @param[in] A \f$m \times n\f$ integer matrix.
 * @param[out] H Unimodular \f$n \times n\f$ matrix with the property above.
 * @param[out] p Rank of \f$A\f$ (number of image columns of H).
 *
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR if any argument
 * is NULL, or RABIN_ERR_INVALID_ARG if an internal check fails, or a
 * memory allocation error.
 *
 * @pre H must be uninitialized or previously cleared.
 *
 * @par Algorithm Reference:
 * Cohen, A Course in Computational Algebraic Number Theory, Alg. 2.7.1.
 *
 * @see rmat_kernel(), rmat_lll()
 */
rabin_err_t rmat_kernel_image(const rmat_t* A, rmat_t* H, u64* p);

/**
 * @brief Search for integer linear dependence relations among real (or
 * complex) numbers.
 *
 * Builds the lattice spanned by the rows
 * \f$(e_1, \dots, e_n, \alpha \lambda)\f$ (and \f$(0,\dots,0,\beta\lambda)\f$
 * for complex input), where \f$\lambda\f$ is `scale`, then LLL-reduces it.
 * Each relation (a short lattice vector with zero identity part) yields an
 * integer vector \f$(c_1,\dots,c_n)\f$ with \f$\sum_j c_j z_j \approx 0\f$.
 * A relation is a genuine dependence when the corresponding residual
 * `err[c]` is small relative to \f$\lambda^{-1}\f$ and the coefficients.
 *
 * @param[in] re Real parts of the \f$n \ge 2\f$ numbers \f$z_j\f$.
 * @param[in] im Imaginary parts, or NULL for real input.
 * @param[in] n Number of numbers.
 * @param[in] scale Approximation factor \f$\lambda > 0\f$.
 * @param[out] REL \f$m \times n\f$ matrix (\f$m = n+1\f$ for real input,
 * \f$m = n+2\f$ for complex input) whose columns are the LLL-reduced
 * lattice basis; a column with zero entries in its last one/two
 * coordinates is a genuine relation with coefficient vector
 * \f$(c_1,\dots,c_n)\f$.
 * @param[out] err Optional array of \f$n\f$ residuals
 * \f$|\sum_j c_j^{(c)} z_j|\f$ for each column of REL, or NULL.
 *
 * @return RABIN_SUCCESS on success, RABIN_ERR_NULL_PTR if re or REL is
 * NULL, RABIN_ERR_INVALID_ARG if n < 2 or the input/scale is not finite or
 * a coefficient would overflow, or a memory allocation error.
 *
 * @pre REL must be uninitialized or previously cleared.
 *
 * @par Algorithm Reference:
 * Cohen, A Course in Computational Algebraic Number Theory, Alg. 2.7.4.
 *
 * @see rmat_algebraic_relation(), rmat_lll()
 */
rabin_err_t rmat_linear_dependence(const double* re, const double* im, u64 n,
                                   double scale, rmat_t* REL, double* err);

/**
 * @brief Find an integer polynomial relation for a real/complex number.
 *
 * Builds the powers \f$z_i = a^i\f$, \f$i = 0,\dots,\mathrm{degree}\f$, of
 * \f$a = \mathrm{re} + i\,\mathrm{im}\f$ and searches for an integer linear
 * dependence among them with rmat_linear_dependence(). A short column
 * \f$c_0,\dots,c_d\f$ of the output gives a candidate polynomial
 * \f$\sum_i c_i x^i\f$ with \f$\sum_i c_i a^i \approx 0\f$; err is the
 * corresponding residual.
 *
 * @param[in] re Real part of \f$a\f$.
 * @param[in] im Imaginary part of \f$a\f$ (zero for real \f$a\f$).
 * @param[in] degree Maximum degree \f$d \le 64\f$ of the relation.
 * @param[in] scale Approximation factor \f$\lambda > 0\f$.
 * @param[out] poly \f$m \times (\mathrm{degree}+1)\f$ matrix
 * (\f$m = \mathrm{degree}+2\f$ for real \f$a\f$, \f$m = \mathrm{degree}+3\f$
 * for complex \f$a\f$); each column is a candidate integer polynomial in
 * natural order (\f$c_0\f$ is the constant term), followed by one/two
 * trailing coordinates that are zero for a genuine relation.
 * @param[out] err Optional array of residuals \f$|\sum_i c_i^{(c)} a^i|\f$.
 *
 * @return RABIN_SUCCESS on success, RABIN_ERR_NULL_PTR if poly is NULL,
 * RABIN_ERR_INVALID_ARG if degree is out of range, the input/scale is not
 * finite, a power of \f$a\f$ overflows, or a memory allocation error.
 *
 * @pre poly must be uninitialized or previously cleared.
 *
 * @par Algorithm Reference:
 * Cohen, A Course in Computational Algebraic Number Theory, Alg. 2.7.5.
 *
 * @see rmat_linear_dependence()
 */
rabin_err_t rmat_algebraic_relation(double re, double im, u64 degree,
                                    double scale, rmat_t* poly, double* err);

/**
 * @brief Number of bytes rmat_to_string() needs for \p A, NUL included.
 *
 * Complexity:
 *   - Time: \f$O(rc \cdot n^2)\f$ decimal conversion, see rz_to_string()
 *   - Auxiliary memory: \f$O(n)\f$ limbs
 *
 * @param[in] A Matrix to measure.
 * @return Required buffer size including the NUL terminator, or 0 if
 * \p A is NULL or conversion fails.
 * @see rmat_to_string()
 */
size_t rmat_to_string_len(const rmat_t* A);

/**
 * @brief Write a matrix as nested bracketed decimal lists
 * \f$[[r_{0,0},\dots],[r_{1,0},\dots],\dots]\f$ into a caller buffer.
 *
 * The output is exactly NUL-terminated and needs rmat_to_string_len()
 * bytes including the terminator; a smaller \p cap leaves \p out
 * untouched and fails.
 *
 * Complexity:
 *   - Time: \f$O(rc \cdot n^2)\f$ decimal conversion, see rz_to_string()
 *   - Auxiliary memory: \f$O(n)\f$ limbs
 *   - Output memory: \f$O(rc)\f$ characters
 *
 * @param[out] out Caller buffer receiving the string.
 * @param[in]  A   Matrix to print.
 * @param[in]  cap Size of \p out in bytes, NUL included.
 * @return RABIN_SUCCESS on success; RABIN_ERR_NULL_PTR if \p out or
 * \p A is NULL; RABIN_ERR_INVALID_ARG if \p cap is too small;
 * RABIN_ERR_OUT_OF_MEMORY on allocation failure.
 * @see rmat_init_val()
 */
rabin_err_t rmat_to_string(const rmat_t* A, char* out, size_t cap);

/**
 * @brief Parse nested bracketed decimal lists into an existing matrix.
 *
 * Elements are signed decimal integers; whitespace is allowed between
 * tokens. \p A must be initialized and the parsed row and column
 * counts must match A->rows and A->cols (RABIN_ERR_MATRIX_DIM
 * otherwise); a failed parse leaves \p A unchanged. The result
 * round-trips exactly with rmat_to_string().
 *
 * Complexity:
 *   - Time: \f$O(rc \cdot n^2)\f$ where \f$n =\f$ the digit count
 *   - Auxiliary memory: \f$O(rc \cdot n)\f$ limbs
 *   - Output memory: \f$O(rc \cdot n)\f$ limbs
 *
 * @param[out] A   Matrix to fill; must be initialized.
 * @param[in]  str String to parse.
 * @return RABIN_SUCCESS on success; RABIN_ERR_NULL_PTR if \p A or
 * \p str is NULL; RABIN_ERR_INVALID_ARG if \p str is malformed;
 * RABIN_ERR_MATRIX_DIM if the parsed shape differs;
 * RABIN_ERR_OUT_OF_MEMORY on allocation failure.
 * @see rmat_to_string()
 */
rabin_err_t rmat_init_val(rmat_t* A, const char* str);

#endif
