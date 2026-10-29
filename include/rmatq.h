#ifndef RMATQ_H
#define RMATQ_H

#include "rabin.h"
#include "rabin_errors.h"
#include "rq.h"

/**
 * @brief A rational matrix with a single common denominator.
 *
 * The value is \f$M[i][j] / den\f$ for every entry. \f$den = 1\f$ for
 * integer matrices; rmatq_normalize() removes any common divisor shared by
 * \f$den\f$ and all entries.
 */
typedef struct rmatq_t {
  rmat_t M;
  rz_t den;
} rmatq_t;

/**
 * @brief Initialize a \f$m \times n\f$ rational matrix of zeros.
 *
 * Sets \f$den = 1\f$.
 *
 * Complexity:
 *   - Time: \f$O(m \cdot n)\f$
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(m \cdot n)\f$ bignums
 *
 * @param[out] A Matrix to initialize.
 * @param[in]  m Number of rows.
 * @param[in]  n Number of columns.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR,
 * RABIN_ERR_OVERFLOW, or RABIN_ERR_OUT_OF_MEMORY.
 */
rabin_err_t rmatq_init(rmatq_t* A, u64 m, u64 n);

/**
 * @brief Free a rational matrix.
 *
 * Complexity:
 *   - Time: \f$O(m \cdot n)\f$
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(1)\f$
 *
 * @param[in,out] A Matrix to free.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR.
 */
rabin_err_t rmatq_clear(rmatq_t* A);

/**
 * @brief Deep-copy one rational matrix into another.
 *
 * R must already be initialized with the same dimensions as A.
 *
 * Complexity:
 *   - Time: \f$O(m \cdot n \cdot d)\f$ where \f$d =\f$ the entry size in
 * limbs
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(m \cdot n)\f$ bignums
 *
 * @param[out] R Destination.
 * @param[in]  A Source.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR /
 * RABIN_ERR_OUT_OF_MEMORY.
 */
rabin_err_t rmatq_copy(rmatq_t* R, const rmatq_t* A);

/**
 * @brief Deep-copy an integer matrix, setting the denominator to 1.
 *
 * R must already be initialized with the same dimensions as A.
 *
 * Complexity:
 *   - Time: \f$O(m \cdot n \cdot d)\f$ where \f$d =\f$ the entry size in
 * limbs
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(m \cdot n)\f$ bignums
 *
 * @param[out] R Destination.
 * @param[in]  A Integer source matrix.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR /
 * RABIN_ERR_OUT_OF_MEMORY.
 */
rabin_err_t rmatq_copyZ(rmatq_t* R, rmat_t* A);

/**
 * @brief Sum of two rational matrices: \f$R = A + B\f$.
 *
 * All entries share the result denominator
 * \f$R.den = \mathrm{lcm}(A.den, B.den)\f$; the operand matrices are
 * scaled to that denominator before the element-wise addition. R must
 * already be initialized with the same dimensions as A and B.
 *
 * Complexity:
 *   - Time: \f$O(m \cdot n \cdot d^2)\f$ where \f$d =\f$ the entry size in
 * limbs
 *   - Auxiliary memory: \f$O(m \cdot n)\f$ bignums (scaled temporaries)
 *   - Output memory: \f$O(m \cdot n)\f$ bignums
 *
 * @param[out] R Result.
 * @param[in]  A First operand.
 * @param[in]  B Second operand.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR,
 * RABIN_ERR_MATRIX_DIM, or RABIN_ERR_OUT_OF_MEMORY.
 */
rabin_err_t rmatq_add(rmatq_t* R, const rmatq_t* A, const rmatq_t* B);

/**
 * @brief Difference of two rational matrices: \f$R = A - B\f$.
 *
 * All entries share the result denominator
 * \f$R.den = \mathrm{lcm}(A.den, B.den)\f$; the operand matrices are
 * scaled to that denominator before the element-wise subtraction. R must
 * already be initialized with the same dimensions as A and B.
 *
 * Complexity:
 *   - Time: \f$O(m \cdot n \cdot d^2)\f$ where \f$d =\f$ the entry size in
 * limbs
 *   - Auxiliary memory: \f$O(m \cdot n)\f$ bignums (scaled temporaries)
 *   - Output memory: \f$O(m \cdot n)\f$ bignums
 *
 * @param[out] R Result.
 * @param[in]  A First operand.
 * @param[in]  B Second operand.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR,
 * RABIN_ERR_MATRIX_DIM, or RABIN_ERR_OUT_OF_MEMORY.
 */
rabin_err_t rmatq_sub(rmatq_t* R, const rmatq_t* A, const rmatq_t* B);

/**
 * @brief Reduce a rational matrix in place.
 *
 * Divides every entry and the common denominator by
 * \f$\gcd(den, M[0][0], M[0][1], \ldots)\f$, removing any common factor
 * shared by the denominator and all entries.
 *
 * Complexity:
 *   - Time: \f$O(m \cdot n \cdot d^2)\f$ where \f$d =\f$ the entry size in
 * limbs
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(1)\f$
 *
 * @param[in,out] M Matrix to reduce.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR /
 * RABIN_ERR_OUT_OF_MEMORY.
 */
rabin_err_t rmatq_normalize(rmatq_t* M);

/**
 * @brief Print a rational matrix to stdout, one row per line.
 *
 * Every entry is printed in the unreduced form \f$M[i][j] / den\f$.
 *
 * Complexity:
 *   - Time: \f$O(m \cdot n \cdot d)\f$ where \f$d =\f$ the entry size in
 * limbs
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(m \cdot n \cdot d)\f$ characters written
 *
 * @param[in] A Matrix to print.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR.
 */
rabin_err_t rmatq_print(const rmatq_t* A);

/**
 * @brief Read one entry as a rational: \f$a = M[r][c] / den\f$.
 *
 * The result is normalized (\f$den > 0\f$, coprime).
 *
 * Complexity:
 *   - Time: \f$O(d^2)\f$ where \f$d =\f$ the entry size in limbs
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(d)\f$ limbs
 *
 * @param[out] a Entry value.
 * @param[in]  A Matrix.
 * @param[in]  r Row index.
 * @param[in]  c Column index.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR,
 * RABIN_ERR_INVALID_ARG, or RABIN_ERR_OUT_OF_MEMORY.
 * @see rmatq_set()
 */
rabin_err_t rmatq_get(rq_t* a, const rmatq_t* A, u64 r, u64 c);

/**
 * @brief Write one entry: \f$A[r][c] = a\f$.
 *
 * If \f$a.den\f$ does not divide the common denominator, every entry of
 * \f$A.M\f$ and \f$A.den\f$ are first scaled by
 * \f$\mathrm{lcm}(A.den, a.den) / A.den\f$, so the exact value of \f$a\f$
 * is representable. \f$A\f$ is normalized afterwards.
 *
 * Complexity:
 *   - Time: \f$O(m \cdot n \cdot d^2)\f$ where \f$m \times n\f$ is the
 * dimension of \f$A\f$ and \f$d =\f$ the entry size in limbs
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(1)\f$
 *
 * @param[in,out] A Matrix.
 * @param[in]     a Entry value.
 * @param[in]     r Row index.
 * @param[in]     c Column index.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR,
 * RABIN_ERR_INVALID_ARG, or RABIN_ERR_OUT_OF_MEMORY.
 * @see rmatq_get()
 */
rabin_err_t rmatq_set(rmatq_t* A, const rq_t* a, u64 r, u64 c);

/**
 * @brief Set one entry to zero.
 *
 * The common denominator is left unchanged.
 *
 * Complexity:
 *   - Time: \f$O(1)\f$
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(1)\f$
 *
 * @param[in,out] A Matrix.
 * @param[in]     r Row index.
 * @param[in]     c Column index.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR,
 * RABIN_ERR_INVALID_ARG, or RABIN_ERR_OUT_OF_MEMORY.
 * @see rmatq_set()
 */
rabin_err_t rmatq_set_zero(rmatq_t* A, u64 r, u64 c);

/**
 * @brief Set all entries of a rational matrix to zero, in place.
 *
 * The common denominator is reset to \f$1\f$.
 *
 * Complexity:
 *   - Time: \f$O(r \cdot c)\f$
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(1)\f$
 *
 * @param[in,out] A Matrix to clear in place.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR.
 * @see rmatq_ones(), rmatq_id()
 */
rabin_err_t rmatq_zero(rmatq_t* A);

/**
 * @brief Set all entries of a rational matrix to one, in place.
 *
 * The common denominator is reset to \f$1\f$.
 *
 * Complexity:
 *   - Time: \f$O(r \cdot c)\f$
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(1)\f$
 *
 * @param[in,out] A Matrix to set in place.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR,
 * or RABIN_ERR_OUT_OF_MEMORY.
 * @see rmatq_zero()
 */
rabin_err_t rmatq_ones(rmatq_t* A);

/**
 * @brief Build the \f$n \times n\f$ rational identity matrix, in place.
 *
 * \f$I\f$ must already be initialized as \f$n \times n\f$ and every
 * entry is overwritten; the common denominator is reset to \f$1\f$.
 *
 * Complexity:
 *   - Time: \f$O(n^2)\f$
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(1)\f$
 *
 * @param[out] I Identity matrix (\f$n \times n\f$, pre-initialized).
 * @param[in]  n Dimension.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR,
 * RABIN_ERR_MATRIX_DIM, or RABIN_ERR_OUT_OF_MEMORY.
 * @see rmat_id(), rmatq_diag()
 */
rabin_err_t rmatq_id(rmatq_t* I, u64 n);

/**
 * @brief Build the diagonal matrix \f$D = \mathrm{diag}(v)\f$ over
 * \f$\mathbb{Q}\f$.
 *
 * \f$D\f$ must already be initialized as \f$|v| \times |v|\f$; the
 * integer entries of \f$v\f$ are stored on common denominator
 * \f$1\f$.
 *
 * Complexity:
 *   - Time: \f$O(n^2 \cdot d)\f$ where \f$n = |v|\f$
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(n^2)\f$ bignums
 *
 * @param[out] D Diagonal result (\f$|v| \times |v|\f$, pre-initialized).
 * @param[in]  v Vector of integer diagonal entries.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR,
 * RABIN_ERR_MATRIX_DIM, or RABIN_ERR_OUT_OF_MEMORY.
 * @see rmat_diag()
 */
rabin_err_t rmatq_diag(rmatq_t* D, const rvec_t* v);

/**
 * @brief Swap rows \f$i\f$ and \f$j\f$, in place.
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
 * @see rmatq_swap_cols()
 */
rabin_err_t rmatq_swap_rows(rmatq_t* A, u64 i, u64 j);

/**
 * @brief Swap columns \f$i\f$ and \f$j\f$, in place.
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
 * @see rmatq_swap_rows()
 */
rabin_err_t rmatq_swap_cols(rmatq_t* A, u64 i, u64 j);

/**
 * @brief Transpose: \f$T = A^t\f$.
 *
 * \f$T\f$ must already be initialized with the transposed dimensions;
 * \f$T\f$ may alias \f$A\f$. The common denominator is copied.
 *
 * Complexity:
 *   - Time: \f$O(m \cdot n \cdot d)\f$
 *   - Auxiliary memory: \f$O(n^2)\f$ only when \f$T = A\f$
 *   - Output memory: \f$O(m \cdot n)\f$ bignums
 *
 * @param[out] T Result matrix.
 * @param[in]  A Source matrix.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR,
 * RABIN_ERR_MATRIX_DIM, or RABIN_ERR_OUT_OF_MEMORY.
 * @see rmat_transpose()
 */
rabin_err_t rmatq_transpose(rmatq_t* T, const rmatq_t* A);

/**
 * @brief Horizontal concatenation: \f$R = [A \,|\, B]\f$.
 *
 * \f$A\f$ and \f$B\f$ must have the same number of rows; \f$R\f$ must
 * already be initialized as \f$m \times (c_1 + c_2)\f$. The result
 * denominator is a common multiple of the input denominators.
 *
 * Complexity:
 *   - Time: \f$O(r \cdot (c_1 + c_2) \cdot M(d))\f$
 *   - Auxiliary memory: \f$O(1)\f$ rationals
 *   - Output memory: \f$O(r \cdot (c_1 + c_2))\f$ bignums
 *
 * @param[out] R Result matrix.
 * @param[in]  A Left block.
 * @param[in]  B Right block.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR,
 * RABIN_ERR_MATRIX_DIM, RABIN_ERR_OVERFLOW, or RABIN_ERR_OUT_OF_MEMORY.
 * @see rmatq_vstack()
 */
rabin_err_t rmatq_hstack(rmatq_t* R, const rmatq_t* A, const rmatq_t* B);

/**
 * @brief Vertical concatenation: \f$R = \begin{bmatrix} A \\ B
 * \end{bmatrix}\f$.
 *
 * \f$A\f$ and \f$B\f$ must have the same number of columns; \f$R\f$
 * must already be initialized as \f$(r_1 + r_2) \times c\f$. The result
 * denominator is a common multiple of the input denominators.
 *
 * Complexity:
 *   - Time: \f$O((r_1 + r_2) \cdot c \cdot M(d))\f$
 *   - Auxiliary memory: \f$O(1)\f$ rationals
 *   - Output memory: \f$O((r_1 + r_2) \cdot c)\f$ bignums
 *
 * @param[out] R Result matrix.
 * @param[in]  A Top block.
 * @param[in]  B Bottom block.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR,
 * RABIN_ERR_MATRIX_DIM, RABIN_ERR_OVERFLOW, or RABIN_ERR_OUT_OF_MEMORY.
 * @see rmatq_hstack()
 */
rabin_err_t rmatq_vstack(rmatq_t* R, const rmatq_t* A, const rmatq_t* B);

/**
 * @brief Scalar multiplication: \f$R = a \cdot A\f$.
 *
 * \f$R\f$ must already be initialized like \f$A\f$ and may alias it.
 *
 * Complexity:
 *   - Time: \f$O(r \cdot c \cdot M(d))\f$
 *   - Auxiliary memory: \f$O(1)\f$ rationals
 *   - Output memory: \f$O(r \cdot c)\f$ bignums
 *
 * @param[out] R Result matrix.
 * @param[in]  A Source matrix.
 * @param[in]  a Scalar factor.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR,
 * RABIN_ERR_MATRIX_DIM, or RABIN_ERR_OUT_OF_MEMORY.
 * @see rmatq_mul()
 */
rabin_err_t rmatq_scalar(rmatq_t* R, const rmatq_t* A, const rq_t* a);

/**
 * @brief Rational matrix power: \f$R = A^e\f$ by binary exponentiation.
 *
 * \f$A\f$ must be square and \f$R\f$ initialized with the same
 * dimensions; \f$e = 0\f$ produces the identity. \f$R\f$ may alias
 * \f$A\f$.
 *
 * Complexity:
 *   - Time: \f$O(n^3 \cdot d^2 \cdot \log e)\f$
 *   - Auxiliary memory: \f$O(n^2)\f$ rationals
 *   - Output memory: \f$O(n^2)\f$ bignums
 *
 * @param[out] R Result matrix (\f$A^0 = I\f$).
 * @param[in]  A Square base matrix.
 * @param[in]  e Exponent.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR,
 * RABIN_ERR_MATRIX_DIM, or RABIN_ERR_OUT_OF_MEMORY.
 * @see rmatq_mul(), rmatq_id()
 */
rabin_err_t rmatq_pow(rmatq_t* R, const rmatq_t* A, u64 e);

/**
 * @brief Equality test, normalize-safe.
 *
 * Compares entries by cross-multiplication
 * \f$A.M[i][j] \cdot B.den = B.M[i][j] \cdot A.den\f$; neither matrix
 * is modified.
 *
 * Complexity:
 *   - Time: \f$O(r \cdot c \cdot M(d))\f$
 *   - Auxiliary memory: \f$O(1)\f$ bignums
 *   - Output memory: \f$O(1)\f$
 *
 * @param[in] A First matrix.
 * @param[in] B Second matrix.
 * @return true if both matrices have equal dimensions and values.
 * @see rmatq_normalize()
 */
bool rmatq_equal(const rmatq_t* A, const rmatq_t* B);

/**
 * @brief Check whether the matrix is square.
 *
 * @param[in] A Matrix to test.
 * @return true if \f$A\f$ is square, false if not or on NULL input.
 * @see rmat_is_square()
 */
bool rmatq_is_square(const rmatq_t* A);

/**
 * @brief Check whether all entries are zero.
 *
 * @param[in] A Matrix to test.
 * @return true if \f$A\f$ is the zero matrix, false if not or on NULL
 * input.
 * @see rmatq_zero()
 */
bool rmatq_is_zero(const rmatq_t* A);

/**
 * @brief Check whether the matrix is the identity.
 *
 * Entries are compared against \f$1\f$ / \f$0\f$ through the common
 * denominator without modifying the matrix.
 *
 * @param[in] A Matrix to test.
 * @return true if \f$A = I\f$, false if not or on NULL / non-square
 * input.
 * @see rmatq_id()
 */
bool rmatq_is_identity(const rmatq_t* A);

/**
 * @brief Check whether \f$A = A^t\f$.
 *
 * @param[in] A Matrix to test.
 * @return true if \f$A\f$ is square and symmetric, false if not or on
 * NULL input.
 * @see rmatq_transpose()
 */
bool rmatq_is_symmetric(const rmatq_t* A);

/**
 * @brief Trace: \f$t = \sum_i A[i][i]\f$.
 *
 * Complexity:
 *   - Time: \f$O(n \cdot M(d))\f$
 *   - Auxiliary memory: \f$O(1)\f$ rationals
 *   - Output memory: \f$O(1)\f$
 *
 * @param[out] t Result trace (caller-allocated).
 * @param[in]  A Square rational matrix.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR,
 * RABIN_ERR_MATRIX_DIM, or RABIN_ERR_OUT_OF_MEMORY.
 * @see rmatq_det()
 */
rabin_err_t rmatq_trace(rq_t* t, const rmatq_t* A);

/**
 * @brief Determinant: \f$d = \det(A)\f$.
 *
 * Uses \f$\det(A) = \det(M) / den^n\f$: the numerator determinant is
 * computed with Bareiss elimination over the integers and rescaled by
 * the \f$n\f$-th power of the common denominator.
 *
 * Complexity:
 *   - Time: as rmat_det_bareiss() on the numerator matrix
 *   - Auxiliary memory: \f$O(n^2)\f$ bignums
 *   - Output memory: \f$O(1)\f$
 *
 * @param[out] d Result determinant (caller-allocated).
 * @param[in]  A Square rational matrix.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR,
 * RABIN_ERR_MATRIX_DIM, or RABIN_ERR_OUT_OF_MEMORY.
 * @par Algorithm Reference:
 * E. G. Bareiss, "Recursive determinants", Numer. Math. 13 (1968).
 * @see rmat_det_bareiss()
 */
rabin_err_t rmatq_det(rq_t* d, const rmatq_t* A);

/**
 * @brief Cholesky decomposition \f$A = C^t C\f$ of a rational matrix.
 *
 * \f$C\f$ is upper triangular and every entry of \f$A\f$ must be
 * representable exactly: each pivot \f$q_{i,i}\f$ of the reduction must be
 * positive and a perfect square in \f$\mathbb{Q}\f$, otherwise the
 * decomposition does not exist over the rationals and the call fails.
 * \f$C\f$ must already be initialized as \f$n \times n\f$ and is left
 * untouched when the decomposition fails.
 *
 * Complexity:
 *   - Time: \f$O(n^3 \cdot d^2)\f$ where \f$n \times n\f$ is the dimension
 * of \f$A\f$ and \f$d =\f$ the entry size in limbs
 *   - Auxiliary memory: \f$O(n^2)\f$ rationals (scratch matrix)
 *   - Output memory: \f$O(n^2)\f$ bignums
 *
 * @param[out] C Upper triangular factor with \f$A = C^t C\f$.
 * @param[in]  A Symmetric positive definite rational matrix.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR,
 * RABIN_ERR_MATRIX_DIM, RABIN_ERR_OVERFLOW, RABIN_ERR_OUT_OF_MEMORY,
 * RABIN_ERR_DIV_BY_ZERO (zero pivot), or RABIN_ERR_INVALID_ARG (negative
 * pivot or a pivot that is not a rational square).
 * @par Algorithm Reference:
 * Algorithm 2.7.6 (Cholesky decomposition), H. Cohen, A Course in
 * Computational Algebraic Number Theory.
 * @see rmatq_get(), rmatq_set()
 */
rabin_err_t rmatq_cholesky(rmatq_t* C, const rmatq_t* A);

/**
 * @brief Product of two rational matrices: \f$R = A \cdot B\f$.
 *
 * \f$R.den = A.den \cdot B.den\f$ before normalization; entries are exact
 * integer combinations \f$\sum_k A.M[i][k] \cdot B.M[k][j]\f$. \f$R\f$ must
 * already be initialized as \f$m \times p\f$ where \f$A\f$ is
 * \f$m \times n\f$ and \f$B\f$ is \f$n \times p\f$; \f$R\f$ may not alias
 * \f$A\f$ or \f$B\f$.
 *
 * Complexity:
 *   - Time: \f$O(m \cdot n \cdot p \cdot d^2)\f$ where \f$d =\f$ the entry
 * size in limbs
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(m \cdot p)\f$ bignums
 *
 * @param[out] R Result matrix.
 * @param[in]  A Left factor.
 * @param[in]  B Right factor.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR,
 * RABIN_ERR_MATRIX_DIM, RABIN_ERR_OVERFLOW, or RABIN_ERR_OUT_OF_MEMORY.
 * @see rmatq_normalize()
 */
rabin_err_t rmatq_mul(rmatq_t* R, const rmatq_t* A, const rmatq_t* B);

/**
 * @brief Exact inverse of a rational matrix: \f$B = A^{-1}\f$.
 *
 * Gauss-Jordan elimination with partial pivoting over \f$\mathbb{Q}\f$;
 * every value stays exact. \f$A\f$ must be square and invertible,
 * otherwise the call fails and \f$B\f$ is left untouched. \f$B\f$ must
 * already be initialized as \f$n \times n\f$; \f$B\f$ may not alias
 * \f$A\f$. The result is verified against the identity before being
 * stored.
 *
 * Complexity:
 *   - Time: \f$O(n^3 \cdot d^2)\f$ where \f$n \times n\f$ is the dimension
 * of \f$A\f$ and \f$d =\f$ the entry size in limbs
 *   - Auxiliary memory: \f$O(n^2)\f$ rationals (augmented system)
 *   - Output memory: \f$O(n^2)\f$ bignums
 *
 * @param[out] B Result inverse matrix.
 * @param[in]  A Square invertible rational matrix.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR,
 * RABIN_ERR_MATRIX_DIM, RABIN_ERR_OVERFLOW, RABIN_ERR_OUT_OF_MEMORY, or
 * RABIN_ERR_DIV_BY_ZERO (singular matrix).
 * @see rmatq_mul()
 */
rabin_err_t rmatq_inv(rmatq_t* B, const rmatq_t* A);

/**
 * @brief Rank of a rational matrix over \f$\mathbb{Q}\f$.
 *
 * The numerator matrix has the same rank as \f$A = M / den\f$, so the
 * rank is computed on the integer numerators.
 *
 * Complexity:
 *   - Time: as rmat_rank() on the numerator matrix
 *   - Auxiliary memory: \f$O(m \cdot n)\f$ bignums (working copy)
 *   - Output memory: \f$O(1)\f$
 *
 * @param[out] rank Result rank.
 * @param[in]  A Rational matrix (\f$m \times n\f$).
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR,
 * or RABIN_ERR_OUT_OF_MEMORY.
 * @see rmat_rank(), rmatq_rref()
 */
rabin_err_t rmatq_rank(u64* rank, const rmatq_t* A);

/**
 * @brief Reduced row echelon form over \f$\mathbb{Q}\f$, in place.
 *
 * Exact Gauss-Jordan elimination: each pivot is scaled to \f$1\f$ and
 * its column is cleared in all other rows; \f$R\f$ is overwritten and
 * must already be initialized with the desired dimensions. The rank
 * out-param is optional (\f$NULL\f$ allowed).
 *
 * Complexity:
 *   - Time: \f$O(m \cdot n \cdot \min(m,n) \cdot M(d))\f$
 *   - Auxiliary memory: \f$O(1)\f$ rationals
 *   - Output memory: in place
 *
 * @param[in,out] R Matrix to reduce; on output it is in RREF.
 * @param[out]    rank Result rank (may be NULL).
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR,
 * RABIN_ERR_DIV_BY_ZERO, or RABIN_ERR_OUT_OF_MEMORY.
 * @par Algorithm Reference:
 * H. Cohen, "A Course in Computational Algebraic Number Theory",
 * Alg. 2.7 (Gauss-Jordan elimination, exact).
 * @see rmat_rref(), rmatq_solve()
 */
rabin_err_t rmatq_rref(rmatq_t* R, u64* rank);

/**
 * @brief LU decomposition with row pivoting over \f$\mathbb{Q}\f$.
 *
 * Produces a permutation array \f$P\f$ (\f$P[i]\f$ is the original row
 * index at final position \f$i\f$), a unit lower triangular \f$L\f$
 * and an upper triangular \f$U\f$ with
 *
 * \f$(P A) = L \cdot U\f$ exactly.
 *
 * A zero pivot leaves \f$U[i][i] = 0\f$ and the corresponding column
 * of \f$L\f$ zero; the factors still multiply to the permuted input,
 * so a zero diagonal entry of \f$U\f$ flags a singular \f$A\f$ and the
 * return code is RABIN_SUCCESS. \f$L\f$, \f$U\f$ must already be
 * initialized as \f$n \times n\f$ and may not alias \f$A\f$
 * simultaneously with each other.
 *
 * Complexity:
 *   - Time: \f$O(n^3 \cdot M(d))\f$
 *   - Auxiliary memory: \f$O(n^2)\f$ rationals
 *   - Output memory: \f$O(n^2)\f$ bignums
 *
 * @param[out] L Unit lower triangular factor.
 * @param[out] U Upper triangular factor.
 * @param[out] P Permutation array of \f$n\f$ u64 slots allocated by
 * the caller.
 * @param[in]  A Square rational matrix.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR,
 * RABIN_ERR_MATRIX_DIM, or RABIN_ERR_OUT_OF_MEMORY.
 * @see rmat_lu(), rmatq_mul()
 */
rabin_err_t rmatq_lu(rmatq_t* L, rmatq_t* U, u64* P, const rmatq_t* A);

/**
 * @brief Solve \f$A X = B\f$ over \f$\mathbb{Q}\f$.
 *
 * \f$A\f$ is \f$n \times n\f$, \f$B\f$ is \f$n \times k\f$, and \f$X\f$
 * must already be initialized as \f$n \times k\f$. The system is
 * solved from the RREF of the augmented matrix \f$[A | B]\f$; the
 * solution is unique when \f$\mathrm{rank}(A) = n\f$.
 *
 * Error semantics: RABIN_ERR_INCONSISTENT when the augmented rank
 * exceeds the rank of \f$A\f$, RABIN_ERR_SINGULAR when \f$A\f$ has a
 * free column but the system is consistent.
 *
 * Complexity:
 *   - Time: \f$O(n^2 (n + k) \cdot M(d))\f$
 *   - Auxiliary memory: \f$O(n (n + k))\f$ rationals
 *   - Output memory: \f$O(n k)\f$ bignums
 *
 * @param[out] X Solution matrix.
 * @param[in]  A Square coefficient matrix.
 * @param[in]  B Right-hand side (\f$B\text{.rows} = A\text{.rows}\f$).
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR,
 * RABIN_ERR_MATRIX_DIM, RABIN_ERR_OUT_OF_MEMORY, RABIN_ERR_SINGULAR,
 * or RABIN_ERR_INCONSISTENT.
 * @see rmatq_rref(), rmat_solve()
 */
rabin_err_t rmatq_solve(rmatq_t* X, const rmatq_t* A, const rmatq_t* B);

/**
 * @brief Shortest non-zero vector of the lattice \f$(\mathbb{Z}^n, Q)\f$.
 *
 * \f$Q(x) = x^t A x\f$ is the positive definite quadratic form defined by
 * the symmetric rational matrix \f$A\f$; the function returns (one of) the
 * non-zero \f$x \in \mathbb{Z}^n\f$ minimizing \f$Q(x)\f$ and the exact
 * value \f$Q(x)\f$. It runs the full Fincke-Pohst pipeline of Cohen
 * 2.7.7: LLL reduction (rmat_lll()) of the dual lattice of
 * \f$(\mathbb{Z}^n, Q)\f$, whose rows are obtained from \f$A^{-1}\f$ via
 * rmat_inv_unimodular() (they coincide with the rows of Cholesky's
 * \f$R^{-1}\f$ but
 * never need square roots), the unimodular \f$U = (H^t)^{-1}\f$,
 * reordering of the columns \f$s_i = Ru_i\f$ of \f$S = RU\f$ by
 * increasing \f$\Vert s_i \Vert^2 = Q(u_i)\f$, the Cholesky coefficients
 * (steps 1-3 of Cohen 2.7.6) of \f$U^t A U\f$, and exact enumeration of
 * all \f$y\f$ with \f$Q'(y) \le \min_i \Vert s_i \Vert^2\f$. Since those
 * are genuine lattice vectors the minimum is attained inside the search
 * box, and all bounds use integer square roots (rz_isqrt()), so the
 * result is exact with no floating point.
 *
 * The number of enumerated nodes is exponential in \f$n\f$; keep
 * \f$n\f$ small.
 *
 * @param[out] x    Shortest vector; must already be initialized with
 *                  rvec_init() to \f$n = A \to M.rows\f$ elements.
 * @param[out] norm Exact norm \f$Q(x)\f$.
 * @param[in]  A    Symmetric positive definite rational Gram matrix.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR,
 * RABIN_ERR_MATRIX_DIM, RABIN_ERR_OVERFLOW, RABIN_ERR_OUT_OF_MEMORY,
 * RABIN_ERR_INVALID_ARG (\f$n = 0\f$, or a Cholesky pivot is not strictly
 * positive), or RABIN_ERR_DIV_BY_ZERO (zero pivot).
 * @par Algorithm Reference:
 * H. Cohen, "A Course in Computational Algebraic Number Theory,"
 * Springer-Verlag, Berlin, 1993, Algorithms 2.7.5 (short vectors),
 * 2.7.6 (Cholesky) and 2.7.7 (Fincke-Pohst); P. Fincke, M. Pohst,
 * "A new procedure for computing local height bounds," Mathematics of
 * Computation, vol. 35, 1980.
 * @see rmatq_cholesky(), rmat_lll()
 */
rabin_err_t rmatq_fincke_pohst(rvec_t* x, rq_t* norm, const rmatq_t* A);

/*
 function rmatq_rref(R, A):
    // Phase 0: Copy input and set up tracking
    rmatq_copy(R, A)
    rows = R.M.rows
    cols = R.M.cols

    pivot_cols = array of size cols
    pivot_rows = array of size rows
    rank = 0
    prev_pivot = 1

    // =========================================================================
    // Phase 1: Forward Pass (Bareiss Elimination to Row Echelon Form)
    // =========================================================================
    r = 0
    for c = 0 to cols - 1:
        // Find pivot row 'p' at or below 'r' with non-zero entry
        p = find_first_row_where(R.M, row >= r, col == c, value != 0)

        if p is not found:
            continue // Column has no pivot, move to next column

        // Swap current row with pivot row if necessary
        if p != r:
            swap_rows(R.M, r, p)

        pivot_rows[rank] = r
        pivot_cols[rank] = c
        pivot = R.M[r, c]

        // Bareiss integer reduction step
        for i = r + 1 to rows - 1:
            for j = c + 1 to cols - 1:
                R.M[i, j] = (pivot * R.M[i, j] - R.M[i, c] * R.M[r, j]) /
 prev_pivot R.M[i, c] = 0

        prev_pivot = pivot
        r = r + 1
        rank = rank + 1

    // =========================================================================
    // Phase 2: Backward Pass (Eliminate entries above pivots)
    // =========================================================================
    for k = rank - 1 down to 0:
        r = pivot_rows[k]
        c = pivot_cols[k]
        P_k = R.M[r, c]

        for i = 0 to r - 1:
            factor = R.M[i, c]
            if factor != 0:
                for j = c to cols - 1:
                    R.M[i, j] = P_k * R.M[i, j] - factor * R.M[r, j]
                R.M[i, c] = 0

    // =========================================================================
    // Phase 3: Construct Rational Rows (Divide each row by its pivot)
    // =========================================================================
    // Temporary 2D grid of rq_t fractions
    Frac[rows][cols] initialized to (0 / 1)

    for k = 0 to rank - 1:
        r = pivot_rows[k]
        c = pivot_cols[k]
        pivot_val = R.M[r, c]

        for j = 0 to cols - 1:
            Frac[r][j].num = R.M[r, j]
            Frac[r][j].den = pivot_val
            rq_normalize(&Frac[r][j]) // Simplifies num/den & ensures den > 0

    // =========================================================================
    // Phase 4: Find Global Denominator and Rebuild rmatq_t
    // =========================================================================
    common_den = 1
    for i = 0 to rows - 1:
        for j = 0 to cols - 1:
            common_den = lcm(common_den, Frac[i][j].den)

    for i = 0 to rows - 1:
        for j = 0 to cols - 1:
            scale = common_den / Frac[i][j].den
            R.M[i, j] = Frac[i][j].num * scale

    R.den = common_den

    // Reduce integer matrix entries and denominator by overall GCD
    rmatq_normalize(R)


        */

#endif
