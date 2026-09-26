#ifndef BIGMATRIX_H
#define BIGMATRIX_H

#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>

#include "bigcore.h"
#include "bigpoly.h"
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
 */
void bigmatrix_init(bigmatrix* M, u64 r, u64 c);

/**
 * @brief Free all storage of a matrix.
 *
 * Complexity:
 *   - Time: \f$O(r \cdot c)\f$
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(1)\f$
 *
 * @param[in,out] A Matrix to free.
 */
void bigmatrix_free(bigmatrix* A);

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
 */
void bigmatrix_print_python(const bigmatrix* A);

/**
 * @brief Copy a matrix: \f$R = A\f$.
 *
 * Let \f$r =\f$ A->r_size, \f$c =\f$ A->c_size.
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
 */
void bigmatrix_copy(bigmatrix* R, const bigmatrix* A);

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
 */
void bigmatrix_get(bignum* R, const bigmatrix* A, u64 r, u64 c);

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
 */
void bigmatrix_set(bigmatrix* A, const bignum* a, u64 r, u64 c);

/**
 * @brief Extract a column into a vector: \f$c = A[:, col]\f$.
 *
 * The vector \f$c\f$ must already be allocated with r_size elements; a
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
 */
void bigmatrix_get_col(bigvector* c, const bigmatrix* A, u64 col);

/**
 * @brief Extract a row into a vector: \f$r = A[row, :]\f$.
 *
 * The vector \f$r\f$ must already be allocated with c_size elements; a
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
 */
void bigmatrix_get_row(bigvector* r, const bigmatrix* A, u64 row);

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
bool bigmatrix_is_square(const bigmatrix* A);

/**
 * @brief Test whether all entries of a matrix are zero.
 *
 * Complexity:
 *   - Time: \f$O(r \cdot c)\f$ where \f$r =\f$ A->r_size, \f$c =\f$ A->c_size,
 *           stopping at the first nonzero entry
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(1)\f$
 *
 * @param[in] A Matrix.
 *
 * @return true If every entry is zero.
 */
bool bigmatrix_is_zero(const bigmatrix* A);

/**
 * @brief Test whether a matrix is the identity matrix.
 *
 * The matrix must be square with 1s on the main diagonal and 0s
 * elsewhere. Non-square matrices are rejected before any entry is
 * inspected.
 *
 * Complexity:
 *   - Time: \f$O(n^2)\f$ where \f$n =\f$ A->r_size, stopping at the first
 *   violation
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(1)\f$
 *
 * @param[in] A Matrix.
 *
 * @return true If A is the identity matrix.
 */
bool bigmatrix_is_identity(const bigmatrix* A);

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
bool bigmatrix_is_diagonal(const bigmatrix* A);

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
bool bigmatrix_is_upper_triangular(const bigmatrix* A);

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
bool bigmatrix_is_lower_triangular(const bigmatrix* A);

/**
 * @brief Test whether a matrix is symmetric (\f$A[i][j] = A[j][i]\f$).
 *
 * The matrix must be square; only the strict upper triangle is walked,
 * comparing each mirror pair once. Non-square matrices are rejected
 * before any entry is inspected.
 *
 * Complexity:
 *   - Time: \f$O(n^2 \cdot n_{b})\f$ where \f$n =\f$ A->r_size and
 *   \f$n_{b} =\f$ the size of the entries in limbs, stopping at the
 *   first mismatching pair
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(1)\f$
 *
 * @param[in] A Matrix.
 *
 * @return true If A is square and A[i][j] == A[j][i] for all (i, j).
 */
bool bigmatrix_is_symmetric(const bigmatrix* A);

void bigmatrix_scalar(bigmatrix* R, const bigmatrix* A, const bignum* a);

void bigmatrix_div_exact_scalar(bigmatrix* R, const bigmatrix* A,
                                const bignum* a);

/**
 * @brief Matrix-vector product: \f$r = A \cdot v\f$.
 *
 * Let \f$r =\f$ A->r_size, \f$c =\f$ A->c_size.
 *
 * Each output element is the dot product of the corresponding row of
 * \f$A\f$ with \f$v\f$. Sizes must match (\f$r\f$ has r_size elements, \f$v\f$
 * has c_size).
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
 */
void bigmatrix_mv(bigvector* r, const bigmatrix* A, bigvector* v);

/**
 * @brief Vector-matrix product: \f$r = v \cdot A\f$.
 *
 * Let \f$r =\f$ A->r_size, \f$c =\f$ A->c_size.
 *
 * Each output element is the dot product of \f$v\f$ with the corresponding
 * column of \f$A\f$. Sizes must match (\f$v\f$ has r_size elements, \f$r\f$ has
 * c_size).
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
 */
void bigmatrix_vm(bigvector* r, const bigmatrix* A, bigvector* v);

/**
 * @brief Hadamard bound on the determinant: \f$r = \prod_i \lVert col_i
 * \rVert\f$.
 *
 * Let \f$n =\f$ A->c_size.
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
 */
void bigmatrix_hadamard(bignum* r, const bigmatrix* A);

/**
 * @brief Component-wise addition of two matrices: \f$R = A + B\f$.
 *
 * Let \f$r =\f$ A->r_size, \f$c =\f$ A->c_size.
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
 */
void bigmatrix_add(bigmatrix* R, const bigmatrix* A, const bigmatrix* B);

void bigmatrix_sub(bigmatrix* R, const bigmatrix* A, const bigmatrix* B);

/**
 * @brief Print a matrix to stdout, one row per line, with right-aligned
 * columns.
 *
 * The width of each column is derived from its widest entry. Entries whose
 * decimal form is wider than 16 characters are truncated to their leading
 * digits followed by an ellipsis, so large matrices stay compact. A single
 * oversized entry (e.g. the last invariant factor of a Smith form sitting
 * in a column of zeros) does not inflate its column: such outliers overflow
 * the column to the right instead. Use bigmatrix_print_full() for
 * untruncated output and bigmatrix_print_tail() for truncation to the
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
 */
void bigmatrix_print(const bigmatrix* A);

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
 */
void bigmatrix_print_full(const bigmatrix* A);

/**
 * @brief Print a matrix to stdout, one row per line, with right-aligned
 * columns.
 *
 * Like bigmatrix_print(), entries wider than 16 characters are truncated,
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
 */
void bigmatrix_print_tail(const bigmatrix* A);

/**
 * @brief Schoolbook matrix multiplication: \f$R = A \cdot B\f$.
 *
 * Let \f$r =\f$ A->r_size, \f$k =\f$ A->c_size, \f$c =\f$ B->c_size.
 *
 * No-op if A->c_size != B->r_size.
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
 */
void bigmatrix_mul(bigmatrix* R, const bigmatrix* A, const bigmatrix* B);

/**
 * @brief Determinant of a square bignum matrix.
 *
 * Let \f$n =\f$ A->c_size.
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
 */
void bigmatrix_det(bignum* d, const bigmatrix* A);

/**
 * @brief Returns the identity matrix
 *
 * @param I
 * @param n
 */
void bigmatrix_id(bigmatrix* I, const u64 n);

void bigmatrix_gcd_all(bignum* g, const bigmatrix* A);

void bigmatrix_swap(bigmatrix* a, bigmatrix* b);

void bigmatrix_trace(bignum* t, const bigmatrix* A);

void bigmatrix_charpoly_adj(bigpoly* p, bigmatrix* J, const bigmatrix* A);

void bigmatrix_trace(bignum* t, const bigmatrix* A);

/**
 * @brief Compute the LLL algorithm on the basis matrix B
 *
 * @param B
 * @param n
 * @param delta
 * @param H
 */
void bigmatrix_LLL(bigmatrix* B, u64 n, double delta, bigmatrix* H);

/**
 * @brief Hermite normal form via the exact (coefficient-blowup prone)
 * column reduction (Cohen's Algorithm 2.4.4).
 *
 * Let \f$m =\f$ A->r_size, \f$n =\f$ A->c_size.
 *
 * Handles zero entries in the pivot positions \f$a_{i,k}\f$ (the extended
 * GCD is safe on a zero argument). This is the reference implementation:
 * use it to cross-check bigmatrix_hermite_mod_d().
 *
 * @param[out] W Result storing the HNF of A.
 * @param[in]  A Input matrix.
 */
void bigmatrix_hermite(bigmatrix* W, const bigmatrix* A);

/**
 * @brief Hermite normal form, Cohen's Algorithm 2.4.5 (exact integer
 * arithmetic).
 *
 * Let \f$m =\f$ A->r_size, \f$n =\f$ A->c_size.
 *
 * Works entirely with integers and Euclidean divisions, but is subject to
 * the coefficient explosion phenomenon; prefer
 * bigmatrix_hermite_mod_d() when a multiple of the module determinant is
 * known.
 *
 * Handles zero entries in the pivot positions, and each Euclidean step uses
 * the Bezout pair normalized to the canonical range of Cohen's "Important
 * Remark" (v = 0 when a_{i,k} | a_{i,j}, |v| <= |a_{i,k} / d| / 2
 * otherwise), which keeps the intermediate coefficients small.
 *
 * @param[out] W Result storing the HNF of A.
 * @param[in]  A Input matrix.
 */
void bigmatrix_hermite_gcd(bigmatrix* W, const bigmatrix* A);

/**
 * @brief Hermite normal form modulo D, Cohen's Algorithm 2.4.8 (essentially
 * due to Domich et al. [DKT]).
 *
 * Let \f$m =\f$ A->r_size, \f$n =\f$ A->c_size.
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
 * \f$D = |\det A|\f$ (e.g. from bigmatrix_det()) always works. A is allowed
 * to have zero entries, including at the trailing pivot \f$a_{m - 1, n - 1}\f$.
 * If \f$m > n\f$, the computation falls back to
 * bigmatrix_hermite_gcd().
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
 */
void bigmatrix_hermite_mod_d(bigmatrix* W, const bigmatrix* A, const bignum* D);

void bigmatrix_smith(bigmatrix* S, const bigmatrix* A);

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
bool bigmatrix_equal(const bigmatrix* A, const bigmatrix* B);

/**
 * @brief Validate that H strictly satisfies the canonical Hermite normal
 * form structure.
 *
 * Let \f$r =\f$ H->r_size, \f$c =\f$ H->c_size.
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
bool bigmatrix_hnf_check_structure(const bigmatrix* H);

/**
 * @brief Verify the unimodular transformation \f$A \cdot U = H\f$.
 *
 * This library's HNF routines act by column operations, so the
 * transformation matrix \f$U\f$ right-multiplies: \f$H = A \cdot U\f$.
 *
 * The product is computed with bigmatrix_mul() into a temporary matrix
 * and compared entrywise with H; the temporary is freed before
 * returning.
 *
 * Complexity:
 *   - Time: \f$O(r \cdot c \cdot n \cdot n_{b}^2)\f$ where \f$r =\f$ A->r_size,
 *   \f$c =\f$ A->c_size, \f$n =\f$ U->c_size, \f$n_{b} =\f$ the size of
 *   the entries in limbs
 *   - Auxiliary memory: \f$O(r \cdot n)\f$ bignums for the product
 *   - Output memory: \f$O(1)\f$
 *
 * @param[in] A Original matrix.
 * @param[in] H Claimed normal form of A.
 * @param[in] U Transformation matrix (expected square, \f$c =\f$ A->c_size).
 *
 * @return true If the dimensions are compatible and A * U == H entrywise.
 */
bool bigmatrix_hnf_check_transformation(const bigmatrix* A, const bigmatrix* H,
                                        const bigmatrix* U);

/**
 * @brief Verify that U is unimodular (integer matrix with
 * \f$|\det U| = 1\f$).
 *
 * Complexity:
 *   - Time: same as bigmatrix_det() for an \f$n \times n\f$ matrix
 *   - Auxiliary memory: as in bigmatrix_det()
 *   - Output memory: \f$O(1)\f$
 *
 * @param[in] U Matrix to validate.
 *
 * @return true If U is square and |det(U)| == 1.
 */
bool bigmatrix_hnf_check_unimodular(const bigmatrix* U);

/**
 * @brief Master HNF verification: structure, transformation, and
 * unimodularity.
 *
 * Runs bigmatrix_hnf_check_structure(),
 * bigmatrix_hnf_check_transformation(), and
 * bigmatrix_hnf_check_unimodular() in that order, returning early on the
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
bool bigmatrix_hnf_verify(const bigmatrix* A, const bigmatrix* H,
                          const bigmatrix* U);
#endif
