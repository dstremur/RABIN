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
typedef struct {
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
rabin_err_t rmatq_copy(rmatq_t* R, rmatq_t* A);

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
