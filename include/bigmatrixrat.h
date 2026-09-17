#ifndef BIGMATRIXRAT_H
#define BIGMATRIXRAT_H

#include "bignum.h"
#include "bigrat.h"

typedef struct {
  bigmatrix M;
  bignum den;
} bigmatrixrat;

void bigmatrixrat_init(bigmatrixrat* A, u64 m, u64 n);
void bigmatrixrat_free(bigmatrixrat* A);

void bigmatrixrat_copy(bigmatrixrat* R, bigmatrixrat* A);
void bigmatrixrat_copyZ(bigmatrixrat* R, bigmatrix* A);

void bigmatrixrat_add(bigmatrixrat* R, const bigmatrixrat* A,
                      const bigmatrixrat* B);
void bigmatrixrat_sub(bigmatrixrat* R, const bigmatrixrat* A,
                      const bigmatrixrat* B);

void bigmatrixrat_normalize(bigmatrixrat* M);
void bigmatrixrat_print(const bigmatrixrat* A);

/*
 function bigmatrixrat_rref(R, A):
    // Phase 0: Copy input and set up tracking
    bigmatrixrat_copy(R, A)
    rows = R.M.r_size
    cols = R.M.c_size

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
    // Temporary 2D grid of bigrat fractions
    Frac[rows][cols] initialized to (0 / 1)

    for k = 0 to rank - 1:
        r = pivot_rows[k]
        c = pivot_cols[k]
        pivot_val = R.M[r, c]

        for j = 0 to cols - 1:
            Frac[r][j].num = R.M[r, j]
            Frac[r][j].den = pivot_val
            br_normalize(&Frac[r][j]) // Simplifies num/den & ensures den > 0

    // =========================================================================
    // Phase 4: Find Global Denominator and Rebuild bigmatrixrat
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
    bigmatrixrat_normalize(R)


        */

#endif
