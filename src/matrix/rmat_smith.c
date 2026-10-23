/*
 * rmat_smith.c
 *
 * Smith normal form for rmat_t.
 *
 * Copyright (C) 2026 Diego Strebel
 * SPDX-License-Identifier: GPL-3.0-or-later
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#include <string.h>

#include "../../include/rmat.h"
#include "../../include/rzmath.h"
rabin_err_t rmat_smith(rmat_t* S, const rmat_t* A)
{
  if (S == NULL || A == NULL) return RABIN_ERR_NULL_PTR;
  if (A->cols != A->rows) return RABIN_ERR_MATRIX_DIM;
  if (S->rows != A->rows || S->cols != A->cols) return RABIN_ERR_MATRIX_DIM;

  // 1. [Initialize i]
  i64 n = A->rows;
  i64 i = n;
  rabin_err_t err = RABIN_SUCCESS;
  rz_t R, u, v, d, q_i, q_j, t1, t2, b;
  rz_t a_abs, a_ii_abs, a_abs_col, a_ii_abs_col, zero;
  rz_init_multi(&R, &u, &v, &d, &q_i, &q_j, &t1, &t2, &b, NULL);
  rz_init_multi(&a_abs, &a_ii_abs, &a_abs_col, &a_ii_abs_col, &zero, NULL);
  TRY(rmat_det(&R, A));
  R.is_neg = false;

  if (n == 1) {
    err = rmat_set(S, &R, 0, 0);
    goto cleanup;
  }

  rmat_t A_work = {0};
  rz_t* B = NULL;
  TRY(rmat_init(&A_work, n, n));
  TRY(rmat_copy(&A_work, A));

  B = malloc(n * sizeof(rz_t));
  if (B == NULL) {
    err = RABIN_ERR_OUT_OF_MEMORY;
    goto cleanup;
  }
  for (u64 x = 0; x < (u64)n; x++) {
    rz_init(&B[x]);
  }

  rz_t* a_i_i = NULL;
  rz_t* a_i_j = NULL;

// 2. [Initialize j for row reduction]
step2:
  i64 j = i;
  i64 c = 0;

// 3. [Check zero]
step3:
  if (j == 1) {
    goto step5;
  }

  j--;
  if (rz_is_zero(RMAT_GET(&A_work, i - 1, j - 1))) {
    goto step3;
  }

  // 4. [Euclidean step]
  a_i_i = RMAT_GET(&A_work, i - 1, i - 1);
  a_i_j = RMAT_GET(&A_work, i - 1, j - 1);

  if ((err = rz_gcd_extended_lehmer(&u, &v, &d, a_i_i, a_i_j)) != RABIN_SUCCESS)
    goto cleanup;

  // use Remark to find u,v minimal
  TRY(rz_copy(&a_abs, RMAT_GET(&A_work, i - 1, j - 1)));
  a_abs.is_neg = false;
  TRY(rz_copy(&a_ii_abs, RMAT_GET(&A_work, i - 1, i - 1)));
  a_ii_abs.is_neg = false;

  // if d == |a_{i,i}|, force u = sign(a_{i,i}), v = 0
  if (rz_cmp(&d, &a_ii_abs) == 0) {
    TRY(rz_set_u64(&u, 1));
    TRY(rz_set_u64(&v, 0));
    if (RMAT_GET(&A_work, i - 1, i - 1)->is_neg) TRY(rz_set_i64(&u, -1));
  }
  // if d == |a_{i,j}|, force u = 0, v = sign(a_{i,j})
  else if (rz_cmp(&d, &a_abs) == 0) {
    TRY(rz_set_u64(&u, 0));
    TRY(rz_set_u64(&v, 1));
    if (RMAT_GET(&A_work, i - 1, j - 1)->is_neg) TRY(rz_set_i64(&v, -1));
  }

  rz_clear(&a_abs);
  rz_clear(&a_ii_abs);

  // B = uA_i + vA_j
  for (u64 x = 0; x < (u64)n; x++) {
    TRY(rz_mul(&B[x], &u, RMAT_GET(&A_work, x, i - 1)));
    TRY(rz_mul(&t1, &v, RMAT_GET(&A_work, x, j - 1)));
    TRY(rz_add(&B[x], &B[x], &t1));
  }

  // A_j = ((a_ii / d)A_j - (a_ij / d)A_i) mod R
  TRY(rz_div(&q_i, RMAT_GET(&A_work, i - 1, i - 1), &d));
  TRY(rz_div(&q_j, RMAT_GET(&A_work, i - 1, j - 1), &d));
  for (u64 x = 0; x < (u64)n; x++) {
    TRY(rz_mul(&t1, &q_i, RMAT_GET(&A_work, x, j - 1)));
    TRY(rz_mul(&t2, &q_j, RMAT_GET(&A_work, x, i - 1)));
    TRY(rz_sub(&t1, &t1, &t2));
    TRY(rz_mod(&t1, &t1, &R));
    TRY(rz_copy(RMAT_GET(&A_work, x, j - 1), &t1));
  }

  // A_i = B mod R
  for (u64 x = 0; x < (u64)n; x++) {
    TRY(rz_mod(RMAT_GET(&A_work, x, i - 1), &B[x], &R));
  }

  goto step3;

// 5. [Initialize j for column reduction]
step5:
  j = i;

// 6. [Check zero]
step6:
  if (j == 1) {
    goto step8;
  }
  j--;
  if (rz_is_zero(RMAT_GET(&A_work, j - 1, i - 1))) {
    goto step6;
  }

  // 7. [Euclidean step]
  a_i_i = RMAT_GET(&A_work, i - 1, i - 1);
  a_i_j = RMAT_GET(&A_work, j - 1, i - 1);

  if ((err = rz_gcd_extended_lehmer(&u, &v, &d, a_i_i, a_i_j)) != RABIN_SUCCESS)
    goto cleanup;

  // use Remark to find u,v minimal safely
  TRY(rz_copy(&a_abs_col, RMAT_GET(&A_work, j - 1, i - 1)));
  a_abs_col.is_neg = false;
  TRY(rz_copy(&a_ii_abs_col, RMAT_GET(&A_work, i - 1, i - 1)));
  a_ii_abs_col.is_neg = false;

  // if d == |a_{i,i}|, force u = sign(a_{i,i}), v = 0
  if (rz_cmp(&d, &a_ii_abs_col) == 0) {
    TRY(rz_set_u64(&u, 1));
    TRY(rz_set_u64(&v, 0));
    if (RMAT_GET(&A_work, i - 1, i - 1)->is_neg) TRY(rz_set_i64(&u, -1));
  }
  // if d == |a_{j,i}|, force u = 0, v = sign(a_{j,i})
  else if (rz_cmp(&d, &a_abs_col) == 0) {
    TRY(rz_set_u64(&u, 0));
    TRY(rz_set_u64(&v, 1));
    if (RMAT_GET(&A_work, j - 1, i - 1)->is_neg) TRY(rz_set_i64(&v, -1));
  }

  rz_clear(&a_abs_col);
  rz_clear(&a_ii_abs_col);

  // B = uA'_i + vA'_j
  for (u64 x = 0; x < (u64)n; x++) {
    TRY(rz_mul(&B[x], &u, RMAT_GET(&A_work, i - 1, x)));
    TRY(rz_mul(&t1, &v, RMAT_GET(&A_work, j - 1, x)));
    TRY(rz_add(&B[x], &B[x], &t1));
  }

  // A'_j = ((a_ii / d)A'_j - (a_ij / d)A'_i) mod R
  TRY(rz_div(&q_i, RMAT_GET(&A_work, i - 1, i - 1), &d));
  TRY(rz_div(&q_j, RMAT_GET(&A_work, j - 1, i - 1), &d));
  for (u64 x = 0; x < (u64)n; x++) {
    TRY(rz_mul(&t1, &q_i, RMAT_GET(&A_work, j - 1, x)));
    TRY(rz_mul(&t2, &q_j, RMAT_GET(&A_work, i - 1, x)));
    TRY(rz_sub(&t1, &t1, &t2));
    TRY(rz_mod(&t1, &t1, &R));
    TRY(rz_copy(RMAT_GET(&A_work, j - 1, x), &t1));
  }

  // A'_i = B mod R
  for (u64 x = 0; x < (u64)n; x++) {
    TRY(rz_mod(RMAT_GET(&A_work, i - 1, x), &B[x], &R));
  }

  c++;
  goto step6;

// 8. [Repeat stage i?]
step8:
  if (c > 0) {
    goto step2;
  }

  // 9. [Check the rest of the matrix]
  TRY(rz_copy(&b, RMAT_GET(&A_work, i - 1, i - 1)));
  for (i64 k_0 = 0; k_0 < i - 1; k_0++) {
    for (i64 l_0 = 0; l_0 < i - 1; l_0++) {
      bool not_divisible = false;
      if (rz_is_zero(&b)) {
        if (!rz_is_zero(RMAT_GET(&A_work, k_0, l_0))) {
          not_divisible = true;
        }
      } else {
        TRY(rz_mod(&t1, RMAT_GET(&A_work, k_0, l_0), &b));
        if (!rz_is_zero(&t1)) {
          not_divisible = true;
        }
      }

      if (not_divisible) {
        for (i64 x = 0; x < n; x++) {
          TRY(rz_add(RMAT_GET(&A_work, i - 1, x), RMAT_GET(&A_work, i - 1, x),
                     RMAT_GET(&A_work, k_0, x)));
        }
        goto step2;
      }
    }
  }

  // 10. [Next stage]
  TRY(rz_gcd_lehmer(&t1, RMAT_GET(&A_work, i - 1, i - 1), &R));
  TRY(rmat_set(&A_work, &t1, i - 1, i - 1));
  TRY(rz_div(&R, &R, &t1));

  if (i == 2) {
    TRY(rz_gcd_lehmer(&t2, RMAT_GET(&A_work, 0, 0), &R));
    TRY(rmat_set(&A_work, &t2, 0, 0));

    // order diagonal entries in proper order
    for (i64 r = 0; r < n; r++) {
      for (i64 c = 0; c < n; c++) {
        if (r == c) {
          TRY(rmat_set(S, RMAT_GET(&A_work, n - 1 - r, n - 1 - r), r, r));
        } else {
          TRY(rmat_set(S, &zero, r, c));
        }
      }
    }
    rz_clear(&zero);
    goto cleanup;
  }

  i--;
  goto step2;

cleanup:
  if (B != NULL) {
    for (i64 x = 0; x < n; x++) {
      rz_clear(&B[x]);
    }
    free(B);
  }
  rz_clear(&a_abs_col);
  rz_clear(&a_ii_abs_col);
  rz_clear(&a_abs);
  rz_clear(&a_ii_abs);
  rz_clear(&zero);
  rmat_clear(&A_work);
  rz_clear(&R);
  rz_clear_multi(&u, &v, &d, &q_i, &q_j, &t1, &t2, &b, NULL);
  return err;
}
