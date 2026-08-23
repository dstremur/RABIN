/*
 * bigmatrix.c
 *
 * Matrix arithmetic over bignums.
 *
 * This file implements matrices with bignum entries: initialization,
 * freeing, copying, element/row/column access, printing,
 * component-wise addition, schoolbook multiplication, matrix-vector
 * products, the Hadamard bound on the determinant, and the
 * determinant (computed via the RNS path; a direct Bareiss
 * implementation is kept below it, currently disabled).
 *
 * A bigmatrix is a row-major dynamic array of bignums with row and
 * column counts.
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

#include "../../include/bigmatrix.h"

#include <stdio.h>

#include "../../include/bigrns.h"
#include "../../include/primes.h"

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
void bigmatrix_init(bigmatrix* M, u64 r, u64 c)
{
  M->c_size = c;
  M->r_size = r;

  M->data = malloc(r * c * sizeof(bignum));

  for (u64 i = 0; i < r * c; i++) {
    bn_init(&M->data[i]);
  }
}

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
void bigmatrix_free(bigmatrix* A)
{
  if (!A || !A->data) return;

  for (u64 i = 0; i < A->r_size * A->c_size; i++) {
    bn_free(&A->data[i]);
  }

  free(A->data);

  A->data = NULL;
}

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
void bigmatrix_print_python(const bigmatrix* A)
{
  printf("[");  // Start outer list
  for (u64 i = 0; i < A->r_size; i++) {
    printf("[");  // Start row
    for (u64 j = 0; j < A->c_size; j++) {
      bn_print(GET(A, i, j));
      if (j < A->c_size - 1) printf(", ");
    }
    printf("]");  // End row
    if (i < A->r_size - 1) printf(",\n ");
  }
  printf("]\n");  // End outer list
}

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
void bigmatrix_copy(bigmatrix* R, bigmatrix* A)
{
  // check if sizes match
  if (A->c_size != R->c_size || A->r_size != R->r_size) return;

  for (u64 i = 0; i < A->c_size * A->r_size; i++) {
    bn_copy(&R->data[i], &A->data[i]);
  }
}

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
void bigmatrix_get(bignum* R, const bigmatrix* A, u64 r, u64 c)
{
  if (r >= A->r_size || c >= A->c_size) return;

  bn_copy(R, GET(A, r, c));
}

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
void bigmatrix_set(bigmatrix* A, const bignum* a, u64 r, u64 c)
{
  if (r >= A->r_size || c >= A->c_size) return;

  bn_copy(GET(A, r, c), a);
}

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
void bigmatrix_get_col(bigvector* c, const bigmatrix* A, u64 col)
{
  if (c->size != A->r_size) printf("Size does not match \n");

  for (u64 i = 0; i < A->r_size; i++) {
    bigvector_set(c, GET(A, i, col), i);
  }
}

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
void bigmatrix_get_row(bigvector* r, const bigmatrix* A, u64 row)
{
  if (r->size != A->c_size) printf("Size does not match \n");

  for (u64 i = 0; i < A->c_size; i++) {
    bigvector_set(r, GET(A, row, i), i);
  }
}

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
void bigmatrix_mv(bigvector* r, const bigmatrix* A, bigvector* v)
{
  if (r->size != A->r_size || v->size != A->c_size) {
    printf("Size does not match \n");
    return;
  }

  bigvector tmp;
  bigvector_init(&tmp, A->c_size);

  bignum dot;
  bn_init(&dot);

  for (u64 i = 0; i < A->r_size; i++) {
    bigmatrix_get_row(&tmp, A, i);
    bigvector_dot(&dot, v, &tmp);
    bigvector_set(r, &dot, i);
  }

  bn_free(&dot);
  bigvector_free(&tmp);
}

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
void bigmatrix_vm(bigvector* r, const bigmatrix* A, bigvector* v)
{
  if (r->size != A->c_size || v->size != A->r_size) {
    printf("Size does not match \n");
    return;
  }

  bigvector tmp;
  bigvector_init(&tmp, A->r_size);

  bignum dot;
  bn_init(&dot);

  for (u64 i = 0; i < A->c_size; i++) {
    bigmatrix_get_col(&tmp, A, i);
    bigvector_dot(&dot, v, &tmp);
    bigvector_set(r, &dot, i);
  }

  bn_free(&dot);
  bigvector_free(&tmp);
}

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
void bigmatrix_hadamard(bignum* r, const bigmatrix* A)
{
  bn_set_u64(r, 1);

  bigvector v;
  bigvector_init(&v, A->r_size);
  bignum tmp;
  bn_init(&tmp);

  for (u64 i = 0; i < A->c_size; i++) {
    bigmatrix_get_col(&v, A, i);

    bigvector_norm(&tmp, &v);
    // add 1 since norm uses isqrt
    bn_add_u64(&tmp, &tmp, 1);
    bn_mul(r, r, &tmp);
  }

  bn_free(&tmp);
  bigvector_free(&v);
}

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
void bigmatrix_add(bigmatrix* R, const bigmatrix* A, const bigmatrix* B)
{
  // check if sizes match
  if (A->c_size != B->c_size || A->r_size != B->r_size) return;

  for (u64 i = 0; i < A->c_size * A->r_size; i++) {
    bn_add(&R->data[i], &A->data[i], &B->data[i]);
  }

  return;
}

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
void bigmatrix_print(const bigmatrix* A)
{
  bignum temp;
  bn_init(&temp);
  for (u64 i = 0; i < A->r_size; i++) {
    for (u64 j = 0; j < A->c_size; j++) {
      bigmatrix_get(&temp, A, i, j);
      bn_print(&temp);
      printf(" ");
    }
    printf("\n");
  }

  bn_free(&temp);
}

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
void bigmatrix_mul(bigmatrix* R, const bigmatrix* A, const bigmatrix* B)
{
  // check if sizes match
  if (A->c_size != B->r_size) return;

  bignum sum, tmp;
  bn_init_multi(&sum, &tmp, NULL);
  for (u64 i = 0; i < A->r_size; i++) {
    for (u64 j = 0; j < B->c_size; j++) {
      bn_set_u64(&sum, 0);

      for (u64 k = 0; k < A->c_size; k++) {
        // tmp = A[i][k] * B[k][j]
        bn_mul(&tmp, &A->data[i * A->c_size + k], &B->data[k * B->c_size + j]);

        bn_add(&sum, &sum, &tmp);
      }

      bn_copy(&R->data[i * B->c_size + j], &sum);
    }
  }
  bn_free(&sum);
  bn_free(&tmp);
}

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
void bigmatrix_det(bignum* d, const bigmatrix* A)
{
  if (A->c_size != A->r_size) {
    printf("Matrix must be square\n");
    return;
  }
  u64 n = A->c_size;

  bigmatrix T;
  bigmatrix_init(&T, A->c_size, A->r_size);
  bigmatrix_copy(&T, A);

  u64 k = rns_estimate_determinant(&T);

  ctx_rns ctx;
  rns_context_init(&ctx, RNS_PRIMES, k);

  bigmatrix_det_rns(d, &T, &ctx);

  rns_context_free(&ctx);
  bigmatrix_free(&T);

  return;

  // broken

  bignum prev;
  bn_init(&prev);

  // Track sign changes
  i64 sign = 1;

  bigmatrix_hadamard(&prev, A);

  bn_alloc(d, prev.size);

  bn_set_u64(&prev, 1);

  for (u64 k = 0; k < n - 1; k++) {
    bignum* pivot = GET(&T, k, k);
    // swap if pivot is zero
    if (bn_is_zero(pivot)) {
      u64 swap = k + 1;
      bool found = false;

      while (swap < n) {
        bignum* t1 = GET(&T, swap, k);
        if (!bn_is_zero(t1)) {
          found = true;
          break;
        }
        swap++;
      }

      if (!found) {
        bn_set_u64(d, 0);
        goto cleanup;
      }

      // swap
      for (u64 j = k; j < n; j++) {
        bn_swap(GET(&T, k, j), GET(&T, swap, j));
      }
      sign *= -1;
      pivot = GET(&T, k, k);
    }

    // #pragma omp parallel
    {
      bignum temp1, temp2, temp3;
      bn_init_multi(&temp1, &temp2, &temp3, NULL);
      // #pragma omp for collapse(2) schedule(static)
      for (u64 i = k + 1; i < n; i++) {
        bignum* t = GET(&T, i, k);
        for (u64 j = k + 1; j < n; j++) {
          // T_ij = (T_ij * T_kk - T_ik * T_kj) / T_kk

          bn_mul(&temp1, GET(&T, i, j), pivot);
          bn_mul(&temp2, t, GET(&T, k, j));
          bn_sub(&temp3, &temp1, &temp2);

          bn_div(GET(&T, i, j), &temp3, &prev);
        }
      }

      bn_free_multi(&temp1, &temp2, &temp3, NULL);
    }

    bn_copy(&prev, pivot);
  }

  // det = M_nn
  bn_copy(d, GET(&T, n - 1, n - 1));
  if (sign == -1) {
    d->is_neg = !d->is_neg;
  }

cleanup:
  bigmatrix_free(&T);
  bn_free_multi(&prev, NULL);
}
