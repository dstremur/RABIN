#include "../../include/bigmatrix.h"

#include <stdio.h>
void bigmatrix_init(bigmatrix* M, u64 r, u64 c)
{
  M->c_size = c;
  M->r_size = r;

  M->data = malloc(r * c * sizeof(bignum));

  for (u64 i = 0; i < r * c; i++) {
    bn_init(&M->data[i]);
  }
}

void bigmatrix_free(bigmatrix* A)
{
  if (!A || !A->data) return;

  for (u64 i = 0; i < A->r_size * A->c_size; i++) {
    bn_free(&A->data[i]);
  }

  free(A->data);

  A->data = NULL;
}

void bigmatrix_copy(bigmatrix* R, bigmatrix* A)
{
  // check if sizes match
  if (A->c_size != R->c_size || A->r_size != R->r_size) return;

  for (u64 i = 0; i < A->c_size * A->r_size; i++) {
    bn_copy(&R->data[i], &A->data[i]);
  }
}

void bigmatrix_get(bignum* R, bigmatrix* A, u64 r, u64 c)
{
  if (r >= A->r_size || c >= A->c_size) return;

  bn_copy(R, GET(A, r, c));
}

void bigmatrix_set(bigmatrix* A, bignum* a, u64 r, u64 c)
{
  if (r >= A->r_size || c >= A->c_size) return;

  bn_copy(GET(A, r, c), a);
}

void bigmatrix_add(bigmatrix* R, bigmatrix* A, bigmatrix* B)
{
  // check if sizes match
  if (A->c_size != B->c_size || A->r_size != B->r_size) return;

  for (u64 i = 0; i < A->c_size * A->r_size; i++) {
    bn_add(&R->data[i], &A->data[i], &B->data[i]);
  }

  return;
}

void bigmatrix_print(bigmatrix* A)
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

void bigmatrix_mul(bigmatrix* R, bigmatrix* A, bigmatrix* B)
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

// Bareiss algorithm
/* optimize:
 use exact division
 better cache locality
 preallocat using hadamards bound
 openmp
 calc mod primes larger than hadamard and then reconstruct using CRT
 Jebelean’s algorithm
 */

void bigmatrix_det(bignum* d, bigmatrix* A)
{
  if (A->c_size != A->r_size) {
    printf("Matrix must be square\n");
    return;
  }
  u64 n = A->c_size;

  bigmatrix T;
  bigmatrix_init(&T, A->c_size, A->r_size);
  bigmatrix_copy(&T, A);

  bignum prev, temp1, temp2, temp3;
  bn_init_multi(&prev, &temp1, &temp2, &temp3, NULL);

  // Track sign changes
  i64 sign = 1;

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

    for (u64 i = k + 1; i < n; i++) {
      for (u64 j = k + 1; j < n; j++) {
        // T_ij = (T_ij * T_kk - T_ik * T_kj) / T_kk

        bn_mul(&temp1, GET(&T, i, j), pivot);
        bn_mul(&temp2, GET(&T, i, k), GET(&T, k, j));
        bn_sub(&temp3, &temp1, &temp2);

        bn_div(GET(&T, i, j), &temp3, &prev);
      }
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
  bn_free_multi(&prev, &temp1, &temp2, &temp3, NULL);
}
