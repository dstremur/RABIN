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

#include "../../include/bigpoly.h"
#include "../../include/bigrns.h"
#include "../../include/primes.h"

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

void bigmatrix_copy(bigmatrix* R, const bigmatrix* A)
{
  // check if sizes match
  if (A->c_size != R->c_size || A->r_size != R->r_size) return;

  for (u64 i = 0; i < A->c_size * A->r_size; i++) {
    bn_copy(&R->data[i], &A->data[i]);
  }
}

void bigmatrix_get(bignum* R, const bigmatrix* A, u64 r, u64 c)
{
  if (r >= A->r_size || c >= A->c_size) return;

  bn_copy(R, GET(A, r, c));
}

void bigmatrix_set(bigmatrix* A, const bignum* a, u64 r, u64 c)
{
  if (r >= A->r_size || c >= A->c_size) return;

  bn_copy(GET(A, r, c), a);
}

void bigmatrix_get_col(bigvector* c, const bigmatrix* A, u64 col)
{
  if (c->size != A->r_size) printf("Size does not match \n");

  for (u64 i = 0; i < A->r_size; i++) {
    bigvector_set(c, GET(A, i, col), i);
  }
}

void bigmatrix_get_row(bigvector* r, const bigmatrix* A, u64 row)
{
  if (r->size != A->c_size) printf("Size does not match \n");

  for (u64 i = 0; i < A->c_size; i++) {
    bigvector_set(r, GET(A, row, i), i);
  }
}

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

void bigmatrix_add(bigmatrix* R, const bigmatrix* A, const bigmatrix* B)
{
  // check if sizes match
  if (A->c_size != B->c_size || A->r_size != B->r_size) return;

  for (u64 i = 0; i < A->c_size * A->r_size; i++) {
    bn_add(&R->data[i], &A->data[i], &B->data[i]);
  }

  return;
}

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

void bigmatrix_mul(bigmatrix* R, const bigmatrix* A, const bigmatrix* B)
{
  // check if sizes match
  if (A->c_size != B->r_size) return;

  bignum sum, tmp;
  bn_init_multi(&sum, &tmp, NULL);
  for (u64 i = 0; i < A->r_size; i++) {
    for (u64 j = 0; j < B->c_size; j++) {
      bn_set_i64(&sum, 0);
      for (u64 k = 0; k < A->c_size; k++) {
        // tmp = A[i][k] * B[k][j]
        bn_mul(&tmp, &A->data[i * A->c_size + k], &B->data[k * B->c_size + j]);

        bn_add(&sum, &sum, &tmp);
      }

      bigmatrix_set(R, &sum, i, j);
    }
  }
  bn_free(&sum);
  bn_free(&tmp);
}

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
}

void bigmatrix_neg(const bigmatrix* A)
{
  for (u64 i = 0; i < A->r_size; i++) {
    for (u64 j = 0; j < A->c_size; j++) {
      bn_neg(GET(A, i, j), GET(A, i, j));
    }
  }
}

void bigmatrix_trace(bignum* t, const bigmatrix* A)
{
  assert(A->c_size == A->r_size);
  u64 n = A->c_size;
  bn_set_u64(t, 0);
  for (u64 i = 0; i < n; i++) {
    bn_add(t, t, GET(A, i, i));
  }
}

void bigmatrix_swap(bigmatrix* a, bigmatrix* b)
{
  bigmatrix t = *a;
  *a = *b;
  *b = t;
}

void bigmatrix_charpoly_adj(bigpoly* p, bigmatrix* J, const bigmatrix* A)
{
  assert(A->c_size == A->r_size);

  // 1. [Initialize]
  u64 n = A->c_size;

  bigmatrix C, T;
  bigmatrix_init(&C, n, n);
  bigmatrix_init(&T, n, n);
  bigmatrix_id(&C, n);

  printf("A: \n");
  bigmatrix_print(A);
  printf("\n");
  printf("C: \n");
  bigmatrix_print(&C);
  printf("\n");

  bignum temp, trace;
  bn_init_multi(&temp, &trace, NULL);
  bignum* a = malloc((n + 1) * sizeof(bignum));
  for (u64 i = 0; i <= n; i++) {
    bn_init(&a[i]);
  }
  bn_set_u64(&a[0], 1);

  // 2. [Finished?]
  for (u64 i = 1; i <= n; i++) {
    // 3. [Compute next a_i and C]
    bigmatrix_mul(&T, A, &C);

    bigmatrix_trace(&trace, &T);

    bn_set_u64(&temp, i);
    // bn_div(&a[i], &trace, &temp);
    bn_divmod_u64(&a[i], &trace, i);

    bn_neg(&a[i], &a[i]);

    if (J && i == n) {
      bigmatrix_copy(J, &C);
      if ((n - 1) & 1) {
        bigmatrix_neg(J);
      }
    }
    bigmatrix_copy(&C, &T);
    for (u64 j = 0; j < n; j++) {
      bn_add(&temp, GET(&C, j, j), &a[i]);
      bigmatrix_set(&C, &temp, j, j);
    }
  }

  bigpoly_alloc(p, n + 1);

  for (u64 i = 0; i <= n; i++) {
    bn_println(&a[i]);
    // bn_init(&p->coeff[n - i]);
    bn_copy(&p->coeff[n - i], &a[i]);
    bn_free(&a[i]);
  }

  p->deg = n;

  free(a);
  bn_free_multi(&temp, &trace, NULL);
  bigmatrix_free(&C);
  bigmatrix_free(&T);
}

void bigmatrix_id(bigmatrix* I, const u64 n)
{
  bignum a, b;
  bn_init_multi(&a, &b, NULL);
  bn_set_u64(&a, 1);
  bn_set_u64(&b, 0);

  for (u64 i = 0; i < n; i++) {
    for (u64 j = 0; j < n; j++) {
      if (i == j) {
        bigmatrix_set(I, &a, i, j);
      } else {
        bigmatrix_set(I, &b, i, j);
      }
    }
  }

  bn_free_multi(&a, &b, NULL);
}

// Computation Number theory p. 86
void bigmatrix_LLL(bigmatrix* B, u64 n, double delta, bigmatrix* H)
{
  // 1 [Initialize]
  u64 k = 2;
  u64 k_max = 1;

  // allocate lists for orth. vectors and squared norms
  bigvector* b_star = malloc(n * sizeof(bigvector));
  bignum* B_vals = malloc(n * sizeof(bignum));
  for (u64 i = 0; i < n; i++) {
    bigvector_init(&b_star[i], B->c_size);
    bn_init(&B_vals[i]);
  }

  // Gram schmidt coefficients
  double* mu = calloc(n * n, sizeof(double));

  bigvector b_k;
  bigvector_init(&b_k, n);
  bigmatrix_get_col(&b_star[0], B, 0);
  bigvector_copy(&b_k, &b_star[0]);
  bigvector_dot(&B_vals[0], &b_k, &b_k);

  bigmatrix_id(H, n);

  bigvector temp_vec;
  bigvector_init(&temp_vec, B->c_size);

  // 2 [Incremental Gram-Schmidt]

  if (k > k_max) {
    k_max = k;
  }
}

// void bigmatrix_LLL(bigmatrix* B, u64 n, double delta, bigmatrix* H)
// {
//   u64 k = 2;
//   u64 k_max = 1;

//   // Allocate Gram-Schmidt orthogonalized vectors and squared norms
//   bigvector* b_star = malloc(n * sizeof(bigvector));
//   bignum* B_vals = malloc(n * sizeof(bignum));
//   for (u64 i = 0; i < n; i++) {
//     bigvector_init(&b_star[i], B->rows);
//     bn_init(&B_vals[i]);
//   }

//   // Allocate table for Gram-Schmidt coefficients mu[k][j] stored as doubles
//   double* mu = calloc(n * n, sizeof(double));

//   // Step 1 [Initialize]
//   bigvector b_k;
//   bigvector_init(&b_k, B->rows);
//   bigmatrix_get_col(&b_star[0], B, 0);
//   bigvector_copy(&b_k, &b_star[0]);
//   bigvector_dot(&B_vals[0], &b_k, &b_k);

//   bigmatrix_id(H, n);

//   bigvector temp_vec;
//   bigvector_init(&temp_vec, B->rows);

//   // Helper lambda or inline logic for RED(k, l)
//   // RED(k, l): makes |mu[k][l]| <= 0.5 by reducing b_k using b_l
//   #auto_red = ^(u64 ki, u64 li) {
//     double mukl = mu[ki * n + li];
//     if (fabs(mukl) > 0.5) {
//       long long x = (long long)round(mukl);
//       // b_k = b_k - x * b_l
//       bigvector_sub_mul(B, ki, li, x); // Assumes a function to subtract
//       scaled column
//       // H_k = H_k - x * H_l
//       bigmatrix_sub_mul(H, ki, li, x);
//       mu[ki * n + li] -= x;
//       for (u64 j = 0; j < li; j++) {
//         mu[ki * n + j] -= x * mu[li * n + j];
//       }
//     }
//   };

//   while (k <= n) {
//     // Step 2 [Incremental Gram-Schmidt]
//     if (k <= k_max) {
//       // Proceed to Step 3
//     } else {
//       k_max = k;
//       bigmatrix_get_col(&b_k, B, k - 1);
//       bigvector_copy(&b_star[k - 1], &b_k);

//       for (u64 j = 1; j < k; j++) {
//         bignum dot_val;
//         bn_init(&dot_val);
//         bigvector_dot(&dot_val, &b_k, &b_star[j - 1]);

//         double dot_d = bn_to_double(&dot_val);
//         double Bj_d = bn_to_double(&B_vals[j - 1]);
//         mu[(k - 1) * n + (j - 1)] = dot_d / Bj_d;

//         // b*_k = b*_k - mu_{k,j} * b*_j (approximate or exact vector
//         subtraction)
//         // Implementation depends on scalar-vector vector subtraction helpers
//         // ...
//       }
//       bigvector_dot(&B_vals[k - 1], &b_star[k - 1], &b_star[k - 1]);

//       if (bn_is_zero(&B_vals[k - 1])) {
//         // Error: vectors did not form a basis (linearly dependent)
//         return;
//       }
//     }

//     // Step 3 [Test LLL condition]
//     // Execute RED(k, k - 1)
//     // (Inline or helper call for RED calculation)
//     double muk_k1 = mu[(k - 1) * n + (k - 2)];
//     if (fabs(muk_k1) > 0.5) {
//       long long x = (long long)round(muk_k1);
//       // Apply reduction to basis column k and transform matrix H
//       mu[(k - 1) * n + (k - 2)] -= x;
//       for (u64 j = 0; j < k - 2; j++) {
//         mu[(k - 1) * n + j] -= x * mu[(k - 2) * n + j];
//       }
//     }

//     bignum Bk, Bk_prev, threshold_term;
//     bn_init(&Bk); bn_init(&Bk_prev); bn_init(&threshold_term);

//     // Check LLL inequality: B_k < (delta - mu_{k,k-1}^2) * B_{k-1}
//     // Using floating point conversion for norm comparison or exact bignum
//     arithmetic double Bk_d = bn_to_double(&B_vals[k - 1]); double Bk_prev_d =
//     bn_to_double(&B_vals[k - 2]); double rhs = (delta - muk_k1 * muk_k1) *
//     Bk_prev_d;

//     if (Bk_d < rhs) {
//       // SWAP(k)
//       // Swap columns k and k-1 in B and H, update Gram-Schmidt data
//       accordingly
//       // ...
//       k = (k > 2) ? k - 1 : 2;
//     } else {
//       for (long long l = (long long)k - 3; l >= 0; l--) {
//         // Execute RED(k, l)
//         double mukl = mu[(k - 1) * n + l];
//         if (fabs(mukl) > 0.5) {
//           long long x = (long long)round(mukl);
//           mu[(k - 1) * n + l] -= x;
//           for (u64 j = 0; j < (u64)l; j++) {
//             mu[(k - 1) * n + j] -= x * mu[l * n + j];
//           }
//         }
//       }
//       k++;
//     }

//     // Step 4 [Finished?] handled by the while (k <= n) loop condition.
//   }

//   // Cleanup allocations
//   for (u64 i = 0; i < n; i++) {
//     bigvector_free(&b_star[i]);
//     bn_free(&B_vals[i]);
//   }
//   free(b_star);
//   free(B_vals);
//   free(mu);
//   bigvector_free(&b_k);
//   bigvector_free(&temp_vec);
// }
