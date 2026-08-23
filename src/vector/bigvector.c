/*
 * bigvector.c
 *
 * Vectors of bignums.
 *
 * This file implements fixed-size and dynamically growing vectors of
 * bignum elements: initialization, freeing, appending, element
 * assignment, component-wise addition and subtraction, dot product,
 * Euclidean norm, and printing.
 *
 * A bigvector is a dynamic array of bignums with a size, a capacity,
 * and a flag marking it as dynamically growable.
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

#include "../../include/bigvector.h"

#include <stdio.h>

/*
 * Initialize a static vector of d zero bignums.
 *
 * Allocates exactly d slots (capacity = size = d) and initializes
 * each to a zero bignum. The vector cannot be grown afterwards.
 *
 * Complexity:
 *   Time: O(d)
 *   Auxiliary memory: O(1)
 *   Output memory: O(d) bignums
 */
void bigvector_init(bigvector* a, u64 d)
{
  a->size = d;
  a->capacity = d;
  a->data = malloc(d * sizeof(bignum));
  a->dynamic = false;

  for (u64 i = 0; i < d; i++) {
    bn_init(&a->data[i]);
  }
}

/*
 * Initialize an empty dynamic vector.
 *
 * The vector starts with size = capacity = 0 and may be grown with
 * bigvector_append().
 *
 * Complexity:
 *   Time: O(1)
 *   Auxiliary memory: O(1)
 *   Output memory: O(1)
 */
void bigvector_init_dynamic(bigvector* a)
{
  bigvector_init(a, 0);
  a->dynamic = true;
}

/*
 * Free all storage of a vector.
 *
 * Complexity:
 *   Time: O(size)
 *   Auxiliary memory: O(1)
 *   Output memory: O(1)
 */
void bigvector_free(bigvector* a)
{
  for (u64 i = 0; i < a->size; i++) {
    bn_free(&a->data[i]);
  }

  free(a->data);
  a->size = 0;
  a->capacity = 0;
  a->dynamic = false;
}

/*
 * Append a copy of a bignum to a dynamic vector.
 *
 * Only allowed on vectors created with bigvector_init_dynamic();
 * otherwise an error is printed and nothing happens. If the capacity
 * is exhausted it is doubled (starting at 4) via realloc and the new
 * slots are initialized to zero bignums.
 *
 * Complexity:
 *   Time: O(n) where n is the size of a in limbs, amortized O(1) for
 *         the growth
 *   Auxiliary memory: O(capacity) during the realloc
 *   Output memory: O(size + 1) bignums
 */
void bigvector_append(bigvector* v, bignum* a)
{
  if (!v->dynamic) {
    perror("no appending on a static vector\n");
    return;
  }

  if (v->size >= v->capacity) {
    u64 new_cap = (v->capacity == 0) ? 4 : v->capacity * 2;
    bignum* new_data = realloc(v->data, new_cap * sizeof(bignum));

    if (!new_data) {
      perror("memory allocation failed\n");
      return;
    }

    v->data = new_data;

    for (u64 i = v->capacity; i < new_cap; i++) {
      bn_init(&v->data[i]);
    }

    v->capacity = new_cap;
  }

  bn_copy(&v->data[v->size], a);
  v->size++;

  return;
}

/*
 * Set element i of a vector to a copy of v.
 *
 * No-op if i is out of range.
 *
 * Complexity:
 *   Time: O(n) where n is the size of v in limbs
 *   Auxiliary memory: O(1)
 *   Output memory: O(n) limbs
 */
void bigvector_set(bigvector* a, bignum* v, u64 i)
{
  if (a->size < i) return;

  bn_copy(&a->data[i], v);
}

/*
 * Component-wise addition of two vectors: r = a + b.
 *
 * Let d = a->size.
 *
 * No-op if the sizes differ. r must have capacity for at least d
 * elements.
 *
 * Complexity:
 *   Time: O(d * n) where n is the size of the elements in limbs
 *   Auxiliary memory: O(1)
 *   Output memory: O(d) bignums
 */
void bigvector_add(bigvector* r, const bigvector* a, const bigvector* b)
{
  if (a->size != b->size) return;

  for (u64 i = 0; i < a->size; i++) {
    bn_add(&r->data[i], &a->data[i], &b->data[i]);
  }
}

/*
 * Component-wise subtraction of two vectors: r = a - b.
 *
 * Let d = a->size.
 *
 * No-op if the sizes differ. r must have capacity for at least d
 * elements.
 *
 * Complexity:
 *   Time: O(d * n) where n is the size of the elements in limbs
 *   Auxiliary memory: O(1)
 *   Output memory: O(d) bignums
 */
void bigvector_sub(bigvector* r, const bigvector* a, const bigvector* b)
{
  if (a->size != b->size) return;

  for (u64 i = 0; i < a->size; i++) {
    bn_sub(&r->data[i], &a->data[i], &b->data[i]);
  }
}

/*
 * Euclidean norm of a vector: r = sqrt(sum_i a_i^2).
 *
 * Let d = a->size.
 *
 * Only defined for static vectors (an error is printed for dynamic
 * ones). Computes the dot product of a with itself and takes the
 * integer square root.
 *
 * Complexity:
 *   Time: O(d * n^2) where n is the size of the elements in limbs,
 *         plus O(n^2 log n) for the square root
 *   Auxiliary memory: O(n) limbs
 *   Output memory: O(n) limbs
 */
void bigvector_norm(bignum* r, const bigvector* a)
{
  if (a->dynamic) {
    perror("Not possible for dynamic arrays\n");
  }
  bignum temp;
  bn_init(&temp);

  bigvector_dot(&temp, a, a);
  bn_isqrt(r, &temp);

  bn_free(&temp);
}

/*
 * Dot product of two vectors: r = sum_i a_i * b_i.
 *
 * Let d = a->size.
 *
 * Only defined for static vectors (an error is printed for dynamic
 * ones); no-op if the sizes differ.
 *
 * Complexity:
 *   Time: O(d * n^2) where n is the size of the elements in limbs
 *   Auxiliary memory: O(n) limbs
 *   Output memory: O(n) limbs
 */
void bigvector_dot(bignum* r, const bigvector* a, const bigvector* b)
{
  if (a->dynamic) {
    perror("Not possible for dynamic arrays\n");
  }
  if (a->size != b->size) return;
  bignum temp;
  bn_init(&temp);

  for (u64 i = 0; i < a->size; i++) {
    bn_mul(&temp, &a->data[i], &b->data[i]);
    bn_add(r, r, &temp);
  }

  bn_free(&temp);
}

/*
 * Print a vector to stdout as [a_0, a_1, ...].
 *
 * Complexity:
 *   Time: O(d * n) where d is the size and n the element size
 *   Auxiliary memory: O(1)
 *   Output memory: O(d * n) characters written
 */
void bigvector_print(bigvector* a)
{
  printf("[");

  for (u64 i = 0; i < a->size; i++) {
    bn_print(&a->data[i]);

    if (i + 1 < a->size) printf(", ");
  }

  printf("]");
}

/*
 * Print a vector to stdout followed by a newline.
 *
 * Complexity:
 *   Time: O(d * n) where d is the size and n the element size
 *   Auxiliary memory: O(1)
 *   Output memory: O(d * n) characters written
 */
void bigvector_println(bigvector* a)
{
  bigvector_print(a);
  printf("\n");
}
