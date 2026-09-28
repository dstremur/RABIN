/*
 * rvec.c
 *
 * Vectors of bignums.
 *
 * This file implements fixed-size and dynamically growing vectors of
 * rz_t elements: initialization, freeing, appending, element
 * assignment, component-wise addition and subtraction, dot product,
 * Euclidean norm, and printing.
 *
 * A rvec_t is a dynamic array of bignums with a size, a capacity,
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

#include "../../include/rvec.h"

#include <stdio.h>

#include "../../include/rmat.h"

rabin_err_t rvec_init(rvec_t* a, u64 d)
{
  if (a == NULL) return RABIN_ERR_NULL_PTR;

  // idempotent init: release any previous allocation
  if (a->data != NULL) {
    rvec_clear(a);
  }

  a->size = d;
  a->capacity = d;
  a->dynamic = false;

  if (d == 0) {
    a->data = NULL;
    return RABIN_SUCCESS;
  }

  a->data = malloc(d * sizeof(rz_t));
  if (a->data == NULL) return RABIN_ERR_OUT_OF_MEMORY;

  for (u64 i = 0; i < d; i++) {
    rz_init(&a->data[i]);
  }
  return RABIN_SUCCESS;
}

rabin_err_t rvec_init_dynamic(rvec_t* a)
{
  if (a == NULL) return RABIN_ERR_NULL_PTR;

  rabin_err_t err = rvec_init(a, 0);
  if (err != RABIN_SUCCESS) return err;

  a->dynamic = true;
  return RABIN_SUCCESS;
}

rabin_err_t rvec_clear(rvec_t* a)
{
  if (a == NULL) return RABIN_ERR_NULL_PTR;

  for (u64 i = 0; i < a->size; i++) {
    rz_clear(&a->data[i]);
  }

  free(a->data);
  a->data = NULL;
  a->size = 0;
  a->capacity = 0;
  a->dynamic = false;
  return RABIN_SUCCESS;
}

rabin_err_t rvec_append(rvec_t* v, const rz_t* a)
{
  if (v == NULL || a == NULL) return RABIN_ERR_NULL_PTR;

  if (!v->dynamic) {
    return RABIN_ERR_INVALID_ARG;
  }

  if (v->size >= v->capacity) {
    u64 new_cap = (v->capacity == 0) ? 4 : v->capacity * 2;
    if (v->capacity != 0 && new_cap < v->capacity)  // doubling overflowed
      return RABIN_ERR_OVERFLOW;
    if (new_cap > SIZE_MAX / sizeof(rz_t)) return RABIN_ERR_OVERFLOW;
    rz_t* new_data = realloc(v->data, new_cap * sizeof(rz_t));

    if (new_data == NULL) {
      return RABIN_ERR_OUT_OF_MEMORY;
    }

    v->data = new_data;

    for (u64 i = v->capacity; i < new_cap; i++) {
      rz_init(&v->data[i]);
    }

    v->capacity = new_cap;
  }

  rabin_err_t err = rz_copy(&v->data[v->size], a);
  if (err != RABIN_SUCCESS) return err;
  v->size++;

  return RABIN_SUCCESS;
}

rabin_err_t rvec_copy(rvec_t* a, const rvec_t* b)
{
  if (a == NULL || b == NULL) return RABIN_ERR_NULL_PTR;

  // check for same size
  if (a->size != b->size) return RABIN_ERR_MATRIX_DIM;

  for (u64 i = 0; i < a->size; i++) {
    rabin_err_t err = rz_copy(&a->data[i], &b->data[i]);
    if (err != RABIN_SUCCESS) return err;
  }
  return RABIN_SUCCESS;
}

rabin_err_t rvec_set(rvec_t* a, const rz_t* v, u64 i)
{
  if (a == NULL || v == NULL) return RABIN_ERR_NULL_PTR;

  if (i >= a->size) return RABIN_ERR_INVALID_ARG;

  return rz_copy(&a->data[i], v);
}

rabin_err_t rvec_add(rvec_t* r, const rvec_t* a, const rvec_t* b)
{
  if (r == NULL || a == NULL || b == NULL) return RABIN_ERR_NULL_PTR;

  if (a->size != b->size) return RABIN_ERR_MATRIX_DIM;
  if (r->size < a->size) return RABIN_ERR_MATRIX_DIM;

  for (u64 i = 0; i < a->size; i++) {
    rabin_err_t err = rz_add(&r->data[i], &a->data[i], &b->data[i]);
    if (err != RABIN_SUCCESS) return err;
  }
  return RABIN_SUCCESS;
}

rabin_err_t rvec_sub(rvec_t* r, const rvec_t* a, const rvec_t* b)
{
  if (r == NULL || a == NULL || b == NULL) return RABIN_ERR_NULL_PTR;

  if (a->size != b->size) return RABIN_ERR_MATRIX_DIM;
  if (r->size < a->size) return RABIN_ERR_MATRIX_DIM;

  for (u64 i = 0; i < a->size; i++) {
    rabin_err_t err = rz_sub(&r->data[i], &a->data[i], &b->data[i]);
    if (err != RABIN_SUCCESS) return err;
  }
  return RABIN_SUCCESS;
}

rabin_err_t rvec_norm(rz_t* r, const rvec_t* a)
{
  if (r == NULL || a == NULL) return RABIN_ERR_NULL_PTR;

  if (a->dynamic) {
    return RABIN_ERR_INVALID_ARG;
  }

  rz_t temp;
  rz_init(&temp);

  rabin_err_t err = rvec_dot(&temp, a, a);
  if (err != RABIN_SUCCESS) goto out;
  if ((err = rz_isqrt(r, &temp)) != RABIN_SUCCESS) goto out;

  err = RABIN_SUCCESS;
out:
  rz_clear(&temp);
  return err;
}

rabin_err_t rvec_neg(rvec_t* r, const rvec_t* a)
{
  if (r == NULL || a == NULL) {
    return RABIN_ERR_NULL_PTR;
  }

  if (r != a) {
    rabin_err_t err = rvec_copy(r, a);
    if (err != RABIN_SUCCESS) return err;
  }

  for (u64 i = 0; i < a->size; i++) {
    rabin_err_t err = rz_neg(&r->data[i], &a->data[i]);
    if (err != RABIN_SUCCESS) return err;
  }
  return RABIN_SUCCESS;
}

rabin_err_t rvec_dot(rz_t* r, const rvec_t* a, const rvec_t* b)
{
  if (r == NULL || a == NULL || b == NULL) return RABIN_ERR_NULL_PTR;

  if (a->dynamic) {
    return RABIN_ERR_INVALID_ARG;
  }
  if (a->size != b->size) return RABIN_ERR_MATRIX_DIM;

  rz_t temp;
  rz_init(&temp);
  rabin_err_t err;

  for (u64 i = 0; i < a->size; i++) {
    if ((err = rz_mul(&temp, &a->data[i], &b->data[i])) != RABIN_SUCCESS)
      goto out;
    if ((err = rz_add(r, r, &temp)) != RABIN_SUCCESS) goto out;
  }

  err = RABIN_SUCCESS;
out:
  rz_clear(&temp);
  return err;
}

rabin_err_t rvec_print(rvec_t* a)
{
  if (a == NULL) return RABIN_ERR_NULL_PTR;

  printf("[");

  for (u64 i = 0; i < a->size; i++) {
    rz_print(&a->data[i]);

    if (i + 1 < a->size) printf(", ");
  }

  printf("]");
  return RABIN_SUCCESS;
}

rabin_err_t rvec_println(rvec_t* a)
{
  if (a == NULL) return RABIN_ERR_NULL_PTR;

  rabin_err_t err = rvec_print(a);
  if (err != RABIN_SUCCESS) return err;
  printf("\n");
  return RABIN_SUCCESS;
}
