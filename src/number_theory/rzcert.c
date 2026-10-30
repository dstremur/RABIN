/*
 * rzcert.c
 *
 * Pocklington primality certificates: construction (rzcert_add_elem),
 * verification (rzcert_verify) and inspection helpers.
 *
 * Copyright (C) 2026 Diego Strebel
 * SPDX-License-Identifier: GPL-3.0-or-later
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 */

#include "../../include/rzcert.h"

#include <stdio.h>
#include <stdlib.h>

#include "../../include/rzmath.h"
#include "../../include/rzprime.h"

rabin_err_t print_pocklington_cert(rzcert_t* cert)
{
  if (cert == NULL) return RABIN_ERR_NULL_PTR;

  printf("Prime N: ");
  rz_print(&cert->N);
  printf("\n");

  if (cert->size == 0) {
    printf("-> Base case (verified via bpsw)\n");
    return RABIN_SUCCESS;
  }

  rabin_err_t err = RABIN_SUCCESS;
  for (u64 i = 0; i < cert->size && err == RABIN_SUCCESS; i++) {
    printf("-> Factor q: ");
    rz_print(&cert->data[i].q);
    printf(" | Base alpha: ");
    rz_print(&cert->data[i].alpha_q);
    printf("\n");

    // Recursively print the certificate for q
    err = print_pocklington_cert(cert->data[i].q_cert);
  }

  return err;
}

rabin_err_t rzcert_clear(rzcert_t* cert)
{
  if (cert == NULL) return RABIN_ERR_NULL_PTR;

  rz_clear(&cert->N);

  if (cert->data) {
    for (u64 i = 0; i < cert->size; i++) {
      rz_clear(&cert->data[i].q);
      rz_clear(&cert->data[i].alpha_q);
      rzcert_clear(cert->data[i].q_cert);
    }
    free(cert->data);
  }
  free(cert);
  return RABIN_SUCCESS;
}

bool rzcert_add_elem(rzcert_t* cert, const rz_t* q, const rz_t* alpha_q,
                     rzcert_t* q_cert)
{
  if (cert == NULL || q == NULL || alpha_q == NULL) return false;

  if (cert->size >= cert->capacity) {
    u64 new_cap = cert->capacity ? cert->capacity * 2 : 4;
    if (new_cap > SIZE_MAX / sizeof(rzcert_elem_t)) return false;

    rzcert_elem_t* nd =
        realloc(cert->data, (size_t)new_cap * sizeof(rzcert_elem_t));
    if (nd == NULL) return false;
    cert->data = nd;
    cert->capacity = new_cap;
  }

  rzcert_elem_t* e = &cert->data[cert->size];
  if (rz_init(&e->q) != RABIN_SUCCESS) return false;
  if (rz_init(&e->alpha_q) != RABIN_SUCCESS) {
    rz_clear(&e->q);
    return false;
  }
  if (rz_copy(&e->q, q) != RABIN_SUCCESS ||
      rz_copy(&e->alpha_q, alpha_q) != RABIN_SUCCESS) {
    rz_clear(&e->q);
    rz_clear(&e->alpha_q);
    return false;
  }

  // Ownership of q_cert transfers to the certificate only on success.
  e->q_cert = q_cert;
  cert->size++;
  return true;
}

bool rzcert_is_base_case(const rzcert_t* cert)
{
  return cert != NULL && cert->size == 0;
}

u64 rzcert_get_depth(const rzcert_t* cert)
{
  if (cert == NULL) return 0;

  u64 max_child = 0;
  if (cert->data != NULL) {
    for (u64 i = 0; i < cert->size; i++) {
      u64 d = rzcert_get_depth(cert->data[i].q_cert);
      if (d > max_child) max_child = d;
    }
  }
  return 1 + max_child;
}

u64 rzcert_get_node_count(const rzcert_t* cert)
{
  if (cert == NULL) return 0;

  u64 count = 1;
  if (cert->data != NULL) {
    for (u64 i = 0; i < cert->size; i++) {
      count += rzcert_get_node_count(cert->data[i].q_cert);
    }
  }
  return count;
}

size_t rzcert_get_byte_size(const rzcert_t* cert)
{
  if (cert == NULL) return 0;

  size_t bytes = sizeof(rzcert_t);
  bytes += (size_t)cert->capacity * sizeof(rzcert_elem_t);
  bytes += (size_t)cert->N.capacity * sizeof(u64);

  if (cert->data != NULL) {
    for (u64 i = 0; i < cert->size; i++) {
      bytes += (size_t)cert->data[i].q.capacity * sizeof(u64);
      bytes += (size_t)cert->data[i].alpha_q.capacity * sizeof(u64);
      bytes += rzcert_get_byte_size(cert->data[i].q_cert);
    }
  }
  return bytes;
}

// Verify one certificate node; see rzcert_verify() for the conditions.
// Recursion terminates because every child certificate proves a factor
// q with q | N - 1, hence q < N: the proved value strictly decreases.
static bool rzcert_verify_node(const rzcert_t* cert, u64 base_case_bound)
{
  if (cert == NULL) return false;

  const rz_t* N = &cert->N;

  if (N->is_neg || rz_is_zero(N)) return false;

  // Leaf: accepted only within the trust bound, via a BPSW test.
  // rz_bpsw() decides small values correctly (including N = 2).
  if (cert->size == 0) {
    return rz_bit_length(N) <= (int)base_case_bound && rz_bpsw(N);
  }

  // Non-leaf N must be odd: no even integer can satisfy the Pocklington
  // conditions below.
  if (N->size == 0 || (N->limbs[0] & 1u) == 0) return false;
  if (cert->data == NULL) return false;

  bool ok = false;
  rz_t n1, one, quo, rem, tmp, y, x, g;
  rz_init_multi(&n1, &one, &quo, &rem, &tmp, &y, &x, &g, NULL);

  if (rz_set_u64(&one, 1) != RABIN_SUCCESS) goto done;
  if (rz_sub(&n1, N, &one) != RABIN_SUCCESS) goto done;

  for (u64 i = 0; i < cert->size; i++) {
    const rzcert_elem_t* e = &cert->data[i];

    // Witness candidates must be positive and at least 2: alpha in {0, 1}
    // trivially satisfies or vacuously fails the conditions below, and a
    // negative alpha is not a valid certificate entry.
    if (e->q.is_neg || rz_is_zero(&e->q) || rz_is_one(&e->q)) goto done;
    if (e->alpha_q.is_neg || rz_is_zero(&e->alpha_q) || rz_is_one(&e->alpha_q))
      goto done;

    // q must divide N - 1 exactly.
    if (rz_divmod(&quo, &rem, &n1, &e->q) != RABIN_SUCCESS) goto done;
    if (!rz_is_zero(&rem)) goto done;

    // Pocklington size condition: q > sqrt(N) - 1, checked exactly as
    // (q + 1)^2 > N to avoid any rounding through a square root.
    if (rz_add_u64(&tmp, &e->q, 1) != RABIN_SUCCESS) goto done;
    if (rz_mul(&tmp, &tmp, &tmp) != RABIN_SUCCESS) goto done;
    if (rz_cmp(&tmp, N) <= 0) goto done;

    // Fermat condition: alpha^(N-1) = 1 mod N.
    if (rz_mod_exp(&y, &e->alpha_q, &n1, N) != RABIN_SUCCESS) goto done;
    if (rz_cmp(&y, &one) != 0) goto done;

    // Order condition: gcd(alpha^((N-1)/q) - 1, N) = 1.
    if (rz_mod_exp(&x, &e->alpha_q, &quo, N) != RABIN_SUCCESS) goto done;
    if (rz_sub(&x, &x, &one) != RABIN_SUCCESS) goto done;
    if (rz_gcd(&g, &x, N) != RABIN_SUCCESS) goto done;
    if (rz_cmp(&g, &one) != 0) goto done;

    // The factor q itself must be proven: child certificate, or BPSW for
    // leaves within the trust bound.
    if (e->q_cert == NULL) {
      if ((u64)rz_bit_length(&e->q) > base_case_bound) goto done;
      if (!rz_bpsw(&e->q)) goto done;
    } else {
      if (rz_cmp(&e->q_cert->N, &e->q) != 0) goto done;
      if (!rzcert_verify_node(e->q_cert, base_case_bound)) goto done;
    }
  }

  ok = true;
done:
  rz_clear_multi(&n1, &one, &quo, &rem, &tmp, &y, &x, &g, NULL);
  return ok;
}

bool rzcert_verify(const rzcert_t* cert, u64 base_case_bound)
{
  return rzcert_verify_node(cert, base_case_bound);
}
