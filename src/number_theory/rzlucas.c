/*
 * rzlucas.c
 *
 * Lucas sequence computation.
 *
 * This file computes terms of the Lucas sequences U_n(P, Q), V_n(P, Q)
 * (and Q^n) using the standard doubling identities:
 *
 *   U_2n   = U_n V_n
 *   V_2n   = V_n^2 - 2 Q^n
 *   U_2n+1 = (P U_2n + V_2n) / 2 - Q^n
 *   V_2n+1 = P U_2n+1 - 2 Q U_2n
 *
 * Both an exact (big integer) version and a modular version (in the
 * Montgomery domain) are provided.
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

#include "../include/rzlucas.h"

#include <assert.h>
#include <ctype.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>

#include "alloca.h"

rabin_err_t rz_lucas_solve(rz_t* u, rz_t* v, const rz_t* p, const rz_t* q,
                           rz_t* qn, const rz_t* n)
{
  if (u == NULL || v == NULL || p == NULL || q == NULL || qn == NULL ||
      n == NULL)
    return RABIN_ERR_NULL_PTR;

  if (rz_is_eq_i64(n, 0)) {
    rabin_err_t err = rz_set_u64(u, 0);
    if (err == RABIN_SUCCESS) err = rz_set_u64(v, 2);
    if (err == RABIN_SUCCESS) err = rz_set_u64(qn, 1);
    return err;
  }

  if (rz_is_eq_i64(n, 1)) {
    rabin_err_t err = rz_set_u64(u, 1);
    if (err == RABIN_SUCCESS) err = rz_copy(v, p);
    if (err == RABIN_SUCCESS) err = rz_copy(qn, q);
    return err;
  }

  rz_t n_half, u_n, v_n, q_n, a, b, tmp;
  if (rz_init_multi(&n_half, &u_n, &v_n, &tmp, &q_n, &a, &b, NULL) !=
      RABIN_SUCCESS)
    return RABIN_ERR_OUT_OF_MEMORY;

  // rz_divmod_u64 returns the remainder; n_half = n / 2
  rz_divmod_u64(&n_half, n, 2);
  rabin_err_t err = RABIN_SUCCESS;

  if ((err = rz_lucas_solve(&u_n, &v_n, p, q, &q_n, &n_half)) != RABIN_SUCCESS)
    goto cleanup;

  // Precompute used values
  if ((err = rz_mul(&a, &u_n, &v_n)) != RABIN_SUCCESS) goto cleanup;
  if ((err = rz_mul(&b, &v_n, &v_n)) != RABIN_SUCCESS) goto cleanup;

  if (rz_is_even(n)) {
    // U_2n = a
    if ((err = rz_copy(u, &a)) != RABIN_SUCCESS) goto cleanup;
    // V_2n = b - 2 * Q^k
    if ((err = rz_copy(v, &b)) != RABIN_SUCCESS) goto cleanup;
    if ((err = rz_sub(v, v, &q_n)) != RABIN_SUCCESS) goto cleanup;
    if ((err = rz_sub(v, v, &q_n)) != RABIN_SUCCESS) goto cleanup;
    // Q^n = (Q^n/2)^2
    if ((err = rz_mul(qn, &q_n, &q_n)) != RABIN_SUCCESS) goto cleanup;
  } else {
    if ((err = rz_mul(&tmp, p, &a)) != RABIN_SUCCESS) goto cleanup;
    if ((err = rz_add(&tmp, &tmp, &b)) != RABIN_SUCCESS) goto cleanup;
    if ((err = rz_rshift1(&tmp)) != RABIN_SUCCESS) goto cleanup;
    if ((err = rz_sub(u, &tmp, &q_n)) != RABIN_SUCCESS) goto cleanup;

    if ((err = rz_mul(v, p, u)) != RABIN_SUCCESS) goto cleanup;
    if ((err = rz_mul(&tmp, q, &a)) != RABIN_SUCCESS) goto cleanup;
    if ((err = rz_sub(v, v, &tmp)) != RABIN_SUCCESS) goto cleanup;
    if ((err = rz_sub(v, v, &tmp)) != RABIN_SUCCESS) goto cleanup;

    if ((err = rz_mul(qn, &q_n, &q_n)) != RABIN_SUCCESS) goto cleanup;
    if ((err = rz_mul(qn, qn, q)) != RABIN_SUCCESS) goto cleanup;
  }

  err = RABIN_SUCCESS;
cleanup:
  rz_clear(&n_half);
  rz_clear(&u_n);
  rz_clear(&v_n);
  rz_clear(&q_n);
  rz_clear(&a);
  rz_clear(&b);
  rz_clear(&tmp);
  return err;
}

rabin_err_t rz_lucas_solve_mod(rz_t* u, rz_t* v, const rz_t* p, const rz_t* q,
                               rz_t* qn, const rz_t* n, const rz_t* m)
{
  if (u == NULL || v == NULL || p == NULL || q == NULL || qn == NULL ||
      n == NULL || m == NULL)
    return RABIN_ERR_NULL_PTR;

  if (rz_is_eq_i64(n, 0)) {
    rabin_err_t err = rz_set_u64(u, 0);
    if (err == RABIN_SUCCESS) err = rz_set_u64(v, 2);
    if (err == RABIN_SUCCESS) err = rz_set_u64(qn, 1);
    return err;
  }

  if (rz_is_eq_i64(n, 1)) {
    rabin_err_t err = rz_set_u64(u, 1);
    if (err == RABIN_SUCCESS) err = rz_copy(v, p);
    if (err == RABIN_SUCCESS) err = rz_copy(qn, q);
    if (err == RABIN_SUCCESS) err = rz_mod(v, v, m);
    if (err == RABIN_SUCCESS) err = rz_mod(qn, qn, m);
    return err;
  }

  rz_mont_ctx ctx;
  rabin_err_t err = rz_mont_ctx_init(&ctx, m);
  if (err != RABIN_SUCCESS) return err;

  rz_t u_bar, v_bar, qn_bar, p_bar, q_bar, one_bar, a, b, tmp;
  if (rz_init_multi(&u_bar, &v_bar, &qn_bar, &q_bar, &p_bar, &one_bar, &a, &b,
                    &tmp, NULL) != RABIN_SUCCESS) {
    rz_mont_ctx_clear(&ctx);
    return RABIN_ERR_OUT_OF_MEMORY;
  }

  // map constants into Montgomery space
  if ((err = rz_mont_in(&p_bar, p, &ctx)) != RABIN_SUCCESS) goto cleanup;
  if ((err = rz_mont_in(&q_bar, q, &ctx)) != RABIN_SUCCESS) goto cleanup;
  if ((err = rz_copy(&one_bar, &ctx.one_mont)) != RABIN_SUCCESS) goto cleanup;

  // initialize u, v, qn
  if ((err = rz_copy(&u_bar, &one_bar)) != RABIN_SUCCESS) goto cleanup;
  if ((err = rz_copy(&v_bar, &p_bar)) != RABIN_SUCCESS) goto cleanup;
  if ((err = rz_copy(&qn_bar, &q_bar)) != RABIN_SUCCESS) goto cleanup;

  u64 len = (u64)rz_bit_length(n);
#pragma GCC unroll 8
  for (i64 i = (i64)len - 2; i >= 0; i--) {
    // Precompute used values
    // a = U_n * V_n = U_2n
    // b = V_n * V_n
    if ((err = rz_mont_mul(&a, &u_bar, &v_bar, &ctx)) != RABIN_SUCCESS)
      goto cleanup;
    if ((err = rz_mont_mul(&b, &v_bar, &v_bar, &ctx)) != RABIN_SUCCESS)
      goto cleanup;

    if (rz_get_bit(n, i) == 0) {
      // U_2n = a
      if ((err = rz_copy(&u_bar, &a)) != RABIN_SUCCESS) goto cleanup;
      // V_2n = b - 2 * Q^k
      if ((err = rz_sub(&v_bar, &b, &qn_bar)) != RABIN_SUCCESS) goto cleanup;
      if (v_bar.is_neg &&
          (err = rz_add_abs(&v_bar, &v_bar, &ctx.n)) != RABIN_SUCCESS)
        goto cleanup;
      if ((err = rz_sub(&v_bar, &v_bar, &qn_bar)) != RABIN_SUCCESS)
        goto cleanup;
      if (v_bar.is_neg &&
          (err = rz_add_abs(&v_bar, &v_bar, &ctx.n)) != RABIN_SUCCESS)
        goto cleanup;

      // Q^n = (Q^n/2)^2
      if ((err = rz_mont_mul(&qn_bar, &qn_bar, &qn_bar, &ctx)) != RABIN_SUCCESS)
        goto cleanup;
    } else {
      // tmp = P * U_2n + V_2n / 2
      if ((err = rz_mont_mul(&tmp, &p_bar, &a, &ctx)) != RABIN_SUCCESS)
        goto cleanup;
      if ((err = rz_add(&tmp, &tmp, &b)) != RABIN_SUCCESS) goto cleanup;

      // fix for modular arethmetic division by 2
      if (!rz_is_even(&tmp)) {
        if ((err = rz_add(&tmp, &tmp, &ctx.n)) != RABIN_SUCCESS) goto cleanup;
      }
      if ((err = rz_rshift1(&tmp)) != RABIN_SUCCESS) goto cleanup;

      if ((err = rz_copy(&u_bar, &tmp)) != RABIN_SUCCESS) goto cleanup;

      // V_2n+1 = P * U_2n+1 - 2 * Q * U_2n
      if ((err = rz_mont_mul(&v_bar, &p_bar, &u_bar, &ctx)) != RABIN_SUCCESS)
        goto cleanup;
      if ((err = rz_mont_mul(&tmp, &q_bar, &a, &ctx)) != RABIN_SUCCESS)
        goto cleanup;

      // Subtract Q*a twice
      for (int j = 0; j < 2; j++) {
        if ((err = rz_sub(&v_bar, &v_bar, &tmp)) != RABIN_SUCCESS) goto cleanup;
        // this needs to be add_abs. I have no idea why
        if (v_bar.is_neg &&
            (err = rz_add_abs(&v_bar, &v_bar, &ctx.n)) != RABIN_SUCCESS)
          goto cleanup;
      }

      // Qn = Qn^2 * Q
      if ((err = rz_mont_mul(&qn_bar, &qn_bar, &qn_bar, &ctx)) != RABIN_SUCCESS)
        goto cleanup;
      if ((err = rz_mont_mul(&qn_bar, &qn_bar, &q_bar, &ctx)) != RABIN_SUCCESS)
        goto cleanup;
    }
  }

  if ((err = rz_mont_out(u, &u_bar, &ctx)) != RABIN_SUCCESS) goto cleanup;
  if ((err = rz_mont_out(v, &v_bar, &ctx)) != RABIN_SUCCESS) goto cleanup;
  if ((err = rz_mont_out(qn, &qn_bar, &ctx)) != RABIN_SUCCESS) goto cleanup;

  err = RABIN_SUCCESS;
cleanup:
  rz_clear(&a);
  rz_clear(&b);
  rz_clear(&tmp);
  rz_clear(&u_bar);
  rz_clear(&v_bar);
  rz_clear(&qn_bar);
  rz_clear(&p_bar);
  rz_clear(&q_bar);
  rz_clear(&one_bar);
  rz_mont_ctx_clear(&ctx);
  return err;
}

rabin_err_t rz_lucas(rz_t* u, rz_t* v, const rz_t* p, const rz_t* q,
                     const rz_t* n)
{
  if (u == NULL || v == NULL || p == NULL || q == NULL || n == NULL)
    return RABIN_ERR_NULL_PTR;

  rz_t qn;

  rabin_err_t err = rz_init(&qn);
  if (err == RABIN_SUCCESS) {
    err = rz_lucas_solve(u, v, p, q, &qn, n);
  }
  rz_clear(&qn);
  return err;
}

rabin_err_t rz_lucas_mod(rz_t* u, rz_t* v, const rz_t* p, const rz_t* q,
                         const rz_t* n, const rz_t* m)
{
  if (u == NULL || v == NULL || p == NULL || q == NULL || n == NULL ||
      m == NULL)
    return RABIN_ERR_NULL_PTR;

  rz_t qn;

  rabin_err_t err = rz_init(&qn);
  if (err == RABIN_SUCCESS) {
    err = rz_lucas_solve_mod(u, v, p, q, &qn, n, m);
  }
  rz_clear(&qn);
  return err;
}
