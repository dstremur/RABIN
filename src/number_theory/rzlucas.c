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
 *   U_2n+1 = (P U_2n + V_n^2) / 2 - Q^n   (equivalently (P U_2n + V_2n)/2)
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
  // sequence index must be a non-negative integer
  if (n->is_neg) return RABIN_ERR_INVALID_ARG;

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

  TRY(rz_lucas_solve(&u_n, &v_n, p, q, &q_n, &n_half));

  // Precompute used values
  TRY(rz_mul(&a, &u_n, &v_n));
  TRY(rz_mul(&b, &v_n, &v_n));

  if (rz_is_even(n)) {
    // U_2n = a
    TRY(rz_copy(u, &a));
    // V_2n = b - 2 * Q^k
    TRY(rz_copy(v, &b));
    TRY(rz_sub(v, v, &q_n));
    TRY(rz_sub(v, v, &q_n));
    // Q^n = (Q^n/2)^2
    TRY(rz_mul(qn, &q_n, &q_n));
  } else {
    TRY(rz_mul(&tmp, p, &a));
    TRY(rz_add(&tmp, &tmp, &b));
    TRY(rz_rshift1(&tmp));
    TRY(rz_sub(u, &tmp, &q_n));

    TRY(rz_mul(v, p, u));
    TRY(rz_mul(&tmp, q, &a));
    TRY(rz_sub(v, v, &tmp));
    TRY(rz_sub(v, v, &tmp));

    TRY(rz_mul(qn, &q_n, &q_n));
    TRY(rz_mul(qn, qn, q));
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
  // sequence index must be a non-negative integer
  if (n->is_neg) return RABIN_ERR_INVALID_ARG;

  /*
   * The modulus must be valid for Montgomery arithmetic.
   * Initialize the context before handling n == 0 or n == 1.
   */
  rz_mont_ctx ctx;
  rabin_err_t err = rz_mont_ctx_init(&ctx, m);
  if (err != RABIN_SUCCESS) return err;

  /*
   * U_0 = 0, V_0 = 2, Q^0 = 1
   * Reduce every output modulo m.
   */
  if (rz_is_eq_i64(n, 0)) {
    err = rz_set_u64(u, 0);
    if (err == RABIN_SUCCESS) err = rz_set_u64(v, 2);
    if (err == RABIN_SUCCESS) err = rz_set_u64(qn, 1);

    if (err == RABIN_SUCCESS) err = rz_mod(u, u, m);
    if (err == RABIN_SUCCESS) err = rz_mod(v, v, m);
    if (err == RABIN_SUCCESS) err = rz_mod(qn, qn, m);

    rz_mont_ctx_clear(&ctx);
    return err;
  }

  /*
   * U_1 = 1, V_1 = P, Q^1 = Q
   * Reduce every output modulo m.
   */
  if (rz_is_eq_i64(n, 1)) {
    err = rz_set_u64(u, 1);
    if (err == RABIN_SUCCESS) err = rz_copy(v, p);
    if (err == RABIN_SUCCESS) err = rz_copy(qn, q);

    if (err == RABIN_SUCCESS) err = rz_mod(u, u, m);
    if (err == RABIN_SUCCESS) err = rz_mod(v, v, m);
    if (err == RABIN_SUCCESS) err = rz_mod(qn, qn, m);

    rz_mont_ctx_clear(&ctx);
    return err;
  }

  rz_t u_bar, v_bar, qn_bar, p_bar, q_bar;
  rz_t one_bar, a, b, tmp;

  err = rz_init_multi(&u_bar, &v_bar, &qn_bar, &q_bar, &p_bar, &one_bar, &a, &b,
                      &tmp, NULL);

  if (err != RABIN_SUCCESS) {
    rz_mont_ctx_clear(&ctx);
    return err;
  }

  /* Map constants into Montgomery space. */
  TRY(rz_mont_in(&p_bar, p, &ctx));
  TRY(rz_mont_in(&q_bar, q, &ctx));
  TRY(rz_copy(&one_bar, &ctx.one_mont));

  /* Initialize U_1, V_1 and Q^1 in Montgomery space. */
  TRY(rz_copy(&u_bar, &one_bar));
  TRY(rz_copy(&v_bar, &p_bar));
  TRY(rz_copy(&qn_bar, &q_bar));

  u64 len = (u64)rz_bit_length(n);

#pragma GCC unroll 8
  for (i64 i = (i64)len - 2; i >= 0; i--) {
    /*
     * a = U_k * V_k = U_(2k)
     * b = V_k^2
     */
    TRY(rz_mont_mul(&a, &u_bar, &v_bar, &ctx));
    TRY(rz_mont_mul(&b, &v_bar, &v_bar, &ctx));

    if (rz_get_bit(n, i) == 0) {
      /*
       * U_(2k) = a
       */
      TRY(rz_copy(&u_bar, &a));

      /*
       * V_(2k) = V_k^2 - 2 * Q^k
       *
       * Both operands are canonical residues. After subtraction,
       * adding m to a negative result gives a canonical residue.
       */
      TRY(rz_sub(&v_bar, &b, &qn_bar));
      if (v_bar.is_neg) TRY(rz_add(&v_bar, &v_bar, &ctx.n));

      TRY(rz_sub(&v_bar, &v_bar, &qn_bar));
      if (v_bar.is_neg) TRY(rz_add(&v_bar, &v_bar, &ctx.n));

      /* Q^(2k) = (Q^k)^2 */
      TRY(rz_mont_mul(&qn_bar, &qn_bar, &qn_bar, &ctx));

    } else {
      /*
       * U_(2k+1) = (P * U_(2k) + V_k^2) / 2 - Q^k
       *
       * First reduce the sum to [0, m). Since m is odd, adding
       * m when the sum is odd makes integer division by two valid
       * modulo m.
       */
      TRY(rz_mont_mul(&tmp, &p_bar, &a, &ctx));
      TRY(rz_add(&tmp, &tmp, &b));
      TRY(rz_mod(&tmp, &tmp, &ctx.n));

      if (!rz_is_even(&tmp)) TRY(rz_add(&tmp, &tmp, &ctx.n));

      TRY(rz_rshift1(&tmp));

      /* Subtract Q^k to obtain U_(2k+1). */
      TRY(rz_sub(&u_bar, &tmp, &qn_bar));
      if (u_bar.is_neg) TRY(rz_add(&u_bar, &u_bar, &ctx.n));

      /*
       * V_(2k+1) = P * U_(2k+1) - 2 * Q * U_(2k)
       */
      TRY(rz_mont_mul(&v_bar, &p_bar, &u_bar, &ctx));
      TRY(rz_mont_mul(&tmp, &q_bar, &a, &ctx));

      /* Subtract Q * U_(2k) twice. */
      for (int j = 0; j < 2; j++) {
        TRY(rz_sub(&v_bar, &v_bar, &tmp));
        if (v_bar.is_neg) TRY(rz_add(&v_bar, &v_bar, &ctx.n));
      }

      /* Q^(2k+1) = (Q^k)^2 * Q */
      TRY(rz_mont_mul(&qn_bar, &qn_bar, &qn_bar, &ctx));
      TRY(rz_mont_mul(&qn_bar, &qn_bar, &q_bar, &ctx));
    }
  }

  /* Convert the results out of Montgomery space. */
  TRY(rz_mont_out(u, &u_bar, &ctx));
  TRY(rz_mont_out(v, &v_bar, &ctx));
  TRY(rz_mont_out(qn, &qn_bar, &ctx));

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
