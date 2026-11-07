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
#include "rzlucas_internal.h"

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
  rz_clear_multi(&n_half, &u_n, &v_n, &q_n, &a, &b, &tmp, NULL);
  return err;
}

/*
 * One left-to-right doubling step (U_k, V_k, Q^k) -> (U_2k, V_2k, Q^2k) or
 * (U_2k+1, V_2k+1, Q^2k+1) in the Montgomery domain. All values are
 * canonical residues in [0, n); a, b, tmp are scratch registers.
 */
static rabin_err_t lucas_step_mont(rz_t* u_bar, rz_t* v_bar, rz_t* qn_bar,
                                   const rz_t* p_bar, const rz_t* q_bar,
                                   rz_t* a, rz_t* b, rz_t* tmp, int bit,
                                   rz_mont_ctx* ctx)
{
  rabin_err_t err = RABIN_SUCCESS;

  /* a = U_k * V_k = U_2k,  b = V_k^2 */
  TRY(rz_mont_mul(a, u_bar, v_bar, ctx));
  TRY(rz_mont_mul(b, v_bar, v_bar, ctx));

  if (bit == 0) {
    /* U_2k = a */
    TRY(rz_copy(u_bar, a));

    /*
     * V_2k = V_k^2 - 2 * Q^k
     *
     * Both operands are canonical residues. After subtraction,
     * adding m to a negative result gives a canonical residue.
     */
    TRY(rz_sub(v_bar, b, qn_bar));
    if (v_bar->is_neg) TRY(rz_add(v_bar, v_bar, &ctx->n));

    TRY(rz_sub(v_bar, v_bar, qn_bar));
    if (v_bar->is_neg) TRY(rz_add(v_bar, v_bar, &ctx->n));

    /* Q^2k = (Q^k)^2 */
    TRY(rz_mont_mul(qn_bar, qn_bar, qn_bar, ctx));

  } else {
    /*
     * U_2k+1 = (P * U_2k + V_k^2) / 2 - Q^k
     *
     * First reduce the sum to [0, m). Since m is odd, adding
     * m when the sum is odd makes integer division by two valid
     * modulo m.
     */
    TRY(rz_mont_mul(tmp, p_bar, a, ctx));
    TRY(rz_add(tmp, tmp, b));
    TRY(rz_mod(tmp, tmp, &ctx->n));

    if (!rz_is_even(tmp)) TRY(rz_add(tmp, tmp, &ctx->n));

    TRY(rz_rshift1(tmp));

    /* Subtract Q^k to obtain U_2k+1. */
    TRY(rz_sub(u_bar, tmp, qn_bar));
    if (u_bar->is_neg) TRY(rz_add(u_bar, u_bar, &ctx->n));

    /*
     * V_2k+1 = P * U_2k+1 - 2 * Q * U_2k
     */
    TRY(rz_mont_mul(v_bar, p_bar, u_bar, ctx));
    TRY(rz_mont_mul(tmp, q_bar, a, ctx));

    /* Subtract Q * U_2k twice. */
    for (int j = 0; j < 2; j++) {
      TRY(rz_sub(v_bar, v_bar, tmp));
      if (v_bar->is_neg) TRY(rz_add(v_bar, v_bar, &ctx->n));
    }

    /* Q^2k+1 = (Q^k)^2 * Q */
    TRY(rz_mont_mul(qn_bar, qn_bar, qn_bar, ctx));
    TRY(rz_mont_mul(qn_bar, qn_bar, q_bar, ctx));
  }

cleanup:
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
    TRY(lucas_step_mont(&u_bar, &v_bar, &qn_bar, &p_bar, &q_bar, &a, &b, &tmp,
                        rz_get_bit(n, i), &ctx));
  }

  /* Convert the results out of Montgomery space. */
  TRY(rz_mont_out(u, &u_bar, &ctx));
  TRY(rz_mont_out(v, &v_bar, &ctx));
  TRY(rz_mont_out(qn, &qn_bar, &ctx));

  err = RABIN_SUCCESS;

cleanup:
  rz_clear_multi(&a, &b, &tmp, &u_bar, &v_bar, &qn_bar, &p_bar, &q_bar,
                 &one_bar, NULL);
  rz_mont_ctx_clear(&ctx);
  return err;
}

/*
 * Strong Lucas test on an existing Montgomery context, never leaving the
 * Montgomery domain. Requires n + 1 = d * 2^s with d >= 1 odd and
 * ctx->n == n. Zero tests run directly on the Montgomery residues: the
 * results of rz_mont_mul() are canonical residues in [0, n) and
 * 0 * R = 0, so U == 0 (mod n) iff u_bar == 0.
 */
rabin_err_t rz_strong_lucas_mont_ctx(bool* res, rz_mont_ctx* ctx, const rz_t* n,
                                     const rz_t* P, const rz_t* Q,
                                     const rz_t* d, u64 s)
{
  if (res == NULL || ctx == NULL || n == NULL || P == NULL || Q == NULL ||
      d == NULL)
    return RABIN_ERR_NULL_PTR;
  // sequence index must be a positive odd integer
  if (d->is_neg || rz_is_zero(d)) return RABIN_ERR_INVALID_ARG;
  // the context must match the modulus under test
  if (ctx->n.size == 0 || rz_cmp(&ctx->n, n) != 0) return RABIN_ERR_INVALID_ARG;

  *res = false;

  rz_t u_bar, v_bar, qn_bar, p_bar, q_bar, one_bar, a, b, tmp;
  rabin_err_t err = rz_init_multi(&u_bar, &v_bar, &qn_bar, &q_bar, &p_bar,
                                  &one_bar, &a, &b, &tmp, NULL);
  if (err != RABIN_SUCCESS) return err;

  /* Map constants into Montgomery space (P and Q may be negative). */
  TRY(rz_mont_in(&p_bar, P, ctx));
  TRY(rz_mont_in(&q_bar, Q, ctx));
  TRY(rz_copy(&one_bar, &ctx->one_mont));

  /* Initialize U_1, V_1 and Q^1 in Montgomery space; d >= 1. */
  TRY(rz_copy(&u_bar, &one_bar));
  TRY(rz_copy(&v_bar, &p_bar));
  TRY(rz_copy(&qn_bar, &q_bar));

  /* Double-and-add up to the odd cofactor d. */
#pragma GCC unroll 8
  for (i64 i = (i64)rz_bit_length(d) - 2; i >= 0; i--) {
    TRY(lucas_step_mont(&u_bar, &v_bar, &qn_bar, &p_bar, &q_bar, &a, &b, &tmp,
                        rz_get_bit(d, i), ctx));
  }

  if (rz_is_zero(&u_bar) || rz_is_zero(&v_bar)) {
    *res = true;
    goto cleanup;
  }

  /*
   * Doubling chain V_(d * 2^r), r = 1..s-1, entered directly from the
   * d-loop without leaving Montgomery space. U is not needed here, so
   * each step costs two Montgomery multiplications.
   */
  for (u64 r = 1; r < s; r++) {
    /* b = V_k^2 */
    TRY(rz_mont_mul(&b, &v_bar, &v_bar, ctx));

    /* V_2k = V_k^2 - 2 * Q^k */
    TRY(rz_sub(&v_bar, &b, &qn_bar));
    if (v_bar.is_neg) TRY(rz_add(&v_bar, &v_bar, &ctx->n));

    TRY(rz_sub(&v_bar, &v_bar, &qn_bar));
    if (v_bar.is_neg) TRY(rz_add(&v_bar, &v_bar, &ctx->n));

    /* Q^2k = (Q^k)^2 */
    TRY(rz_mont_mul(&qn_bar, &qn_bar, &qn_bar, ctx));

    if (rz_is_zero(&v_bar)) {
      *res = true;
      goto cleanup;
    }
  }

cleanup:
  rz_clear_multi(&a, &b, &tmp, &u_bar, &v_bar, &qn_bar, &p_bar, &q_bar,
                 &one_bar, NULL);
  return err;
}

rabin_err_t rz_strong_lucas_solve_mod(bool* res, const rz_t* n, const rz_t* P,
                                      const rz_t* Q, const rz_t* d, u64 s)
{
  if (res == NULL || n == NULL || P == NULL || Q == NULL || d == NULL)
    return RABIN_ERR_NULL_PTR;

  rz_mont_ctx ctx = {0};
  rabin_err_t err = rz_mont_ctx_init(&ctx, n);
  if (err != RABIN_SUCCESS) return err;

  err = rz_strong_lucas_mont_ctx(res, &ctx, n, P, Q, d, s);
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
