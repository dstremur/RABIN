/*
 * biglucas.c
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

#include <assert.h>
#include <ctype.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>

#include "../include/bignum.h"
#include "alloca.h"

/*
 * Compute (U_n, V_n, Q^n) exactly, by recursion on the binary
 * expansion of n.
 *
 * Let n_l = n->size, measured in 64-bit limbs.
 *
 * Base cases: n = 0 gives (U, V, Q^n) = (0, 2, 1); n = 1 gives
 * (1, P, Q). The recursion halves n and applies the doubling
 * identities above, so the sequence of values is built from the
 * halved result.
 *
 * Complexity:
 *   Time: O(n_l^2 log n) - O(log n) levels of recursion, each with a
 *         few multiplications of up to n_l-limb values
 *   Auxiliary memory: O(n_l) limbs per recursion level, O(n_l log n)
 *                     total on the stack of temporaries
 *   Output memory: O(n_l) limbs per result
 */
void bn_lucas_solve(bignum* u, bignum* v, const bignum* p, const bignum* q,
                    bignum* qn, const bignum* n)
{
  if (bn_is_eq_i64(n, 0)) {
    bn_set_u64(u, 0);
    bn_set_u64(v, 2);
    bn_set_u64(qn, 1);
    return;
  }

  if (bn_is_eq_i64(n, 1)) {
    bn_set_u64(u, 1);
    bn_copy(v, p);
    bn_copy(qn, q);
    return;
  }

  bignum n_half, u_n, v_n, q_n, a, b, tmp;
  bn_init_multi(&n_half, &u_n, &v_n, &tmp, &q_n, &a, &b, NULL);

  bn_divmod_u64(&n_half, n, 2);

  bn_lucas_solve(&u_n, &v_n, p, q, &q_n, &n_half);

  // Precompute used values
  bn_mul(&a, &u_n, &v_n);
  bn_mul(&b, &v_n, &v_n);

  if (bn_is_even(n)) {
    // U_2n = a
    bn_copy(u, &a);
    // V_2n = b - 2 * Q^k
    bn_copy(v, &b);
    bn_sub(v, v, &q_n);
    bn_sub(v, v, &q_n);
    // Q^n = (Q^n/2)^2
    bn_mul(qn, &q_n, &q_n);
  } else {
    bn_mul(&tmp, p, &a);
    bn_add(&tmp, &tmp, &b);
    bn_rshift1(&tmp);
    bn_sub(u, &tmp, &q_n);

    bn_mul(v, p, u);
    bn_mul(&tmp, q, &a);
    bn_sub(v, v, &tmp);
    bn_sub(v, v, &tmp);

    bn_mul(qn, &q_n, &q_n);
    bn_mul(qn, qn, q);
  }

  bn_free(&n_half);
  bn_free(&u_n);
  bn_free(&v_n);
  bn_free(&q_n);
  bn_free(&a);
  bn_free(&b);
  bn_free(&tmp);
}

/*
 * Compute (U_n mod m, V_n mod m, Q^n mod m) iteratively over the bits
 * of n, in the Montgomery domain.
 *
 * Let n_l = m->size, measured in 64-bit limbs.
 *
 * Same doubling identities as bn_lucas_solve(), but applied
 * left-to-right over the bits of n (starting from (U_1, V_1, Q^1) =
 * (1, P, Q)) with all multiplications as Montgomery multiplications
 * modulo m. Subtractions that would go negative are brought back into
 * [0, m) by adding m. Division by 2 is done by adding m when the
 * value is odd, then shifting.
 *
 * The results are converted out of the Montgomery domain before
 * return.
 *
 * Complexity:
 *   Time: O(n_l^3) for the context initialization, then O(log n *
 *         n_l^2) for the bit loop
 *   Auxiliary memory: O(n_l) limbs for the context and temporaries
 *   Output memory: O(n_l) limbs per result
 */
void bn_lucas_solve_mod(bignum* u, bignum* v, const bignum* p, const bignum* q,
                        bignum* qn, const bignum* n, const bignum* m)
{
  if (bn_is_eq_i64(n, 0)) {
    bn_set_u64(u, 0);
    bn_set_u64(v, 2);
    bn_set_u64(qn, 1);
    return;
  }

  if (bn_is_eq_i64(n, 1)) {
    bn_set_u64(u, 1);
    bn_copy(v, p);
    bn_copy(qn, q);

    bn_mod(v, v, m);
    bn_mod(qn, qn, m);
    return;
  }

  bn_mont_ctx ctx;
  bn_mont_ctx_init(&ctx, m);
  bignum u_bar, v_bar, qn_bar, p_bar, q_bar, one_bar, a, b, tmp;
  bn_init_multi(&u_bar, &v_bar, &qn_bar, &q_bar, &p_bar, &one_bar, &a, &b, &tmp,
                NULL);

  // map constants into Montgomery space
  bn_mont_in(&p_bar, p, &ctx);
  bn_mont_in(&q_bar, q, &ctx);
  bn_copy(&one_bar, &ctx.one_mont);

  // initialize u, v, qn
  bn_copy(&u_bar, &one_bar);
  bn_copy(&v_bar, &p_bar);
  bn_copy(&qn_bar, &q_bar);

  u64 len = bn_bit_length(n);
#pragma GCC unroll 8
  for (i64 i = len - 2; i >= 0; i--) {
    // Precompute used values
    // a = U_n * V_n = U_2n
    // b = V_n * V_n
    bn_mont_mul(&a, &u_bar, &v_bar, &ctx);
    bn_mont_mul(&b, &v_bar, &v_bar, &ctx);

    if (bn_get_bit(n, i) == 0) {
      // U_2n = a
      bn_copy(&u_bar, &a);
      // V_2n = b - 2 * Q^k
      bn_sub(&v_bar, &b, &qn_bar);
      if (v_bar.is_neg) bn_add_abs(&v_bar, &v_bar, &ctx.n);
      bn_sub(&v_bar, &v_bar, &qn_bar);
      if (v_bar.is_neg) bn_add_abs(&v_bar, &v_bar, &ctx.n);

      // Q^n = (Q^n/2)^2
      bn_mont_mul(&qn_bar, &qn_bar, &qn_bar, &ctx);
    } else {
      // tmp = P * U_2n + V_2n / 2
      bn_mont_mul(&tmp, &p_bar, &a, &ctx);
      bn_add(&tmp, &tmp, &b);

      // fix for modular arethmetic division by 2
      if (!bn_is_even(&tmp)) {
        bn_add(&tmp, &tmp, &ctx.n);
      }
      bn_rshift1(&tmp);

      bn_copy(&u_bar, &tmp);

      // V_2n+1 = P * U_2n+1 - 2 * Q * U_2n
      bn_mont_mul(&v_bar, &p_bar, &u_bar, &ctx);
      bn_mont_mul(&tmp, &q_bar, &a, &ctx);

      // Subtract Q*a twice
      for (int j = 0; j < 2; j++) {
        bn_sub(&v_bar, &v_bar, &tmp);
        // this needs to be add_abs. I have no idea why
        if (v_bar.is_neg) bn_add_abs(&v_bar, &v_bar, &ctx.n);
      }

      // Qn = Qn^2 * Q
      bn_mont_mul(&qn_bar, &qn_bar, &qn_bar, &ctx);
      bn_mont_mul(&qn_bar, &qn_bar, &q_bar, &ctx);
    }
  }

  bn_mont_out(u, &u_bar, &ctx);
  bn_mont_out(v, &v_bar, &ctx);
  bn_mont_out(qn, &qn_bar, &ctx);

  bn_free(&a);
  bn_free(&b);
  bn_free(&tmp);
  bn_free(&u_bar);
  bn_free(&v_bar);
  bn_free(&qn_bar);
  bn_free(&p_bar);
  bn_free(&q_bar);
  bn_free(&one_bar);
  bn_mont_ctx_free(&ctx);
}

/*
 * Compute the n-th terms of the Lucas sequences U_n(P, Q) and
 * V_n(P, Q) exactly.
 *
 * Let n_l = n->size, measured in 64-bit limbs.
 *
 * Thin wrapper around bn_lucas_solve() that discards Q^n.
 *
 * Complexity:
 *   Time: O(n_l^2 log n), see bn_lucas_solve()
 *   Auxiliary memory: O(n_l log n) limbs
 *   Output memory: O(n_l) limbs per result
 */
void bn_lucas(bignum* u, bignum* v, const bignum* p, const bignum* q,
              const bignum* n)
{
  bignum qn;

  bn_init(&qn);
  bn_lucas_solve(u, v, p, q, &qn, n);
  bn_free(&qn);
}

/*
 * Compute U_n(P, Q) mod m and V_n(P, Q) mod m.
 *
 * Let n_l = m->size, measured in 64-bit limbs.
 *
 * Thin wrapper around bn_lucas_solve_mod() that discards Q^n mod m.
 *
 * Complexity:
 *   Time: O(n_l^3) for the context initialization, then O(log n *
 *         n_l^2), see bn_lucas_solve_mod()
 *   Auxiliary memory: O(n_l) limbs
 *   Output memory: O(n_l) limbs per result
 */
void bn_lucas_mod(bignum* u, bignum* v, const bignum* p, const bignum* q,
                  const bignum* n, const bignum* m)
{
  bignum qn;

  bn_init(&qn);
  bn_lucas_solve_mod(u, v, p, q, &qn, n, m);
  bn_free(&qn);
}
