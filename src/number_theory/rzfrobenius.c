/*
 * rzfrobenius.c
 *
 * Frobenius pseudo prime test
 *
 * This file implements the quadratic frobenius primality test
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

#include <fcntl.h>
#include <limits.h>
#include <math.h>
#include <stdarg.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "../../include/rns_primes.h"
#include "../../include/rz.h"
#include "../../include/rzcert.h"
#include "../../include/rzfactor.h"
#include "../../include/rzlucas.h"
#include "../../include/rzmath.h"
#include "../../include/rzprime.h"
#include "../../include/rzrabin.h"
#include "../../include/rzrand.h"
#include "../../include/u64.h"

bool qft_select_parameters(const rz_t* n, rz_t* out_b, rz_t* out_c,
                           rz_t* out_delta)
{
  rz_t b, c, delta, neg_c, temp;
  if (rz_init_multi(&b, &c, &delta, &neg_c, &temp, NULL) != RABIN_SUCCESS)
    return false;

  bool found = false;
  for (int b_val = 1; b_val <= 1000; b_val++) {
    rz_set_i64(&b, b_val);
    for (int c_val = 1; c_val <= 1000; c_val++) {
      rz_set_i64(&c, c_val);

      // 1. Calculate delta = b^2 + 4c
      rz_mul(&delta, &b, &b);
      rz_mul_i64(&temp, &c, 4);
      rz_add(&delta, &delta, &temp);

      // 2. Calculate -c for the second Jacobi check
      rz_copy(&neg_c, &c);
      neg_c.is_neg = !neg_c.is_neg;

      // Check Grantham's conditions: ((b^2+4c)/n) == -1 and ((-c)/n) == 1
      if (rz_jacobi(&delta, n) == -1 && rz_jacobi(&neg_c, n) == 1) {
        rz_copy(out_b, &b);
        rz_copy(out_c, &c);
        rz_copy(out_delta, &delta);
        found = true;
        goto cleanup1;
      }
    }
  }

cleanup1:
  rz_clear_multi(&b, &c, &delta, &neg_c, &temp, NULL);
  return found;
}

/*
We use
x^k = U_k x + (V_k - bU_k) / 2 mod n
Prove later
*/

rabin_err_t rz_qft_core(rz_t* n, rz_t* b, rz_t* c, rz_t* delta, u64 B)
{
  // 1) Trial division by primes <= min{B, sqrt(n)}
  if (!rz_trialdiv(n, B)) {
    return RABIN_ERR_COMPOSITE;
  }

  // 2) Check for perfect square
  if (rz_is_square(NULL, n)) {
    return RABIN_ERR_COMPOSITE;
  }

  rabin_err_t result = RABIN_SUCCESS;

  // Compute q_param = -c mod n
  rz_t q_param, two_c, two, minus_two;
  rz_init_multi(&q_param, &two_c, &two, &minus_two, NULL);

  rz_neg(&q_param, c);
  rz_mod(&q_param, &q_param, n);

  if (q_param.is_neg) {
    rz_add(&q_param, &q_param, n);
  }

  // Prepare 2 mod n and -2 mod n
  rz_set_i64(&two, 2);
  rz_mod(&two, &two, n);

  rz_set_i64(&minus_two, -2);
  rz_mod(&minus_two, &minus_two, n);

  if (minus_two.is_neg) {
    rz_add(&minus_two, &minus_two, n);
  }

  // 3) x^{(n+1)/2} mod (n, x^2 - bx - c)
  rz_t idx, u, v, q;
  rz_init_multi(&idx, &u, &v, &q, NULL);

  // idx = (n + 1) / 2
  rz_add_u64(&idx, n, 1);
  rz_rshift1(&idx);

  // Solve Lucas sequence
  if (rz_lucas_solve_mod(&u, &v, b, &q_param, &q, &idx, n) != RABIN_SUCCESS) {
    result = RABIN_ERR_COMPOSITE;
    goto cleanup2;
  }

  // Check whether x^{(n+1)/2} is in Z/nZ
  // This is equivalent to U_{(n+1)/2} = 0 mod n
  if (!rz_is_zero(&u)) {
    result = RABIN_ERR_COMPOSITE;
    goto cleanup2;
  }

  // 4) x^{n+1} mod (n, x^2 - bx - c) = -c

  // idx = n + 1
  rz_add_u64(&idx, n, 1);

  // compute Lucas sequence
  if (rz_lucas_solve_mod(&u, &v, b, &q_param, &q, &idx, n) != RABIN_SUCCESS) {
    result = RABIN_ERR_COMPOSITE;
    goto cleanup2;
  }

  // Check U_{n+1} = 0 and V_{n+1} = -2c mod n
  rz_mul_i64(&two_c, c, -2);
  rz_mod(&two_c, &two_c, n);

  if (two_c.is_neg) {
    rz_add(&two_c, &two_c, n);
  }

  if (!rz_is_zero(&u) || rz_cmp(&v, &two_c) != 0) {
    result = RABIN_ERR_COMPOSITE;
    goto cleanup2;
  }

  // 5)
  rz_t n_sq, s, exp, pow_2;
  rz_init_multi(&n_sq, &s, &exp, &pow_2, NULL);

  rz_set_i64(&s, 1);

  rz_sqr(&n_sq, n);
  rz_sub(&n_sq, &n_sq, &s);

  // n^2 - 1 = 2^r * s, where s is odd
  u64 r = rz_cnt_trailing_zeros(&n_sq);
  rz_rshift(&s, &n_sq, r);

  bool pass = false;

  // First check whether x^s = 1
  // x^s = U_s x + (V_s - bU_s) / 2 => U_s = 0 and V_s = 2
  if (rz_lucas_solve_mod(&u, &v, b, &q_param, &q, &s, n) != RABIN_SUCCESS) {
    result = RABIN_ERR_COMPOSITE;
    goto cleanup3;
  }

  // x^s = 1 iff U_s = 0 and V_s = 2 mod n
  if (rz_is_zero(&u) && rz_cmp(&v, &two) == 0) {
    pass = true;
  }

  // Otherwise check whether x^(2^j * s) = -1
  // for some 0 <= j <= r - 2
  for (u64 j = 0; !pass && j <= r - 2; j++) {
    rz_set_u64(&pow_2, 1);
    rz_lshift(&pow_2, &pow_2, j);

    rz_mul(&exp, &pow_2, &s);

    if (rz_lucas_solve_mod(&u, &v, b, &q_param, &q, &exp, n) != RABIN_SUCCESS) {
      result = RABIN_ERR_COMPOSITE;
      goto cleanup3;
    }

    // x^(2^j * s) = -1 iff
    // U_(2^j * s) = 0 and V_(2^j * s) = -2 mod n
    if (rz_is_zero(&u) && rz_cmp(&v, &minus_two) == 0) {
      pass = true;
    }
  }

  if (!pass) {
    result = RABIN_ERR_COMPOSITE;
  } else {
    result = RABIN_SUCCESS;
  }

cleanup3:
  rz_clear_multi(&n_sq, &s, &exp, &pow_2, NULL);

cleanup2:
  rz_clear_multi(&idx, &u, &v, &q, &q_param, &two_c, &two, &minus_two, NULL);

  return result;
}

rabin_err_t rz_qft_deterministic(rz_t* n)
{
  rz_t b, c, delta;
  if (rz_init_multi(&b, &c, &delta, NULL) != RABIN_SUCCESS)
    return RABIN_ERR_OUT_OF_MEMORY;

  // parameter selection
  if (!qft_select_parameters(n, &b, &c, &delta)) {
    rz_clear_multi(&b, &c, &delta, NULL);
    return RABIN_ERR_INCONSISTENT;
  }

  rabin_err_t result = rz_qft_core(n, &b, &c, &delta, 50000);

  rz_clear_multi(&b, &c, &delta, NULL);

  return result;
}

bool rz_is_prime_qft(rz_t* n)
{
  if (rz_qft_deterministic(n) == RABIN_SUCCESS) {
    return true;
  } else {
    return false;
  }
}
