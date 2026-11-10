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
        goto cleanup;
      }
    }
  }

cleanup:
  rz_clear_multi(&b, &c, &delta, &neg_c, &temp, NULL);
  return found;
}

/*
We use
x^k = U_k x + (V_k - bU_k) / 2 mod n
*/

rabin_err_t rz_qft(rz_t* n)
{
  rz_t b, c, delta;
  if (rz_init_multi(&b, &c, &delta, NULL) != RABIN_SUCCESS)
    return RABIN_ERR_OUT_OF_MEMORY;

  // parameter selection
  if (!qft_select_parameters(n, &b, &c, &delta)) {
    rz_clear_multi(&b, &c, &delta, NULL);
    return RABIN_ERR_INCONSISTENT;
  }

  rabin_err_t result = RABIN_SUCCESS;

  // 1) Trial division by primes <= min{B, sqrt(n)}

  if (!rz_trialdiv(n, 50000)) {
    result = RABIN_ERR_COMPOSITE;
    goto cleanup;
  }

  // 2) Check for perfect square
  if (rz_is_square(NULL, n)) {
    result = RABIN_ERR_COMPOSITE;
    goto cleanup;
  }

  // Pre-calculate Q = -c mod n using rz_neg
  rz_t neg_c, two_c_neg;
  rz_init_multi(&neg_c, &two_c_neg, NULL);

  // Set neg_c = -c and normalize modulo n
  rz_neg(&neg_c, &c);
  rz_mod(&neg_c, &neg_c, n);
  if (neg_c.is_neg) {
    rz_add(&neg_c, &neg_c, n);
  }

  // 3) x^{n+1} mod (n, x^2 -bx - c)
  rz_t idx, u, v, q;
  rz_init_multi(&idx, &u, &v, &q, NULL);

  // idx = (n + 1) / 2
  rz_add_u64(&idx, n, 1);
  rz_rshift1(&idx);

  // solve lucas sequence
  if (rz_lucas_solve_mod(&u, &v, &b, &neg_c, &q, &idx, n) != RABIN_SUCCESS) {
    result = RABIN_ERR_COMPOSITE;
    goto cleanstep3;
  }

  // 4. Check if x^{(n+1)/2} is in Z/nZ -> u_res (U_k) =0 mod n
  if (!rz_is_eq_i64(&u, 0)) {
    result = RABIN_ERR_COMPOSITE;  // Failed Step 3!
    goto cleanstep3;
  }

cleanstep3:
  rz_clear_multi(&idx, &u, &v, &q, NULL);
  if (result == RABIN_ERR_COMPOSITE) goto cleanup;

cleanup:
  rz_clear_multi(&b, &c, &delta, NULL);
  return result;
}