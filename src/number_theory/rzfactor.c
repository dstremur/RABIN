/*
 * rzfactor.c
 *
 * Integer factorization routines.
 *
 * This file implements complete factorization of a rz_t using a
 * combination of trial division (against a precomputed prime table),
 * Pollard's rho algorithm, and Pollard's p - 1 method, with BPSW
 * primality testing to detect prime factors.
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

#include "../../include/rabin.h"
#include "../../include/rns_primes.h"
#include "../../include/rvec.h"
#include "stdio.h"

bool rz_trialdiv(const rz_t* n, u64 g)
{
  if (n == NULL) return false;

  u64 i = 0;

  while (i < 50000 && SMALL_PRIMES[i] <= g) {
    if (rz_mod_u64(n, SMALL_PRIMES[i]) == 0) {
      return false;
    }
    i++;
  }

  return true;
}

bool rz_trialdiv_factor(rz_t* f, const rz_t* n, u64 g)
{
  if (f == NULL || n == NULL) return false;

  u64 i = 0;

  while (i < 70000 && SMALL_PRIMES[i] < g) {
    if (rz_mod_u64(n, SMALL_PRIMES[i]) == 0) {
      if (rz_set_u64(f, SMALL_PRIMES[i]) != RABIN_SUCCESS) return false;
      return false;
    }
    i++;
  }

  return true;
}

bool rz_pollard_rho_inner(rz_t* f, const rz_t* n)
{
  if (f == NULL || n == NULL) return false;
  if (rz_is_even(n)) {
    return rz_set_u64(f, 2) == RABIN_SUCCESS;
  }

  rz_t a, b, c, d, tmp;
  rz_init_multi(&a, &b, &c, &d, &tmp, NULL);

  if (rz_set_u64(&a, 2) != RABIN_SUCCESS) goto out;
  if (rz_set_u64(&b, 2) != RABIN_SUCCESS) goto out;
  if (rz_gen_random_range(&c, &RZ_TWO, n) != RABIN_SUCCESS) goto out;

  for (u64 i = 0; i < 1000000; i++) {
    if (rz_mul(&a, &a, &a) != RABIN_SUCCESS) goto out;
    if (rz_add(&a, &a, &c) != RABIN_SUCCESS) goto out;
    if (rz_mod(&a, &a, n) != RABIN_SUCCESS) goto out;

    if (rz_mul(&b, &b, &b) != RABIN_SUCCESS) goto out;
    if (rz_add(&b, &b, &c) != RABIN_SUCCESS) goto out;
    if (rz_mod(&b, &b, n) != RABIN_SUCCESS) goto out;

    if (rz_mul(&b, &b, &b) != RABIN_SUCCESS) goto out;
    if (rz_add(&b, &b, &c) != RABIN_SUCCESS) goto out;
    if (rz_mod(&b, &b, n) != RABIN_SUCCESS) goto out;

    if (rz_cmp(&a, &b) >= 0) {
      if (rz_sub(&tmp, &a, &b) != RABIN_SUCCESS) goto out;
    } else {
      if (rz_sub(&tmp, &b, &a) != RABIN_SUCCESS) goto out;
    }

    if (rz_gcd(&d, &tmp, n) != RABIN_SUCCESS) goto out;

    if (rz_set_u64(&tmp, 1) != RABIN_SUCCESS) goto out;
    if (rz_cmp(&tmp, &d) == -1 && rz_cmp(&d, n) == -1) {
      if (rz_copy(f, &d) != RABIN_SUCCESS) goto out;
      rz_clear_multi(&a, &b, &c, &d, &tmp, NULL);
      return true;
    }

    if (rz_cmp(&d, n) == 0) {
      goto out;
    }
  }
out:
  rz_clear_multi(&a, &b, &c, &d, &tmp, NULL);
  return false;
}

bool rz_pollard_rho(rz_t* f, const rz_t* n)
{
  if (f == NULL || n == NULL) return false;
  if (rz_is_even(n)) {
    return rz_set_u64(f, 2) == RABIN_SUCCESS;
  }

  // try several times
  for (u64 i = 0; i < 64; i++) {
    if (rz_pollard_rho_inner(f, n)) {
      return true;
    }
  }

  return false;
}

/**
 * @brief Append f to the factor vector unless it equals the last factor
 * already stored.
 *
 * Keeps consecutive factors distinct so repeated factors are recorded
 * once per run; the comparison is a deep value comparison (rz_cmp),
 * never a raw struct comparison.
 *
 * Complexity:
 *   Time: O(n) where n is the size of f in limbs
 *   Auxiliary memory: O(1)
 *   Output memory: O(n) limbs appended
 *
 * @param[in,out] v Factor vector to append to.
 * @param[in]     f Factor to append.
 */
static rabin_err_t rvec_append_distinct(rvec_t* v, const rz_t* f)
{
  // 1. Compare values deeply using rz_cmp (never compare raw struct memory
  // blocks)
  if (v->size > 0) {
    if (rz_cmp(&v->data[v->size - 1], f) == 0) {
      return RABIN_SUCCESS;  // Element matches the last factor found, exit
                             // out to stay distinct
    }
  }

  // 2. Hand off to rvec_append(), which performs its own rz_copy
  return rvec_append(v, f);
}

rabin_err_t rz_factorize(rvec_t* v, rz_t* n)
{
  if (v == NULL || n == NULL) return RABIN_ERR_NULL_PTR;

  if (rz_cmp(n, &RZ_ONE) == 0 || rz_is_zero(n) || n->is_neg) {
    return RABIN_SUCCESS;
  }

  if (rz_is_eq_i64(n, 2) || rz_is_eq_i64(n, 3)) {
    return rvec_append_distinct(v, n);
  }

  if (rz_bpsw(n)) {
    return rvec_append_distinct(v, n);
  }

  rz_t tmp, f, n1, rem;
  rz_init_multi(&tmp, &f, &n1, &rem, NULL);
  rabin_err_t err = rz_copy(&tmp, n);
  if (err != RABIN_SUCCESS) goto cleanup;

  // even case
  if (rz_is_even(&tmp)) {
    if ((err = rz_set_u64(&f, 2)) != RABIN_SUCCESS) goto cleanup;
    if ((err = rvec_append_distinct(v, &f)) != RABIN_SUCCESS) goto cleanup;

    while (rz_is_even(&tmp)) {
      if ((err = rz_div(&tmp, &tmp, &f)) != RABIN_SUCCESS) goto cleanup;
    }

    err = rz_factorize(v, &tmp);
    goto cleanup;
  }

  // first rz_trialdiv
  if (!rz_trialdiv_factor(&f, &tmp, 50000)) {
    if ((err = rvec_append_distinct(v, &f)) != RABIN_SUCCESS) goto cleanup;

    if ((err = rz_div(&tmp, &tmp, &f)) != RABIN_SUCCESS) goto cleanup;
    if ((err = rz_mod(&rem, &tmp, &f)) != RABIN_SUCCESS) goto cleanup;
    while (rz_is_zero(&rem)) {
      if ((err = rz_div(&tmp, &tmp, &f)) != RABIN_SUCCESS) goto cleanup;
      if ((err = rz_mod(&rem, &tmp, &f)) != RABIN_SUCCESS) goto cleanup;
    }

    err = rz_factorize(v, &tmp);
    goto cleanup;
  }

  bool res = false;

  // pollard rho
  for (u64 i = 0; i < 16 && !res; i++) {
    if (rz_pollard_rho(&f, &tmp)) {
      res = true;
    }
  }

  // p - 1
  if (!res) {
    if (rz_pollard_p_minus_one(&f, &tmp)) {
      res = true;
    }
  }

  if (!res) {
    // no algorithm managed to split this composite
    err = RABIN_ERR_INVALID_ARG;
    goto cleanup;
  }

  // rz_println(&f);

  if (rz_trialdiv(&f, 10000) && (!rz_is_zero(&f))) {
    if ((err = rvec_append_distinct(v, &f)) != RABIN_SUCCESS) goto cleanup;

    if ((err = rz_div(&tmp, &tmp, &f)) != RABIN_SUCCESS) goto cleanup;
    if ((err = rz_mod(&rem, &tmp, &f)) != RABIN_SUCCESS) goto cleanup;
    while (rz_is_zero(&rem)) {
      if ((err = rz_div(&tmp, &tmp, &f)) != RABIN_SUCCESS) goto cleanup;
      if ((err = rz_mod(&rem, &tmp, &f)) != RABIN_SUCCESS) goto cleanup;
    }

    err = rz_factorize(v, &tmp);
    goto cleanup;
  } else {
    if ((err = rz_div(&n1, &tmp, &f)) != RABIN_SUCCESS) goto cleanup;
    if ((err = rz_factorize(v, &f)) != RABIN_SUCCESS) goto cleanup;
    err = rz_factorize(v, &n1);
    goto cleanup;
  }

cleanup:
  rz_clear_multi(&tmp, &f, &n1, &rem, NULL);
  return err;
}

bool rz_pollard_p_minus_one_stage_1(rz_t* f, const rz_t* n, u64 B,
                                    u64 iterations)
{
  if (f == NULL || n == NULL) return false;

  rz_t M, a, tmp, n_min1, a_min_1, q, ln_q, ln_n, l;
  rz_init_multi(&M, &a, &tmp, &n_min1, &a_min_1, &q, &ln_q, &ln_n, &l, NULL);

  bool res = false;
  if (rz_sub(&n_min1, n, &RZ_ONE) != RABIN_SUCCESS) goto cleanup;

  for (u64 attempt = 0; attempt < iterations; attempt++) {
    // set a to random number in 2 <= a <= n - 1
    if (rz_gen_random_range(&a, &RZ_TWO, &n_min1) != RABIN_SUCCESS)
      goto cleanup;

    if (rz_gcd(&tmp, &a, n) != RABIN_SUCCESS) goto cleanup;

    // if gcd(a,n) >= 2 found a factor
    if (rz_cmp(&tmp, &RZ_ONE) > 0 && rz_cmp(&tmp, n) < 0) {
      if (rz_copy(f, &tmp) != RABIN_SUCCESS) goto cleanup;
      res = true;
      goto cleanup;
    }

    int i = 0;

    // go through all primes upto B
    while (SMALL_PRIMES[i] <= B) {
      u64 q_val = SMALL_PRIMES[i];
      u64 q_pow = q_val;

      while (q_pow <= B / q_val) {
        q_pow *= q_val;
      }

      if (rz_set_u64(&tmp, q_pow) != RABIN_SUCCESS) goto cleanup;
      if (rz_mod_exp(&a, &a, &tmp, n) != RABIN_SUCCESS) goto cleanup;

      i++;
    }

    // nmin1 = a - 1
    if (rz_sub(&a_min_1, &a, &RZ_ONE) != RABIN_SUCCESS) goto cleanup;

    if (rz_gcd(&tmp, &a_min_1, n) != RABIN_SUCCESS) goto cleanup;

    if (rz_cmp(&tmp, &RZ_ONE) > 0 && rz_cmp(&tmp, n) < 0) {
      if (rz_copy(f, &tmp) != RABIN_SUCCESS) goto cleanup;
      res = true;
      goto cleanup;
    }
  }

cleanup:
  rz_clear_multi(&M, &a, &tmp, &n_min1, &a_min_1, &q, &ln_q, &ln_n, &l, NULL);
  return res;
}

bool rz_pollard_p_minus_one(rz_t* f, const rz_t* n)
{
  if (rz_pollard_p_minus_one_stage_1(f, n, 10000, 1000)) {
    return true;
  }

  if (rz_pollard_p_minus_one_stage_1(f, n, 10000, 100)) {
    return true;
  }

  return rz_pollard_p_minus_one_stage_1(f, n, 70000, 100);
}
