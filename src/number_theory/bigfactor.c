/*
 * bigfactor.c
 *
 * Integer factorization routines.
 *
 * This file implements complete factorization of a bignum using a
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

#include "../../include/bignum.h"
#include "../../include/bigvector.h"
#include "../../include/primes.h"
#include "stdio.h"

/**
 * @brief Trial-divide n against the prime table up to the bound g.
 *
 * Returns false if any prime p < g (up to 50000 table entries) divides
 * n, true if no such divisor was found.
 *
 * Complexity:
 *   Time: O(min(50000, pi(g)) * n) - one O(n) modular reduction per
 *         prime
 *   Auxiliary memory: O(1)
 *   Output memory: O(1)
 *
 * @param[in,out] n Number to trial-divide (reduced in place).
 * @param[in]     g Upper bound on the trial primes.
 *
 * @return true  If no prime p < g divides n.
 * @return false If a prime p < g divides n.
 */
bool trialdiv(bignum* n, u64 g)
{
  u64 i = 0;

  while (i < 50000 && primes[i] < g) {
    if (bn_mod_u64(n, primes[i]) == 0) {
      return false;
    }
    i++;
  }

  return true;
}

/**
 * @brief Trial-divide n against the prime table up to the bound g, storing
 * the first divisor found in f.
 *
 * Returns false and sets f to the smallest prime p < g (up to 70000
 * table entries) that divides n. Returns true if no such divisor was
 * found (f is left unchanged).
 *
 * Complexity:
 *   Time: O(min(70000, pi(g)) * n)
 *   Auxiliary memory: O(1)
 *   Output memory: O(1) limbs for f
 *
 * @param[out]    f Result storing the first divisor found.
 * @param[in,out] n Number to trial-divide (reduced in place).
 * @param[in]     g Upper bound on the trial primes.
 *
 * @return true  If no prime p < g divides n (f left unchanged).
 * @return false If a prime p < g divides n (f set to it).
 */
bool trialdiv_factor(bignum* f, bignum* n, u64 g)
{
  u64 i = 0;

  while (i < 70000 && primes[i] < g) {
    if (bn_mod_u64(n, primes[i]) == 0) {
      bn_set_u64(f, primes[i]);
      return false;
    }
    i++;
  }

  return true;
}

/**
 * @brief One attempt at Pollard's rho with Floyd's cycle detection.
 *
 * Let n = n->size, measured in 64-bit limbs.
 *
 * Iterates x -> x^2 + c (mod n) with two pointers (one step and two
 * steps per round) and computes gcd(|x - y|, n) each round. A
 * nontrivial gcd (1 < d < n) is a factor of n.
 *
 * Returns true and stores the factor in f on success. Returns false if
 * the cycle closes without a factor (d == n) or after 1000000
 * iterations.
 *
 * Complexity:
 *   Time: O(n^{1/4} * n^2) expected - about n^{1/4} iterations, each
 *         costing a few n-limb multiplications and a gcd
 *   Auxiliary memory: O(n) limbs for temporaries
 *   Output memory: O(n) limbs for f
 *
 * @param[out] f Result storing the factor found.
 * @param[in]  n Number to factor.
 *
 * @return true  On success (a nontrivial factor is stored in f).
 * @return false If the cycle closes without a factor or after 1000000
 *               iterations.
 */
bool bn_pollard_rho_inner(bignum* f, const bignum* n)
{
  if (bn_is_even(n)) {
    bn_set_u64(f, 2);
    return true;
  }

  bignum a, b, c, d, tmp;
  bn_init_multi(&a, &b, &c, &d, &tmp, NULL);

  bn_set_u64(&a, 2);
  bn_set_u64(&b, 2);
  bn_gen_random_range(&c, &BN_TWO, n);

  for (u64 i = 0; i < 1000000; i++) {
    bn_mul(&a, &a, &a);
    bn_add(&a, &a, &c);
    bn_mod(&a, &a, n);

    bn_mul(&b, &b, &b);
    bn_add(&b, &b, &c);
    bn_mod(&b, &b, n);

    bn_mul(&b, &b, &b);
    bn_add(&b, &b, &c);
    bn_mod(&b, &b, n);

    if (bn_cmp(&a, &b) >= 0) {
      bn_sub(&tmp, &a, &b);
    } else {
      bn_sub(&tmp, &b, &a);
    }

    bn_gcd(&d, &tmp, n);

    bn_set_u64(&tmp, 1);
    if (bn_cmp(&tmp, &d) == -1 && bn_cmp(&d, n) == -1) {
      bn_copy(f, &d);
      bn_free_multi(&a, &b, &c, &d, &tmp, NULL);
      return true;
    }

    if (bn_cmp(&d, n) == 0) {
      bn_free_multi(&a, &b, &c, &d, &tmp, NULL);
      return false;
    }
  }
  bn_free_multi(&a, &b, &c, &d, &tmp, NULL);

  return false;
}

/**
 * @brief Pollard's rho factorization of n.
 *
 * Let n = n->size, measured in 64-bit limbs.
 *
 * Tries bn_pollard_rho_inner() up to 64 times with fresh random
 * constants. Returns true and stores a nontrivial factor in f on
 * success, false if all attempts failed.
 *
 * Complexity:
 *   Time: O(n^{1/4} * n^2) expected per attempt, up to 64 attempts
 *   Auxiliary memory: O(n) limbs
 *   Output memory: O(n) limbs for f
 *
 * @param[out] f Result storing the factor found.
 * @param[in]  n Number to factor.
 *
 * @return true  On success (a nontrivial factor is stored in f).
 * @return false If all 64 attempts failed.
 */
bool bn_pollard_rho(bignum* f, const bignum* n)
{
  if (bn_is_even(n)) {
    bn_set_u64(f, 2);
    return true;
  }

  // try several times
  for (u64 i = 0; i < 64; i++) {
    if (bn_pollard_rho_inner(f, n)) {
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
 * once per run; the comparison is a deep value comparison (bn_cmp),
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
static void bigvector_append_distinct(bigvector* v, const bignum* f)
{
  // 1. Compare values deeply using bn_cmp (never compare raw struct memory
  // blocks)
  if (v->size > 0) {
    if (bn_cmp(&v->data[v->size - 1], f) == 0) {
      return;  // Element matches the last factor found, exit out to stay
               // distinct
    }
  }

  // 2. Safely hand off to your append function, which executes its own safe
  // bn_copy
  bigvector_append(v, (bignum*)f);
}

/**
 * @brief Completely factorize n, appending the factors to the vector v.
 *
 * Let n = n->size, measured in 64-bit limbs.
 *
 * Strategy, applied recursively:
 *
 *   1. stop at 0 or 1,
 *   2. append 2 or 3 directly,
 *   3. append n if BPSW says it is prime,
 *   4. peel off all factors of 2,
 *   5. trial division up to 50000 (peel off all powers of the factor),
 *   6. Pollard's rho (up to 16 attempts), then Pollard p - 1,
 *   7. if the found factor f is prime (trial division up to 10000),
 *      peel off all powers of f and recurse on the cofactor; otherwise
 *      recurse on both f and the cofactor.
 *
 * n is modified (it is a working copy); v receives the factors in no
 * particular order, with consecutive duplicates suppressed.
 *
 * Complexity:
 *   Time: subexponential in practice; dominated by the Pollard rho /
 *         p - 1 attempts on the hardest composite branch
 *   Auxiliary memory: O(n) limbs for temporaries
 *   Output memory: O(n) limbs for the factor list
 *
 * @param[out]    v Vector receiving the factors.
 * @param[in,out] n Number to factorize (modified in place).
 */
void bn_factorize(bigvector* v, bignum* n)
{
  if (bn_cmp(n, &BN_ONE) == 0 || bn_is_zero(n)) {
    return;
  }

  if (bn_is_eq_i64(n, 2) || bn_is_eq_i64(n, 3)) {
    bigvector_append_distinct(v, n);
    return;
  }

  if (bn_bpsw(n)) {
    bigvector_append_distinct(v, n);
    return;
  }

  bignum tmp, f, n1, rem;
  bn_init_multi(&tmp, &f, &n1, &rem, NULL);
  bn_copy(&tmp, n);

  // even case
  if (bn_is_even(&tmp)) {
    bn_set_u64(&f, 2);
    bigvector_append_distinct(v, &f);

    while (bn_is_even(&tmp)) {
      bn_div(&tmp, &tmp, &f);
    }

    bn_factorize(v, &tmp);
    goto cleanup;
  }

  // first trialdiv
  if (!trialdiv_factor(&f, &tmp, 50000)) {
    bigvector_append_distinct(v, &f);

    bn_div(&tmp, &tmp, &f);
    bn_mod(&rem, &tmp, &f);
    while (bn_is_zero(&rem)) {
      bn_div(&tmp, &tmp, &f);
      bn_mod(&rem, &tmp, &f);
    }

    bn_factorize(v, &tmp);
    goto cleanup;
  }

  bool res = false;

  // pollard rho
  for (u64 i = 0; i < 16 && !res; i++) {
    if (bn_pollard_rho(&f, &tmp)) {
      res = true;
    }
  }

  // p - 1
  if (!res) {
    printf("trying p-1 \n");
    if (bn_pollard_p_minus_one(&f, &tmp)) {
      res = true;
    }
  }

  if (!res) {
    printf("failed to split composite branch completely!\n");
    goto cleanup;
  }

  // bn_println(&f);

  if (trialdiv(&f, 10000) && (&f)) {
    bigvector_append_distinct(v, &f);

    bn_div(&tmp, &tmp, &f);
    bn_mod(&rem, &tmp, &f);
    while (bn_is_zero(&rem)) {
      bn_div(&tmp, &tmp, &f);
      bn_mod(&rem, &tmp, &f);
    }

    bn_factorize(v, &tmp);
  } else {
    bn_div(&n1, &tmp, &f);
    bn_factorize(v, &f);
    bn_factorize(v, &n1);
  }

cleanup:
  bn_free_multi(&tmp, &f, &n1, &rem, NULL);
}

/**
 * @brief One stage of Pollard's p - 1 method.
 *
 * Let n = n->size, measured in 64-bit limbs.
 *
 * Repeatedly (up to `iterations` times) picks a random base a in
 * [2, n-1] and computes a^M mod n where M is the product of the
 * highest prime powers q^e <= B for all primes q <= B. If p - 1 is
 * B-smooth for some prime factor p of n, then a^M = 1 mod p, so
 * gcd(a^M - 1, n) reveals p.
 *
 * Returns true and stores the factor in f on success, false otherwise.
 *
 * Complexity:
 *   Time: O(iterations * pi(B) * n^2) - one modular exponentiation per
 *         prime <= B per iteration
 *   Auxiliary memory: O(n) limbs for temporaries
 *   Output memory: O(n) limbs for f
 *
 * @param[out]         f          Result storing the factor found.
 * @param[in]          n          Number to factor.
 * @param[in]          B          Smoothness bound.
 * @param[in]          iterations Number of random-base attempts.
 *
 * @return true  On success (a nontrivial factor is stored in f).
 * @return false If no factor was found.
 */
bool bn_pollard_p_minus_one_stage_1(bignum* f, bignum* n, u64 B, u64 iterations)
{
  bignum M, a, tmp, n_min1, a_min_1, q, ln_q, ln_n, l;
  bn_init_multi(&M, &a, &tmp, &n_min1, &a_min_1, &q, &ln_q, &ln_n, &l, NULL);

  bool res = false;
  bn_sub(&n_min1, n, &BN_ONE);

  for (u64 attempt = 0; attempt < iterations; attempt++) {
    // set a to random number in 2 <= a <= a - 1
    bn_gen_random_range(&a, &BN_TWO, &n_min1);

    bn_gcd(&tmp, &a, n);

    // if gcd(a,n) >= 2 found a factor
    if (bn_cmp(&tmp, &BN_ONE) > 0 && bn_cmp(&tmp, n) < 0) {
      bn_copy(f, &tmp);
      res = true;
      goto cleanup;
    }

    int i = 0;

    // go through all primes upto B
    while (primes[i] <= B) {
      u64 q_val = primes[i];
      u64 q_pow = q_val;

      while (q_pow <= B / q_val) {
        q_pow *= q_val;
      }

      bn_set_u64(&tmp, q_pow);
      bn_mod_exp(&a, &a, &tmp, n);

      i++;
    }

    // nmin1 = a - 1
    bn_sub(&a_min_1, &a, &BN_ONE);

    bn_gcd(&tmp, &a_min_1, n);

    if (bn_cmp(&tmp, &BN_ONE) > 0 && bn_cmp(&tmp, n) < 0) {
      bn_copy(f, &tmp);
      res = true;
      goto cleanup;
    }
  }

cleanup:
  bn_free_multi(&M, &a, &tmp, &n_min1, &a_min_1, &q, &ln_q, &ln_n, &l, NULL);
  return res;
}

/**
 * @brief Pollard's p - 1 factorization of n, with three escalating stages.
 *
 * Let n = n->size, measured in 64-bit limbs.
 *
 * Runs bn_pollard_p_minus_one_stage_1() with increasing smoothness
 * bounds and iteration counts:
 *
 *   1. B = 10000, 1000 rounds
 *   2. B = 10000, 100 rounds
 *   3. B = 70000, 100 rounds
 *
 * Returns true and stores a nontrivial factor in f as soon as any
 * stage succeeds.
 *
 * Complexity:
 *   Time: sum of the three stages, O(pi(B) * n^2) per iteration
 *   Auxiliary memory: O(n) limbs
 *   Output memory: O(n) limbs for f
 *
 * @param[out] f Result storing the factor found.
 * @param[in]  n Number to factor.
 *
 * @return true  On success (a nontrivial factor is stored in f).
 * @return false If all stages failed.
 */
bool bn_pollard_p_minus_one(bignum* f, bignum* n)
{
  if (bn_pollard_p_minus_one_stage_1(f, n, 10000, 1000)) {
    return true;
  }

  if (bn_pollard_p_minus_one_stage_1(f, n, 10000, 100)) {
    return true;
  }

  return bn_pollard_p_minus_one_stage_1(f, n, 70000, 100);
}
