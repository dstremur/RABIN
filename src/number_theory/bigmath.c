/*
 * bigmath.c
 *
 * Classical number-theory routines.
 *
 * This file implements the Jacobi symbol, the Tonelli-Shanks square
 * root algorithm modulo a prime, and the binary (Stein) GCD.
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

#include <ctype.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>

#include "../include/bignum.h"

/**
 * @brief Compute the Jacobi symbol (a / m).
 *
 * Let n = m->size, measured in 64-bit limbs.
 *
 * Returns 1 if a is a quadratic residue mod m, -1 if it is a
 * nonresidue, and 0 if gcd(a, m) > 1. m must be positive and odd;
 * otherwise 0 is returned and a message is printed.
 *
 * Uses the binary algorithm: repeatedly strip factors of 2 from a
 * (flipping the sign when both the stripped power and m mod 8 are
 * odd), apply the mutual-reciprocity rule when both a and m are 3
 * mod 4, and reduce m mod a.
 *
 * Complexity:
 *   Time: O(n^2) - O(n) iterations of O(n) modular reductions
 *   Auxiliary memory: O(n) limbs for temporaries
 *   Output memory: O(1)
 *
 * @param[in] a Numerator of the Jacobi symbol.
 * @param[in] m Denominator (must be positive and odd).
 *
 * @return 1 If a is a quadratic residue mod m, -1 if a nonresidue, 0
 *           if gcd(a, m) > 1 or m is not positive and odd.
 */
i64 bn_jacobi(const bignum* a, const bignum* m)
{
  if (bn_is_zero(m) || bn_is_even(m)) {
    printf("m must be positive and odd\n");
    return 0;
  }

  bignum A, M, R, one;
  bn_init_multi(&A, &M, &R, &one, NULL);

  bn_set_u64(&one, 1);
  bn_copy(&M, m);

  bn_mod(&A, a, &M);
  if (a->is_neg && !bn_is_zero(&A)) {
    bn_add(&A, &M, &A);
  }

  int t = 0;

  while (!bn_is_zero(&A)) {
    u64 z = bn_cnt_trailing_zeros(&A);

    if (z > 0) {
      bn_rshift(&A, &A, z);

      if ((z & 1) == 1) {
        u64 m_val = M.size > 0 ? M.limbs[0] : 0;
        u64 m_mod8 = m_val & 7;
        if (m_mod8 == 3 || m_mod8 == 5) {
          t ^= 1;
        }
      }
    }

    u64 a_val = A.size > 0 ? A.limbs[0] : 0;
    u64 m_val = M.size > 0 ? M.limbs[0] : 0;

    if ((a_val & 3) == 3 && (m_val & 3) == 3) {
      t ^= 1;
    }
    bn_mod(&R, &M, &A);
    bn_copy(&M, &A);
    bn_copy(&A, &R);
  }

  i64 res;

  if (bn_cmp(&M, &one) != 0) {
    res = 0;
  } else {
    res = (t % 2 == 0) ? 1 : -1;
  }

  bn_free(&A);
  bn_free(&M);
  bn_free(&R);
  bn_free(&one);

  return res;
}

/**
 * @brief Tonelli-Shanks: square root of n modulo the prime p.
 *
 * Let n_l = p->size, measured in 64-bit limbs.
 *
 * Inputs:
 * p, a prime
 * n, an element of Z / p Z such that solutions to the congruence
 * r^2 = n exist; when this is so we say that n is a quadratic
 * residue mod p.
 *
 * Computes r with r^2 = n (mod p) and stores it in r. The algorithm
 * factors p - 1 = Q * 2^S with Q odd, finds a quadratic nonresidue z,
 * and iteratively refines the candidate root R until t = n^Q * c^2
 * reaches 1.
 *
 * If n is zero, r is set to 0. If n is not a quadratic residue mod p
 * (Jacobi symbol != 1), a message is printed and r is left unchanged.
 *
 * Complexity:
 *   Time: O(n_l^2 * S^2) worst case - O(S) iterations, each with O(S)
 *         modular squarings, plus O(n_l^2) for the initial
 *         exponentiations and the nonresidue search
 *   Auxiliary memory: O(n_l) limbs for temporaries
 *   Output memory: O(n_l) limbs
 *
 * @param[out] r Result storing a square root of n mod p (if it exists).
 * @param[in]  n Value whose square root is computed (a quadratic
 *               residue mod p).
 * @param[in]  p Prime modulus.
 */
void tonelli_shanks(bignum* r, const bignum* n, const bignum* p)
{
  if (bn_is_zero(n)) {
    bn_set_u64(r, 0);
    return;
  }

  if (bn_jacobi(n, p) != 1) {
    printf("No square roots exist\n");
    return;
  }

  bignum p_minus_one, one, Q, z, M, c, t, R, exp2, tmp2, b2;
  bignum i, tmp, b, b_exp, j;

  bn_init_multi(&p_minus_one, &one, &Q, &z, &M, &c, &t, &R, &exp2, &tmp2, &b2,
                NULL);
  bn_init_multi(&i, &tmp, &b, &b_exp, &j, NULL);

  bn_set_i64(&one, 1);
  bn_copy(&p_minus_one, p);
  bn_sub(&p_minus_one, &p_minus_one, &one);
  bn_copy(&Q, &p_minus_one);

  u64 S = 0;
  while (bn_is_even(&Q)) {
    bn_rshift1(&Q);
    S++;
  }

  // now p - 1 = Q2^S

  bn_set_u64(&z, 2);
  while (bn_cmp(&z, p) < 0) {
    if (bn_jacobi(&z, p) == -1) {
      break;
    }
    bn_add(&z, &z, &one);
  }

  bn_set_u64(&M, S);
  bn_mod_exp(&c, &z, &Q, p);
  bn_mod_exp(&t, n, &Q, p);

  bn_copy(&exp2, &Q);
  bn_add_u64(&exp2, &exp2, 1);
  bn_rshift1(&exp2);

  bn_mod_exp(&R, n, &exp2, p);

  while (1) {
    if (bn_is_zero(&t)) {
      bn_set_u64(r, 0);
      break;
    }

    if (bn_cmp(&t, &one) == 0) {
      bn_copy(r, &R);
      break;
    }

    bn_set_u64(&i, 0);
    bn_copy(&tmp, &t);

    while (!bn_is_eq_i64(&tmp, 1) && bn_cmp(&i, &M) < 0) {
      bn_mul(&tmp2, &tmp, &tmp);
      bn_mod(&tmp, &tmp2, p);
      bn_add(&i, &i, &one);
    }

    if (bn_cmp(&i, &M) == 0) {
      printf("No quadratic residue\n");
      break;
    }

    bn_copy(&b_exp, &M);
    bn_sub(&b_exp, &b_exp, &i);
    bn_sub(&b_exp, &b_exp, &one);

    bn_copy(&b, &c);
    bn_set_u64(&j, 0);

    while (bn_cmp(&j, &b_exp) < 0) {
      bn_mul(&b2, &b, &b);
      bn_mod(&b, &b2, p);
      bn_add(&j, &j, &one);
    }

    bn_copy(&M, &i);

    bn_mul(&c, &b, &b);
    bn_mod(&c, &c, p);

    bn_mul(&t, &t, &c);
    bn_mod(&t, &t, p);

    bn_mul(&R, &R, &b);
    bn_mod(&R, &R, p);
  }

  bn_free(&p_minus_one);
  bn_free(&one);
  bn_free(&Q);
  bn_free(&z);
  bn_free(&M);
  bn_free(&c);
  bn_free(&t);
  bn_free(&R);
  bn_free(&exp2);
  bn_free(&i);
  bn_free(&tmp);
  bn_free(&b);
  bn_free(&b_exp);
  bn_free(&j);
  bn_free(&tmp2);
  bn_free(&b2);
}

/**
 * @brief Greatest common divisor via the binary (Stein) algorithm:
 * d = gcd(a, b).
 *
 * Let n = max(a->size, b->size), measured in 64-bit limbs.
 *
 * Strips the common power of two from both operands, then repeatedly
 * removes factors of 2 and replaces the larger operand by half the
 * difference, until one operand is zero. The common power of two is
 * restored at the end. The result is always nonnegative.
 *
 * Complexity:
 *   Time: O(n^2) - O(n) iterations of O(n) subtractions/shifts
 *   Auxiliary memory: O(n) limbs for temporaries
 *   Output memory: O(n) limbs
 *
 * @param[out] d Result storing gcd(a, b).
 * @param[in]  a First operand.
 * @param[in]  b Second operand.
 */
void bn_gcd(bignum* d, const bignum* a, const bignum* b)
{
  if (bn_is_zero(a)) {
    bn_copy(d, b);
    return;
  }
  if (bn_is_zero(b)) {
    bn_copy(d, a);
    return;
  }

  bn_set_u64(d, 1);
  bignum t, tmp_a, tmp_b;
  bn_init_multi(&t, &tmp_a, &tmp_b, NULL);
  bn_copy(&tmp_a, a);
  bn_copy(&tmp_b, b);

  u64 shifts = 0;

  while (bn_is_even(&tmp_a) && bn_is_even(&tmp_b)) {
    bn_rshift1(&tmp_a);
    bn_rshift1(&tmp_b);
    shifts++;
  }

  while (!bn_is_zero(&tmp_a)) {
    while (bn_is_even(&tmp_a)) {
      bn_rshift1(&tmp_a);
    }
    while (bn_is_even(&tmp_b)) {
      bn_rshift1(&tmp_b);
    }

    if (bn_cmp(&tmp_a, &tmp_b) < 0) {
      bn_sub(&t, &tmp_b, &tmp_a);
      bn_rshift1(&t);
      bn_copy(&tmp_b, &t);
    } else {
      bn_sub(&t, &tmp_a, &tmp_b);
      bn_rshift1(&t);
      bn_copy(&tmp_a, &t);
    }
  }

  bn_lshift(d, &tmp_b, shifts);

  bn_free(&tmp_a);
  bn_free(&tmp_b);
  bn_free(&t);
}
