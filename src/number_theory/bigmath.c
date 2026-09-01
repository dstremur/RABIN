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
 * @brief  Computes the greatest common divisor (GCD) of two bignums.
 *
 * @details Uses the classic Euclidean algorithm,
 *          @c gcd(a, b) == gcd(b, a % b), iterating until the remainder
 *          is zero; the last non-zero value is the GCD. The loop rotates
 *          three internal buffers by pointer swap, so no temporary
 *          bignums are allocated per iteration.
 *
 *          Zero handling follows the standard conventions:
 *          - @p a is zero  -> result is @p b
 *          - @p b is zero  -> result is @p a
 *          - both zero     -> result is zero, i.e. `gcd(0, 0) == 0`
 *
 *          The result is the principal (non-negative) GCD; the sign of the
 *          inputs is ignored, so `bn_gcd(d, &a, &b) == bn_gcd(d, &a, &nb)`.
 *
 * @param[out] d  Receives `gcd(a, b)`. Must be initialized before the call.
 *                May alias @p a or @p b — both operands are copied into
 *                temporaries before any work is done — so
 *                `bn_gcd(&a, &a, &b)` is valid.
 * @param[in]  a  First operand. Not modified. Must not be @c NULL.
 * @param[in]  b  Second operand. Not modified. Must not be @c NULL.
 *
 * @pre @p d, @p a and @p b are initialized (e.g. via bn_init() /
 *      bn_init_multi()) and have valid size fields.
 *
 * @note The three temporaries created here are freed before returning;
 *       the caller only has to manage the lifetime of @p d.
 * @note bn_mod() is only called with a non-zero divisor, which is
 *       guaranteed by the `while (!bn_is_zero(v))` loop condition.
 *
 * @par Memory
 * Allocates 3 temporary bignums (each up to `max(size(a), size(b))`),
 * peak extra memory ≈ 3 operands. Fails silently on allocation error
 * if the underlying allocator does, so check @p d if that matters to you.
 *
 * @par Complexity
 * O(log min(a, b)) modulo operations; each modulo is O(n·m) word
 * divisions for n- and m-word operands.
 *
 * @warning <b>Not constant-time.</b> The number and shape of divisions
 *          depend on the operand values, which leaks information through
 *          execution time. Do not use on secret inputs (e.g. private keys
 *          in `invmod`/key-derivation paths) without a constant-time
 *          variant such as binary GCD.
 *
 * @par Example
 * @code
 * bignum a, b, g;
 * bn_init_multi(&a, &b, &g, NULL);
 *
 * bn_set_str(&a, "1071", 10);   // 1071 = 3 * 3 * 7 * 17
 * bn_set_str(&b, "462",  10);   //  462 = 2 * 3 * 7 * 11
 * bn_gcd(&g, &a, &b);           // g == 21
 *
 * bn_gcd(&a, &a, &b);           // also fine: a == 21, b unchanged
 *
 * bn_free_multi(&a, &b, &g, NULL);
 * @endcode
 *
 * @see bn_mod(), bn_copy(), bn_is_zero()
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

  bignum tmp_a, tmp_b, rem;
  bn_init_multi(&tmp_a, &tmp_b, &rem, NULL);
  bn_copy(&tmp_a, a);
  bn_copy(&tmp_b, b);

  bignum* u = &tmp_a;
  bignum* v = &tmp_b;
  bignum* r = &rem;

  while (!bn_is_zero(v)) {
    // r = u % v
    bn_mod(r, u, v);

    bignum* temp = u;
    u = v;
    v = r;
    r = temp;
  }

  // The GCD is left in 'u'
  bn_copy(d, u);

  bn_free(&tmp_a);
  bn_free(&tmp_b);
  bn_free(&rem);
}
