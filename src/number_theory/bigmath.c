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

#include "../../include/bignum.h"

// TODO: make faster
// Comp. Number Theory p.29: Algorithm 1.4.10
i64 bn_kronecker(const bignum* a, const bignum* b)
{
  static const short tab2[8] = {0, 1, 0, -1, 0, -1, 0, 1};

  // 1. [Test b == 0]
  if (bn_is_zero(b)) {
    bignum a_abs, one;
    bn_init_multi(&a_abs, &one);
    bn_set_u64(&one, 1);
    bn_copy(&a_abs, a);
    a_abs.is_neg = false;
    i64 r = (bn_cmp(&a_abs, &one) == 0) ? 1 : 0;
    bn_free_multi(&a_abs, &one);
    return r;
  }

  // 2. [Remove 2's from b]
  if (bn_is_even(a) && bn_is_even(b)) {
    return 0;
  }

  bignum a_work, b_work, r;
  bn_init_multi(&a_work, &b_work, &r, NULL);
  bn_copy(&a_work, a);
  bn_copy(&b_work, b);
  u64 v = bn_cnt_trailing_zeros(b);
  bn_rshift(&b_work, b, v);

  i64 k = 1;
  if ((v & 1) == 1) {
    // k = (-1)^((a^2 - 1) / 8)
    k = tab2[a_work.limbs[0] & 7];
  }

  if (b_work.is_neg) {
    b_work.is_neg = !b_work.is_neg;

    if (a_work.is_neg) {
      k = -k;
    }
  }

  if (a_work.is_neg) {
    a_work.is_neg = false;
    if (b_work.limbs[0] & 2) {
      k = -k;
    }
  }

  // 3. [Finished?]
  for (;;) {
    if (bn_is_zero(&a_work)) {
      bignum one;
      bn_init(&one);
      bn_set_u64(&one, 1);
      if (bn_cmp(&b_work, &one) == 1) {
        bn_free(&one);
        k = 0;
        goto end;
      }
      if (bn_cmp(&b_work, &one) == 0) {
        bn_free(&one);
        goto end;
      }
    }

    u64 v = bn_cnt_trailing_zeros(&a_work);
    bn_rshift(&a_work, &a_work, v);
    if ((v & 1) == 1) {
      k *= tab2[b_work.limbs[0] & 7];
    }

    // 4. [Apply reprocity]
    if ((a_work.limbs[0] & 2) & (b_work.limbs[0] & 2)) k = -k;

    bn_copy(&r, &a_work);
    r.is_neg = false;
    bn_mod(&a_work, &b_work, &r);
    bn_copy(&b_work, &r);
  }

end:
  bn_free_multi(&a_work, &b_work, &r, NULL);
  return k;
}

i64 bn_jacobi(const bignum* a, const bignum* b) { return bn_kronecker(a, b); }

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

void bn_gcd(bignum* d, const bignum* a, const bignum* b)
{
  if (bn_is_zero(a)) {
    bn_copy(d, b);
    d->is_neg = false;
    return;
  }
  if (bn_is_zero(b)) {
    bn_copy(d, a);
    d->is_neg = false;
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
  d->is_neg = false;

  bn_free(&tmp_a);
  bn_free(&tmp_b);
  bn_free(&rem);
}

// Computational number theory p.15
void bn_gcd_binary(bignum* d, const bignum* a, const bignum* b)
{
  bignum a_work, b_work, r;
  bn_init_multi(&a_work, &b_work, &r, NULL);
  bn_copy(&a_work, a);
  bn_copy(&b_work, b);

  // 1 [Reduce size]
  if (bn_cmp(&a_work, &b_work) == -1) {
    bn_swap(&a_work, &b_work);
  }

  if (bn_is_zero(&b_work)) {
    bn_copy(d, &a_work);
    bn_free_multi(&a_work, &b_work, &r, NULL);
    return;
  }

  bn_mod(&r, &a_work, &b_work);
  bn_copy(&a_work, &b_work);
  bn_copy(&b_work, &r);

  // 2. [Compute power of 2] & 3
  if (bn_is_zero(&b_work)) {
    bn_copy(d, &a_work);
    bn_free_multi(&a_work, &b_work, &r, NULL);
    return;
  }

  u64 k_a = bn_cnt_trailing_zeros(&a_work);
  u64 k_b = bn_cnt_trailing_zeros(&b_work);
  u64 k = MIN(k_a, k_b);

  bn_rshift(&a_work, &a_work, k_a);
  bn_rshift(&b_work, &b_work, k_b);

  // 4 & 5. [Subtract and Loop]
  while (true) {
    int cmp = bn_cmp(&a_work, &b_work);

    if (cmp == 0) {
      // If t = 0, output 2^k * a and terminate
      bn_lshift(d, &a_work, k);
      break;
    }

    if (cmp > 0) {
      // If a > b, t > 0. Set t = a - b, remove powers of 2, set a = t
      bn_sub(&r, &a_work, &b_work);
      bn_rshift(&a_work, &r, bn_cnt_trailing_zeros(&r));
    } else {
      // If a < b, t < 0. Set -t = b - a, remove powers of 2, set b = -t
      bn_sub(&r, &b_work, &a_work);
      bn_rshift(&b_work, &r, bn_cnt_trailing_zeros(&r));
    }
  }
}

void bn_gcd_lehmer(bignum* d, const bignum* a, const bignum* b)
{
  i64 a_hat, b_hat, A, B, C, D, T, q;
  bignum a_work, b_work, t, r, p1, p2;

  bn_init_multi(&a_work, &b_work, &t, &r, &p1, &p2, NULL);
  bn_copy(&a_work, a);
  bn_copy(&b_work, b);

  while (!bn_is_zero(&b_work)) {
    A = 1;
    B = 0;
    C = 0;
    D = 1;

    if (b_work.size > 1) {
      a_hat = a_work.limbs[a_work.size - 1];

      if (b_work.size == a_work.size) {
        b_hat = b_work.limbs[a_work.size - 1];
      } else {
        b_hat = 0;
      }

      // 1. [Test Quotient]
      while (true) {
        __int128 num1 = (__int128)(u64)a_hat + A;
        __int128 den1 = (__int128)(u64)b_hat + C;
        __int128 num2 = (__int128)(u64)a_hat + B;
        __int128 den2 = (__int128)(u64)b_hat + D;

        if (den1 == 0 || den2 == 0) break;

        q = (i64)(num1 / den1);
        if (q != (i64)(num2 / den2)) break;

        // 2. [Euclidian Step]
        T = A - q * C;
        A = C;
        C = T;
        T = B - q * D;
        B = D;
        D = T;
        T = a_hat - q * b_hat;
        a_hat = b_hat;
        b_hat = T;
      }
    }

    // 3. [Multi-precision step]
    if (B == 0) {
      // Replaced divmod with mod since quotient is not needed
      bn_mod(&t, &a_work, &b_work);
      bn_copy(&a_work, &b_work);
      bn_copy(&b_work, &t);
    } else {
      bn_mul_i64(&p1, &a_work, A);
      bn_mul_i64(&p2, &b_work, B);
      bn_add(&t, &p1, &p2);

      bn_mul_i64(&p1, &a_work, C);
      bn_mul_i64(&p2, &b_work, D);
      bn_add(&r, &p1, &p2);

      bn_copy(&a_work, &t);
      bn_copy(&b_work, &r);
    }
  }

  bn_copy(d, &a_work);
  bn_free_multi(&a_work, &b_work, &t, &r, &p1, &p2, NULL);
}

// Computational number theory p.16
void bn_gcd_extended(bignum* u, bignum* v, bignum* d, const bignum* a,
                     const bignum* b)
{
  bignum v_1, v_3, t_1, t_3, temp;
  bn_init_multi(&v_1, &v_3, &t_1, &t_3, &temp, NULL);

  bn_set_u64(u, 1);
  bn_copy(d, a);

  if (bn_is_zero(b)) {
    bn_set_u64(v, 0);
    // normalize: gcd is non-negative (GMP convention)
    if (d->is_neg) {
      d->is_neg = false;
      u->is_neg = !u->is_neg;
    }
    bn_free_multi(&v_1, &v_3, &t_1, &t_3, &temp, NULL);
    return;
  }

  bn_set_u64(&v_1, 0);
  bn_copy(&v_3, b);

  while (!bn_is_zero(&v_3)) {
    bn_divmod(&temp, &t_3, d, &v_3);

    // t_1 = u - qv_1
    bn_mul(&temp, &temp, &v_1);
    bn_copy(&t_1, u);
    bn_sub(&t_1, &t_1, &temp);

    bn_copy(u, &v_1);
    bn_copy(d, &v_3);
    bn_copy(&v_1, &t_1);
    bn_copy(&v_3, &t_3);
  }

  // normalize: gcd is non-negative (GMP convention); flip d and u
  // together so the Bezout identity u*a + v*b == d is preserved
  if (d->is_neg) {
    d->is_neg = false;
    u->is_neg = !u->is_neg;
  }

  // v = (d - au) / b
  bn_mul(&temp, a, u);
  bn_copy(v, d);
  bn_sub(v, v, &temp);
  bn_div(v, v, b);
  bn_free_multi(&v_1, &v_3, &t_1, &t_3, &temp, NULL);
  return;
}

void bn_gcd_extended_lehmer(bignum* u, bignum* v, bignum* d, const bignum* a,
                            const bignum* b)
{
  i64 a_hat, b_hat, A, B, C, D, T, q;
  bignum a_work, b_work, t, r, v_1, Q, p1, p2;

  bn_init_multi(&a_work, &b_work, &t, &r, &v_1, &Q, &p1, &p2, NULL);
  bn_copy(&a_work, a);
  bn_copy(&b_work, b);

  // the Lehmer loop reasons about top limbs as positive magnitudes, so it
  // must run on |a|, |b|; the sign of a is re-applied to u at the end
  a_work.is_neg = false;
  b_work.is_neg = false;

  // 1. [Initialize]
  bn_set_u64(u, 1);
  bn_set_u64(&v_1, 0);

  // 2. [Finished?]
  while (!bn_is_zero(&b_work)) {
    A = 1;
    B = 0;
    C = 0;
    D = 1;

    if (b_work.size > 1) {
      a_hat = a_work.limbs[a_work.size - 1];

      if (b_work.size == a_work.size) {
        b_hat = b_work.limbs[a_work.size - 1];
      } else {
        b_hat = 0;
      }

      // 3. [Test Quotient]
      while (true) {
        __int128 num1 = (__int128)(u64)a_hat + A;
        __int128 den1 = (__int128)(u64)b_hat + C;
        __int128 num2 = (__int128)(u64)a_hat + B;
        __int128 den2 = (__int128)(u64)b_hat + D;

        if (den1 == 0 || den2 == 0) break;

        q = (i64)(num1 / den1);
        if (q != (i64)(num2 / den2)) break;

        // 4. [Euclidian Step]
        T = A - q * C;
        A = C;
        C = T;
        T = B - q * D;
        B = D;
        D = T;
        T = a_hat - q * b_hat;
        a_hat = b_hat;
        b_hat = T;
      }
    }

    // 5. [Multi-precision step]
    if (B == 0) {
      bn_divmod(&Q, &t, &a_work, &b_work);
      bn_copy(&a_work, &b_work);
      bn_copy(&b_work, &t);

      // t = u - Q * v_1
      bn_mul(&t, &Q, &v_1);
      bn_sub(&t, u, &t);
      bn_copy(u, &v_1);
      bn_copy(&v_1, &t);
    } else {
      bn_mul_i64(&p1, &a_work, A);
      bn_mul_i64(&p2, &b_work, B);
      bn_add(&t, &p1, &p2);

      bn_mul_i64(&p1, &a_work, C);
      bn_mul_i64(&p2, &b_work, D);
      bn_add(&r, &p1, &p2);

      bn_copy(&a_work, &t);
      bn_copy(&b_work, &r);

      // Cofactor update:
      // u_new = A*u + B*v_1,  v1_new = C*u + D*v_1
      bn_mul_i64(&p1, u, A);
      bn_mul_i64(&p2, &v_1, B);
      bn_add(&t, &p1, &p2);

      bn_mul_i64(&p1, u, C);
      bn_mul_i64(&p2, &v_1, D);
      bn_add(&r, &p1, &p2);

      bn_copy(u, &t);
      bn_copy(&v_1, &r);
    }
  }

  // Set final GCD: d = a_work (>= 0, the loop ran on magnitudes)
  bn_copy(d, &a_work);

  // u is the Bezout cofactor of |a|; re-apply the sign of a so that
  // u*a + v*b == d holds for the original (signed) inputs
  if (a->is_neg) {
    u->is_neg = !u->is_neg;
  }

  // Calculate final cofactor v = (d - a * u) / b
  bn_mul(&t, a, u);
  bn_sub(&t, d, &t);
  bn_divmod(v, &r, &t, b);

  // Free local bignum temporaries
  bn_free_multi(&a_work, &b_work, &t, &r, &v_1, &Q, &p1, &p2, NULL);
}

// Implement modified version later
void bn_cornacchia(bignum* x, bignum* y, const bignum* p, const bignum* d)
{
  bignum d_work, x_0, p_half, a, b, l, b_sq, t, r, c;
  bn_init_multi(&d_work, &x_0, &p_half, &a, &b, &l, &b_sq, &t, &r, &c, NULL);

  bn_copy(&d_work, d);
  // d > 0
  d_work.is_neg = true;
  i64 k = bn_jacobi(&d_work, p);
  if (k == -1) {
    printf("Equation has no solutions \n");
    goto cleanup;
  }

  // 2. [Compute square root]
  // move -d into Z_p
  bn_add(&d_work, &d_work, p);
  tonelli_shanks(&x_0, &d_work, p);

  printf("sqrt(-d) mod p = ");
  bn_print(&x_0);
  printf("\n");

  bn_rshift(&p_half, p, 1);

  if (bn_cmp(&x_0, &p_half) <= 0) {
    bn_sub(&x_0, p, &x_0);
  }

  bn_copy(&a, p);
  bn_copy(&b, &x_0);
  bn_isqrt(&l, p);

  // 3. [Euclidean algorithm]

  while (bn_cmp(&b, &l) == 1) {
    bn_mod(&r, &a, &b);
    bn_copy(&a, &b);
    bn_copy(&b, &r);
  }

  // 4. [Test solution]
  bn_sqr(&b_sq, &b);
  bn_sub(&t, p, &b_sq);

  // check d | p - b^2
  bn_divmod(&c, &r, &t, d);
  if (!bn_is_zero(&r)) {
    printf("Equation has no solutions \n");
    goto cleanup;
  }

  // check if (p - b^2) / d is not a square
  bn_isqrt(y, &c);
  bn_sqr(&b_sq, y);
  if (bn_cmp(&b_sq, &c) != 0) {
    printf("Equation has no solutions \n");
    goto cleanup;
  }

  bn_copy(x, &b);

cleanup:
  bn_free_multi(&d_work, &x_0, &p_half, &a, &b, &l, &b_sq, &t, &r, &c, NULL);

  return;
}