/*
 * rzmath.c
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

#include "../../include/rabin.h"

// TODO: make faster
// Comp. Number Theory p.29: Algorithm 1.4.10
i64 rz_kronecker(const rz_t* a, const rz_t* b)
{
  static const short tab2[8] = {0, 1, 0, -1, 0, -1, 0, 1};

  // the i64 return type cannot carry an error code; invalid inputs and
  // allocation failures both yield 0 (the "no symbol defined" value)
  if (a == NULL || b == NULL) return 0;

  // 1. [Test b == 0]
  if (rz_is_zero(b)) {
    rz_t a_abs, one;
    rz_init_multi(&a_abs, &one, NULL);
    if (rz_set_u64(&one, 1) != RABIN_SUCCESS) return 0;
    if (rz_copy(&a_abs, a) != RABIN_SUCCESS) {
      rz_clear_multi(&a_abs, &one, NULL);
      return 0;
    }
    a_abs.is_neg = false;
    i64 r = (rz_cmp(&a_abs, &one) == 0) ? 1 : 0;
    rz_clear_multi(&a_abs, &one, NULL);
    return r;
  }

  // 2. [Remove 2's from b]
  if (rz_is_even(a) && rz_is_even(b)) {
    return 0;
  }

  rz_t a_work, b_work, r;
  rz_init_multi(&a_work, &b_work, &r, NULL);
  if (rz_copy(&a_work, a) != RABIN_SUCCESS) goto fail;
  if (rz_copy(&b_work, b) != RABIN_SUCCESS) goto fail;
  u64 v = rz_cnt_trailing_zeros(b);
  if (rz_rshift(&b_work, b, v) != RABIN_SUCCESS) goto fail;

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
    if (rz_is_zero(&a_work)) {
      rz_t one;
      rz_init(&one);
      if (rz_set_u64(&one, 1) != RABIN_SUCCESS) {
        rz_clear(&one);
        k = 0;
        goto end;
      }
      if (rz_cmp(&b_work, &one) == 1) {
        rz_clear(&one);
        k = 0;
        goto end;
      }
      if (rz_cmp(&b_work, &one) == 0) {
        rz_clear(&one);
        goto end;
      }
    }

    u64 v = rz_cnt_trailing_zeros(&a_work);
    if (rz_rshift(&a_work, &a_work, v) != RABIN_SUCCESS) {
      k = 0;
      goto end;
    }
    if ((v & 1) == 1) {
      k *= tab2[b_work.limbs[0] & 7];
    }

    // 4. [Apply reprocity]
    if ((a_work.limbs[0] & 2) & (b_work.limbs[0] & 2)) k = -k;

    if (rz_copy(&r, &a_work) != RABIN_SUCCESS) {
      k = 0;
      goto end;
    }
    r.is_neg = false;
    if (rz_mod(&a_work, &b_work, &r) != RABIN_SUCCESS) {
      k = 0;
      goto end;
    }
    if (rz_copy(&b_work, &r) != RABIN_SUCCESS) {
      k = 0;
      goto end;
    }
  }

fail:
  k = 0;
end:
  rz_clear_multi(&a_work, &b_work, &r, NULL);
  return k;
}

i64 rz_jacobi(const rz_t* a, const rz_t* b) { return rz_kronecker(a, b); }

rabin_err_t rz_tonelli_shanks(rz_t* r, const rz_t* n, const rz_t* p)
{
  if (r == NULL || n == NULL || p == NULL) return RABIN_ERR_NULL_PTR;

  if (rz_is_zero(n)) {
    return rz_set_u64(r, 0);
  }

  if (rz_jacobi(n, p) != 1) {
    return RABIN_ERR_INVALID_ARG;
  }

  rz_t p_minus_one, one, Q, z, M, c, t, R, exp2, tmp2, b2;
  rz_t i, tmp, b, b_exp, j;

  rz_init_multi(&p_minus_one, &one, &Q, &z, &M, &c, &t, &R, &exp2, &tmp2, &b2,
                NULL);
  rz_init_multi(&i, &tmp, &b, &b_exp, &j, NULL);

  rabin_err_t err = rz_set_i64(&one, 1);
  if (err != RABIN_SUCCESS) goto cleanup;
  if ((err = rz_copy(&p_minus_one, p)) != RABIN_SUCCESS) goto cleanup;
  if ((err = rz_sub(&p_minus_one, &p_minus_one, &one)) != RABIN_SUCCESS)
    goto cleanup;
  if ((err = rz_copy(&Q, &p_minus_one)) != RABIN_SUCCESS) goto cleanup;

  u64 S = 0;
  while (rz_is_even(&Q)) {
    if ((err = rz_rshift1(&Q)) != RABIN_SUCCESS) goto cleanup;
    S++;
  }

  // now p - 1 = Q2^S

  if ((err = rz_set_u64(&z, 2)) != RABIN_SUCCESS) goto cleanup;
  while (rz_cmp(&z, p) < 0) {
    if (rz_jacobi(&z, p) == -1) {
      break;
    }
    if ((err = rz_add(&z, &z, &one)) != RABIN_SUCCESS) goto cleanup;
  }

  if ((err = rz_set_u64(&M, S)) != RABIN_SUCCESS) goto cleanup;
  if ((err = rz_mod_exp(&c, &z, &Q, p)) != RABIN_SUCCESS) goto cleanup;
  if ((err = rz_mod_exp(&t, n, &Q, p)) != RABIN_SUCCESS) goto cleanup;

  if ((err = rz_copy(&exp2, &Q)) != RABIN_SUCCESS) goto cleanup;
  if ((err = rz_add_u64(&exp2, &exp2, 1)) != RABIN_SUCCESS) goto cleanup;
  if ((err = rz_rshift1(&exp2)) != RABIN_SUCCESS) goto cleanup;

  if ((err = rz_mod_exp(&R, n, &exp2, p)) != RABIN_SUCCESS) goto cleanup;

  while (1) {
    if (rz_is_zero(&t)) {
      if ((err = rz_set_u64(r, 0)) != RABIN_SUCCESS) goto cleanup;
      err = RABIN_SUCCESS;
      break;
    }

    if (rz_cmp(&t, &one) == 0) {
      if ((err = rz_copy(r, &R)) != RABIN_SUCCESS) goto cleanup;
      err = RABIN_SUCCESS;
      break;
    }

    if ((err = rz_set_u64(&i, 0)) != RABIN_SUCCESS) goto cleanup;
    if ((err = rz_copy(&tmp, &t)) != RABIN_SUCCESS) goto cleanup;

    while (!rz_is_eq_i64(&tmp, 1) && rz_cmp(&i, &M) < 0) {
      if ((err = rz_mul(&tmp2, &tmp, &tmp)) != RABIN_SUCCESS) goto cleanup;
      if ((err = rz_mod(&tmp, &tmp2, p)) != RABIN_SUCCESS) goto cleanup;
      if ((err = rz_add(&i, &i, &one)) != RABIN_SUCCESS) goto cleanup;
    }

    if (rz_cmp(&i, &M) == 0) {
      err = RABIN_ERR_INVALID_ARG;
      goto cleanup;
    }

    if ((err = rz_copy(&b_exp, &M)) != RABIN_SUCCESS) goto cleanup;
    if ((err = rz_sub(&b_exp, &b_exp, &i)) != RABIN_SUCCESS) goto cleanup;
    if ((err = rz_sub(&b_exp, &b_exp, &one)) != RABIN_SUCCESS) goto cleanup;

    if ((err = rz_copy(&b, &c)) != RABIN_SUCCESS) goto cleanup;
    if ((err = rz_set_u64(&j, 0)) != RABIN_SUCCESS) goto cleanup;

    while (rz_cmp(&j, &b_exp) < 0) {
      if ((err = rz_mul(&b2, &b, &b)) != RABIN_SUCCESS) goto cleanup;
      if ((err = rz_mod(&b, &b2, p)) != RABIN_SUCCESS) goto cleanup;
      if ((err = rz_add(&j, &j, &one)) != RABIN_SUCCESS) goto cleanup;
    }

    if ((err = rz_copy(&M, &i)) != RABIN_SUCCESS) goto cleanup;

    if ((err = rz_mul(&c, &b, &b)) != RABIN_SUCCESS) goto cleanup;
    if ((err = rz_mod(&c, &c, p)) != RABIN_SUCCESS) goto cleanup;

    if ((err = rz_mul(&t, &t, &c)) != RABIN_SUCCESS) goto cleanup;
    if ((err = rz_mod(&t, &t, p)) != RABIN_SUCCESS) goto cleanup;

    if ((err = rz_mul(&R, &R, &b)) != RABIN_SUCCESS) goto cleanup;
    if ((err = rz_mod(&R, &R, p)) != RABIN_SUCCESS) goto cleanup;
  }

cleanup:
  rz_clear_multi(&p_minus_one, &one, &Q, &z, &M, &c, &t, &R, &exp2, &tmp2, &b2,
                 NULL);
  rz_clear_multi(&i, &tmp, &b, &b_exp, &j, NULL);
  return err;
}

rabin_err_t rz_gcd(rz_t* d, const rz_t* a, const rz_t* b)
{
  if (d == NULL || a == NULL || b == NULL) return RABIN_ERR_NULL_PTR;

  if (rz_is_zero(a)) {
    rabin_err_t err = rz_copy(d, b);
    if (err != RABIN_SUCCESS) return err;
    d->is_neg = false;
    return RABIN_SUCCESS;
  }
  if (rz_is_zero(b)) {
    rabin_err_t err = rz_copy(d, a);
    if (err != RABIN_SUCCESS) return err;
    d->is_neg = false;
    return RABIN_SUCCESS;
  }

  if (a->size >= 16 || b->size >= 16) {
    return rz_gcd_lehmer(d, a, b);
  }

  rz_t tmp_a, tmp_b, rem;
  rz_init_multi(&tmp_a, &tmp_b, &rem, NULL);
  rabin_err_t err = rz_copy(&tmp_a, a);
  if (err != RABIN_SUCCESS) goto out;
  if ((err = rz_copy(&tmp_b, b)) != RABIN_SUCCESS) goto out;

  rz_t* u = &tmp_a;
  rz_t* v = &tmp_b;
  rz_t* r = &rem;

  while (!rz_is_zero(v)) {
    // r = u % v
    if ((err = rz_mod(r, u, v)) != RABIN_SUCCESS) goto out;

    rz_t* temp = u;
    u = v;
    v = r;
    r = temp;
  }

  // The GCD is left in 'u'
  if ((err = rz_copy(d, u)) != RABIN_SUCCESS) goto out;
  d->is_neg = false;

  err = RABIN_SUCCESS;
out:
  rz_clear(&tmp_a);
  rz_clear(&tmp_b);
  rz_clear(&rem);
  return err;
}

rabin_err_t rz_lcm(rz_t* l, const rz_t* a, const rz_t* b)
{
  if (l == NULL || a == NULL || b == NULL) return RABIN_ERR_NULL_PTR;

  rz_t d;
  rz_init(&d);

  rabin_err_t err = rz_gcd(&d, a, b);
  if (err != RABIN_SUCCESS) goto out;
  if ((err = rz_mul(l, a, b)) != RABIN_SUCCESS) goto out;
  if ((err = rz_div(l, l, &d)) != RABIN_SUCCESS) goto out;

  err = RABIN_SUCCESS;
out:
  rz_clear(&d);
  return err;
}

// Computational number theory p.15
rabin_err_t rz_gcd_binary(rz_t* d, const rz_t* a, const rz_t* b)
{
  if (d == NULL || a == NULL || b == NULL) return RABIN_ERR_NULL_PTR;

  rz_t a_work, b_work, r;
  rz_init_multi(&a_work, &b_work, &r, NULL);
  rabin_err_t err = rz_copy(&a_work, a);
  if (err != RABIN_SUCCESS) goto out;
  if ((err = rz_copy(&b_work, b)) != RABIN_SUCCESS) goto out;

  // 1 [Reduce size]
  if (rz_cmp(&a_work, &b_work) == -1) {
    if ((err = rz_swap(&a_work, &b_work)) != RABIN_SUCCESS) goto out;
  }

  if (rz_is_zero(&b_work)) {
    if ((err = rz_copy(d, &a_work)) != RABIN_SUCCESS) goto out;
    d->is_neg = false;
    err = RABIN_SUCCESS;
    goto out;
  }

  if ((err = rz_mod(&r, &a_work, &b_work)) != RABIN_SUCCESS) goto out;
  if ((err = rz_copy(&a_work, &b_work)) != RABIN_SUCCESS) goto out;
  if ((err = rz_copy(&b_work, &r)) != RABIN_SUCCESS) goto out;

  // 2. [Compute power of 2] & 3
  if (rz_is_zero(&b_work)) {
    if ((err = rz_copy(d, &a_work)) != RABIN_SUCCESS) goto out;
    d->is_neg = false;
    err = RABIN_SUCCESS;
    goto out;
  }

  u64 k_a = rz_cnt_trailing_zeros(&a_work);
  u64 k_b = rz_cnt_trailing_zeros(&b_work);
  u64 k = RZ_MIN(k_a, k_b);

  if ((err = rz_rshift(&a_work, &a_work, k_a)) != RABIN_SUCCESS) goto out;
  if ((err = rz_rshift(&b_work, &b_work, k_b)) != RABIN_SUCCESS) goto out;

  // 4 & 5. [Subtract and Loop]
  while (true) {
    int cmp = rz_cmp(&a_work, &b_work);

    if (cmp == 0) {
      // If t = 0, output 2^k * a and terminate
      if ((err = rz_lshift(d, &a_work, k)) != RABIN_SUCCESS) goto out;
      d->is_neg = false;
      err = RABIN_SUCCESS;
      break;
    }

    if (cmp > 0) {
      // If a > b, t > 0. Set t = a - b, remove powers of 2, set a = t
      if ((err = rz_sub(&r, &a_work, &b_work)) != RABIN_SUCCESS) goto out;
      if ((err = rz_rshift(&a_work, &r, rz_cnt_trailing_zeros(&r))) !=
          RABIN_SUCCESS)
        goto out;
    } else {
      // If a < b, t < 0. Set -t = b - a, remove powers of 2, set b = -t
      if ((err = rz_sub(&r, &b_work, &a_work)) != RABIN_SUCCESS) goto out;
      if ((err = rz_rshift(&b_work, &r, rz_cnt_trailing_zeros(&r))) !=
          RABIN_SUCCESS)
        goto out;
    }
  }

out:
  rz_clear_multi(&a_work, &b_work, &r, NULL);
  return err;
}

rabin_err_t rz_gcd_lehmer(rz_t* d, const rz_t* a, const rz_t* b)
{
  if (d == NULL || a == NULL || b == NULL) return RABIN_ERR_NULL_PTR;

  i64 a_hat, b_hat, A, B, C, D, T, q;
  rz_t a_work, b_work, t, r, p1, p2;

  // 1. [Initialize]
  rz_init_multi(&a_work, &b_work, &t, &r, &p1, &p2, NULL);
  rabin_err_t err = rz_copy(&a_work, a);
  if (err != RABIN_SUCCESS) goto out;
  if ((err = rz_copy(&b_work, b)) != RABIN_SUCCESS) goto out;
  a_work.is_neg = false;
  b_work.is_neg = false;

  if (rz_cmp(&a_work, &b_work) < 0) {
    if ((err = rz_swap(&a_work, &b_work)) != RABIN_SUCCESS) goto out;
  }

  while (!rz_is_zero(&b_work)) {
    if (a_work.size - b_work.size > 1) {
      if ((err = rz_mod(&t, &a_work, &b_work)) != RABIN_SUCCESS) goto out;
      if ((err = rz_copy(&a_work, &b_work)) != RABIN_SUCCESS) goto out;
      if ((err = rz_copy(&b_work, &t)) != RABIN_SUCCESS) goto out;
      continue;
    }

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

      // 2. [Test Quotient]
      while (true) {
        __int128 num1 = (__int128)(u64)a_hat + A;
        __int128 den1 = (__int128)(u64)b_hat + C;
        __int128 num2 = (__int128)(u64)a_hat + B;
        __int128 den2 = (__int128)(u64)b_hat + D;

        if (den1 == 0 || den2 == 0) break;

        q = (i64)(num1 / den1);

        if (q == 0) break;

        if (q != (i64)(num2 / den2)) break;

        // 3. [Euclidian Step]
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

    // 4. [Multi-precision step]
    if (B == 0) {
      if ((err = rz_mod(&t, &a_work, &b_work)) != RABIN_SUCCESS) goto out;
      if ((err = rz_copy(&a_work, &b_work)) != RABIN_SUCCESS) goto out;
      if ((err = rz_copy(&b_work, &t)) != RABIN_SUCCESS) goto out;
    } else {
      if (A > 0) {
        // A > 0, B <= 0, C <= 0, D > 0
        if ((err = rz_mul_i64(&p1, &a_work, (u64)A)) != RABIN_SUCCESS) goto out;
        if ((err = rz_mul_i64(&p2, &b_work, (u64)(-B))) != RABIN_SUCCESS)
          goto out;
        if ((err = rz_sub(&t, &p1, &p2)) != RABIN_SUCCESS)
          goto out;  // t = A*a - |B|*b

        if ((err = rz_mul_i64(&p1, &a_work, (u64)(-C))) != RABIN_SUCCESS)
          goto out;
        if ((err = rz_mul_i64(&p2, &b_work, (u64)D)) != RABIN_SUCCESS) goto out;
        if ((err = rz_sub(&r, &p2, &p1)) != RABIN_SUCCESS)
          goto out;  // r = D*b - |C|*a
      } else {
        // A <= 0, B > 0, C > 0, D <= 0
        if ((err = rz_mul_i64(&p1, &a_work, (u64)(-A))) != RABIN_SUCCESS)
          goto out;
        if ((err = rz_mul_i64(&p2, &b_work, (u64)B)) != RABIN_SUCCESS) goto out;
        if ((err = rz_sub(&t, &p2, &p1)) != RABIN_SUCCESS)
          goto out;  // t = B*b - |A|*a

        if ((err = rz_mul_i64(&p1, &a_work, (u64)C)) != RABIN_SUCCESS) goto out;
        if ((err = rz_mul_i64(&p2, &b_work, (u64)(-D))) != RABIN_SUCCESS)
          goto out;
        if ((err = rz_sub(&r, &p1, &p2)) != RABIN_SUCCESS)
          goto out;  // r = C*a - |D|*b
      }

      if ((err = rz_copy(&a_work, &t)) != RABIN_SUCCESS) goto out;
      if ((err = rz_copy(&b_work, &r)) != RABIN_SUCCESS) goto out;
    }
  }

  if ((err = rz_copy(d, &a_work)) != RABIN_SUCCESS) goto out;
  d->is_neg = false;
  err = RABIN_SUCCESS;
out:
  rz_clear_multi(&a_work, &b_work, &t, &r, &p1, &p2, NULL);
  return err;
}
// Computational number theory p.16
rabin_err_t rz_gcd_extended(rz_t* u, rz_t* v, rz_t* d, const rz_t* a,
                            const rz_t* b)
{
  if (u == NULL || v == NULL || d == NULL || a == NULL || b == NULL)
    return RABIN_ERR_NULL_PTR;

  rz_t v_1, v_3, t_1, t_3, temp;
  rz_init_multi(&v_1, &v_3, &t_1, &t_3, &temp, NULL);

  rabin_err_t err = rz_set_u64(u, 1);
  if (err != RABIN_SUCCESS) goto out;
  if ((err = rz_copy(d, a)) != RABIN_SUCCESS) goto out;

  if (rz_is_zero(b)) {
    if ((err = rz_set_u64(v, 0)) != RABIN_SUCCESS) goto out;
    // normalize: gcd is non-negative (GMP convention)
    if (d->is_neg) {
      d->is_neg = false;
      u->is_neg = !u->is_neg;
    }
    err = RABIN_SUCCESS;
    goto out;
  }

  if ((err = rz_set_u64(&v_1, 0)) != RABIN_SUCCESS) goto out;
  if ((err = rz_copy(&v_3, b)) != RABIN_SUCCESS) goto out;

  while (!rz_is_zero(&v_3)) {
    if ((err = rz_divmod(&temp, &t_3, d, &v_3)) != RABIN_SUCCESS) goto out;

    // t_1 = u - qv_1
    if ((err = rz_mul(&temp, &temp, &v_1)) != RABIN_SUCCESS) goto out;
    if ((err = rz_copy(&t_1, u)) != RABIN_SUCCESS) goto out;
    if ((err = rz_sub(&t_1, &t_1, &temp)) != RABIN_SUCCESS) goto out;

    if ((err = rz_copy(u, &v_1)) != RABIN_SUCCESS) goto out;
    if ((err = rz_copy(d, &v_3)) != RABIN_SUCCESS) goto out;
    if ((err = rz_copy(&v_1, &t_1)) != RABIN_SUCCESS) goto out;
    if ((err = rz_copy(&v_3, &t_3)) != RABIN_SUCCESS) goto out;
  }

  // normalize: gcd is non-negative (GMP convention); flip d and u
  // together so the Bezout identity u*a + v*b == d is preserved
  if (d->is_neg) {
    d->is_neg = false;
    u->is_neg = !u->is_neg;
  }

  // v = (d - au) / b
  if ((err = rz_mul(&temp, a, u)) != RABIN_SUCCESS) goto out;
  if ((err = rz_copy(v, d)) != RABIN_SUCCESS) goto out;
  if ((err = rz_sub(v, v, &temp)) != RABIN_SUCCESS) goto out;
  if ((err = rz_div(v, v, b)) != RABIN_SUCCESS) goto out;

  err = RABIN_SUCCESS;
out:
  rz_clear_multi(&v_1, &v_3, &t_1, &t_3, &temp, NULL);
  return err;
}

rabin_err_t rz_gcd_extended_lehmer(rz_t* u, rz_t* v, rz_t* d, const rz_t* a,
                                   const rz_t* b)
{
  if (u == NULL || v == NULL || d == NULL || a == NULL || b == NULL)
    return RABIN_ERR_NULL_PTR;

  i64 a_hat, b_hat, A, B, C, D, T, q;
  rz_t a_work, b_work, t, r, v_1, Q, p1, p2;

  rz_init_multi(&a_work, &b_work, &t, &r, &v_1, &Q, &p1, &p2, NULL);
  rabin_err_t err = rz_copy(&a_work, a);
  if (err != RABIN_SUCCESS) goto out;
  if ((err = rz_copy(&b_work, b)) != RABIN_SUCCESS) goto out;

  // the Lehmer loop reasons about top limbs as positive magnitudes, so it
  // must run on |a|, |b|; the sign of a is re-applied to u at the end
  a_work.is_neg = false;
  b_work.is_neg = false;

  // zero-argument guards: the loop below reads the top limb of a_work,
  // which does not exist for a size-0 rz_t
  if (rz_is_zero(a)) {
    // u = 0, v = sign(b), d = |b|
    if ((err = rz_set_u64(u, 0)) != RABIN_SUCCESS) goto out;
    if (b->is_neg) {
      if ((err = rz_set_i64(v, -1)) != RABIN_SUCCESS) goto out;
    } else {
      if ((err = rz_set_u64(v, 1)) != RABIN_SUCCESS) goto out;
    }
    if ((err = rz_copy(d, b)) != RABIN_SUCCESS) goto out;
    d->is_neg = false;
    err = RABIN_SUCCESS;
    goto out;
  }
  if (rz_is_zero(b)) {
    // u = sign(a), v = 0, d = |a|
    if (a->is_neg) {
      if ((err = rz_set_i64(u, -1)) != RABIN_SUCCESS) goto out;
    } else {
      if ((err = rz_set_u64(u, 1)) != RABIN_SUCCESS) goto out;
    }
    if ((err = rz_set_u64(v, 0)) != RABIN_SUCCESS) goto out;
    if ((err = rz_copy(d, a)) != RABIN_SUCCESS) goto out;
    d->is_neg = false;
    err = RABIN_SUCCESS;
    goto out;
  }

  // 1. [Initialize]
  if ((err = rz_set_u64(u, 1)) != RABIN_SUCCESS) goto out;
  if ((err = rz_set_u64(&v_1, 0)) != RABIN_SUCCESS) goto out;

  // 2. [Finished?]
  while (!rz_is_zero(&b_work)) {
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
      if ((err = rz_divmod(&Q, &t, &a_work, &b_work)) != RABIN_SUCCESS)
        goto out;
      if ((err = rz_copy(&a_work, &b_work)) != RABIN_SUCCESS) goto out;
      if ((err = rz_copy(&b_work, &t)) != RABIN_SUCCESS) goto out;

      // t = u - Q * v_1
      if ((err = rz_mul(&t, &Q, &v_1)) != RABIN_SUCCESS) goto out;
      if ((err = rz_sub(&t, u, &t)) != RABIN_SUCCESS) goto out;
      if ((err = rz_copy(u, &v_1)) != RABIN_SUCCESS) goto out;
      if ((err = rz_copy(&v_1, &t)) != RABIN_SUCCESS) goto out;
    } else {
      if ((err = rz_mul_i64(&p1, &a_work, A)) != RABIN_SUCCESS) goto out;
      if ((err = rz_mul_i64(&p2, &b_work, B)) != RABIN_SUCCESS) goto out;
      if ((err = rz_add(&t, &p1, &p2)) != RABIN_SUCCESS) goto out;

      if ((err = rz_mul_i64(&p1, &a_work, C)) != RABIN_SUCCESS) goto out;
      if ((err = rz_mul_i64(&p2, &b_work, D)) != RABIN_SUCCESS) goto out;
      if ((err = rz_add(&r, &p1, &p2)) != RABIN_SUCCESS) goto out;

      if ((err = rz_copy(&a_work, &t)) != RABIN_SUCCESS) goto out;
      if ((err = rz_copy(&b_work, &r)) != RABIN_SUCCESS) goto out;

      // Cofactor update:
      // u_new = A*u + B*v_1,  v1_new = C*u + D*v_1
      if ((err = rz_mul_i64(&p1, u, A)) != RABIN_SUCCESS) goto out;
      if ((err = rz_mul_i64(&p2, &v_1, B)) != RABIN_SUCCESS) goto out;
      if ((err = rz_add(&t, &p1, &p2)) != RABIN_SUCCESS) goto out;

      if ((err = rz_mul_i64(&p1, u, C)) != RABIN_SUCCESS) goto out;
      if ((err = rz_mul_i64(&p2, &v_1, D)) != RABIN_SUCCESS) goto out;
      if ((err = rz_add(&r, &p1, &p2)) != RABIN_SUCCESS) goto out;

      if ((err = rz_copy(u, &t)) != RABIN_SUCCESS) goto out;
      if ((err = rz_copy(&v_1, &r)) != RABIN_SUCCESS) goto out;
    }
  }

  // Set final GCD: d = a_work (>= 0, the loop ran on magnitudes)
  if ((err = rz_copy(d, &a_work)) != RABIN_SUCCESS) goto out;
  d->is_neg = false;

  // u is the Bezout cofactor of |a|; re-apply the sign of a so that
  // u*a + v*b == d holds for the original (signed) inputs
  if (a->is_neg) {
    u->is_neg = !u->is_neg;
  }

  // Calculate final cofactor v = (d - a * u) / b
  if ((err = rz_mul(&t, a, u)) != RABIN_SUCCESS) goto out;
  if ((err = rz_sub(&t, d, &t)) != RABIN_SUCCESS) goto out;
  if ((err = rz_divmod(v, &r, &t, b)) != RABIN_SUCCESS) goto out;

  err = RABIN_SUCCESS;
out:
  // Free local rz_t temporaries
  rz_clear_multi(&a_work, &b_work, &t, &r, &v_1, &Q, &p1, &p2, NULL);
  return err;
}

// Implement modified version later
rabin_err_t rz_cornacchia(rz_t* x, rz_t* y, const rz_t* p, const rz_t* d)
{
  if (x == NULL || y == NULL || p == NULL || d == NULL)
    return RABIN_ERR_NULL_PTR;

  rz_t d_work, x_0, p_half, a, b, l, b_sq, t, r, c;
  rz_init_multi(&d_work, &x_0, &p_half, &a, &b, &l, &b_sq, &t, &r, &c, NULL);

  rabin_err_t err = rz_copy(&d_work, d);
  if (err != RABIN_SUCCESS) goto cleanup;

  // d > 0
  d_work.is_neg = true;
  i64 k = rz_jacobi(&d_work, p);
  if (k == -1) {
    err = RABIN_ERR_INVALID_ARG;
    goto cleanup;
  }

  // 2. [Compute square root]
  // move -d into Z_p
  if ((err = rz_add(&d_work, &d_work, p)) != RABIN_SUCCESS) goto cleanup;
  if ((err = rz_tonelli_shanks(&x_0, &d_work, p)) != RABIN_SUCCESS)
    goto cleanup;

  if ((err = rz_rshift(&p_half, p, 1)) != RABIN_SUCCESS) goto cleanup;

  if (rz_cmp(&x_0, &p_half) <= 0) {
    if ((err = rz_sub(&x_0, p, &x_0)) != RABIN_SUCCESS) goto cleanup;
  }

  if ((err = rz_copy(&a, p)) != RABIN_SUCCESS) goto cleanup;
  if ((err = rz_copy(&b, &x_0)) != RABIN_SUCCESS) goto cleanup;
  if ((err = rz_isqrt(&l, p)) != RABIN_SUCCESS) goto cleanup;

  // 3. [Euclidean algorithm]

  while (rz_cmp(&b, &l) == 1) {
    if ((err = rz_mod(&r, &a, &b)) != RABIN_SUCCESS) goto cleanup;
    if ((err = rz_copy(&a, &b)) != RABIN_SUCCESS) goto cleanup;
    if ((err = rz_copy(&b, &r)) != RABIN_SUCCESS) goto cleanup;
  }

  // 4. [Test solution]
  if ((err = rz_sqr(&b_sq, &b)) != RABIN_SUCCESS) goto cleanup;
  if ((err = rz_sub(&t, p, &b_sq)) != RABIN_SUCCESS) goto cleanup;

  // check d | p - b^2
  if ((err = rz_divmod(&c, &r, &t, d)) != RABIN_SUCCESS) goto cleanup;
  if (!rz_is_zero(&r)) {
    err = RABIN_ERR_INVALID_ARG;
    goto cleanup;
  }

  // check if (p - b^2) / d is not a square
  if ((err = rz_isqrt(y, &c)) != RABIN_SUCCESS) goto cleanup;
  if ((err = rz_sqr(&b_sq, y)) != RABIN_SUCCESS) goto cleanup;
  if (rz_cmp(&b_sq, &c) != 0) {
    err = RABIN_ERR_INVALID_ARG;
    goto cleanup;
  }

  if ((err = rz_copy(x, &b)) != RABIN_SUCCESS) goto cleanup;

  err = RABIN_SUCCESS;
cleanup:
  rz_clear_multi(&d_work, &x_0, &p_half, &a, &b, &l, &b_sq, &t, &r, &c, NULL);

  return err;
}
