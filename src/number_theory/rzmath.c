/*
 * rzmath.c
 *
 * Classical number-theory routines.
 *
 * This file implements the Jacobi and Kronecker symbols, the
 * Tonelli-Shanks square root algorithm modulo a prime, binary (Stein)
 * and Lehmer GCD variants, and Cornacchia's algorithm.
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
// True when |a| == 1; the i64 symbol helpers need this test in several
// places and cannot signal allocation failures, so a failed copy maps
// to false (the caller then yields the 0 "undefined" result).
static bool rz_abs_is_one(const rz_t* a)
{
  rz_t a_abs, one;
  rz_init_multi(&a_abs, &one, NULL);
  if (rz_set_u64(&one, 1) != RABIN_SUCCESS) {
    rz_clear_multi(&a_abs, &one, NULL);
    return false;
  }
  if (rz_copy(&a_abs, a) != RABIN_SUCCESS) {
    rz_clear_multi(&a_abs, &one, NULL);
    return false;
  }
  a_abs.is_neg = false;
  bool res = (rz_cmp(&a_abs, &one) == 0);
  rz_clear_multi(&a_abs, &one, NULL);
  return res;
}

i64 rz_kronecker(const rz_t* a, const rz_t* b)
{
  static const short tab2[8] = {0, 1, 0, -1, 0, -1, 0, 1};

  // the i64 return type cannot carry an error code; invalid inputs and
  // allocation failures both yield 0 (the "no symbol defined" value)
  if (a == NULL || b == NULL) return 0;

  // 1. [Test b == 0] the Kronecker symbol (a/0) is 1 for a == +-1
  // and 0 otherwise
  if (rz_is_zero(b)) {
    return (rz_abs_is_one(a)) ? 1 : 0;
  }

  // 1a. [Test a == 0] the Kronecker symbol (0/b) is 1 for b == +-1
  // and 0 otherwise
  if (rz_is_zero(a)) {
    return (rz_abs_is_one(b)) ? 1 : 0;
  }

  // 2. [Remove 2's from b]
  if (rz_is_even(a) && rz_is_even(b)) {
    return 0;
  }

  rz_t a_work, b_work, r;
  rz_init_multi(&a_work, &b_work, &r, NULL);
  i64 k = 1;
  if (rz_copy(&a_work, a) != RABIN_SUCCESS) {
    k = 0;
    goto cleanup;
  }
  if (rz_copy(&b_work, b) != RABIN_SUCCESS) {
    k = 0;
    goto cleanup;
  }
  u64 v = rz_cnt_trailing_zeros(b);
  if (rz_rshift(&b_work, b, v) != RABIN_SUCCESS) {
    k = 0;
    goto cleanup;
  }

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
      // gcd has collapsed to |b|; the symbol is k only when b == 1
      rz_t one;
      rz_init(&one);
      if (rz_set_u64(&one, 1) != RABIN_SUCCESS) {
        rz_clear(&one);
        k = 0;
        goto cleanup;
      }
      k = (rz_cmp(&b_work, &one) == 0) ? k : 0;
      rz_clear(&one);
      break;
    }

    u64 v = rz_cnt_trailing_zeros(&a_work);
    if (rz_rshift(&a_work, &a_work, v) != RABIN_SUCCESS) {
      k = 0;
      goto cleanup;
    }
    if ((v & 1) == 1) {
      k *= tab2[b_work.limbs[0] & 7];
    }

    // 4. [Apply reprocity]
    if ((a_work.limbs[0] & 2) & (b_work.limbs[0] & 2)) k = -k;

    if (rz_copy(&r, &a_work) != RABIN_SUCCESS) {
      k = 0;
      goto cleanup;
    }
    r.is_neg = false;
    if (rz_mod(&a_work, &b_work, &r) != RABIN_SUCCESS) {
      k = 0;
      goto cleanup;
    }
    if (rz_copy(&b_work, &r) != RABIN_SUCCESS) {
      k = 0;
      goto cleanup;
    }
  }

cleanup:
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
  TRY(rz_copy(&p_minus_one, p));
  TRY(rz_sub(&p_minus_one, &p_minus_one, &one));
  TRY(rz_copy(&Q, &p_minus_one));

  u64 S = 0;
  while (rz_is_even(&Q)) {
    TRY(rz_rshift1(&Q));
    S++;
  }

  // now p - 1 = Q2^S

  TRY(rz_set_u64(&z, 2));
  while (rz_cmp(&z, p) < 0) {
    if (rz_jacobi(&z, p) == -1) {
      break;
    }
    TRY(rz_add(&z, &z, &one));
  }

  TRY(rz_set_u64(&M, S));
  TRY(rz_mod_exp(&c, &z, &Q, p));
  TRY(rz_mod_exp(&t, n, &Q, p));

  TRY(rz_copy(&exp2, &Q));
  TRY(rz_add_u64(&exp2, &exp2, 1));
  TRY(rz_rshift1(&exp2));

  TRY(rz_mod_exp(&R, n, &exp2, p));

  while (1) {
    if (rz_is_zero(&t)) {
      TRY(rz_set_u64(r, 0));
      err = RABIN_SUCCESS;
      break;
    }

    if (rz_cmp(&t, &one) == 0) {
      TRY(rz_copy(r, &R));
      err = RABIN_SUCCESS;
      break;
    }

    TRY(rz_set_u64(&i, 0));
    TRY(rz_copy(&tmp, &t));

    while (!rz_is_eq_i64(&tmp, 1) && rz_cmp(&i, &M) < 0) {
      TRY(rz_mul(&tmp2, &tmp, &tmp));
      TRY(rz_mod(&tmp, &tmp2, p));
      TRY(rz_add(&i, &i, &one));
    }

    if (rz_cmp(&i, &M) == 0) {
      err = RABIN_ERR_INVALID_ARG;
      goto cleanup;
    }

    TRY(rz_copy(&b_exp, &M));
    TRY(rz_sub(&b_exp, &b_exp, &i));
    TRY(rz_sub(&b_exp, &b_exp, &one));

    TRY(rz_copy(&b, &c));
    TRY(rz_set_u64(&j, 0));

    while (rz_cmp(&j, &b_exp) < 0) {
      TRY(rz_mul(&b2, &b, &b));
      TRY(rz_mod(&b, &b2, p));
      TRY(rz_add(&j, &j, &one));
    }

    TRY(rz_copy(&M, &i));

    TRY(rz_mul(&c, &b, &b));
    TRY(rz_mod(&c, &c, p));

    TRY(rz_mul(&t, &t, &c));
    TRY(rz_mod(&t, &t, p));

    TRY(rz_mul(&R, &R, &b));
    TRY(rz_mod(&R, &R, p));
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
  if (err != RABIN_SUCCESS) goto cleanup;
  TRY(rz_copy(&tmp_b, b));

  rz_t* u = &tmp_a;
  rz_t* v = &tmp_b;
  rz_t* r = &rem;

  while (!rz_is_zero(v)) {
    // r = u % v
    TRY(rz_mod(r, u, v));

    rz_t* temp = u;
    u = v;
    v = r;
    r = temp;
  }

  // The GCD is left in 'u'
  TRY(rz_copy(d, u));
  d->is_neg = false;

  err = RABIN_SUCCESS;
cleanup:
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
  if (err != RABIN_SUCCESS) goto cleanup;
  TRY(rz_mul(l, a, b));
  TRY(rz_div(l, l, &d));
  // the LCM is the principal (non-negative) value
  l->is_neg = false;

  err = RABIN_SUCCESS;
cleanup:
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
  if (err != RABIN_SUCCESS) goto cleanup;
  TRY(rz_copy(&b_work, b));

  // 1 [Reduce size]
  if (rz_cmp(&a_work, &b_work) == -1) {
    TRY(rz_swap(&a_work, &b_work));
  }

  if (rz_is_zero(&b_work)) {
    TRY(rz_copy(d, &a_work));
    d->is_neg = false;
    err = RABIN_SUCCESS;
    goto cleanup;
  }

  TRY(rz_mod(&r, &a_work, &b_work));
  TRY(rz_copy(&a_work, &b_work));
  TRY(rz_copy(&b_work, &r));

  // 2. [Compute power of 2] & 3
  if (rz_is_zero(&b_work)) {
    TRY(rz_copy(d, &a_work));
    d->is_neg = false;
    err = RABIN_SUCCESS;
    goto cleanup;
  }

  u64 k_a = rz_cnt_trailing_zeros(&a_work);
  u64 k_b = rz_cnt_trailing_zeros(&b_work);
  u64 k = RZ_MIN(k_a, k_b);

  TRY(rz_rshift(&a_work, &a_work, k_a));
  TRY(rz_rshift(&b_work, &b_work, k_b));

  // 4 & 5. [Subtract and Loop]
  while (true) {
    int cmp = rz_cmp(&a_work, &b_work);

    if (cmp == 0) {
      // If t = 0, output 2^k * a and terminate
      TRY(rz_lshift(d, &a_work, k));
      d->is_neg = false;
      err = RABIN_SUCCESS;
      break;
    }

    if (cmp > 0) {
      // If a > b, t > 0. Set t = a - b, remove powers of 2, set a = t
      TRY(rz_sub(&r, &a_work, &b_work));
      TRY(rz_rshift(&a_work, &r, rz_cnt_trailing_zeros(&r)));
    } else {
      // If a < b, t < 0. Set -t = b - a, remove powers of 2, set b = -t
      TRY(rz_sub(&r, &b_work, &a_work));
      TRY(rz_rshift(&b_work, &r, rz_cnt_trailing_zeros(&r)));
    }
  }

cleanup:
  rz_clear_multi(&a_work, &b_work, &r, NULL);
  return err;
}

// Runs the Lehmer [Test Quotient] inner loop shared by rz_gcd_lehmer()
// and rz_gcd_extended_lehmer(). Accumulates the 2x2 transformation
// matrix [[A, B], [C, D]] from the top limbs until no further Euclidian
// step is justified, returning true only if a non-identity, in-range
// matrix was accumulated; false means the caller must fall back to a
// plain remainder step.
static bool lehmer_quotient_loop(const rz_t* a_work, const rz_t* b_work, i64* A,
                                 i64* B, i64* C, i64* D)
{
  *A = 1;
  *B = 0;
  *C = 0;
  *D = 1;

  if (b_work->size <= 1) return false;

  // Truncate the top limbs by one bit so the tracked values always lie
  // in [0, 2^63) on a common scale; this is what keeps the Euclidian
  // step products inside a signed 64-bit word (Knuth, TAOCP II,
  // 4.5.2 Alg. B)
  i64 a_hat = (i64)(a_work->limbs[a_work->size - 1] >> 1);
  i64 b_hat = (b_work->size == a_work->size)
                  ? (i64)(b_work->limbs[b_work->size - 1] >> 1)
                  : 0;

  while (true) {
    // Only accept steps that keep both transformed values positive
    __int128 den1 = (__int128)b_hat + *C;
    __int128 den2 = (__int128)b_hat + *D;

    if (den1 <= 0 || den2 <= 0) break;

    __int128 q1 = ((__int128)a_hat + *A) / den1;
    __int128 q2 = ((__int128)a_hat + *B) / den2;

    if (q1 != q2 || q1 <= 0) break;

    // Euclidian step, computed in 128 bits
    __int128 nA = (__int128)*A - q1 * *C;
    __int128 nB = (__int128)*B - q1 * *D;
    __int128 nh = (__int128)a_hat - q1 * b_hat;

    // Discard a step whose entries leave the signed 64-bit range; the
    // caller falls back to a plain remainder step instead
    if (nA > (__int128)INT64_MAX || nA < (__int128)INT64_MIN ||
        nB > (__int128)INT64_MAX || nB < (__int128)INT64_MIN ||
        nh > (__int128)INT64_MAX || nh < (__int128)INT64_MIN) {
      break;
    }

    *A = *C;
    *C = (i64)nA;
    *B = *D;
    *D = (i64)nB;
    a_hat = b_hat;
    b_hat = (i64)nh;
  }

  // B == 0 means no Euclidian step was accepted
  return *B != 0;
}

rabin_err_t rz_gcd_lehmer(rz_t* d, const rz_t* a, const rz_t* b)
{
  if (d == NULL || a == NULL || b == NULL) return RABIN_ERR_NULL_PTR;

  rz_t a_work, b_work, t, r, p1, p2;

  // 1. [Initialize]
  rz_init_multi(&a_work, &b_work, &t, &r, &p1, &p2, NULL);
  rabin_err_t err = rz_copy(&a_work, a);
  if (err != RABIN_SUCCESS) goto cleanup;
  TRY(rz_copy(&b_work, b));
  a_work.is_neg = false;
  b_work.is_neg = false;

  if (rz_cmp(&a_work, &b_work) < 0) {
    TRY(rz_swap(&a_work, &b_work));
  }

  while (!rz_is_zero(&b_work)) {
    if (a_work.size - b_work.size > 1) {
      TRY(rz_mod(&t, &a_work, &b_work));
      TRY(rz_copy(&a_work, &b_work));
      TRY(rz_copy(&b_work, &t));
      continue;
    }

    // 2. [Test Quotient]
    i64 A, B, C, D;
    bool matrix_ok = lehmer_quotient_loop(&a_work, &b_work, &A, &B, &C, &D);

    // 3. [Multi-precision step]
    if (matrix_ok) {
      // Apply the accumulated matrix exactly:
      // a_new = A*a + B*b, b_new = C*a + D*b (rz_mul_i64 is signed)
      TRY(rz_mul_i64(&p1, &a_work, A));
      TRY(rz_mul_i64(&p2, &b_work, B));
      TRY(rz_add(&t, &p1, &p2));

      TRY(rz_mul_i64(&p1, &a_work, C));
      TRY(rz_mul_i64(&p2, &b_work, D));
      TRY(rz_add(&r, &p1, &p2));

      // The quotient test only sees truncated limbs, so the transformed
      // pair could still be invalid; require 0 <= b_new <= a_new with b
      // strictly decreasing, otherwise discard the matrix
      if (t.is_neg || r.is_neg || rz_cmp(&t, &r) < 0 ||
          rz_cmp(&r, &b_work) >= 0) {
        matrix_ok = false;
      }
    }

    if (!matrix_ok) {
      TRY(rz_copy(&t, &b_work));
      TRY(rz_mod(&r, &a_work, &b_work));
    }

    TRY(rz_copy(&a_work, &t));
    TRY(rz_copy(&b_work, &r));
  }

  TRY(rz_copy(d, &a_work));
  d->is_neg = false;
  err = RABIN_SUCCESS;
cleanup:
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
  if (err != RABIN_SUCCESS) goto cleanup;
  TRY(rz_copy(d, a));

  if (rz_is_zero(b)) {
    TRY(rz_set_u64(v, 0));
    // normalize: gcd is non-negative (GMP convention)
    if (d->is_neg) {
      d->is_neg = false;
      u->is_neg = !u->is_neg;
    }
    err = RABIN_SUCCESS;
    goto cleanup;
  }

  TRY(rz_set_u64(&v_1, 0));
  TRY(rz_copy(&v_3, b));

  while (!rz_is_zero(&v_3)) {
    TRY(rz_divmod(&temp, &t_3, d, &v_3));

    // t_1 = u - qv_1
    TRY(rz_mul(&temp, &temp, &v_1));
    TRY(rz_copy(&t_1, u));
    TRY(rz_sub(&t_1, &t_1, &temp));

    TRY(rz_copy(u, &v_1));
    TRY(rz_copy(d, &v_3));
    TRY(rz_copy(&v_1, &t_1));
    TRY(rz_copy(&v_3, &t_3));
  }

  // normalize: gcd is non-negative (GMP convention); flip d and u
  // together so the Bezout identity u*a + v*b == d is preserved
  if (d->is_neg) {
    d->is_neg = false;
    u->is_neg = !u->is_neg;
  }

  // v = (d - au) / b
  TRY(rz_mul(&temp, a, u));
  TRY(rz_copy(v, d));
  TRY(rz_sub(v, v, &temp));
  TRY(rz_div(v, v, b));

  err = RABIN_SUCCESS;
cleanup:
  rz_clear_multi(&v_1, &v_3, &t_1, &t_3, &temp, NULL);
  return err;
}

rabin_err_t rz_gcd_extended_lehmer(rz_t* u, rz_t* v, rz_t* d, const rz_t* a,
                                   const rz_t* b)
{
  if (u == NULL || v == NULL || d == NULL || a == NULL || b == NULL)
    return RABIN_ERR_NULL_PTR;

  rz_t a_work, b_work, t, r, v_1, Q, p1, p2;

  rz_init_multi(&a_work, &b_work, &t, &r, &v_1, &Q, &p1, &p2, NULL);
  rabin_err_t err = rz_copy(&a_work, a);
  if (err != RABIN_SUCCESS) goto cleanup;
  TRY(rz_copy(&b_work, b));

  // the Lehmer loop reasons about top limbs as positive magnitudes, so it
  // must run on |a|, |b|; the sign of a is re-applied to u at the end
  a_work.is_neg = false;
  b_work.is_neg = false;

  // zero-argument guards: the loop below reads the top limb of a_work,
  // which does not exist for a size-0 rz_t
  if (rz_is_zero(a)) {
    // u = 0, v = sign(b), d = |b|
    TRY(rz_set_u64(u, 0));
    if (b->is_neg) {
      TRY(rz_set_i64(v, -1));
    } else {
      TRY(rz_set_u64(v, 1));
    }
    TRY(rz_copy(d, b));
    d->is_neg = false;
    err = RABIN_SUCCESS;
    goto cleanup;
  }
  if (rz_is_zero(b)) {
    // u = sign(a), v = 0, d = |a|
    if (a->is_neg) {
      TRY(rz_set_i64(u, -1));
    } else {
      TRY(rz_set_u64(u, 1));
    }
    TRY(rz_set_u64(v, 0));
    TRY(rz_copy(d, a));
    d->is_neg = false;
    err = RABIN_SUCCESS;
    goto cleanup;
  }

  // 1. [Initialize]
  TRY(rz_set_u64(u, 1));
  TRY(rz_set_u64(&v_1, 0));

  // 2. [Finished?]
  while (!rz_is_zero(&b_work)) {
    // 3. [Test Quotient]
    i64 A, B, C, D;
    bool matrix_ok = lehmer_quotient_loop(&a_work, &b_work, &A, &B, &C, &D);

    // 4. [Multi-precision step]
    if (matrix_ok) {
      // Apply the accumulated matrix exactly:
      // a_new = A*a + B*b, b_new = C*a + D*b (rz_mul_i64 is signed)
      TRY(rz_mul_i64(&p1, &a_work, A));
      TRY(rz_mul_i64(&p2, &b_work, B));
      TRY(rz_add(&t, &p1, &p2));

      TRY(rz_mul_i64(&p1, &a_work, C));
      TRY(rz_mul_i64(&p2, &b_work, D));
      TRY(rz_add(&r, &p1, &p2));

      // The quotient test only sees truncated limbs, so the transformed
      // pair could still be invalid; require 0 <= b_new <= a_new with b
      // strictly decreasing, otherwise discard the matrix
      if (t.is_neg || r.is_neg || rz_cmp(&t, &r) < 0 ||
          rz_cmp(&r, &b_work) >= 0) {
        matrix_ok = false;
      }
    }

    if (!matrix_ok) {
      TRY(rz_divmod(&Q, &t, &a_work, &b_work));
      TRY(rz_copy(&a_work, &b_work));
      TRY(rz_copy(&b_work, &t));

      // t = u - Q * v_1
      TRY(rz_mul(&t, &Q, &v_1));
      TRY(rz_sub(&t, u, &t));
      TRY(rz_copy(u, &v_1));
      TRY(rz_copy(&v_1, &t));
    } else {
      TRY(rz_copy(&a_work, &t));
      TRY(rz_copy(&b_work, &r));

      // Cofactor update:
      // u_new = A*u + B*v_1,  v1_new = C*u + D*v_1
      TRY(rz_mul_i64(&p1, u, A));
      TRY(rz_mul_i64(&p2, &v_1, B));
      TRY(rz_add(&t, &p1, &p2));

      TRY(rz_mul_i64(&p1, u, C));
      TRY(rz_mul_i64(&p2, &v_1, D));
      TRY(rz_add(&r, &p1, &p2));

      TRY(rz_copy(u, &t));
      TRY(rz_copy(&v_1, &r));
    }
  }

  // Set final GCD: d = a_work (>= 0, the loop ran on magnitudes)
  TRY(rz_copy(d, &a_work));
  d->is_neg = false;

  // u is the Bezout cofactor of |a|; re-apply the sign of a so that
  // u*a + v*b == d holds for the original (signed) inputs
  if (a->is_neg) {
    u->is_neg = !u->is_neg;
  }

  // Calculate final cofactor v = (d - a * u) / b
  TRY(rz_mul(&t, a, u));
  TRY(rz_sub(&t, d, &t));
  TRY(rz_divmod(v, &r, &t, b));

  err = RABIN_SUCCESS;
cleanup:
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

  // cornacchia needs 0 < d < p; the loop below works with -d mod p
  if (d_work.is_neg || rz_is_zero(&d_work) || rz_cmp(&d_work, p) >= 0) {
    err = RABIN_ERR_INVALID_ARG;
    goto cleanup;
  }

  // d_work = -d; ( -d | p ) != -1 is required for a square root to exist
  d_work.is_neg = true;
  i64 k = rz_jacobi(&d_work, p);
  if (k == -1) {
    err = RABIN_ERR_INVALID_ARG;
    goto cleanup;
  }

  // 2. [Compute square root]
  // move -d into Z_p
  TRY(rz_add(&d_work, &d_work, p));
  TRY(rz_tonelli_shanks(&x_0, &d_work, p));

  TRY(rz_rshift(&p_half, p, 1));

  if (rz_cmp(&x_0, &p_half) <= 0) {
    TRY(rz_sub(&x_0, p, &x_0));
  }

  TRY(rz_copy(&a, p));
  TRY(rz_copy(&b, &x_0));
  TRY(rz_isqrt(&l, p));

  // 3. [Euclidean algorithm]

  while (rz_cmp(&b, &l) == 1) {
    TRY(rz_mod(&r, &a, &b));
    TRY(rz_copy(&a, &b));
    TRY(rz_copy(&b, &r));
  }

  // 4. [Test solution]
  TRY(rz_sqr(&b_sq, &b));
  TRY(rz_sub(&t, p, &b_sq));

  // check d | p - b^2
  TRY(rz_divmod(&c, &r, &t, d));
  if (!rz_is_zero(&r)) {
    err = RABIN_ERR_INVALID_ARG;
    goto cleanup;
  }

  // check if (p - b^2) / d is not a square
  TRY(rz_isqrt(y, &c));
  TRY(rz_sqr(&b_sq, y));
  if (rz_cmp(&b_sq, &c) != 0) {
    err = RABIN_ERR_INVALID_ARG;
    goto cleanup;
  }

  TRY(rz_copy(x, &b));

  err = RABIN_SUCCESS;
cleanup:
  rz_clear_multi(&d_work, &x_0, &p_half, &a, &b, &l, &b_sq, &t, &r, &c, NULL);

  return err;
}
