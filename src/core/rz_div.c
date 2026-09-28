/*
 * rz_div.c
 *
 * rz_t division.
 *
 * Three algorithms are implemented:
 *
 * 1. rz_divmod / rz_div / rz_mod
 *    Knuth's Algorithm D (TAoCP 4.3.1) with the Moeller-Granlund
 *    "invariant divisor" quotient-digit computation
 *    (M. Moeller, P. Granlund, "Improved Division by Invariant Integers",
 *    ACM Trans. on Computer Systems, 2011).
 *
 *      - single-limb divisor:  one rz_div2by1 step per digit   (O(n))
 *      - multi-limb divisor:   one rz_div3by2 step per digit   (O(n*m))
 *
 * 2. rz_div_exact
 *    Jebelean's exact division. Only valid when b divides a, but about
 *    2x faster than a real division because it needs neither
 *    normalisation nor quotient correction.
 *
 * 3. rz_newton_div
 *    Reciprocal-based division: Newton-iterate x -> 2^P / d, then
 *    q = (a * x) >> P. Useful when many divisions by the same divisor
 *    are expected (e.g. in primality tests).
 *
 * Division is TRUNCATED (C semantics): q truncates toward 0 and
 * sign(r) == sign(a).
 *
 * Assumptions about the rz_t API used here:
 *     struct rz_t { u64 *limbs; u64 size; bool is_neg; ... };
 *     rz_init / rz_clear / rz_copy / rz_set_u64 / rz_trim / rz_is_zero
 *     rz_alloc(x, n)  -> guarantees capacity >= n limbs, preserves contents
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
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include "../../include/rabin.h"
#include "../../include/rzlimb.h"
#include "rz_internal.h"

#define RZ_UNLIKELY(x) __builtin_expect(!!(x), 0)

/* limbs of scratch we are willing to put on the stack (4 KiB) */
#define RZ_DIV_STACK_LIMBS 512

/**
 * @brief Number of Newton iterations needed to refine the reciprocal of a
 * d_bits-bit divisor.
 *
 * Each iteration roughly doubles the number of correct bits, so about
 * 2/3 * log2(d) iterations are enough.
 *
 * Complexity:
 *   Time: O(1)
 *   Auxiliary memory: O(1)
 *   Output memory: O(1)
 *
 * @param[in] bits Bit length of the divisor d.
 *
 * @return The number of Newton iterations to perform.
 */
static u64 estimate_iterations(u64 bits)
{
  if (bits <= 64) return 2;

  u64 log2 = 64 - __builtin_clzll(bits);
  u64 its = (log2 * 2) / 3;

  return its + 1;
}

rabin_err_t rz_newton_div(rz_t* q, const rz_t* a, const rz_t* d)
{
  if (q == NULL || a == NULL || d == NULL) return RABIN_ERR_NULL_PTR;
  if (rz_is_zero(d)) return RABIN_ERR_DIV_BY_ZERO;
  if (rz_cmp(a, d) < 0) {
    return rz_set_u64(q, 0);
  }

  rz_t x, two_p, tmp, error, r;
  rz_init_multi(&x, &two_p, &tmp, &error, &r, NULL);

  rabin_err_t err;

  /* work with P bits of precision: bit length of a plus a safety margin */
  u64 P = (u64)rz_bit_length(a) + 32;

  /* seed the iteration with one real division: x = 2^P / d */
  if ((err = rz_set_u64(&x, 1)) != RABIN_SUCCESS) goto out;
  if ((err = rz_lshift(&x, &x, P)) != RABIN_SUCCESS) goto out;
  if ((err = rz_div(&x, &x, d)) != RABIN_SUCCESS) goto out;

  /* two_p = 2^{P} */
  if ((err = rz_set_u64(&two_p, 1)) != RABIN_SUCCESS) goto out;
  if ((err = rz_lshift(&two_p, &two_p, P)) != RABIN_SUCCESS) goto out;

  /*
   * Newton iteration for f(x) = 1/x - d, i.e.
   *
   *     x <- x + (x * (2^P - d*x)) >> P
   *
   * Each step roughly doubles the number of correct bits of x.
   */
  for (u64 i = 0; i < estimate_iterations((u64)rz_bit_length(d)); i++) {
    if ((err = rz_mul(&tmp, d, &x)) != RABIN_SUCCESS)
      goto out; /* tmp = d * x */
    if ((err = rz_sub(&error, &two_p, &tmp)) != RABIN_SUCCESS) goto out;
    if (rz_is_zero(&error)) break; /* reciprocal is exact */

    if ((err = rz_mul(&tmp, &x, &error)) != RABIN_SUCCESS) goto out;
    if ((err = rz_rshift(&tmp, &tmp, P)) != RABIN_SUCCESS) goto out;
    if ((err = rz_add(&x, &x, &tmp)) != RABIN_SUCCESS) goto out;
  }

  /* q = (a * x) >> P, then fix up the (rare) off-by-one */
  if ((err = rz_mul(q, a, &x)) != RABIN_SUCCESS) goto out;
  if ((err = rz_rshift(q, q, P)) != RABIN_SUCCESS) goto out;

  if ((err = rz_mul(&tmp, q, d)) != RABIN_SUCCESS) goto out;
  if ((err = rz_sub(&r, a, &tmp)) != RABIN_SUCCESS) goto out;

  while (rz_cmp(&r, d) >= 0) {
    if ((err = rz_add_u64(q, q, 1)) != RABIN_SUCCESS) goto out;
    if ((err = rz_sub(&r, &r, d)) != RABIN_SUCCESS) goto out;
  }

  err = RABIN_SUCCESS;
out:
  rz_clear_multi(&x, &tmp, &two_p, &r, &error, NULL);
  return err;
}

/**
 * @brief v = floor((2^128 - 1) / d) - 2^64 ; requires d >= 2^63.
 *
 * Reciprocal used by the 2/1 division: for n1 < d it yields
 * q = floor((n1 n0) / d) with an error of at most 1.
 *
 * Complexity:
 *   Time: O(1)
 *   Auxiliary memory: O(1)
 *   Output memory: O(1)
 *
 * @param[in] d Normalised divisor limb (d >= 2^63).
 *
 * @return The reciprocal v = floor((2^128 - 1) / d) - 2^64.
 */
static inline u64 rz_invert_limb(u64 d)
{
  return (u64)(((((u128)~d) << 64) | ~(u64)0) / d);
}

/**
 * @brief Reciprocal for the 3/2 division; requires d1 >= 2^63.
 *
 * With v = rz_invert_3by2(d1, d0), q = floor((n2 n1 n0) / (d1 d0))
 * is obtained from a couple of multiplications with an error of at
 * most 1.
 *
 * Complexity:
 *   Time: O(1)
 *   Auxiliary memory: O(1)
 *   Output memory: O(1)
 *
 * @param[in] d1 High limb of the normalised divisor (d1 >= 2^63).
 * @param[in] d0 Low limb of the normalised divisor.
 *
 * @return The 3/2 reciprocal v.
 */
static u64 rz_invert_3by2(u64 d1, u64 d0)
{
  u64 v = rz_invert_limb(d1);
  u64 p = d1 * v + d0; /* mod 2^64 */
  if (p < d0) {
    v--;
    u64 mask = -(u64)(p >= d1);
    p -= d1;
    v += mask;
    p -= mask & d1;
  }
  u128 t = (u128)v * d0;
  u64 t1 = (u64)(t >> 64), t0 = (u64)t;
  p += t1;
  if (p < t1) {
    v--;
    if (RZ_UNLIKELY(p >= d1))
      if (p > d1 || t0 >= d0) v--;
  }
  return v;
}

/**
 * @brief {n1,n0} / d, requires d >= 2^63 and n1 < d.
 *
 * Returns the quotient digit and stores the remainder in *rp.
 *
 * The reciprocal dinv yields a quotient estimate that is off by at
 * most 1; the two corrections below bring it to the exact value.
 *
 * Complexity:
 *   Time: O(1)
 *   Auxiliary memory: O(1)
 *   Output memory: O(1)
 *
 * @param[out] rp   Receives the remainder.
 * @param[in]  n1   High limb of the dividend window.
 * @param[in]  n0   Low limb of the dividend window.
 * @param[in]  d    Normalised divisor (d >= 2^63).
 * @param[in]  dinv Reciprocal from rz_invert_limb(d).
 *
 * @return The exact quotient digit.
 */
static inline u64 rz_div2by1(u64* rp, u64 n1, u64 n0, u64 d, u64 dinv)
{
  u128 t = (u128)n1 * dinv;
  u64 q1 = (u64)(t >> 64), q0 = (u64)t;

  u128 lo = (u128)q0 + n0;
  q0 = (u64)lo;
  q1 += (n1 + 1) + (u64)(lo >> 64);

  u64 r = n0 - q1 * d;
  u64 mask = -(u64)(r > q0);
  q1 += mask;
  r += mask & d;
  if (RZ_UNLIKELY(r >= d)) {
    q1++;
    r -= d;
  }
  *rp = r;
  return q1;
}

/**
 * @brief {n2,n1,n0} / {d1,d0}.  Requires d1 >= 2^63 and {n2,n1} < {d1,d0}.
 *
 * Returns the quotient digit and stores the 128-bit remainder
 * r1:r0 = {n2,n1,n0} - q*{d1,d0} in *r1p:*r0p.
 *
 * Complexity:
 *   Time: O(1)
 *   Auxiliary memory: O(1)
 *   Output memory: O(1)
 *
 * @param[out] r1p  Receives the high limb of the 128-bit remainder.
 * @param[out] r0p  Receives the low limb of the 128-bit remainder.
 * @param[in]  n2   Top limb of the dividend window.
 * @param[in]  n1   Middle limb of the dividend window.
 * @param[in]  n0   Low limb of the dividend window.
 * @param[in]  d1   High limb of the normalised divisor (d1 >= 2^63).
 * @param[in]  d0   Low limb of the normalised divisor.
 * @param[in]  dinv Reciprocal from rz_invert_3by2(d1, d0).
 *
 * @return The exact quotient digit.
 */
static inline u64 rz_div3by2(u64* r1p, u64* r0p, u64 n2, u64 n1, u64 n0, u64 d1,
                             u64 d0, u64 dinv)
{
  u128 t = (u128)n2 * dinv;
  u64 q1 = (u64)(t >> 64), q0 = (u64)t;

  u128 lo = (u128)q0 + n1; /* {q1,q0} += {n2,n1} */
  q0 = (u64)lo;
  q1 += n2 + (u64)(lo >> 64);

  const u128 D = ((u128)d1 << 64) | d0;

  u64 r1 = n1 - d1 * q1;
  u128 r = ((((u128)r1) << 64) | n0) - D - (u128)d0 * q1;
  q1++;

  u64 mask = -(u64)((u64)(r >> 64) >= q0); /* branch-free, ~50/50 case */
  q1 += mask;
  r += ((u128)(mask & d1) << 64) | (u64)(mask & d0);

  if (RZ_UNLIKELY(r >= D)) {
    q1++;
    r -= D;
  } /* almost never taken */

  *r1p = (u64)(r >> 64);
  *r0p = (u64)r;
  return q1;
}

/**
 * @brief qp[0..un-1] = up/d, returns up % d.  qp may be NULL.
 *
 * Let n = un, measured in 64-bit limbs.
 *
 * The divisor is first normalised (dn = d << clz(d) >= 2^63) so that
 * rz_invert_limb / rz_div2by1 can be used. The digits are then
 * computed from the top down, one rz_div2by1 step per digit.
 *
 * Complexity:
 *   Time: O(n)
 *   Auxiliary memory: O(1)
 *   Output memory: O(n) limbs for qp
 *
 * @param[out] qp  Quotient limbs (may be NULL).
 * @param[in]  up  Dividend limbs.
 * @param[in]  un  Number of dividend limbs.
 * @param[in]  d   Single-limb divisor.
 *
 * @return The remainder up % d.
 */
static u64 limbs_divrem_1(u64* qp, const u64* up, u64 un, u64 d)
{
  unsigned s = (unsigned)__builtin_clzll(d);
  u64 dn = d << s;
  u64 dinv = rz_invert_limb(dn);
  u64 r;

  if (s == 0) {
    /*
     * Divisor already normalised. The top quotient digit is 0 or 1;
     * peel it off first so the 2/1 precondition r < d holds.
     */
    if (up[un - 1] >= dn) {
      r = up[un - 1] - dn;
      if (qp) qp[un - 1] = 1;
    } else {
      r = up[un - 1];
      if (qp) qp[un - 1] = 0;
    }
    for (i64 i = (i64)un - 2; i >= 0; i--) {
      u64 q = rz_div2by1(&r, r, up[i], dn, dinv);
      if (qp) qp[i] = q;
    }
    return r;
  }

  /*
   * Normalised divisor dn = d << s. Work from the top down, feeding
   * each digit the s bits that spilled over from the limb above.
   */
  u64 prev = up[un - 1];
  r = prev >> (64 - s);
  for (i64 i = (i64)un - 1; i > 0; i--) {
    u64 cur = up[i - 1];
    u64 q = rz_div2by1(&r, r, (prev << s) | (cur >> (64 - s)), dn, dinv);
    if (qp) qp[i] = q;
    prev = cur;
  }
  {
    u64 q = rz_div2by1(&r, r, prev << s, dn, dinv);
    if (qp) qp[0] = q;
  }
  return r >> s; /* undo the normalisation of the remainder */
}

/**
 * @brief qp[0 .. un-vn]  = up/vp          (un >= vn >= 2, vp[vn-1] != 0)
 *
 *      rp[0 .. vn-1]   = up%vp          (rp may be NULL)
 *      scratch: un + 1 + (normalisation ? vn : 0) limbs
 *
 * How it works:
 *   1. Normalise: shift up and vp left by s = clz(vp[vn-1]) so that the
 *      top limb of the divisor is >= 2^63. This bounds the quotient
 *      digit error of the 3/2 estimate to 1.
 *   2. For each digit position j (from the top down):
 *        - estimate the digit qhat with rz_div3by2 from the top three
 *          limbs of the current window and the top two of the divisor
 *        - subtract qhat * vp from the window
 *        - if the subtraction underflowed, qhat was one too large:
 *          decrement it and add vp back
 *   3. Shift the remainder right by s to undo the normalisation.
 *
 * Complexity:
 *   Time: O((un - vn + 1) * vn) - one O(vn) submul per quotient digit
 *   Auxiliary memory: O(un + vn) limbs of scratch
 *   Output memory: O(un - vn + 1) limbs for qp, O(vn) for rp
 *
 * @param[out] qp      Quotient limbs (un - vn + 1 of them).
 * @param[out] rp      Remainder limbs (vn of them; may be NULL).
 * @param[in]  up      Dividend limbs.
 * @param[in]  un      Number of dividend limbs.
 * @param[in]  vp      Divisor limbs (vp[vn-1] != 0).
 * @param[in]  vn      Number of divisor limbs (>= 2).
 * @param[in]  scratch Scratch area (un + 1 + (s ? vn : 0) limbs).
 */
static void rz_divmod_limbs(u64* qp, u64* rp, const u64* up, u64 un,
                            const u64* vp, u64 vn, u64* scratch)
{
  unsigned s = (unsigned)__builtin_clzll(vp[vn - 1]);
  u64* nu = scratch; /* normalised dividend, un+1 limbs */
  const u64* nv;

  if (s) {
    u64* tv = scratch + un + 1;
    limbs_lshift(tv, vp, vn, s); /* high bits are zero by definition */
    nv = tv;
    nu[un] = limbs_lshift(nu, up, un, s);
  } else {
    nv = vp; /* no copy needed */
    memcpy(nu, up, un * sizeof(u64));
    nu[un] = 0;
  }

  const u64 d1 = nv[vn - 1];
  const u64 d0 = nv[vn - 2];
  const u64 dinv = rz_invert_3by2(d1, d0);
  const u128 D = ((u128)d1 << 64) | d0;

  for (i64 j = (i64)(un - vn); j >= 0; j--) {
    u64 n2 = nu[j + vn], n1 = nu[j + vn - 1], n0 = nu[j + vn - 2];
    u64 qhat;

    if (RZ_UNLIKELY(n2 == d1 && n1 == d0)) {
      /* the only case where the 3/2 quotient would not fit in a limb;
       * one can prove qhat == 2^64-1 exactly here. */
      qhat = ~(u64)0;
      u64 cy = limbs_submul_1(nu + j, nv, vn, qhat);
      (void)cy; /* == n2 */
      assert(cy == n2);
      nu[j + vn] = 0;
    } else {
      u64 r1, r0;
      qhat = rz_div3by2(&r1, &r0, n2, n1, n0, d1, d0, dinv);

      /*
       * Subtract qhat * v[0 .. vn-3] from the window. The top two limbs
       * of the product are already accounted for in r1:r0, so the new
       * top two limbs of the window are
       *
       *     rem = r1:r0 - cy        (cy = borrow out of the low limbs)
       *
       * If rem underflows, qhat was one too large: add v back.
       */
      u64 cy = limbs_submul_1(nu + j, nv, vn - 2, qhat);
      u128 rem = ((u128)r1 << 64) | r0;

      if (RZ_UNLIKELY(rem < (u128)cy)) { /* qhat was 1 too large: add back */
        rem -= cy;                       /* wraps, on purpose */
        qhat--;
        rem += limbs_add_n(nu + j, nv, vn - 2);
        rem += D; /* carry out of 128 bits cancels */
      } else {
        rem -= cy;
      }
      nu[j + vn - 2] = (u64)rem;
      nu[j + vn - 1] = (u64)(rem >> 64);
      nu[j + vn] = 0;
    }
    qp[j] = qhat;
  }

  if (rp) limbs_rshift(rp, nu, vn, s); /* undo normalisation */
}

/**
 * @brief Number of significant limbs of a rz_t (skips trailing zeros).
 *
 * Complexity:
 *   Time: O(number of trailing zero limbs)
 *   Auxiliary memory: O(1)
 *   Output memory: O(1)
 *
 * @param[in] x Bignum to measure.
 *
 * @return The number of significant (non-trailing-zero) limbs.
 */
static inline u64 rz_nsize(const rz_t* x)
{
  return limbs_norm(x->limbs, x->size);
}

rabin_err_t rz_divmod(rz_t* q, rz_t* r, const rz_t* a, const rz_t* b)
{
  if (a == NULL || b == NULL) return RABIN_ERR_NULL_PTR;
  assert(q == NULL || q != r);

  u64 an = rz_nsize(a);
  u64 vn = rz_nsize(b);
  if (RZ_UNLIKELY(vn == 0)) {
    return RABIN_ERR_DIV_BY_ZERO;
  }

  const bool qneg = a->is_neg ^ b->is_neg;
  const bool rneg = a->is_neg;

  /* |a| < |b|  ->  q = 0, r = a   (compare magnitudes, never rz_cmp!) */
  if (an < vn || (an == vn && limbs_cmp(a->limbs, b->limbs, an) < 0)) {
    rabin_err_t err = RABIN_SUCCESS;
    if (r) {
      if ((err = rz_copy(r, a)) != RABIN_SUCCESS) return err;
      rz_trim(r);
      r->is_neg = rneg && !rz_is_zero(r);
    }
    if (q && (err = rz_set_u64(q, 0)) != RABIN_SUCCESS) return err;
    return RABIN_SUCCESS;
  }

  const u64 qn = an - vn + 1;

  /* aliasing: work into temporaries only when we really have to */
  rz_t tq, tr;
  rz_t *Q = q, *R = r;
  bool alq = q && (q == a || q == b);
  bool alr = r && (r == a || r == b);
  if (alq) {
    rz_init(&tq);
    Q = &tq;
  }
  if (alr) {
    rz_init(&tr);
    R = &tr;
  }

  rabin_err_t err = RABIN_SUCCESS;
  if (Q && (err = rz_alloc(Q, qn)) != RABIN_SUCCESS) goto out;
  if (R && (err = rz_alloc(R, vn)) != RABIN_SUCCESS) goto out;

  if (vn == 1) {
    u64 rem = limbs_divrem_1(Q ? Q->limbs : NULL, a->limbs, an, b->limbs[0]);
    if (R) {
      R->limbs[0] = rem;
      R->size = 1;
    }
  } else {
    unsigned s = (unsigned)__builtin_clzll(b->limbs[vn - 1]);
    u64 need = an + 1 + (s ? vn : 0) + (Q ? 0 : qn);
    u64 stackbuf[RZ_DIV_STACK_LIMBS];
    u64* scratch;
    bool from_arena = false;
    if (need <= RZ_DIV_STACK_LIMBS) {
      scratch = stackbuf;
    } else {
      scratch = rz_scratch_get(need);
      if (scratch == NULL) {
        err = RABIN_ERR_OUT_OF_MEMORY;
        goto out;
      }
      from_arena = true;
    }
    u64* qp = Q ? Q->limbs : scratch + an + 1 + (s ? vn : 0);

    rz_divmod_limbs(qp, R ? R->limbs : NULL, a->limbs, an, b->limbs, vn,
                    scratch);

    if (from_arena) rz_scratch_release();
    if (R) R->size = vn;
  }

  if (Q) {
    Q->size = qn;
    rz_trim(Q);
    Q->is_neg = qneg && !rz_is_zero(Q);
  }
  if (R) {
    rz_trim(R);
    R->is_neg = rneg && !rz_is_zero(R);
  }

  if (alq) {
    if ((err = rz_copy(q, &tq)) != RABIN_SUCCESS) goto out;
  }
  if (alr) {
    if ((err = rz_copy(r, &tr)) != RABIN_SUCCESS) goto out;
  }

out:
  if (alq) rz_clear(&tq);
  if (alr) rz_clear(&tr);
  return err;
}

rabin_err_t rz_div_euclid(rz_t* q, const rz_t* a, const rz_t* b)
{
  if (q == NULL || a == NULL || b == NULL) return RABIN_ERR_NULL_PTR;

  rz_t zero, one, r;

  rz_init(&zero);
  rz_init(&one);
  rz_init(&r);

  rabin_err_t err = rz_set_u64(&zero, 0);
  if (err == RABIN_SUCCESS) err = rz_set_u64(&one, 1);

  if (rz_is_zero(b)) {
    /* division by zero */
    err = RABIN_ERR_DIV_BY_ZERO;
    goto cleanup;
  }
  if (err != RABIN_SUCCESS) goto cleanup;

  /*
   * rz_divmod gives truncated division:
   *
   *     a = b*q + r
   *
   * with |r| < |b|.
   */
  if ((err = rz_divmod(q, &r, a, b)) != RABIN_SUCCESS) goto cleanup;

  /*
   * Convert to Euclidean division:
   *
   *     0 <= r < |b|
   */
  if (rz_cmp(&r, &zero) < 0) {
    if (rz_cmp(b, &zero) > 0) {
      /* q <- q - 1, r <- r + b */
      err = rz_sub(q, q, &one);
    } else {
      /* q <- q + 1, r <- r - b */
      err = rz_add(q, q, &one);
    }
    if (err != RABIN_SUCCESS) goto cleanup;
  }

cleanup:
  rz_clear(&zero);
  rz_clear(&one);
  rz_clear(&r);
  return err;
}

rabin_err_t rz_div(rz_t* q, const rz_t* a, const rz_t* b)
{
  return rz_divmod(q, NULL, a, b);
}

rabin_err_t rz_mod(rz_t* r, const rz_t* a, const rz_t* b)
{
  return rz_divmod(NULL, r, a, b);
}

rabin_err_t rz_mod_pos(rz_t* r, const rz_t* a, const rz_t* b)
{
  if (r == NULL) return RABIN_ERR_NULL_PTR;

  rabin_err_t err = rz_mod(r, a, b);
  if (err != RABIN_SUCCESS) return err;

  if (r->is_neg) {
    if (b->is_neg) {
      err = rz_sub(r, r, b);
    } else {
      err = rz_add(r, r, b);
    }
  }
  return err;
}

rabin_err_t rz_div_exact(rz_t* q, const rz_t* a, const rz_t* b)
{
  if (q == NULL || a == NULL || b == NULL) return RABIN_ERR_NULL_PTR;

  u64 an = rz_nsize(a), rz_size = rz_nsize(b);
  if (RZ_UNLIKELY(rz_size == 0)) {
    return RABIN_ERR_DIV_BY_ZERO;
  }
  if (an == 0) {
    return rz_set_u64(q, 0);
  }

  /* strip the power of two from b (and the same amount from a) */
  u64 kl = 0;
  while (b->limbs[kl] == 0) kl++;
  unsigned kb = (unsigned)__builtin_ctzll(b->limbs[kl]);
  assert(an > kl);

  u64 dl = rz_size - kl, al = an - kl;
  u64* buf = rz_scratch_get(dl + al);
  if (buf == NULL) return RABIN_ERR_OUT_OF_MEMORY;
  u64 *D = buf, *A = buf + dl;
  limbs_rshift(D, b->limbs + kl, dl, kb);
  limbs_rshift(A, a->limbs + kl, al, kb);
  dl = limbs_norm(D, dl);
  al = limbs_norm(A, al);

  if (al < dl) {
    rabin_err_t err = rz_set_u64(q, 0);
    rz_scratch_release();
    return err;
  }

  u64 qn = al - dl + 1;
  rz_t tq;
  rz_t* Q = q;
  bool al_ = (q == a || q == b);
  if (al_) {
    rz_init(&tq);
    Q = &tq;
  }
  rabin_err_t err = rz_alloc(Q, qn);
  if (err != RABIN_SUCCESS) {
    if (al_) rz_clear(&tq);
    rz_scratch_release();
    return err;
  }
  Q->size = qn;

  /* D[0] is odd; mod_inverse_u64 returns -D[0]^{-1}, so negate it to get
   * the true inverse with D[0]*dinv == 1 (mod 2^64) */
  const u64 dinv = -mod_inverse_u64(D[0]);

  for (u64 i = 0; i < qn; i++) {
    u64 qi = A[i] * dinv; /* kills limb i of A */
    Q->limbs[i] = qi;
    u64 cy = limbs_submul_1(A + i, D, dl, qi);
    if (i + dl < al) limbs_sub_1(A + i + dl, al - i - dl, cy);
  }

  rz_trim(Q);
  Q->is_neg = (a->is_neg ^ b->is_neg) && !rz_is_zero(Q);
  if (al_) {
    if ((err = rz_copy(q, &tq)) != RABIN_SUCCESS) {
      rz_clear(&tq);
      rz_scratch_release();
      return err;
    }
    rz_clear(&tq);
  }
  rz_scratch_release();
  return RABIN_SUCCESS;
}
