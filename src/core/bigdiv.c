/*
 * bigdiv.c
 *
 * bignum division.
 *
 * Three algorithms are implemented:
 *
 * 1. bn_divmod / bn_div / bn_mod
 *    Knuth's Algorithm D (TAoCP 4.3.1) with the Moeller-Granlund
 *    "invariant divisor" quotient-digit computation
 *    (M. Moeller, P. Granlund, "Improved Division by Invariant Integers",
 *    ACM Trans. on Computer Systems, 2011).
 *
 *      - single-limb divisor:  one bn_div2by1 step per digit   (O(n))
 *      - multi-limb divisor:   one bn_div3by2 step per digit   (O(n*m))
 *
 * 2. bn_div_exact
 *    Jebelean's exact division. Only valid when b divides a, but about
 *    2x faster than a real division because it needs neither
 *    normalisation nor quotient correction.
 *
 * 3. bn_newton_div
 *    Reciprocal-based division: Newton-iterate x -> 2^P / d, then
 *    q = (a * x) >> P. Useful when many divisions by the same divisor
 *    are expected (e.g. in primality tests).
 *
 * Division is TRUNCATED (C semantics): q truncates toward 0 and
 * sign(r) == sign(a).
 *
 * Assumptions about the bignum API used here:
 *     struct bignum { u64 *limbs; u64 size; bool is_neg; ... };
 *     bn_init / bn_free / bn_copy / bn_set_u64 / bn_trim / bn_is_zero
 *     bn_alloc(x, n)  -> guarantees capacity >= n limbs, preserves contents
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

#include "../../include/bighelper.h"
#include "../../include/bignum.h"

#define BN_UNLIKELY(x) __builtin_expect(!!(x), 0)

/* limbs of scratch we are willing to put on the stack (4 KiB) */
#define BN_DIV_STACK_LIMBS 512

/*===========================================================================
 *  1. Newton division
 *=========================================================================*/

/*
 * Number of Newton iterations needed to refine the reciprocal of a
 * d_bits-bit divisor. Each iteration roughly doubles the number of
 * correct bits, so about 2/3 * log2(d) iterations are enough.
 *
 * Complexity:
 *   Time: O(1)
 *   Auxiliary memory: O(1)
 *   Output memory: O(1)
 */
static u64 estimate_iterations(u64 bits)
{
  if (bits <= 64) return 2;

  u64 log2 = 64 - __builtin_clzll(bits);
  u64 its = (log2 * 2) / 3;

  return its + 1;
}

/*
 * Divide a by d using a Newton-iterated reciprocal: q = a / d.
 *
 * Let n = a->size and m = d->size, measured in 64-bit limbs.
 *
 * Works with P = bit_length(a) + 32 bits of precision. The reciprocal
 * x = 2^P / d is seeded with one real division, then refined with the
 * Newton iteration
 *
 *   x <- x + (x * (2^P - d*x)) >> P
 *
 * which roughly doubles the correct bits each step. The quotient is
 * then q = (a * x) >> P, followed by a rare off-by-one fix-up using
 * the remainder r = a - q*d.
 *
 * If d is zero, q is left unchanged. If |a| < |d|, q is set to 0.
 *
 * Complexity:
 *   Time: O(n^2) - O(log m) Newton iterations of n-limb multiplications
 *         plus one seeding division
 *   Auxiliary memory: O(n) limbs for temporaries
 *   Output memory: O(n) limbs
 */
void bn_newton_div(bignum* q, const bignum* a, const bignum* d)
{
  if (bn_is_zero(d)) return;
  if (bn_cmp(a, d) < 0) {
    bn_set_u64(q, 0);
    return;
  }

  bignum x, two_p, tmp, error;
  bn_init_multi(&x, &two_p, &tmp, &error, NULL);

  /* work with P bits of precision: bit length of a plus a safety margin */
  u64 P = (u64)bn_bit_length(a) + 32;

  /* seed the iteration with one real division: x = 2^P / d */
  bn_set_u64(&x, 1);
  bn_lshift(&x, &x, P);
  bn_div(&x, &x, d);

  /* two_p = 2^{P} */
  bn_set_u64(&two_p, 1);
  bn_lshift(&two_p, &two_p, P);

  /*
   * Newton iteration for f(x) = 1/x - d, i.e.
   *
   *     x <- x + (x * (2^P - d*x)) >> P
   *
   * Each step roughly doubles the number of correct bits of x.
   */
  for (u64 i = 0; i < estimate_iterations((u64)bn_bit_length(d)); i++) {
    bn_mul(&tmp, d, &x);           /* tmp   = d * x              */
    bn_sub(&error, &two_p, &tmp);  /* error = 2^P - d * x        */
    if (bn_is_zero(&error)) break; /* reciprocal is exact */

    bn_mul(&tmp, &x, &error); /* tmp = x * error */
    bn_rshift(&tmp, &tmp, P);
    bn_add(&x, &x, &tmp);
  }

  /* q = (a * x) >> P, then fix up the (rare) off-by-one */
  bn_mul(q, a, &x);
  bn_rshift(q, q, P);

  bignum r;
  bn_init(&r);
  bn_mul(&tmp, q, d);
  bn_sub(&r, a, &tmp);

  while (bn_cmp(&r, d) >= 0) {
    bn_add_u64(q, q, 1);
    bn_sub(&r, &r, d);
  }

  bn_free_multi(&x, &tmp, &two_p, &r, &error, NULL);
}

/*===========================================================================
 *  2. Invariant divisor reciprocals  (Moeller & Granlund, TC 2011)
 *
 *  Instead of dividing to estimate each quotient digit, a fixed
 *  "reciprocal" of the top limbs of the divisor is precomputed once,
 *  turning the per-digit estimate into a few multiplications.
 *=========================================================================*/

/*
 * v = floor((2^128 - 1) / d) - 2^64 ; requires d >= 2^63.
 *
 * Reciprocal used by the 2/1 division: for n1 < d it yields
 * q = floor((n1 n0) / d) with an error of at most 1.
 *
 * Complexity:
 *   Time: O(1)
 *   Auxiliary memory: O(1)
 *   Output memory: O(1)
 */
static inline u64 bn_invert_limb(u64 d)
{
#if defined(__x86_64__) && defined(__GNUC__)
  /* (2^64*(2^64-1-d) + (2^64-1)) / d  fits in one limb because d >= 2^63 */
  u64 q, r;
  __asm__("divq %4" : "=a"(q), "=d"(r) : "0"(~(u64)0), "1"(~d), "r"(d));
  return q;
#else
  return (u64)(((((u128)~d) << 64) | ~(u64)0) / d);
#endif
}

/*
 * Reciprocal for the 3/2 division; requires d1 >= 2^63.
 * With v = bn_invert_3by2(d1, d0), q = floor((n2 n1 n0) / (d1 d0))
 * is obtained from a couple of multiplications with an error of at
 * most 1.
 *
 * Complexity:
 *   Time: O(1)
 *   Auxiliary memory: O(1)
 *   Output memory: O(1)
 */
static u64 bn_invert_3by2(u64 d1, u64 d0)
{
  u64 v = bn_invert_limb(d1);
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
    if (BN_UNLIKELY(p >= d1))
      if (p > d1 || t0 >= d0) v--;
  }
  return v;
}

/*===========================================================================
 *  3. Quotient digit computation
 *=========================================================================*/

/*
 * {n1,n0} / d, requires d >= 2^63 and n1 < d.
 * Returns the quotient digit and stores the remainder in *rp.
 *
 * The reciprocal dinv yields a quotient estimate that is off by at
 * most 1; the two corrections below bring it to the exact value.
 *
 * Complexity:
 *   Time: O(1)
 *   Auxiliary memory: O(1)
 *   Output memory: O(1)
 */
static inline u64 bn_div2by1(u64* rp, u64 n1, u64 n0, u64 d, u64 dinv)
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
  if (BN_UNLIKELY(r >= d)) {
    q1++;
    r -= d;
  }
  *rp = r;
  return q1;
}

/*
 * {n2,n1,n0} / {d1,d0}.  Requires d1 >= 2^63 and {n2,n1} < {d1,d0}.
 * Returns the quotient digit and stores the 128-bit remainder
 * r1:r0 = {n2,n1,n0} - q*{d1,d0} in *r1p:*r0p.
 *
 * Complexity:
 *   Time: O(1)
 *   Auxiliary memory: O(1)
 *   Output memory: O(1)
 */
static inline u64 bn_div3by2(u64* r1p, u64* r0p, u64 n2, u64 n1, u64 n0, u64 d1,
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

  if (BN_UNLIKELY(r >= D)) {
    q1++;
    r -= D;
  } /* almost never taken */

  *r1p = (u64)(r >> 64);
  *r0p = (u64)r;
  return q1;
}

/*===========================================================================
 *  4. Single-limb divisor
 *=========================================================================*/

/*
 * qp[0..un-1] = up/d, returns up % d.  qp may be NULL.
 *
 * Let n = un, measured in 64-bit limbs.
 *
 * The divisor is first normalised (dn = d << clz(d) >= 2^63) so that
 * bn_invert_limb / bn_div2by1 can be used. The digits are then
 * computed from the top down, one bn_div2by1 step per digit.
 *
 * Complexity:
 *   Time: O(n)
 *   Auxiliary memory: O(1)
 *   Output memory: O(n) limbs for qp
 */
static u64 limbs_divrem_1(u64* qp, const u64* up, u64 un, u64 d)
{
  unsigned s = (unsigned)__builtin_clzll(d);
  u64 dn = d << s;
  u64 dinv = bn_invert_limb(dn);
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
      u64 q = bn_div2by1(&r, r, up[i], dn, dinv);
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
    u64 q = bn_div2by1(&r, r, (prev << s) | (cur >> (64 - s)), dn, dinv);
    if (qp) qp[i] = q;
    prev = cur;
  }
  {
    u64 q = bn_div2by1(&r, r, prev << s, dn, dinv);
    if (qp) qp[0] = q;
  }
  return r >> s; /* undo the normalisation of the remainder */
}

/*===========================================================================
 *  5. Algorithm D core (multi-limb divisor)
 *=========================================================================*/

/*
 *      qp[0 .. un-vn]  = up/vp          (un >= vn >= 2, vp[vn-1] != 0)
 *      rp[0 .. vn-1]   = up%vp          (rp may be NULL)
 *      scratch: un + 1 + (normalisation ? vn : 0) limbs
 *
 * How it works:
 *   1. Normalise: shift up and vp left by s = clz(vp[vn-1]) so that the
 *      top limb of the divisor is >= 2^63. This bounds the quotient
 *      digit error of the 3/2 estimate to 1.
 *   2. For each digit position j (from the top down):
 *        - estimate the digit qhat with bn_div3by2 from the top three
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
 */
static void bn_divmod_limbs(u64* qp, u64* rp, const u64* up, u64 un,
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
  const u64 dinv = bn_invert_3by2(d1, d0);
  const u128 D = ((u128)d1 << 64) | d0;

  for (i64 j = (i64)(un - vn); j >= 0; j--) {
    u64 n2 = nu[j + vn], n1 = nu[j + vn - 1], n0 = nu[j + vn - 2];
    u64 qhat;

    if (BN_UNLIKELY(n2 == d1 && n1 == d0)) {
      /* the only case where the 3/2 quotient would not fit in a limb;
       * one can prove qhat == 2^64-1 exactly here. */
      qhat = ~(u64)0;
      u64 cy = limbs_submul_1(nu + j, nv, vn, qhat);
      (void)cy; /* == n2 */
      assert(cy == n2);
      nu[j + vn] = 0;
    } else {
      u64 r1, r0;
      qhat = bn_div3by2(&r1, &r0, n2, n1, n0, d1, d0, dinv);

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

      if (BN_UNLIKELY(rem < (u128)cy)) { /* qhat was 1 too large: add back */
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

/*===========================================================================
 *  6. Public API
 *=========================================================================*/

/*
 * Number of significant limbs of a bignum (skips trailing zeros).
 *
 * Complexity:
 *   Time: O(number of trailing zero limbs)
 *   Auxiliary memory: O(1)
 *   Output memory: O(1)
 */
static inline u64 bn_nsize(const bignum* x)
{
  return limbs_norm(x->limbs, x->size);
}

/*
 * q = a / b, r = a % b  (either result may be NULL).
 *
 * Let n = a->size and m = b->size, measured in 64-bit limbs.
 *
 * Truncated division (C semantics): q truncates toward 0 and
 * sign(r) == sign(a).
 *
 * If |a| < |b| the result is q = 0, r = a. Otherwise the magnitude
 * division dispatches to the single-limb path (limbs_divrem_1) or to
 * Knuth's Algorithm D (bn_divmod_limbs). q and r may each alias a or
 * b (handled via temporaries); q and r must not alias each other.
 *
 * Scratch for the multi-limb path is taken from the stack when it fits
 * in BN_DIV_STACK_LIMBS limbs, otherwise from the heap.
 *
 * Complexity:
 *   Time: O(n) for a single-limb divisor, O((n - m + 1) * m) for a
 *         multi-limb divisor
 *   Auxiliary memory: O(n + m) limbs of scratch (stack or heap), plus
 *                     O(n) if q or r aliases an operand
 *   Output memory: O(n - m + 1) limbs for q, O(m) for r
 */
void bn_divmod(bignum* q, bignum* r, const bignum* a, const bignum* b)
{
  assert(q == NULL || q != r);

  u64 an = bn_nsize(a);
  u64 vn = bn_nsize(b);
  if (BN_UNLIKELY(vn == 0)) {
    assert(!"bn_divmod: division by zero");
    return;
  }

  const bool qneg = a->is_neg ^ b->is_neg;
  const bool rneg = a->is_neg;

  /* |a| < |b|  ->  q = 0, r = a   (compare magnitudes, never bn_cmp!) */
  if (an < vn || (an == vn && limbs_cmp(a->limbs, b->limbs, an) < 0)) {
    if (r) {
      bn_copy(r, a);
      bn_trim(r);
      r->is_neg = rneg && !bn_is_zero(r);
    }
    if (q) bn_set_u64(q, 0);
    return;
  }

  const u64 qn = an - vn + 1;

  /* aliasing: work into temporaries only when we really have to */
  bignum tq, tr;
  bignum *Q = q, *R = r;
  bool alq = q && (q == a || q == b);
  bool alr = r && (r == a || r == b);
  if (alq) {
    bn_init(&tq);
    Q = &tq;
  }
  if (alr) {
    bn_init(&tr);
    R = &tr;
  }
  if (Q) bn_alloc(Q, qn);
  if (R) bn_alloc(R, vn);

  if (vn == 1) {
    u64 rem = limbs_divrem_1(Q ? Q->limbs : NULL, a->limbs, an, b->limbs[0]);
    if (R) {
      R->limbs[0] = rem;
      R->size = 1;
    }
  } else {
    unsigned s = (unsigned)__builtin_clzll(b->limbs[vn - 1]);
    u64 need = an + 1 + (s ? vn : 0) + (Q ? 0 : qn);
    u64 stackbuf[BN_DIV_STACK_LIMBS];
    u64* scratch = (need <= BN_DIV_STACK_LIMBS)
                       ? stackbuf
                       : (u64*)malloc(need * sizeof(u64));
    u64* qp = Q ? Q->limbs : scratch + an + 1 + (s ? vn : 0);

    bn_divmod_limbs(qp, R ? R->limbs : NULL, a->limbs, an, b->limbs, vn,
                    scratch);

    if (scratch != stackbuf) free(scratch);
    if (R) R->size = vn;
  }

  if (Q) {
    Q->size = qn;
    bn_trim(Q);
    Q->is_neg = qneg && !bn_is_zero(Q);
  }
  if (R) {
    bn_trim(R);
    R->is_neg = rneg && !bn_is_zero(R);
  }

  if (alq) {
    bn_copy(q, &tq);
    bn_free(&tq);
  }
  if (alr) {
    bn_copy(r, &tr);
    bn_free(&tr);
  }
}

/*
 * q = a / b (truncated division, C semantics).
 *
 * Let n = a->size and m = b->size, measured in 64-bit limbs.
 *
 * Thin wrapper around bn_divmod() that discards the remainder.
 *
 * Complexity:
 *   Time: O(n) for a single-limb divisor, O((n - m + 1) * m) otherwise
 *   Auxiliary memory: O(n + m) limbs
 *   Output memory: O(n - m + 1) limbs
 */
void bn_div(bignum* q, const bignum* a, const bignum* b)
{
  bn_divmod(q, NULL, a, b);
}

/*
 * r = a % b (truncated remainder, C semantics: sign(r) == sign(a)).
 *
 * Let n = a->size and m = b->size, measured in 64-bit limbs.
 *
 * Thin wrapper around bn_divmod() that discards the quotient.
 *
 * Complexity:
 *   Time: O(n) for a single-limb divisor, O((n - m + 1) * m) otherwise
 *   Auxiliary memory: O(n + m) limbs
 *   Output memory: O(m) limbs
 */
void bn_mod(bignum* r, const bignum* a, const bignum* b)
{
  bn_divmod(NULL, r, a, b);
}

/*===========================================================================
 *  7. Exact division (Jebelean)
 *
 *  Only valid when b | a, but ~2x faster than a real division and no
 *  normalisation / correction steps at all.
 *
 *  Idea: strip the common power of two from a and b, then compute the
 *  quotient from the bottom up. Since b[0] is odd, each quotient limb
 *  qi = Ai * b[0]^{-1} mod 2^64 is forced (it must cancel limb i of
 *  the running remainder), so no estimation or correction is needed.
 *=========================================================================*/

/*
 * q = a / b, assuming b divides a exactly (Jebelean's exact division).
 *
 * Let n = a->size and m = b->size, measured in 64-bit limbs.
 *
 * Strips the common power of two from a and b (so the low limb of the
 * divisor becomes odd), then computes the quotient limbs from the
 * bottom up: with D[0] odd, each quotient limb qi = A[i] * D[0]^{-1}
 * (mod 2^64) is forced, because it must cancel limb i of the running
 * remainder. No estimation or correction steps are needed, which makes
 * this about 2x faster than a real division.
 *
 * The behaviour is undefined if b does not divide a. If b is zero, q
 * is left unchanged (asserts in debug builds). If a is zero, q is set
 * to 0. q may alias a or b.
 *
 * Complexity:
 *   Time: O(n * m) - one O(m) submul per quotient limb
 *   Auxiliary memory: O(n + m) limbs for the shifted copies
 *   Output memory: O(n - m + 1) limbs
 */
void bn_div_exact(bignum* q, const bignum* a, const bignum* b)
{
  u64 an = bn_nsize(a), bn_size = bn_nsize(b);
  if (BN_UNLIKELY(bn_size == 0)) {
    assert(!"bn_div_exact: /0");
    return;
  }
  if (an == 0) {
    bn_set_u64(q, 0);
    return;
  }

  /* strip the power of two from b (and the same amount from a) */
  u64 kl = 0;
  while (b->limbs[kl] == 0) kl++;
  unsigned kb = (unsigned)__builtin_ctzll(b->limbs[kl]);
  assert(an > kl);

  u64 dl = bn_size - kl, al = an - kl;
  u64* buf = (u64*)malloc((dl + al) * sizeof(u64));
  u64 *D = buf, *A = buf + dl;
  limbs_rshift(D, b->limbs + kl, dl, kb);
  limbs_rshift(A, a->limbs + kl, al, kb);
  dl = limbs_norm(D, dl);
  al = limbs_norm(A, al);

  if (al < dl) {
    bn_set_u64(q, 0);
    free(buf);
    return;
  }

  u64 qn = al - dl + 1;
  bignum tq;
  bignum* Q = q;
  bool al_ = (q == a || q == b);
  if (al_) {
    bn_init(&tq);
    Q = &tq;
  }
  bn_alloc(Q, qn);
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

  bn_trim(Q);
  Q->is_neg = (a->is_neg ^ b->is_neg) && !bn_is_zero(Q);
  if (al_) {
    bn_copy(q, &tq);
    bn_free(&tq);
  }
  free(buf);
}
