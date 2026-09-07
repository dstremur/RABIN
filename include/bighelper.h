#ifndef BIGHELPER_H
#define BIGHELPER_H

/*===========================================================================
 *  bighelper.h
 *
 *  Low-level helpers that operate directly on raw limb arrays (u64*),
 *  in contrast to the bignum API in bignum.h.
 *
 *  Nothing in here allocates memory, trims, or handles signs: every
 *  function works on plain u64 buffers of a given length, so callers
 *  are responsible for sizing the buffers correctly.
 *
 *  Layout:
 *    - carry / borrow primitives        (adc64, sbb64)
 *    - normalisation & comparison       (limbs_norm, limbs_cmp)
 *    - bit shifts                       (limbs_lshift, limbs_rshift)
 *    - add / subtract                   (limbs_add_n, limbs_sub_1,
 *                                        limbs_submul_1, limbs_add_raw)
 *    - multiplication                   (limbs_mul_school,
 *                                        limbs_mul_karatsuba,
 *                                        limbs_sqr_karatsuba)
 *    - assembly kernels                 (bn_add_inner, bn_sub_inner,
 *                                        bn_mul_add_inner,
 *                                        bn_sub_mul_digit)
 *    - single-limb helpers              (mod_inverse_u64)
 *===========================================================================*/

#include <stddef.h>
#include <string.h>

#include "bigcore.h"

#if defined(__x86_64__) || defined(_M_X64)
#include <immintrin.h>
#endif

typedef unsigned __int128 u128;

/* limb count below which multiplication falls back to schoolbook */
#define KARATSUBA_LIMIT 128

/*---------------------------------------------------------------------------
 *  carry / borrow primitives
 *-------------------------------------------------------------------------*/

/**
 * @brief o = c + a + b, returns the carry out (0 or 1).
 *
 * @param[in]  c Carry in.
 * @param[in]  a First operand.
 * @param[in]  b Second operand.
 * @param[out] o Result.
 *
 * @return The carry out (0 or 1).
 */
static inline unsigned char adc64(unsigned char c, u64 a, u64 b, u64* o)
{
#if defined(__x86_64__) || defined(_M_X64)
  return _addcarry_u64(c, a, b, (unsigned long long*)o);
#else
  u128 t = (u128)a + b + c;
  *o = (u64)t;
  return (unsigned char)(t >> 64);
#endif
}

/**
 * @brief o = c - a - b, returns the borrow out (0 or 1).
 *
 * @param[in]  c Borrow in.
 * @param[in]  a First operand.
 * @param[in]  b Second operand.
 * @param[out] o Result.
 *
 * @return The borrow out (0 or 1).
 */
static inline unsigned char sbb64(unsigned char c, u64 a, u64 b, u64* o)
{
#if defined(__x86_64__) || defined(_M_X64)
  return _subborrow_u64(c, a, b, (unsigned long long*)o);
#else
  u128 t = (u128)a - b - c;
  *o = (u64)t;
  return (unsigned char)((t >> 64) & 1);
#endif
}

/*---------------------------------------------------------------------------
 *  normalisation & comparison
 *-------------------------------------------------------------------------*/

/**
 * @brief number of significant limbs (skips trailing zeros).
 *
 * @param[in] p Limb array.
 * @param[in] n Number of limbs.
 *
 * @return The number of significant limbs.
 */
static inline u64 limbs_norm(const u64* p, u64 n)
{
  while (n && p[n - 1] == 0) n--;
  return n;
}

/**
 * @brief unsigned magnitude compare: 1 if a > b, -1 if a < b, 0 if equal.
 *
 * @param[in] a First operand limbs.
 * @param[in] b Second operand limbs.
 * @param[in] n Number of limbs.
 *
 * @return 1 If a > b, -1 if a < b, 0 if equal.
 */
static inline int limbs_cmp(const u64* a, const u64* b, u64 n)
{
  while (n--)
    if (a[n] != b[n]) return a[n] > b[n] ? 1 : -1;
  return 0;
}

/*---------------------------------------------------------------------------
 *  bit shifts
 *-------------------------------------------------------------------------*/

/**
 * @brief rp = up << s, 0 < s < 64, returns the bits shifted out of the top.
 *
 * @param[out] rp   Result buffer.
 * @param[in]  up   Source limbs.
 * @param[in]  n    Number of limbs.
 * @param[in]  s    Shift amount (0 < s < 64).
 *
 * @return The bits shifted out of the top.
 */
static inline u64 limbs_lshift(u64* rp, const u64* up, u64 n, unsigned s)
{
  u64 high = up[n - 1] >> (64 - s);
  for (u64 i = n - 1; i > 0; i--)
    rp[i] = (up[i] << s) | (up[i - 1] >> (64 - s));
  rp[0] = up[0] << s;
  return high;
}

/**
 * @brief rp = up >> s, 0 <= s < 64 (rp may alias up).
 *
 * @param[out] rp   Result buffer (may alias up).
 * @param[in]  up   Source limbs.
 * @param[in]  n    Number of limbs.
 * @param[in]  s    Shift amount (0 <= s < 64).
 */
static inline void limbs_rshift(u64* rp, const u64* up, u64 n, unsigned s)
{
  if (s == 0) {
    if (rp != up) memmove(rp, up, n * sizeof(u64));
    return;
  }
  for (u64 i = 0; i + 1 < n; i++)
    rp[i] = (up[i] >> s) | (up[i + 1] << (64 - s));
  rp[n - 1] = up[n - 1] >> s;
}

/*---------------------------------------------------------------------------
 *  add / subtract
 *-------------------------------------------------------------------------*/

/**
 * @brief rp += vp (n limbs), returns the carry out.
 *
 * @param[in,out] rp Accumulator (modified in place).
 * @param[in]     vp Value to add.
 * @param[in]     n  Number of limbs.
 *
 * @return The carry out.
 */
static inline u64 limbs_add_n(u64* rp, const u64* vp, u64 n)
{
  unsigned char c = 0;
  for (u64 i = 0; i < n; i++) c = adc64(c, rp[i], vp[i], &rp[i]);
  return c;
}

/**
 * @brief p -= v, with borrow propagation over n limbs.
 *
 * @param[in,out] p Limb array (modified in place).
 * @param[in]     n Number of limbs.
 * @param[in]     v Value to subtract.
 */
static inline void limbs_sub_1(u64* p, u64 n, u64 v)
{
  for (u64 i = 0; i < n && v; i++) {
    u64 x = p[i];
    p[i] = x - v;
    v = (x < v);
  }
}

/**
 * @brief rp -= vp * q  (n limbs), returns the borrow out of the top.
 *
 * The multiply carry and the subtract borrow are kept in two independent
 * carry chains so the CPU can overlap them; on ADX hardware gcc/clang
 * lower this to mulx + adcx + sbb.
 *
 * @param[in,out] rp Accumulator (modified in place).
 * @param[in]     vp Multiplicand limbs.
 * @param[in]     n  Number of limbs.
 * @param[in]     q  Single-limb multiplier.
 *
 * @return The borrow out of the top.
 */
static inline u64 limbs_submul_1(u64* rp, const u64* vp, u64 n, u64 q)
{
  u64 mcy = 0;          /* carry out of the multiply, always <= 2^64-2 */
  unsigned char bw = 0; /* borrow of the subtract chain               */
  u64 i = 0;

  for (; i + 4 <= n; i += 4) { /* unrolled: shortens the loop  */
    u128 p0 = (u128)vp[i + 0] * q + mcy;
    u64 l0 = (u64)p0;
    mcy = (u64)(p0 >> 64);
    u128 p1 = (u128)vp[i + 1] * q + mcy;
    u64 l1 = (u64)p1;
    mcy = (u64)(p1 >> 64);
    u128 p2 = (u128)vp[i + 2] * q + mcy;
    u64 l2 = (u64)p2;
    mcy = (u64)(p2 >> 64);
    u128 p3 = (u128)vp[i + 3] * q + mcy;
    u64 l3 = (u64)p3;
    mcy = (u64)(p3 >> 64);
    bw = sbb64(bw, rp[i + 0], l0, &rp[i + 0]);
    bw = sbb64(bw, rp[i + 1], l1, &rp[i + 1]);
    bw = sbb64(bw, rp[i + 2], l2, &rp[i + 2]);
    bw = sbb64(bw, rp[i + 3], l3, &rp[i + 3]);
  }
  for (; i < n; i++) {
    u128 p = (u128)vp[i] * q + mcy;
    u64 l = (u64)p;
    mcy = (u64)(p >> 64);
    bw = sbb64(bw, rp[i], l, &rp[i]);
  }
  return mcy + bw; /* provably cannot overflow: mcy <= 2^64-2 */
}

/**
 * @brief r = a + b (r must not alias a or b), returns the number of limbs
 * in the result (max(a_len, b_len) or +1 on carry).
 *
 * @param[out] r     Result buffer (must not alias a or b).
 * @param[in]  a     First operand limbs.
 * @param[in]  a_len Number of limbs in a.
 * @param[in]  b     Second operand limbs.
 * @param[in]  b_len Number of limbs in b.
 *
 * @return The number of limbs in the result.
 */
u64 limbs_add_raw(u64* r, const u64* a, u64 a_len, const u64* b, u64 b_len);

/*---------------------------------------------------------------------------
 *  multiplication
 *-------------------------------------------------------------------------*/

/**
 * @brief r = a * b, r must have room for a_size + b_size limbs.
 *
 * @param[out] r      Result buffer (needs a_size + b_size limbs).
 * @param[in]  a      First operand limbs.
 * @param[in]  a_size Number of limbs in a.
 * @param[in]  b      Second operand limbs.
 * @param[in]  b_size Number of limbs in b.
 */
void limbs_mul_school(u64* r, const u64* a, u64 a_size, const u64* b,
                      u64 b_size);

/**
 * @brief r = a * b (both n limbs), r must have room for 2n limbs.
 *
 * @param[out] r       Result buffer (needs 2n limbs).
 * @param[in]  a       First operand limbs.
 * @param[in]  b       Second operand limbs.
 * @param[in]  n       Number of limbs in each operand.
 * @param[in]  scratch Scratch buffer for intermediate results.
 */
void limbs_mul_karatsuba(u64* r, const u64* a, const u64* b, u64 n,
                         u64* scratch);

/**
 * @brief r = a * a (a_len limbs), r must have room for 2 * a_len limbs.
 *
 * @param[out] r       Result buffer (needs 2 * a_len limbs).
 * @param[in]  a       Operand limbs.
 * @param[in]  a_len   Number of limbs in a.
 * @param[in]  scratch Scratch buffer for intermediate results.
 */
void limbs_sqr_karatsuba(u64* r, const u64* a, u64 a_len, u64* scratch);

/*---------------------------------------------------------------------------
 *  assembly kernels (src/arch)
 *-------------------------------------------------------------------------*/

/**
 * @brief r = a + b (a_size >= b_size), returns the carry out.
 *
 * @param[out] r      Result buffer.
 * @param[in]  a      First operand limbs (the longer one).
 * @param[in]  a_size Number of limbs in a.
 * @param[in]  b      Second operand limbs.
 * @param[in]  b_size Number of limbs in b.
 *
 * @return The carry out.
 */
u64 bn_add_inner(u64* r, const u64* a, u64 a_size, const u64* b, u64 b_size);

/**
 * @brief r = a - b (a_size >= b_size), returns the borrow out.
 *
 * @param[out] r      Result buffer.
 * @param[in]  a      First operand limbs (the longer one).
 * @param[in]  a_size Number of limbs in a.
 * @param[in]  b      Second operand limbs.
 * @param[in]  b_size Number of limbs in b.
 *
 * @return The borrow out.
 */
u64 bn_sub_inner(u64* r, const u64* a, u64 a_size, const u64* b, u64 b_size);

/**
 * @brief r += a_limb * b (len limbs), returns the carry out.
 *
 * @param[in,out] r      Accumulator (modified in place).
 * @param[in]     b      Limb array.
 * @param[in]     a_limb Single-limb multiplier.
 * @param[in]     len    Number of limbs in b.
 *
 * @return The carry out.
 */
u64 bn_mul_add_inner(u64* r, const u64* b, u64 a_limb, u64 len);

/**
 * @brief r -= a_limb * b (len limbs), returns the borrow out.
 *
 * @param[in,out] r      Accumulator (modified in place).
 * @param[in]     b      Limb array.
 * @param[in]     a_limb Single-limb multiplier.
 * @param[in]     len    Number of limbs in b.
 *
 * @return The borrow out.
 */
u64 bn_sub_mul_digit(u64* r, const u64* b, u64 a_limb, u64 len);

/*---------------------------------------------------------------------------
 *  single-limb helpers
 *-------------------------------------------------------------------------*/

/**
 * @brief Newton iteration: returns -n^{-1} mod 2^64 (n must be odd).
 *
 * Used by Montgomery reduction and by Jebelean exact division.
 *
 * @param[in] n Odd modulus.
 *
 * @return -n^{-1} mod 2^64.
 */
static inline u64 mod_inverse_u64(u64 n)
{
  u64 inv = 1;
  for (int i = 0; i < 6; i++) {
    inv *= (2 - n * inv);
  }
  return -inv;
}

#endif
