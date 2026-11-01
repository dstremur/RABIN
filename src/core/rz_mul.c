/*
 * rz_mul.c
 *
 * rz_t multiplication routines.
 *
 * This file implements signed rz_t multiplication with a schoolbook
 * and a Karatsuba path (dispatched by size), squaring, and an NTT-based
 * fast multiplication path for very large operands. The raw-limb
 * kernels live in bighelper.c.
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

#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>

#include "../../include/precomp.h"
#include "../../include/rabin.h"
#include "../../include/rzlimb.h"
#include "../../include/u64.h"
#include "rz_internal.h"

/* Fitted cost-model gate for the NTT kernels (benchmark on x86-64):
 * the padded Karatsuba fallback costs ~ 2.1e-9 * max_len^1.585 s and
 * depends only on max_len (both operands are zero-padded), while the
 * NTT costs c * N * log2(N) with N the transform length in 32-bit
 * chunks and c ~ 5.3e-9 for a product (two forward + one inverse
 * transform per prime) or ~ 3.8e-9 for a square (one forward + one
 * inverse per prime). Asymmetry therefore enters only through N:
 * lopsided pairs shrink the transform and FAVOR the NTT, which is why
 * the old 4:1 asymmetry guard is replaced by this model. Callers keep
 * the RZ_NTT_LIMIT floor and the RZ_KARATSUBA_LIMIT fallback gate. */
static bool rz_ntt_beats_karat(u64 na, u64 nb, bool squaring)
{
  u64 mx = RZ_MAX(na, nb);
  /* nextpow2 of the chunk convolution length 2*(na+nb)-1 */
  double N = ldexp(1.0, ceil(log2(2.0 * ((double)na + (double)nb))));
  double kar = 2.1e-9 * pow((double)mx, 1.585);
  double c = squaring ? 1.7e-9 : 2.4e-9;
  return c * N * log2(N) < kar;
}

/* Two-prime RNS exactness. q1 = 5*2^55+1 (golden) and q2 = 69*2^55+1
   both support transform lengths up to 2^54, and
   q1*q2 ~ 2^118.4 > 2^64 * min(na, nb) for every coefficient sum
   c_i < 2^64 * min(na, nb) reachable under the transform cap
   min(na, nb) <= nz <= 2^54. So with 32-bit chunks the CRT merge
   below reconstructs every chunk exactly. */
#define RZ_RNS_Q1_INV_Q2 2291769260378161154ULL /* q1^{-1} mod q2 */

/* CRT merge of one convolution chunk: r1 < q1 < q2 and r2 < q2 are the
   canonical residues of the same chunk; returns the exact value
   c < q1*q2 with c = r1 (mod q1), c = r2 (mod q2). inv_mont is
   u64_mont_in(RZ_RNS_Q1_INV_Q2, &c2->mctx). */
static inline u128 rz_rns_crt_chunk(u64 r1, u64 r2, u64 q1, u64 q2,
                                    u64 inv_mont, const u64_mont_ctx_t* m2)
{
  /* r1 < q1 < q2, so r2 + q2 - r1 lies in (0, 2*q2) with q2 < 2^62:
     no overflow, one masked subtract lands in [0, q2). */
  u64 d = r2 + q2 - r1;
  d -= q2 & -(u64)(d >= q2);
  u64 t = u64_mont_mul(d, inv_mont, m2); /* t = d * q1^{-1} mod q2 */
  return (u128)r1 + (u128)q1 * t;        /* exact: c < q1*q2 < 2^119 */
}

/* Carry-propagate and pack the CRT-reconstructed 32-bit chunks
 * (r_arr: residues mod q2, r1s: residues mod q1) into r->limbs,
 * two chunks per limb. nz is odd, so the operand bound gives at most
 * ceil(nz/2) = (nz+1)/2 limbs: the true product needs no more.
 * Sets r->size; leaves r->is_neg to the caller. */
static void rz_rns_pack(rz_t* r, const u64* r_arr, const u64* r1s, u64 nz,
                        u64 q1, u64 q2, u64 inv_mont, const u64_mont_ctx_t* m2)
{
  u128 carry = 0;
  u64 last = 0, nlimb = 0;
  for (u64 j = 0; 2 * j < nz || carry; j++) {
    u128 t0 = carry;
    if (2 * j < nz)
      t0 += rz_rns_crt_chunk(r1s[2 * j], r_arr[2 * j], q1, q2, inv_mont, m2);
    u64 lo = (u64)t0 & 0xFFFFFFFFu;
    carry = t0 >> 32;

    u128 t1 = carry;
    if (2 * j + 1 < nz)
      t1 += rz_rns_crt_chunk(r1s[2 * j + 1], r_arr[2 * j + 1], q1, q2, inv_mont,
                             m2);
    u64 hi = (u64)t1 & 0xFFFFFFFFu;
    carry = t1 >> 32;

    u64 L = lo | (hi << 32);
    if (L != 0) last = nlimb;
    r->limbs[nlimb++] = L;
  }
  r->size = last + 1;
}

/* Unpack a magnitude into 32-bit chunks in arr (2 per limb) and zero
   the transform padding [2*size, ntt_size). The inverse transform
   clobbers the unpacked buffer (it doubles as the output), so the
   second prime's pass unpacks again. */
static void rz_ntt_unpack_w32(const rz_t* a, u64* arr, u64 ntt_size)
{
  u64 na = a->size * 2;
  for (u64 i = 0; i < a->size; i++) {
    arr[2 * i + 0] = a->limbs[i] & 0xFFFFFFFFu;
    arr[2 * i + 1] = a->limbs[i] >> 32;
  }
  if (na < ntt_size) memset(arr + na, 0, (ntt_size - na) * sizeof(u64));
}

/* Exact two-prime (W=32 chunk) NTT multiplication of magnitudes.
 * a, b must be nonzero; r must not alias a or b (dispatch handles it).
 * Does NOT set r->is_neg (caller does).
 * Five scratch arrays: the inverse transform reuses the dead a_arr
 * region, and r1s keeps the first prime's canonical residues (nz <=
 * ntt_size entries). The inverse clobbers a_arr/b_arr, so the second
 * prime's pass unpacks the operands again. */
static rabin_err_t rz_mul_ntt_u64x2(rz_t* r, const rz_t* a, const rz_t* b)
{
  u64 na, nb, nz, ntt_size, k, i;
  const u64_ntt_ctx_t *c1, *c2;
  u64 *buf, *a_arr, *b_arr, *a_hat, *b_hat, *r_arr, *r1s;
  rabin_err_t err;

  if (a->size > SIZE_MAX / 2 || b->size > SIZE_MAX / 2)
    return RABIN_ERR_OVERFLOW;

  na = a->size * 2; /* 32-bit chunks */
  nb = b->size * 2;
  if (na > UINT64_MAX - (nb - 1)) return RABIN_ERR_OVERFLOW; /* nz overflow */
  nz = na + nb - 1; /* true product length in chunks */

  ntt_size = 1;
  k = 0;
  while (k < 54 && ntt_size < nz) {
    ntt_size <<= 1;
    k++;
  }

  /* Both primes cover k <= 54, and the chunk sums stay below q1*q2
     there (see RZ_RNS_Q1_INV_Q2 note): only size can fail, which
     RABIN_ERR_OVERFLOW reports so the caller falls back to Karatsuba. */
  if (ntt_size < nz) return RABIN_ERR_OVERFLOW; /* transform size k > 54 */
  if (ntt_size > SIZE_MAX / (5 * sizeof(u64))) return RABIN_ERR_OVERFLOW;

  c1 = u64_ntt_ctx_golden_cached(k);
  c2 = u64_ntt_ctx_proth2_cached(k);
  if (c1 == NULL || c2 == NULL) return RABIN_ERR_OUT_OF_MEMORY;
  u64 q1 = c1->q, q2 = c2->q;
  u64 inv_mont = u64_mont_in(RZ_RNS_Q1_INV_Q2, &c2->mctx);

  buf = rz_scratch_get(5 * ntt_size);
  if (buf == NULL) return RABIN_ERR_OUT_OF_MEMORY;
  a_arr = buf;
  b_arr = buf + ntt_size;
  a_hat = buf + 2 * ntt_size;
  b_hat = buf + 3 * ntt_size;
  r1s = buf + 4 * ntt_size;
  r_arr = a_arr; /* dead after the forward transforms */

  /* Unpack: 2 x 32-bit chunks per 64-bit limb. */
  rz_ntt_unpack_w32(a, a_arr, ntt_size);
  rz_ntt_unpack_w32(b, b_arr, ntt_size);

  /* Prime q1: residues land canonical in r_arr, kept in r1s. */
  u64_ntt_cyclic_forward(a_hat, a_arr, c1);
  u64_ntt_cyclic_forward(b_hat, b_arr, c1);
  for (i = 0; i < ntt_size; i++)
    a_hat[i] = u64_mont_mul(a_hat[i], b_hat[i], &c1->mctx);
  u64_ntt_cyclic_inverse_montgomery_in(r_arr, a_hat, c1);
  memcpy(r1s, r_arr, nz * sizeof(u64));

  /* Prime q2: r_arr (aliased with a_arr) was overwritten, re-unpack. */
  rz_ntt_unpack_w32(a, a_arr, ntt_size);
  rz_ntt_unpack_w32(b, b_arr, ntt_size);
  u64_ntt_cyclic_forward(a_hat, a_arr, c2);
  u64_ntt_cyclic_forward(b_hat, b_arr, c2);
  for (i = 0; i < ntt_size; i++)
    a_hat[i] = u64_mont_mul(a_hat[i], b_hat[i], &c2->mctx);
  u64_ntt_cyclic_inverse_montgomery_in(r_arr, a_hat, c2);

  rz_rns_pack(r, r_arr, r1s, nz, q1, q2, inv_mont, &c2->mctx);

  err = RABIN_SUCCESS;
  rz_scratch_release();
  return err;
}

/* Exact two-prime (W=32 chunk) NTT squaring of a magnitude.
 * a must be nonzero; r must not alias a (dispatch handles it).
 * One forward transform per prime; three scratch arrays suffice:
 * a_hat squared pointwise, inverse transform into the dead a_arr. */
static rabin_err_t rz_sqr_ntt_u64x2(rz_t* r, const rz_t* a)
{
  u64 na, nz, ntt_size, k, i;
  const u64_ntt_ctx_t *c1, *c2;
  u64 *buf, *a_arr, *a_hat, *r1s;
  rabin_err_t err;

  if (a->size > SIZE_MAX / 2) return RABIN_ERR_OVERFLOW;

  na = a->size * 2;                                          /* 32-bit chunks */
  if (na > UINT64_MAX - (na - 1)) return RABIN_ERR_OVERFLOW; /* nz overflow */
  nz = 2 * na - 1; /* true square length in chunks */

  ntt_size = 1;
  k = 0;
  while (k < 54 && ntt_size < nz) {
    ntt_size <<= 1;
    k++;
  }

  if (ntt_size < nz) return RABIN_ERR_OVERFLOW; /* transform size k > 54 */
  if (ntt_size > SIZE_MAX / (3 * sizeof(u64))) return RABIN_ERR_OVERFLOW;

  c1 = u64_ntt_ctx_golden_cached(k);
  c2 = u64_ntt_ctx_proth2_cached(k);
  if (c1 == NULL || c2 == NULL) return RABIN_ERR_OUT_OF_MEMORY;
  u64 q1 = c1->q, q2 = c2->q;
  u64 inv_mont = u64_mont_in(RZ_RNS_Q1_INV_Q2, &c2->mctx);

  buf = rz_scratch_get(3 * ntt_size);
  if (buf == NULL) return RABIN_ERR_OUT_OF_MEMORY;
  a_arr = buf;
  a_hat = buf + ntt_size;
  r1s = buf + 2 * ntt_size;

  /* Unpack: 2 x 32-bit chunks per 64-bit limb. */
  rz_ntt_unpack_w32(a, a_arr, ntt_size);

  /* Prime q1 (inverse lands in the dead a_arr), save residues. */
  u64_ntt_cyclic_forward(a_hat, a_arr, c1);
  for (i = 0; i < ntt_size; i++)
    a_hat[i] = u64_mont_mul(a_hat[i], a_hat[i], &c1->mctx);
  u64_ntt_cyclic_inverse_montgomery_in(a_arr, a_hat, c1);
  memcpy(r1s, a_arr, nz * sizeof(u64));

  /* Prime q2: a_arr was overwritten by the inverse, re-unpack. */
  rz_ntt_unpack_w32(a, a_arr, ntt_size);
  u64_ntt_cyclic_forward(a_hat, a_arr, c2);
  for (i = 0; i < ntt_size; i++)
    a_hat[i] = u64_mont_mul(a_hat[i], a_hat[i], &c2->mctx);
  u64_ntt_cyclic_inverse_montgomery_in(a_arr, a_hat, c2);

  rz_rns_pack(r, a_arr, r1s, nz, q1, q2, inv_mont, &c2->mctx);

  err = RABIN_SUCCESS;
  rz_scratch_release();
  return err;
}

rabin_err_t rz_sqr(rz_t* r, const rz_t* a)
{
  if (r == NULL || a == NULL) return RABIN_ERR_NULL_PTR;

  if (a->size == 0) {
    r->size = 0;
    if (r->capacity > 0) r->limbs[0] = 0;
    return RABIN_SUCCESS;
  }

  if (r == a) {
    rz_t tmp;
    rz_init(&tmp);
    rabin_err_t err = rz_sqr(&tmp, a);
    if (err != RABIN_SUCCESS) {
      rz_clear(&tmp);
      return err;
    }
    err = rz_copy(r, &tmp);
    rz_clear(&tmp);
    return err;
  }

  rabin_err_t err = rz_alloc(r, 2 * a->size);
  if (err != RABIN_SUCCESS) return err;

  if (a->size >= RZ_NTT_LIMIT && rz_ntt_beats_karat(a->size, a->size, true)) {
    rabin_err_t nerr = rz_sqr_ntt_u64x2(r, a);
    if (nerr == RABIN_SUCCESS) {
      r->is_neg = false;
      rz_trim(r);
      return RABIN_SUCCESS;
    }
    if (nerr != RABIN_ERR_OVERFLOW && nerr != RABIN_ERR_OUT_OF_MEMORY)
      return nerr;
    /* RABIN_ERR_OVERFLOW (past exactness bound) or OOM: fall through to
       Karatsuba squaring. */
  }

  u64 scratch_size = 8 * a->size + 8;  // conservative heuristic
  u64* scratch = rz_scratch_get(scratch_size);
  if (scratch == NULL) return RABIN_ERR_OUT_OF_MEMORY;

  limbs_sqr_karatsuba(r->limbs, a->limbs, a->size, scratch);

  rz_scratch_release();

  r->size = 2 * a->size;
  r->is_neg = false;
  rz_trim(r);
  return RABIN_SUCCESS;
}

rabin_err_t rz_mul(rz_t* r, const rz_t* a, const rz_t* b)
{
  if (r == NULL || a == NULL || b == NULL) return RABIN_ERR_NULL_PTR;

  if (a->size == 0 || b->size == 0) {
    r->size = 0;
    if (r->capacity > 0) r->limbs[0] = 0;
    return RABIN_SUCCESS;
  }

  if (r == a || r == b) {
    rz_t tmp;
    rz_init(&tmp);
    rabin_err_t err = rz_mul(&tmp, a, b);
    if (err != RABIN_SUCCESS) {
      rz_clear(&tmp);
      return err;
    }
    err = rz_copy(r, &tmp);
    rz_clear(&tmp);
    return err;
  }

  if (a == b) {
    return rz_sqr(r, a);
  }

  u64 max_len = RZ_MAX(a->size, b->size);
  rabin_err_t err = rz_alloc(r, 2 * max_len);
  if (err != RABIN_SUCCESS) return err;

  // use the u64-NTT path for large pairs where the fitted cost model
  // predicts the transform beats padded Karatsuba (asymmetric pairs
  // often shrink N enough for the NTT to win; both operands must be
  // Karatsuba-sized, otherwise the schoolbook fallback is cheaper)
  if (max_len >= RZ_NTT_LIMIT &&
      RZ_MIN(a->size, b->size) >= RZ_KARATSUBA_LIMIT &&
      rz_ntt_beats_karat(a->size, b->size, false)) {
    rabin_err_t nerr = rz_mul_ntt_u64x2(r, a, b);
    if (nerr == RABIN_SUCCESS) {
      /* the kernel already set r->size; do NOT overwrite it */
      r->is_neg = a->is_neg ^ b->is_neg;
      rz_trim(r);
      return RABIN_SUCCESS;
    }
    if (nerr != RABIN_ERR_OVERFLOW && nerr != RABIN_ERR_OUT_OF_MEMORY)
      return nerr; /* real error */
    /* RABIN_ERR_OVERFLOW (past exactness bound) or OOM: fall through to
       Karatsuba/schoolbook. */
  }

  // use karasuba if both inputs are larger then the limit
  if (a->size >= RZ_KARATSUBA_LIMIT && b->size >= RZ_KARATSUBA_LIMIT) {
    // One arena block: pad_a, pad_b, then the Karatsuba scratch
    uint64_t* buf = rz_scratch_get(10 * max_len);
    if (buf == NULL) return RABIN_ERR_OUT_OF_MEMORY;
    uint64_t* pad_a = buf;
    uint64_t* pad_b = buf + max_len;
    uint64_t* scratch = buf + 2 * max_len;

    // Normalize asymmetric arrays with zero-padding up front
    memset(pad_a, 0, max_len * sizeof(uint64_t));
    memset(pad_b, 0, max_len * sizeof(uint64_t));
    memcpy(pad_a, a->limbs, a->size * sizeof(uint64_t));
    memcpy(pad_b, b->limbs, b->size * sizeof(uint64_t));

    limbs_mul_karatsuba(r->limbs, pad_a, pad_b, max_len, scratch);

    rz_scratch_release();
  } else {
    limbs_mul_school(r->limbs, a->limbs, a->size, b->limbs, b->size);
  }

  r->size = a->size + b->size;
  r->is_neg = a->is_neg ^ b->is_neg;
  rz_trim(r);
  return RABIN_SUCCESS;
}

rabin_err_t rz_mul_school(rz_t* r, const rz_t* a, const rz_t* b)
{
  if (r == NULL || a == NULL || b == NULL) return RABIN_ERR_NULL_PTR;

  u64 max = a->size + b->size;
  rabin_err_t err = rz_alloc(r, max);
  if (err != RABIN_SUCCESS) return err;

  memset(r->limbs, 0, max * sizeof(u64));
  r->size = max;
  r->is_neg = a->is_neg ^ b->is_neg;

  for (u64 i = 0; i < a->size; i++) {
    if (a->limbs[i] == 0) continue;

    rz_mul_add_inner(&r->limbs[i], b->limbs, a->limbs[i], b->size);
  }

  rz_trim(r);
  return RABIN_SUCCESS;
}

rabin_err_t rz_mul_fast(rz_t* res, const rz_t* a, const rz_t* b)
{
  if (res == NULL || a == NULL || b == NULL) return RABIN_ERR_NULL_PTR;

  if (a->size == 0 || b->size == 0) {
    res->size = 0;
    if (res->capacity > 0) res->limbs[0] = 0;
    return RABIN_SUCCESS;
  }

  if (res == a || res == b) {
    rz_t tmp;
    rz_init(&tmp);
    rabin_err_t err = rz_mul_fast(&tmp, a, b);
    if (err != RABIN_SUCCESS) {
      rz_clear(&tmp);
      return err;
    }
    err = rz_copy(res, &tmp);
    rz_clear(&tmp);
    return err;
  }

  rabin_err_t err = rz_alloc(res, 2 * RZ_MAX(a->size, b->size));
  if (err != RABIN_SUCCESS) return err;

  err = rz_mul_ntt_u64x2(res, a, b);
  if (err != RABIN_SUCCESS) return err;

  res->is_neg = a->is_neg ^ b->is_neg;
  rz_trim(res);
  return RABIN_SUCCESS;
}

rabin_err_t rz_mul_i64(rz_t* r, const rz_t* a, const i64 c)
{
  if (r == NULL || a == NULL) return RABIN_ERR_NULL_PTR;

  // a zero factor or zero multiplier gives a clean non-negative zero;
  // the loop below would leave r->size stale or flag a negative zero
  if (c == 0 || rz_is_zero(a)) {
    return rz_set_u64(r, 0);
  }

  rabin_err_t err = rz_alloc(r, a->size + 1);
  if (err != RABIN_SUCCESS) return err;
  r->size = a->size;

  u64 carry = 0;

  // In two's complement, 0 - (u64)c safely computes |c| for all negative c,
  // including INT64_MIN.
  u64 abs_c = (c < 0) ? (0 - (u64)c) : (u64)c;

  for (u64 i = 0; i < a->size; i++) {
    u128 prod = (u128)a->limbs[i] * abs_c + carry;
    r->limbs[i] = (u64)prod;
    carry = (u64)(prod >> 64);
  }

  if (carry) {
    r->limbs[r->size++] = carry;
  }

  bool neg = (c < 0);
  r->is_neg = a->is_neg ^ neg;
  return RABIN_SUCCESS;
}

bool rz_is_square(rz_t* q, const rz_t* a)
{
  if (a == NULL) return false;

  // 0 is a square; the limb tests below would read limbs[0] of a size-0
  // zero
  if (rz_is_zero(a)) {
    return rz_set_u64(q, 0) == RABIN_SUCCESS;
  }

  // 1. [Test 64]
  u64 t = a->limbs[0] & 63;

  if (q64[t] == 0) return false;

  u64 r = 0;
  for (int i = a->size - 1; i >= 0; --i)
    r = (r * 16 + (a->limbs[i] % 45045)) % 45045;

  // 2. [Test 63]
  if (q63[r % 63] == 0) return false;
  // 3. [Test 65]
  if (q65[r % 65] == 0) return false;
  // 4. [Test 11]
  if (q11[r % 11] == 0) return false;

  rz_t q_1, q_2;
  rz_init_multi(&q_1, &q_2, NULL);

  rabin_err_t err = rz_isqrt(&q_1, a);
  if (err != RABIN_SUCCESS) {
    rz_clear_multi(&q_1, &q_2, NULL);
    return false;
  }
  err = rz_sqr(&q_2, &q_1);
  if (err != RABIN_SUCCESS) {
    rz_clear_multi(&q_1, &q_2, NULL);
    return false;
  }
  if (rz_cmp(&q_2, a) != 0) {
    rz_clear_multi(&q_1, &q_2, NULL);
    return false;
  }

  if (q) {
    rz_copy(q, &q_1);
  }

  rz_clear_multi(&q_1, &q_2, NULL);
  return true;
}

bool rz_is_prime_power(rz_t* p, const rz_t* n)
{
  if (n == NULL) return false;

  rz_t a, b, p_0, temp, r;
  rz_init_multi(&a, &b, &p_0, &temp, &r, NULL);
  if (rz_set_u64(&a, 1) != RABIN_SUCCESS) {
    rz_clear_multi(&a, &b, &p_0, &temp, &r, NULL);
    return false;
  }

  bool verdict = false;

  // 2. [Compute GCD]
  while (true) {
    if (rz_add_u64(&a, &a, 1) != RABIN_SUCCESS) break;
    if (rz_mod_exp(&b, &a, n, n) != RABIN_SUCCESS) break;

    rz_sub(&temp, &b, &a);
    temp.is_neg = false;

    if (rz_gcd(&p_0, &temp, n) != RABIN_SUCCESS) break;

    // 3. [Finished?]
    if (rz_is_eq_i64(&p_0, 1)) {
      verdict = false;
      break;
    }

    if (!rz_bpsw(&p_0)) {
      continue;
    }

    // 4. [Final Test]
    if (rz_copy(&temp, n) != RABIN_SUCCESS) break;
    // repeatadly divide n by p
    while (true) {
      if (rz_divmod(&b, &r, &temp, &p_0) != RABIN_SUCCESS) break;
      if (!rz_is_zero(&r)) {
        break;
      }
      if (rz_copy(&temp, &b) != RABIN_SUCCESS) break;
    }

    if (rz_is_eq_i64(&temp, 1)) {
      verdict = true;
      if (p) {
        rz_copy(p, &p_0);
      }
    } else {
      verdict = false;
    }

    break;
  }

  rz_clear_multi(&a, &b, &p_0, &temp, &r, NULL);
  return verdict;
}
