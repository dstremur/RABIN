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

#include <stddef.h>
#include <stdio.h>
#include <string.h>

#include "../../include/precomp.h"
#include "../../include/rabin.h"
#include "../../include/rzlimb.h"
#include "../../include/u64.h"
#include "rz_internal.h"

/* Exact single-prime NTT multiplication of magnitudes.
 * a, b must be nonzero; r must not alias a or b (dispatch handles it).
 * Does NOT set r->is_neg (caller does). */
static rabin_err_t rz_mul_ntt_u64(rz_t* r, const rz_t* a, const rz_t* b)
{
  u64 na, nb, nz, ntt_size, k, i;
  const u64_ntt_ctx_t* ctx;
  u64 *buf, *a_arr, *b_arr, *a_hat, *b_hat, *r_arr;
  u64 carry;
  rabin_err_t err;

  if (a->size > SIZE_MAX / 4 || b->size > SIZE_MAX / 4)
    return RABIN_ERR_OVERFLOW;

  na = a->size * 4; /* 16-bit chunks */
  nb = b->size * 4;
  if (na > UINT64_MAX - (nb - 1)) return RABIN_ERR_OVERFLOW; /* nz overflow */
  nz = na + nb - 1; /* true product length in chunks */

  ntt_size = 1;
  k = 0;
  while (k < 54 && ntt_size < nz) {
    ntt_size <<= 1;
    k++;
  }

  /* Exactness + transform-size guards: RABIN_ERR_OVERFLOW signals the
     caller to fall back to Karatsuba. */
  if (RZ_MIN(na, nb) > RZ_NTT_U64_MAX_CHUNKS) return RABIN_ERR_OVERFLOW;
  if (ntt_size < nz) return RABIN_ERR_OVERFLOW; /* transform size k > 54 */
  if (ntt_size > SIZE_MAX / (5 * sizeof(u64))) return RABIN_ERR_OVERFLOW;

  ctx = u64_ntt_ctx_golden_cached(k);
  if (ctx == NULL) return RABIN_ERR_OUT_OF_MEMORY;

  buf = rz_scratch_get(5 * ntt_size);
  if (buf == NULL) return RABIN_ERR_OUT_OF_MEMORY;
  memset(buf, 0, 5 * ntt_size * sizeof(u64)); /* scratch is NOT zeroed */
  a_arr = buf;
  b_arr = buf + ntt_size;
  a_hat = buf + 2 * ntt_size;
  b_hat = buf + 3 * ntt_size;
  r_arr = buf + 4 * ntt_size;

  /* Unpack: 4 x 16-bit chunks per 64-bit limb (rest stays zero). */
  for (i = 0; i < a->size; i++) {
    u64 L = a->limbs[i];
    a_arr[4 * i + 0] = L & 0xFFFFu;
    a_arr[4 * i + 1] = (L >> 16) & 0xFFFFu;
    a_arr[4 * i + 2] = (L >> 32) & 0xFFFFu;
    a_arr[4 * i + 3] = L >> 48;
  }
  for (i = 0; i < b->size; i++) {
    u64 L = b->limbs[i];
    b_arr[4 * i + 0] = L & 0xFFFFu;
    b_arr[4 * i + 1] = (L >> 16) & 0xFFFFu;
    b_arr[4 * i + 2] = (L >> 32) & 0xFFFFu;
    b_arr[4 * i + 3] = L >> 48;
  }

  /* The proven u64-NTT convolution chain. */
  u64_ntt_cyclic_forward(a_hat, a_arr, ctx);
  u64_ntt_cyclic_forward(b_hat, b_arr, ctx);
  for (i = 0; i < ntt_size; i++)
    a_hat[i] = u64_mont_mul(a_hat[i], b_hat[i], &ctx->mctx);
  u64_ntt_cyclic_inverse_montgomery_in(r_arr, a_hat, ctx);

  /* Carry propagation, base 2^16.
     Since a*b < 2^(16*na + 16*nb), every carried chunk at index >= nz + 1
     is zero; the tail chunks are kept locally (tail[4] bounds the loose
     per-step carry analysis, which allows at most indices nz .. nz+3). */
  carry = 0;
  u64 last = 0;
  u64 tail[4] = {0, 0, 0, 0}; /* chunk values at indices nz .. nz+3 */
  i = 0;
  while (i < nz || carry) {
    u64 total = (i < nz) ? r_arr[i] + carry : carry;
    u64 low = total & 0xFFFFu;
    carry = total >> 16;
    if (i < nz)
      r_arr[i] = low;
    else
      tail[i - nz] = low;
    if (low || carry) last = i; /* last nonzero chunk (a,b != 0 => last >= 0) */
    i++;
  }

  /* Pack: 4 x 16-bit chunks per limb, exactly (last+4)/4 limbs.
     Since a*b < 2^(16*(na+nb)), n_limbs <= (na+nb)/4 = a->size + b->size
     <= 2*max_len, so this never overruns r. */
  u64 n_limbs = (last + 4) / 4;
  for (i = 0; i < n_limbs; i++) {
    u64 L = 0;
    for (u64 j = 0; j < 4; j++) {
      u64 ci = 4 * i + j;
      u64 c = (ci < nz) ? r_arr[ci] : tail[ci - nz];
      L |= c << (16 * j);
    }
    r->limbs[i] = L;
  }
  r->size = n_limbs;

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

  if (a->size >= RZ_NTT_LIMIT) {
    rabin_err_t nerr = rz_mul_ntt_u64(r, a, a);
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

  // use the u64-NTT path for large, reasonably balanced operands;
  // lopsided pairs are cheaper with Karatsuba/schoolbook
  if (max_len >= RZ_NTT_LIMIT && RZ_MIN(a->size, b->size) >= max_len / 4) {
    rabin_err_t nerr = rz_mul_ntt_u64(r, a, b);
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

  err = rz_mul_ntt_u64(res, a, b);
  if (err != RABIN_SUCCESS) return err;

  res->is_neg = a->is_neg ^ b->is_neg;
  rz_trim(res);
  return RABIN_SUCCESS;
}

rabin_err_t rz_mul_i64(rz_t* r, const rz_t* a, const i64 c)
{
  if (r == NULL || a == NULL) return RABIN_ERR_NULL_PTR;

  u64 carry = 0;
  r->size = a->size;

  rabin_err_t err = rz_alloc(r, a->size + 1);
  if (err != RABIN_SUCCESS) return err;

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

  // 1. [Test 64]
  u64 t = a->limbs[0] & 63;

  if (q64[t] == 0) return false;

  u64 r = 0;
  for (int i = a->size - 1; i >= 0; --i)
    r = (r * 16 + (a->limbs[i] % 45045)) % 45045;

  // 2. [Test 63
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