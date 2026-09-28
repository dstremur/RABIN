/**
 * rz_logic.c
 *
 * Logic functions
 *
 * Implements the standard logic functions, such as and, or, xor, ...
 *
 * Copyright (C) 2026 Diego Strebel
 *
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

#include "../../include/rabin.h"
#include "../../include/rzlogic.h"
#include "rz_internal.h"

#if defined(__AVX512F__)
#include <immintrin.h>
#endif

static void rz_and_scalar(u64* res, const u64* a, const u64* b, u64 size);

static void rz_or_scalar(u64* res, const u64* a, const u64* b, u64 size);

static void rz_xor_scalar(u64* res, const u64* a, const u64* b, u64 size);

static void rz_not_scalar(u64* res, const u64* a, u64 size);

#if defined(__AVX512F__)

static void rz_and_avx512(u64* res, const u64* a, const u64* b, u64 size);

static void rz_or_avx512(u64* res, const u64* a, const u64* b, u64 size);

static void rz_xor_avx512(u64* res, const u64* a, const u64* b, u64 size);

static void rz_not_avx512(u64* res, const u64* a, u64 size);

static inline int rz_has_avx512(void)
{
  return __builtin_cpu_supports("avx512f") != 0;
}

#else

static inline int rz_has_avx512(void) { return 0; }

#endif

bool rz_supports_avx512(void) { return rz_has_avx512() != 0; }

rabin_err_t rz_and(rz_t* r, const rz_t* a, const rz_t* m)
{
  if (r == NULL || a == NULL || m == NULL) return RABIN_ERR_NULL_PTR;

  // aliasing
  if (r == a || r == m) {
    rz_t tmp;
    rz_init(&tmp);
    rabin_err_t err = rz_and(&tmp, a, m);
    if (err != RABIN_SUCCESS) {
      rz_clear(&tmp);
      return err;
    }
    err = rz_copy(r, &tmp);
    rz_clear(&tmp);
    return err;
  }

  rabin_err_t err = rz_copy(r, a);
  if (err != RABIN_SUCCESS) return err;

  u64 min = RZ_MIN(r->size, m->size);

  if (rz_has_avx512()) {
#if defined(__AVX512F__)
    rz_and_avx512(r->limbs, a->limbs, m->limbs, min);
#else
    rz_and_scalar(r->limbs, a->limbs, m->limbs, min);
#endif
  } else {
    rz_and_scalar(r->limbs, a->limbs, m->limbs, min);
  }

  /* Zero the remaining limbs. */
  for (u64 i = min; i < r->size; i++) {
    r->limbs[i] = 0;
  }

  rz_trim(r);
  return RABIN_SUCCESS;
}

rabin_err_t rz_or(rz_t* r, const rz_t* a, const rz_t* m)
{
  if (r == NULL || a == NULL || m == NULL) return RABIN_ERR_NULL_PTR;

  // aliasing
  if (r == a || r == m) {
    rz_t tmp;
    rz_init(&tmp);
    rabin_err_t err = rz_or(&tmp, a, m);
    if (err != RABIN_SUCCESS) {
      rz_clear(&tmp);
      return err;
    }
    err = rz_copy(r, &tmp);
    rz_clear(&tmp);
    return err;
  }

  const rz_t* larger;
  const rz_t* smaller;

  if (a->size > m->size) {
    larger = a;
    smaller = m;
  } else {
    larger = m;
    smaller = a;
  }

  rabin_err_t err = rz_copy(r, larger);
  if (err != RABIN_SUCCESS) return err;

  if (rz_has_avx512()) {
#if defined(__AVX512F__)
    rz_or_avx512(r->limbs, r->limbs, smaller->limbs, smaller->size);
#else
    rz_or_scalar(r->limbs, r->limbs, smaller->limbs, smaller->size);
#endif
  } else {
    rz_or_scalar(r->limbs, r->limbs, smaller->limbs, smaller->size);
  }

  rz_trim(r);
  return RABIN_SUCCESS;
}

rabin_err_t rz_xor(rz_t* r, const rz_t* a, const rz_t* m)
{
  if (r == NULL || a == NULL || m == NULL) return RABIN_ERR_NULL_PTR;

  // aliasing
  if (r == a || r == m) {
    rz_t tmp;
    rz_init(&tmp);
    rabin_err_t err = rz_xor(&tmp, a, m);
    if (err != RABIN_SUCCESS) {
      rz_clear(&tmp);
      return err;
    }
    err = rz_copy(r, &tmp);
    rz_clear(&tmp);
    return err;
  }

  const rz_t* larger;
  const rz_t* smaller;

  if (a->size > m->size) {
    larger = a;
    smaller = m;
  } else {
    larger = m;
    smaller = a;
  }

  rabin_err_t err = rz_copy(r, larger);
  if (err != RABIN_SUCCESS) return err;

  if (rz_has_avx512()) {
#if defined(__AVX512F__)
    rz_xor_avx512(r->limbs, r->limbs, smaller->limbs, smaller->size);
#else
    rz_xor_scalar(r->limbs, r->limbs, smaller->limbs, smaller->size);
#endif
  } else {
    rz_xor_scalar(r->limbs, r->limbs, smaller->limbs, smaller->size);
  }

  rz_trim(r);
  return RABIN_SUCCESS;
}

rabin_err_t rz_not(rz_t* r, const rz_t* a)
{
  if (r == NULL || a == NULL) return RABIN_ERR_NULL_PTR;

  if (a->size == 0) {
    return rz_alloc(r, 0);
  }

  rabin_err_t err = rz_copy(r, a);
  if (err != RABIN_SUCCESS) return err;

  if (rz_has_avx512()) {
#if defined(__AVX512F__)
    rz_not_avx512(r->limbs, r->limbs, r->size);
#else
    rz_not_scalar(r->limbs, r->limbs, r->size);
#endif
  } else {
    rz_not_scalar(r->limbs, r->limbs, r->size);
  }

  rz_trim(r);
  return RABIN_SUCCESS;
}

static void rz_and_scalar(u64* res, const u64* a, const u64* b, u64 size)
{
  for (u64 i = 0; i < size; i++) {
    res[i] = a[i] & b[i];
  }
}

static void rz_or_scalar(u64* res, const u64* a, const u64* b, u64 size)
{
  for (u64 i = 0; i < size; i++) {
    res[i] = a[i] | b[i];
  }
}

static void rz_xor_scalar(u64* res, const u64* a, const u64* b, u64 size)
{
  for (u64 i = 0; i < size; i++) {
    res[i] = a[i] ^ b[i];
  }
}

static void rz_not_scalar(u64* res, const u64* a, u64 size)
{
  for (u64 i = 0; i < size; i++) {
    res[i] = ~a[i];
  }
}

#if defined(__AVX512F__)

static void rz_and_avx512(u64* res, const u64* a, const u64* b, u64 size)
{
  u64 i = 0;

  for (; i + 8 <= size; i += 8) {
    __m512i vec_a = _mm512_loadu_si512((const void*)&a[i]);

    __m512i vec_b = _mm512_loadu_si512((const void*)&b[i]);

    __m512i vec_res = _mm512_and_si512(vec_a, vec_b);

    _mm512_storeu_si512((void*)&res[i], vec_res);
  }

  for (; i < size; i++) {
    res[i] = a[i] & b[i];
  }
}

static void rz_or_avx512(u64* res, const u64* a, const u64* b, u64 size)
{
  u64 i = 0;

  for (; i + 8 <= size; i += 8) {
    __m512i vec_a = _mm512_loadu_si512((const void*)&a[i]);

    __m512i vec_b = _mm512_loadu_si512((const void*)&b[i]);

    __m512i vec_res = _mm512_or_si512(vec_a, vec_b);

    _mm512_storeu_si512((void*)&res[i], vec_res);
  }

  for (; i < size; i++) {
    res[i] = a[i] | b[i];
  }
}

static void rz_xor_avx512(u64* res, const u64* a, const u64* b, u64 size)
{
  u64 i = 0;

  for (; i + 8 <= size; i += 8) {
    __m512i vec_a = _mm512_loadu_si512((const void*)&a[i]);

    __m512i vec_b = _mm512_loadu_si512((const void*)&b[i]);

    __m512i vec_res = _mm512_xor_si512(vec_a, vec_b);

    _mm512_storeu_si512((void*)&res[i], vec_res);
  }

  for (; i < size; i++) {
    res[i] = a[i] ^ b[i];
  }
}

static void rz_not_avx512(u64* res, const u64* a, u64 size)
{
  u64 i = 0;

  const __m512i ones = _mm512_set1_epi64(-1LL);

  for (; i + 8 <= size; i += 8) {
    __m512i vec_a = _mm512_loadu_si512((const void*)&a[i]);

    __m512i vec_res = _mm512_xor_si512(vec_a, ones);

    _mm512_storeu_si512((void*)&res[i], vec_res);
  }

  for (; i < size; i++) {
    res[i] = ~a[i];
  }
}

#endif