/**
 * biglogic.c
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

#include "../../include/biglogic.h"

#include "../../include/bignum.h"

#if defined(__AVX512F__)
#include <immintrin.h>
#endif

static void bn_and_scalar(u64* res, const u64* a, const u64* b, u64 size);

static void bn_or_scalar(u64* res, const u64* a, const u64* b, u64 size);

static void bn_xor_scalar(u64* res, const u64* a, const u64* b, u64 size);

static void bn_not_scalar(u64* res, const u64* a, u64 size);

#if defined(__AVX512F__)

static void bn_and_avx512(u64* res, const u64* a, const u64* b, u64 size);

static void bn_or_avx512(u64* res, const u64* a, const u64* b, u64 size);

static void bn_xor_avx512(u64* res, const u64* a, const u64* b, u64 size);

static void bn_not_avx512(u64* res, const u64* a, u64 size);

static inline int bn_has_avx512(void)
{
  return __builtin_cpu_supports("avx512f") != 0;
}

#else

static inline int bn_has_avx512(void) { return 0; }

#endif

bool bn_supports_avx512(void) { return bn_has_avx512() != 0; }

void bn_and(bignum* r, const bignum* a, const bignum* m)
{
  // aliasing
  if (r == a || r == m) {
    bignum tmp;
    bn_init(&tmp);
    bn_and(&tmp, a, m);
    bn_copy(r, &tmp);
    bn_free(&tmp);
    return;
  }

  bn_copy(r, a);

  u64 min = MIN(r->size, m->size);

  if (bn_has_avx512()) {
#if defined(__AVX512F__)
    bn_and_avx512(r->limbs, a->limbs, m->limbs, min);
#else
    bn_and_scalar(r->limbs, a->limbs, m->limbs, min);
#endif
  } else {
    bn_and_scalar(r->limbs, a->limbs, m->limbs, min);
  }

  /* Zero the remaining limbs. */
  for (u64 i = min; i < r->size; i++) {
    r->limbs[i] = 0;
  }

  bn_trim(r);
}

void bn_or(bignum* r, const bignum* a, const bignum* m)
{
  // aliasing
  if (r == a || r == m) {
    bignum tmp;
    bn_init(&tmp);
    bn_or(&tmp, a, m);
    bn_copy(r, &tmp);
    bn_free(&tmp);
    return;
  }

  const bignum* larger;
  const bignum* smaller;

  if (a->size > m->size) {
    larger = a;
    smaller = m;
  } else {
    larger = m;
    smaller = a;
  }

  bn_copy(r, larger);

  if (bn_has_avx512()) {
#if defined(__AVX512F__)
    bn_or_avx512(r->limbs, r->limbs, smaller->limbs, smaller->size);
#else
    bn_or_scalar(r->limbs, r->limbs, smaller->limbs, smaller->size);
#endif
  } else {
    bn_or_scalar(r->limbs, r->limbs, smaller->limbs, smaller->size);
  }

  bn_trim(r);
}

void bn_xor(bignum* r, const bignum* a, const bignum* m)
{
  // aliasing
  if (r == a || r == m) {
    bignum tmp;
    bn_init(&tmp);
    bn_xor(&tmp, a, m);
    bn_copy(r, &tmp);
    bn_free(&tmp);
    return;
  }

  const bignum* larger;
  const bignum* smaller;

  if (a->size > m->size) {
    larger = a;
    smaller = m;
  } else {
    larger = m;
    smaller = a;
  }

  bn_copy(r, larger);

  if (bn_has_avx512()) {
#if defined(__AVX512F__)
    bn_xor_avx512(r->limbs, r->limbs, smaller->limbs, smaller->size);
#else
    bn_xor_scalar(r->limbs, r->limbs, smaller->limbs, smaller->size);
#endif
  } else {
    bn_xor_scalar(r->limbs, r->limbs, smaller->limbs, smaller->size);
  }

  bn_trim(r);
}

void bn_not(bignum* r, const bignum* a)
{
  if (a->size == 0) {
    bn_alloc(r, 0);
    return;
  }

  bn_copy(r, a);

  if (bn_has_avx512()) {
#if defined(__AVX512F__)
    bn_not_avx512(r->limbs, r->limbs, r->size);
#else
    bn_not_scalar(r->limbs, r->limbs, r->size);
#endif
  } else {
    bn_not_scalar(r->limbs, r->limbs, r->size);
  }

  bn_trim(r);
}

static void bn_and_scalar(u64* res, const u64* a, const u64* b, u64 size)
{
  for (u64 i = 0; i < size; i++) {
    res[i] = a[i] & b[i];
  }
}

static void bn_or_scalar(u64* res, const u64* a, const u64* b, u64 size)
{
  for (u64 i = 0; i < size; i++) {
    res[i] = a[i] | b[i];
  }
}

static void bn_xor_scalar(u64* res, const u64* a, const u64* b, u64 size)
{
  for (u64 i = 0; i < size; i++) {
    res[i] = a[i] ^ b[i];
  }
}

static void bn_not_scalar(u64* res, const u64* a, u64 size)
{
  for (u64 i = 0; i < size; i++) {
    res[i] = ~a[i];
  }
}

#if defined(__AVX512F__)

static void bn_and_avx512(u64* res, const u64* a, const u64* b, u64 size)
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

static void bn_or_avx512(u64* res, const u64* a, const u64* b, u64 size)
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

static void bn_xor_avx512(u64* res, const u64* a, const u64* b, u64 size)
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

static void bn_not_avx512(u64* res, const u64* a, u64 size)
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