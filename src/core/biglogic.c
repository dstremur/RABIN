/*
 * biglogic.c
 *
 * Logic functions
 *
 * Implements the standard logic functions, such as and, or, xor, ...
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

#include "../../include/biglogic.h"

#include "../../include/bignum.h"
#include "immintrin.h"

static void bn_and_avx(u64* res, const u64* a, const u64* b, const u64 size);
static void bn_or_avx(u64* res, const u64* a, const u64* b, const u64 size);
static void bn_xor_avx(u64* res, const u64* a, const u64* b, const u64 size);
static void bn_not_avx(u64* res, const u64* a, const u64 size);

void bn_and(bignum* r, const bignum* a, const bignum* m)
{
  bn_copy(r, a);

  u64 min = MIN(r->size, m->size);

  bn_and_avx(r->limbs, a->limbs, m->limbs, min);

  // zero the rest
  for (u64 i = min; i < r->size; i++) {
    r->limbs[i] = 0;
  }

  bn_trim(r);
}

void bn_or(bignum* r, const bignum* a, const bignum* m)
{
  const bignum* larger = (a->size > m->size) ? a : m;
  const bignum* smaller = (a->size > m->size) ? m : a;

  bn_copy(r, larger);

  bn_or_avx(r->limbs, r->limbs, smaller->limbs, smaller->size);

  bn_trim(r);
}

void bn_xor(bignum* r, const bignum* a, const bignum* m)
{
  const bignum* larger = (a->size > m->size) ? a : m;
  const bignum* smaller = (a->size > m->size) ? m : a;

  bn_copy(r, larger);

  bn_xor_avx(r->limbs, r->limbs, smaller->limbs, smaller->size);

  bn_trim(r);
}

void bn_not(bignum* r, const bignum* a)
{
  if (a->size == 0) {
    bn_alloc(r, 0);
    return;
  }

  bn_copy(r, a);
  bn_not_avx(r->limbs, r->limbs, r->size);
  bn_trim(r);
}

void bn_and_avx(u64* res, const u64* a, const u64* b, const u64 size)
{
  u64 i = 0;

  // load 8 limbs
  for (; i + 8 <= size; i += 8) {
    __m512i vec_a = _mm512_loadu_si512((const __m512i*)&a[i]);
    __m512i vec_b = _mm512_loadu_si512((const __m512i*)&b[i]);

    __m512i vec_res = _mm512_and_si512(vec_a, vec_b);

    _mm512_storeu_si512((__m512i*)&res[i], vec_res);
  }

  for (; i < size; i++) {
    res[i] = a[i] & b[i];
  }
}

void bn_or_avx(u64* res, const u64* a, const u64* b, const u64 size)
{
  u64 i = 0;

  // load 8 limbs
  for (; i + 8 <= size; i += 8) {
    __m512i vec_a = _mm512_loadu_si512((const __m512i*)&a[i]);
    __m512i vec_b = _mm512_loadu_si512((const __m512i*)&b[i]);

    __m512i vec_res = _mm512_or_si512(vec_a, vec_b);

    _mm512_storeu_si512((__m512i*)&res[i], vec_res);
  }

  for (; i < size; i++) {
    res[i] = a[i] | b[i];
  }
}

void bn_xor_avx(u64* res, const u64* a, const u64* b, const u64 size)
{
  u64 i = 0;

  // load 8 limbs
  for (; i + 8 <= size; i += 8) {
    __m512i vec_a = _mm512_loadu_si512((const __m512i*)&a[i]);
    __m512i vec_b = _mm512_loadu_si512((const __m512i*)&b[i]);

    __m512i vec_res = _mm512_xor_si512(vec_a, vec_b);

    _mm512_storeu_si512((__m512i*)&res[i], vec_res);
  }

  for (; i < size; i++) {
    res[i] = a[i] ^ b[i];
  }
}

void bn_not_avx(u64* res, const u64* a, const u64 size)
{
  u64 i = 0;

  __m512i ones = _mm512_set1_epi64(-1LL);

  // load 8 limbs
  for (; i + 8 <= size; i += 8) {
    __m512i vec_a = _mm512_loadu_si512((const __m512i*)&a[i]);
    __m512i vec_res = _mm512_xor_si512(vec_a, ones);
    _mm512_storeu_si512((__m512i*)&res[i], vec_res);
  }

  for (; i < size; i++) {
    res[i] = ~(a[i]);
  }
}
