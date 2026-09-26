/*
 * bignum.c
 *
 * Core bignum type management.
 *
 * This file implements the fundamental operations on the bignum type:
 * initialization, allocation, copying, freeing, trimming, comparison,
 * bit access, and decimal string conversion (parsing and printing).
 *
 * A bignum is a dynamic array of 64-bit limbs in base 2^64, little
 * endian (limb 0 is least significant), with a size, a capacity, and a
 * sign flag.
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

#include "../../include/bignum.h"

#include <ctype.h>
#include <fcntl.h>
#include <stdarg.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

#define BASE_10_19 10000000000000000000ULL

void bn_init(bignum* r)
{
  r->limbs = NULL;
  r->size = 0;
  r->capacity = 0;
  r->is_neg = false;
}

void bn_init_multi(bignum* r, ...)
{
  if (r == NULL) return;

  bn_init(r);

  va_list arg;
  va_start(arg, r);

  bignum* next;
  while ((next = va_arg(arg, bignum*)) != NULL) {
    bn_init(next);
  }

  va_end(arg);
}

void bn_swap(bignum* a, bignum* b)
{
  if (a == b) return;

  bignum temp = *a;
  *a = *b;
  *b = temp;
}

int bn_is_even(const bignum* n)
{
  if (n->size == 0 || n->limbs == NULL) {
    return 1;
  }

  return (n->limbs[0] & 1) == 0;
}

bool bn_is_zero(const bignum* a)
{
  return (a->size == 0 || (a->size == 1 && a->limbs[0] == 0));
}

bool bn_is_eq_i64(const bignum* n, i64 a)
{
  if (n->size != 1) return false;
  return n->limbs[0] == (uint64_t)a;
}

bool bn_is_one(const bignum* a)
{
  // size == 1 guards the limb access; is_neg is rejected so that -1
  // (single limb 1, sign bit set) does not count as one
  return !a->is_neg && a->size == 1 && a->limbs[0] == 1;
}

bool bn_alloc(bignum* r, u64 capacity)
{
  if (capacity <= r->capacity) return true;

  // Grow memory exponentially
  u64 new_cap = r->capacity * 2;
  if (new_cap < capacity) new_cap = capacity;

  u64* new_limbs = realloc(r->limbs, new_cap * sizeof(u64));
  if (!new_limbs) return false;

  memset(new_limbs + r->capacity, 0, (new_cap - r->capacity) * sizeof(u64));

  r->limbs = new_limbs;
  r->capacity = new_cap;
  return true;
}

void bn_init_val(bignum* n, const char* str)
{
  bn_free(n);
  if (!str) return;

  // Start with a value of 0
  bn_set_u64(n, 0);

  n->is_neg = (str[0] == '-');
  const char* s = n->is_neg ? str + 1 : str;

  for (size_t i = 0; s[i]; i++) {
    if (!isdigit(s[i])) continue;

    int digit = s[i] - '0';

    // go digits by digit propagating carry
    u64 carry = digit;
    for (u64 j = 0; j < n->size; j++) {
      unsigned __int128 prod = (unsigned __int128)n->limbs[j] * 10 + carry;
      n->limbs[j] = (u64)prod;
      carry = (u64)(prod >> 64);
    }

    // If there's still a carry after the last limb, grow the bignum
    while (carry) {
      bn_alloc(n, n->size + 1);
      n->limbs[n->size++] = carry % 0xFFFFFFFFFFFFFFFFULL;
      n->limbs[n->size - 1] = carry;
      carry = 0;
    }
  }
  bn_trim(n);
}

char* bn_to_string(const bignum* n)
{
  if (n->size == 0 || (n->size == 1 && n->limbs[0] == 0)) {
    return strdup("0");
  }

  bignum tmp, q;
  bn_init(&tmp);
  bn_init(&q);
  bn_copy(&tmp, n);

  uint64_t* parts = NULL;
  size_t parts_count = 0;
  uint64_t base_10_19 = 10000000000000000000ULL;  // 10^19

  // Extract 19-digit chunks safely
  while (!(tmp.size == 1 && tmp.limbs[0] == 0)) {
    uint64_t rem = bn_divmod_u64(&q, &tmp, base_10_19);

    parts = realloc(parts, (parts_count + 1) * sizeof(uint64_t));
    if (!parts) {
      bn_free(&tmp);
      bn_free(&q);
      return NULL;
    }
    parts[parts_count++] = rem;

    // Use bn_copy instead of raw struct assignment to prevent pointer aliasing
    bn_copy(&tmp, &q);
  }

  bn_free(&q);
  bn_free(&tmp);

  if (parts_count == 0) {
    free(parts);
    return strdup("0");
  }

  u64 buffer_size = (parts_count * 19) + 2;
  char* result = (char*)malloc(buffer_size);
  if (!result) {
    free(parts);
    return NULL;
  }

  char* ptr = result;

  if (n->is_neg) {
    *ptr++ = '-';
  }

  // 1. Print the most significant chunk (no leading zeros)
  ptr += sprintf(ptr, "%llu", (unsigned long long)parts[parts_count - 1]);

  // 2. Print all other chunks (MUST have 19 digits, pad with zeros)
  for (int64_t i = (int64_t)parts_count - 2; i >= 0; i--) {
    ptr += sprintf(ptr, "%019llu", (unsigned long long)parts[i]);
  }

  free(parts);
  return result;
}

void bn_print(const bignum* n)
{
  if (n->size == 0 || (n->size == 1 && n->limbs[0] == 0)) {
    printf("0");
    return;
  }

  if (n->is_neg) printf("-");

  bignum tmp;
  bn_init(&tmp);
  bn_copy(&tmp, n);

  uint64_t* parts = NULL;
  size_t parts_count = 0;

  // Extract 19-digit chunks
  while (!(tmp.size == 1 && tmp.limbs[0] == 0)) {
    bignum q;
    bn_init(&q);
    // Ensure bn_divmod_u64 is correctly updating 'q' and returning 'rem'
    uint64_t rem = bn_divmod_u64(&q, &tmp, BASE_10_19);

    parts = realloc(parts, (parts_count + 1) * sizeof(uint64_t));
    parts[parts_count++] = rem;

    bn_free(&tmp);
    tmp = q;
  }

  // 1. Print the most significant chunk (no leading zeros)
  printf("%llu", (unsigned long long)parts[parts_count - 1]);

  // 2. Print all other chunks (MUST have 19 digits, pad with zeros)
  for (int64_t i = (int64_t)parts_count - 2; i >= 0; i--) {
    printf("%019llu", (unsigned long long)parts[i]);
  }

  free(parts);
  bn_free(&tmp);
}

void bn_println(const bignum* n)
{
  if (n->size == 0 || (n->size == 1 && n->limbs[0] == 0)) {
    printf("0\n");
    return;
  }
  bn_print(n);

  printf("\n");
}

void bn_free(bignum* r)
{
  if (!r) {
    return;
  }

  if (r->limbs) {
    free(r->limbs);
    r->limbs = NULL;
  }

  r->size = 0;
  r->capacity = 0;
}

void bn_free_multi(bignum* r, ...)
{
  if (r == NULL) return;

  va_list arg;
  va_start(arg, r);

  bignum* next;
  while ((next = va_arg(arg, bignum*)) != NULL) {
    bn_free(next);
  }

  va_end(arg);

  bn_free(r);
}

void bn_trim(bignum* r)
{
  while (r->size > 1 && r->limbs[r->size - 1] == 0) {
    r->size--;
  }
}

int bn_get_bit(const bignum* a, int i)
{
  u64 limb = i / 64;
  u64 offset = i % 64;

  if (limb >= a->size) return 0;

  return (a->limbs[limb] >> offset) & 1;
}

int bn_cmp(const bignum* a, const bignum* b)
{
  if (!a->is_neg && b->is_neg) return 1;

  if (a->is_neg && !b->is_neg) return -1;

  int cmp = 0;

  if (a->size != b->size)
    cmp = (a->size < b->size) ? -1 : 1;
  else {
    for (i64 i = a->size - 1; i >= 0; i--) {
      if (a->limbs[i] != b->limbs[i]) {
        cmp = (a->limbs[i] < b->limbs[i]) ? -1 : 1;
        break;
      }
    }
  }
  if (a->is_neg && b->is_neg) {
    cmp = -cmp;
  }

  return cmp;
}

int bn_cmp_abs(const bignum* a, const bignum* b)
{
  if (a->size > b->size) return 1;
  if (a->size < b->size) return -1;

  for (i64 i = a->size - 1; i >= 0; i--) {
    if (a->limbs[i] > b->limbs[i]) return 1;
    if (a->limbs[i] < b->limbs[i]) return -1;
  }
  return 0;
}

void bn_copy(bignum* r, const bignum* a)
{
  if (r == a) return;

  if (a->size == 0) {
    r->size = 0;
    return;
  }

  bn_alloc(r, a->size);
  r->is_neg = a->is_neg;
  r->size = a->size;
  memcpy(r->limbs, a->limbs, r->size * sizeof(u64));
}

void bn_set_u64(bignum* n, uint64_t val)
{
  bn_free(n);  // free any existing limbs
  bn_alloc(n, 1);
  n->limbs[0] = val;
  n->size = 1;
  n->is_neg = false;
}

void bn_set_i64(bignum* n, int64_t val)
{
  if (val >= 0) {
    bn_set_u64(n, (u64)val);
  } else {
    // In two's complement, 0 - (u64)val safely computes |val| for all negative
    // values, including INT64_MIN (9223372036854775808U).
    u64 abs_val = 0 - (u64)val;
    bn_set_u64(n, abs_val);
    n->is_neg = true;
  }
}

void bn_set_bit(bignum* a, int i)
{
  u64 limb = i / 64;
  u64 offset = i % 64;

  if (limb >= a->size) {
    bn_alloc(a, limb + 1);
    a->size = limb + 1;
  }

  a->limbs[limb] |= ((u64)1 << offset);
}

void bn_clear_bit(bignum* a, int i)
{
  u64 limb = i / 64;
  u64 offset = i % 64;

  if (limb < a->size) {
    a->limbs[limb] &= ~((u64)1 << offset);
  }
}

int bn_bit_length(const bignum* a)
{
  if (a->size == 0) return 0;

  u64 top = a->limbs[a->size - 1];
  int bits = 0;

  while (top) {
    top >>= 1;
    bits++;
  }

  return (a->size - 1) * 64 + bits;
}

/**
 * @brief Count the trailing zero bits of a u64 (64 for a zero word).
 *
 * @param[in] val 64-bit value to count trailing zeros of.
 *
 * @return The number of trailing zero bits (64 if val is 0).
 */
static inline u64 count_trailing_zeros_u64(u64 val)
{
  if (val == 0) return 64;
  return (u64)__builtin_ctzll(val);
}

u64 bn_cnt_trailing_zeros(const bignum* a)
{
  if (bn_is_zero(a)) return 0;

  u64 zeros = 0;
  u64 i = 0;

  while (i < a->size && a->limbs[i] == 0) {
    zeros += 64;
    i++;
  }

  if (i == a->size) {
    return 0;
  }

  zeros += count_trailing_zeros_u64(a->limbs[i]);

  return zeros;
}
