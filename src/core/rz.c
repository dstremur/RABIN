/*
 * rz.c
 *
 * Core rz_t type management.
 *
 * This file implements the fundamental operations on the rz_t type:
 * initialization, allocation, copying, freeing, trimming, comparison,
 * bit access, and decimal string conversion (parsing and printing).
 *
 * A rz_t is a dynamic array of 64-bit limbs in base 2^64, little
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

#include <ctype.h>
#include <fcntl.h>
#include <stdarg.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

#include "../../include/rabin.h"

#define RZ_BASE_10_19 10000000000000000000ULL

rabin_err_t rz_init(rz_t* r)
{
  if (r == NULL) return RABIN_ERR_NULL_PTR;

  r->limbs = NULL;
  r->size = 0;
  r->capacity = 0;
  r->is_neg = false;
  return RABIN_SUCCESS;
}

rabin_err_t rz_init_multi(rz_t* r, ...)
{
  if (r == NULL) return RABIN_ERR_NULL_PTR;

  rabin_err_t err = rz_init(r);
  if (err != RABIN_SUCCESS) return err;

  va_list arg;
  va_start(arg, r);

  rz_t* next;
  while ((next = va_arg(arg, rz_t*)) != NULL) {
    err = rz_init(next);
    if (err != RABIN_SUCCESS) {
      va_end(arg);
      return err;
    }
  }

  va_end(arg);
  return RABIN_SUCCESS;
}

rabin_err_t rz_swap(rz_t* a, rz_t* b)
{
  if (a == NULL || b == NULL) return RABIN_ERR_NULL_PTR;
  if (a == b) return RABIN_SUCCESS;

  rz_t temp = *a;
  *a = *b;
  *b = temp;
  return RABIN_SUCCESS;
}

int rz_is_even(const rz_t* n)
{
  if (n == NULL) return 1;
  if (n->size == 0 || n->limbs == NULL) {
    return 1;
  }

  return (n->limbs[0] & 1) == 0;
}

bool rz_is_zero(const rz_t* a)
{
  if (a == NULL) return false;
  return (a->size == 0 || (a->size == 1 && a->limbs[0] == 0));
}

bool rz_is_eq_i64(const rz_t* n, i64 a)
{
  if (n == NULL) return false;
  if (n->size != 1) return false;
  return n->limbs[0] == (uint64_t)a;
}

bool rz_is_one(const rz_t* a)
{
  if (a == NULL) return false;
  // size == 1 guards the limb access; is_neg is rejected so that -1
  // (single limb 1, sign bit set) does not count as one
  return !a->is_neg && a->size == 1 && a->limbs[0] == 1;
}

#include "rz_internal.h"

rabin_err_t rz_alloc(rz_t* r, u64 capacity)
{
  if (r == NULL) return RABIN_ERR_NULL_PTR;
  if (capacity <= r->capacity) return RABIN_SUCCESS;

  // Grow memory exponentially
  u64 new_cap = r->capacity * 2;
  if (new_cap < capacity) new_cap = capacity;

  if (new_cap > SIZE_MAX / sizeof(u64)) return RABIN_ERR_OUT_OF_MEMORY;

  u64* new_limbs = realloc(r->limbs, new_cap * sizeof(u64));
  if (!new_limbs) return RABIN_ERR_OUT_OF_MEMORY;

  memset(new_limbs + r->capacity, 0, (new_cap - r->capacity) * sizeof(u64));

  r->limbs = new_limbs;
  r->capacity = new_cap;
  return RABIN_SUCCESS;
}

rabin_err_t rz_init_val(rz_t* n, const char* str)
{
  if (n == NULL || str == NULL) return RABIN_ERR_NULL_PTR;

  // Start with a value of 0
  rabin_err_t err = rz_set_u64(n, 0);
  if (err != RABIN_SUCCESS) return err;

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

    // If there's still a carry after the last limb, grow the rz_t
    while (carry) {
      err = rz_alloc(n, n->size + 1);
      if (err != RABIN_SUCCESS) return err;
      n->limbs[n->size++] = carry % 0xFFFFFFFFFFFFFFFFULL;
      n->limbs[n->size - 1] = carry;
      carry = 0;
    }
  }
  rz_trim(n);
  return RABIN_SUCCESS;
}

char* rz_to_string(const rz_t* n)
{
  if (n == NULL) return NULL;

  if (n->size == 0 || (n->size == 1 && n->limbs[0] == 0)) {
    return strdup("0");
  }

  rz_t tmp, q;
  rz_init(&tmp);
  rz_init(&q);
  if (rz_copy(&tmp, n) != RABIN_SUCCESS) {
    rz_clear(&tmp);
    rz_clear(&q);
    return NULL;
  }

  uint64_t* parts = NULL;
  size_t parts_count = 0;
  uint64_t base_10_19 = 10000000000000000000ULL;  // 10^19

  // Extract 19-digit chunks safely
  while (!(tmp.size == 1 && tmp.limbs[0] == 0)) {
    uint64_t rem = rz_divmod_u64(&q, &tmp, base_10_19);

    uint64_t* grown = realloc(parts, (parts_count + 1) * sizeof(uint64_t));
    if (!grown) {
      rz_clear(&tmp);
      rz_clear(&q);
      free(parts);
      return NULL;
    }
    parts = grown;
    parts[parts_count++] = rem;

    // Use rz_copy instead of raw struct assignment to prevent pointer aliasing
    rz_copy(&tmp, &q);
  }

  rz_clear(&q);
  rz_clear(&tmp);

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

rabin_err_t rz_print(const rz_t* n)
{
  if (n == NULL) return RABIN_ERR_NULL_PTR;

  if (n->size == 0 || (n->size == 1 && n->limbs[0] == 0)) {
    printf("0");
    return RABIN_SUCCESS;
  }

  if (n->is_neg) printf("-");

  rz_t tmp;
  rz_init(&tmp);
  if (rz_copy(&tmp, n) != RABIN_SUCCESS) {
    rz_clear(&tmp);
    return RABIN_ERR_OUT_OF_MEMORY;
  }

  uint64_t* parts = NULL;
  size_t parts_count = 0;

  // Extract 19-digit chunks
  while (!(tmp.size == 1 && tmp.limbs[0] == 0)) {
    rz_t q;
    rz_init(&q);
    // Ensure rz_divmod_u64 is correctly updating 'q' and returning 'rem'
    uint64_t rem = rz_divmod_u64(&q, &tmp, RZ_BASE_10_19);

    uint64_t* grown = realloc(parts, (parts_count + 1) * sizeof(uint64_t));
    if (!grown) {
      rz_clear(&q);
      rz_clear(&tmp);
      free(parts);
      return RABIN_ERR_OUT_OF_MEMORY;
    }
    parts = grown;
    parts[parts_count++] = rem;

    rz_clear(&tmp);
    tmp = q;
  }

  // 1. Print the most significant chunk (no leading zeros)
  printf("%llu", (unsigned long long)parts[parts_count - 1]);

  // 2. Print all other chunks (MUST have 19 digits, pad with zeros)
  for (int64_t i = (int64_t)parts_count - 2; i >= 0; i--) {
    printf("%019llu", (unsigned long long)parts[i]);
  }

  free(parts);
  rz_clear(&tmp);
  return RABIN_SUCCESS;
}

rabin_err_t rz_println(const rz_t* n)
{
  if (n == NULL) return RABIN_ERR_NULL_PTR;

  if (n->size == 0 || (n->size == 1 && n->limbs[0] == 0)) {
    printf("0\n");
    return RABIN_SUCCESS;
  }

  rabin_err_t err = rz_print(n);
  if (err != RABIN_SUCCESS) return err;

  printf("\n");
  return RABIN_SUCCESS;
}

rabin_err_t rz_clear(rz_t* r)
{
  if (r == NULL) return RABIN_ERR_NULL_PTR;

  if (r->limbs) {
    free(r->limbs);
    r->limbs = NULL;
  }

  r->size = 0;
  r->capacity = 0;
  return RABIN_SUCCESS;
}

rabin_err_t rz_clear_multi(rz_t* r, ...)
{
  if (r == NULL) return RABIN_ERR_NULL_PTR;

  va_list arg;
  va_start(arg, r);

  rz_t* next;
  while ((next = va_arg(arg, rz_t*)) != NULL) {
    rabin_err_t err = rz_clear(next);
    if (err != RABIN_SUCCESS) {
      va_end(arg);
      return err;
    }
  }

  va_end(arg);

  return rz_clear(r);
}

rabin_err_t rz_trim(rz_t* r)
{
  if (r == NULL) return RABIN_ERR_NULL_PTR;

  while (r->size > 1 && r->limbs[r->size - 1] == 0) {
    r->size--;
  }
  return RABIN_SUCCESS;
}

int rz_get_bit(const rz_t* a, int i)
{
  if (a == NULL || i < 0) return 0;

  u64 limb = i / 64;
  u64 offset = i % 64;

  if (limb >= a->size) return 0;

  return (a->limbs[limb] >> offset) & 1;
}

int rz_cmp(const rz_t* a, const rz_t* b)
{
  if (a == NULL || b == NULL) return 0;

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

int rz_cmp_abs(const rz_t* a, const rz_t* b)
{
  if (a == NULL || b == NULL) return 0;
  if (a->size > b->size) return 1;
  if (a->size < b->size) return -1;

  for (i64 i = a->size - 1; i >= 0; i--) {
    if (a->limbs[i] > b->limbs[i]) return 1;
    if (a->limbs[i] < b->limbs[i]) return -1;
  }
  return 0;
}

rabin_err_t rz_copy(rz_t* r, const rz_t* a)
{
  if (r == NULL || a == NULL) return RABIN_ERR_NULL_PTR;
  if (r == a) return RABIN_SUCCESS;

  if (a->size == 0) {
    r->size = 0;
    return RABIN_SUCCESS;
  }

  rabin_err_t err = rz_alloc(r, a->size);
  if (err != RABIN_SUCCESS) return err;

  r->is_neg = a->is_neg;
  r->size = a->size;
  memcpy(r->limbs, a->limbs, r->size * sizeof(u64));
  return RABIN_SUCCESS;
}

rabin_err_t rz_set_u64(rz_t* n, uint64_t val)
{
  if (n == NULL) return RABIN_ERR_NULL_PTR;

  rabin_err_t err = rz_clear(n);  // free any existing limbs
  if (err != RABIN_SUCCESS) return err;

  err = rz_alloc(n, 1);
  if (err != RABIN_SUCCESS) return err;

  n->limbs[0] = val;
  n->size = 1;
  n->is_neg = false;
  return RABIN_SUCCESS;
}

rabin_err_t rz_set_i64(rz_t* n, int64_t val)
{
  if (n == NULL) return RABIN_ERR_NULL_PTR;

  rabin_err_t err;
  if (val >= 0) {
    err = rz_set_u64(n, (u64)val);
  } else {
    // In two's complement, 0 - (u64)val safely computes |val| for all negative
    // values, including INT64_MIN (9223372036854775808U).
    u64 abs_val = 0 - (u64)val;
    err = rz_set_u64(n, abs_val);
    if (err != RABIN_SUCCESS) return err;
    n->is_neg = true;
  }
  return err;
}

rabin_err_t rz_set_bit(rz_t* a, int i)
{
  if (a == NULL) return RABIN_ERR_NULL_PTR;
  if (i < 0) return RABIN_ERR_INVALID_ARG;

  u64 limb = i / 64;
  u64 offset = i % 64;

  if (limb >= a->size) {
    rabin_err_t err = rz_alloc(a, limb + 1);
    if (err != RABIN_SUCCESS) return err;
    a->size = limb + 1;
  }

  a->limbs[limb] |= ((u64)1 << offset);
  return RABIN_SUCCESS;
}

rabin_err_t rz_clear_bit(rz_t* a, int i)
{
  if (a == NULL) return RABIN_ERR_NULL_PTR;
  if (i < 0) return RABIN_ERR_INVALID_ARG;

  u64 limb = i / 64;
  u64 offset = i % 64;

  if (limb < a->size) {
    a->limbs[limb] &= ~((u64)1 << offset);
  }
  return RABIN_SUCCESS;
}

int rz_bit_length(const rz_t* a)
{
  if (a == NULL) return 0;
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

u64 rz_cnt_trailing_zeros(const rz_t* a)
{
  if (a == NULL || rz_is_zero(a)) return 0;

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
