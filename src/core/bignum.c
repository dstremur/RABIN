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

/*
 * Initialize a bignum to zero, releasing no memory (r must not have
 * been allocated yet).
 *
 * Complexity:
 *   Time: O(1)
 *   Auxiliary memory: O(1)
 *   Output memory: O(1)
 */
void bn_init(bignum* r)
{
  r->limbs = NULL;
  r->size = 0;
  r->capacity = 0;
  r->is_neg = false;
}

/*
 * Initialize multiple bignums to zero.
 *
 * Takes a NULL-terminated list of bignum pointers; the first argument
 * is the first bignum and the variadic arguments continue the list.
 * IMPORTANT: the list must be terminated with NULL.
 *
 * Complexity:
 *   Time: O(k) for k bignums
 *   Auxiliary memory: O(1)
 *   Output memory: O(1)
 */
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

/*
 * Swap the contents of two bignums by exchanging their structs.
 *
 * Only the struct fields (pointers, sizes, sign) are swapped, so this
 * is O(1) regardless of size.
 *
 * Complexity:
 *   Time: O(1)
 *   Auxiliary memory: O(1)
 *   Output memory: O(1)
 */
void bn_swap(bignum* a, bignum* b)
{
  if (a == b) return;

  bignum temp = *a;
  *a = *b;
  *b = temp;
}

/*
 * Test whether a bignum is even.
 *
 * Returns 1 if n is even (including zero), 0 otherwise. Only the least
 * significant bit of the lowest limb is inspected.
 *
 * Complexity:
 *   Time: O(1)
 *   Auxiliary memory: O(1)
 *   Output memory: O(1)
 */
int bn_is_even(const bignum* n)
{
  if (n->size == 0 || n->limbs == NULL) {
    return 1;
  }

  return (n->limbs[0] & 1) == 0;
}

/*
 * Test whether a bignum is zero.
 *
 * Complexity:
 *   Time: O(1)
 *   Auxiliary memory: O(1)
 *   Output memory: O(1)
 */
bool bn_is_zero(const bignum* a)
{
  return (a->size == 0 || (a->size == 1 && a->limbs[0] == 0));
}

/*
 * Test whether a bignum equals a 64-bit signed integer.
 *
 * Returns true iff n has exactly one limb and that limb equals a.
 *
 * Complexity:
 *   Time: O(1)
 *   Auxiliary memory: O(1)
 *   Output memory: O(1)
 */
bool bn_is_eq_i64(const bignum* n, i64 a)
{
  if (n->size != 1) return false;
  return n->limbs[0] == (uint64_t)a;
}

/*
 * Grow the limb storage of r to at least capacity limbs.
 *
 * If r->capacity is already sufficient, nothing happens. Otherwise the
 * storage is reallocated with exponential growth (at least double the
 * old capacity, or the requested capacity, whichever is larger) and
 * the newly added limbs are zeroed. Existing limb contents are
 * preserved.
 *
 * Returns true on success, false if the realloc fails (r is left
 * unchanged).
 *
 * Complexity:
 *   Time: O(capacity) for the zero-fill and the realloc copy
 *   Auxiliary memory: O(capacity) during the realloc
 *   Output memory: O(capacity) limbs
 */
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

/*
 * Parse a base-10 decimal string into a bignum.
 *
 * Let d = number of decimal digits in str.
 *
 * Frees any previous contents of n and sets it to the value of str,
 * which may start with a '-' sign. Non-digit characters are skipped.
 * Digits are processed one at a time with a multiply-by-10-and-add
 * carry loop over the limbs.
 *
 * Complexity:
 *   Time: O(d * n) where n is the number of limbs, i.e. O(d^2 / 64)
 *   Auxiliary memory: O(1)
 *   Output memory: O(n) limbs
 */
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

/*
 * Convert a bignum to a base-10 decimal string.
 *
 * Let n = a->size, measured in 64-bit limbs.
 *
 * Returns a newly allocated string (caller must free it) containing
 * the decimal representation of n, with a leading '-' for negative
 * values. The magnitude is reduced repeatedly by 10^19 (which fits in
 * a u64), collecting 19-digit chunks; the most significant chunk is
 * printed without padding and the rest with zero padding.
 *
 * Returns NULL on allocation failure.
 *
 * Complexity:
 *   Time: O(n^2) - O(n) divisions by a 64-bit divisor, each O(n)
 *   Auxiliary memory: O(n) limbs for temporaries, O(n) for the chunk
 *                    array and the output string
 *   Output memory: O(n) characters
 */
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

/*
 * Print a bignum to stdout in base 10.
 *
 * Let n = a->size, measured in 64-bit limbs.
 *
 * Same chunking strategy as bn_to_string(): the magnitude is reduced
 * repeatedly by 10^19 and the chunks are printed most-significant
 * first, with a leading '-' for negative values.
 *
 * Complexity:
 *   Time: O(n^2) - O(n) divisions by a 64-bit divisor, each O(n)
 *   Auxiliary memory: O(n) limbs for temporaries, O(n) for the chunk
 *                    array
 *   Output memory: O(n) characters written
 */
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

/*
 * Print a bignum to stdout in base 10 followed by a newline.
 *
 * Complexity:
 *   Time: O(n^2), see bn_print()
 *   Auxiliary memory: O(n)
 *   Output memory: O(n) characters written
 */
void bn_println(const bignum* n)
{
  if (n->size == 0 || (n->size == 1 && n->limbs[0] == 0)) {
    printf("0\n");
    return;
  }
  bn_print(n);

  printf("\n");
}

/*
 * Free the limb storage of a bignum and reset it to the zero state.
 *
 * Safe to call on a NULL pointer or on an already-freed bignum.
 *
 * Complexity:
 *   Time: O(1)
 *   Auxiliary memory: O(1)
 *   Output memory: O(1)
 */
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

/*
 * Free multiple bignums.
 *
 * Takes a NULL-terminated list of bignum pointers; the first argument
 * is the first bignum and the variadic arguments continue the list.
 * IMPORTANT: the list must be terminated with NULL.
 *
 * Complexity:
 *   Time: O(k) for k bignums
 *   Auxiliary memory: O(1)
 *   Output memory: O(1)
 */
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

/*
 * Remove trailing zero limbs from a bignum.
 *
 * Reduces r->size so that the most significant limb is nonzero (or
 * size is 1). The capacity is left unchanged.
 *
 * Complexity:
 *   Time: O(number of trimmed limbs)
 *   Auxiliary memory: O(1)
 *   Output memory: O(1)
 */
void bn_trim(bignum* r)
{
  while (r->size > 1 && r->limbs[r->size - 1] == 0) {
    r->size--;
  }
}

/*
 * Return the value of the i-th bit of a (0 or 1).
 *
 * Bit 0 is the least significant bit. Bits beyond the current size
 * read as 0.
 *
 * Complexity:
 *   Time: O(1)
 *   Auxiliary memory: O(1)
 *   Output memory: O(1)
 */
int bn_get_bit(const bignum* a, int i)
{
  u64 limb = i / 64;
  u64 offset = i % 64;

  if (limb >= a->size) return 0;

  return (a->limbs[limb] >> offset) & 1;
}

/*
 * Compare two signed bignums.
 *
 * Returns 1 if a > b, -1 if a < b, 0 if a == b.
 *
 * Signs are handled first (a positive number is greater than a
 * negative one); for equal signs the magnitudes are compared, most
 * significant limb first, and the result is negated when both
 * operands are negative.
 *
 * Complexity:
 *   Time: O(n) worst case, where n = max(a->size, b->size)
 *   Auxiliary memory: O(1)
 *   Output memory: O(1)
 */
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

/*
 * Compare the absolute values (magnitudes) of two bignums.
 *
 * Returns 1 if |a| > |b|, -1 if |a| < |b|, 0 if |a| == |b|. Signs are
 * ignored. The comparison is by size first, then most significant
 * limb first.
 *
 * Complexity:
 *   Time: O(n) worst case, where n = max(a->size, b->size)
 *   Auxiliary memory: O(1)
 *   Output memory: O(1)
 */
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

/*
 * Deep copy a bignum: r = a.
 *
 * Copies the limb storage and the sign. r may alias a (a no-op in
 * that case). If a is zero-sized, r is reset to size 0.
 *
 * Complexity:
 *   Time: O(n) where n = a->size
 *   Auxiliary memory: O(1)
 *   Output memory: O(n) limbs
 */
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

/*
 * Set a bignum to an unsigned 64-bit value.
 *
 * Frees any previous contents and stores val in a single limb with a
 * positive sign.
 *
 * Complexity:
 *   Time: O(1)
 *   Auxiliary memory: O(1)
 *   Output memory: O(1) limbs
 */
void bn_set_u64(bignum* n, uint64_t val)
{
  bn_free(n);  // free any existing limbs
  bn_alloc(n, 1);
  n->limbs[0] = val;
  n->size = 1;
  n->is_neg = false;
}

/*
 * Set a bignum to a signed 64-bit value.
 *
 * For negative val the magnitude is computed as -val (with overflow
 * protection for INT64_MIN) and the sign flag is set.
 *
 * Complexity:
 *   Time: O(1)
 *   Auxiliary memory: O(1)
 *   Output memory: O(1) limbs
 */
void bn_set_i64(bignum* n, int64_t val)
{
  if (val >= 0) {
    bn_set_u64(n, (u64)val);
  } else {
    u64 abs = (u64) - (val - 1) + 1;
    bn_set_u64(n, abs);
    n->is_neg = true;
  }
}

/*
 * Set the i-th bit of a to 1.
 *
 * Bit 0 is the least significant bit. If the bit lies beyond the
 * current size, the bignum is grown (with zeroed limbs) to reach it.
 *
 * Complexity:
 *   Time: O(1) amortized (O(limb) for the zero-fill on growth)
 *   Auxiliary memory: O(1)
 *   Output memory: O(1) limbs, possibly grown
 */
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

/*
 * Clear the i-th bit of a to 0.
 *
 * Bit 0 is the least significant bit. Bits beyond the current size
 * are already 0 and are left untouched.
 *
 * Complexity:
 *   Time: O(1)
 *   Auxiliary memory: O(1)
 *   Output memory: O(1)
 */
void bn_clear_bit(bignum* a, int i)
{
  u64 limb = i / 64;
  u64 offset = i % 64;

  if (limb < a->size) {
    a->limbs[limb] &= ~((u64)1 << offset);
  }
}

/*
 * Return the bit length of a: the index of the highest set bit plus
 * one (0 for zero).
 *
 * Only the most significant limb is inspected.
 *
 * Complexity:
 *   Time: O(1)
 *   Auxiliary memory: O(1)
 *   Output memory: O(1)
 */
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

/*
 * Count the trailing zero bits of a u64 (64 for a zero word).
 */
static inline u64 count_trailing_zeros_u64(u64 val)
{
  if (val == 0) return 64;
  return (u64)__builtin_ctzll(val);
}

/*
 * Return the number of trailing zero bits of a.
 *
 * Counts whole zero limbs (64 bits each) plus the trailing zeros of
 * the first nonzero limb. Returns 0 for a == 0.
 *
 * Complexity:
 *   Time: O(number of zero limbs)
 *   Auxiliary memory: O(1)
 *   Output memory: O(1)
 */
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
