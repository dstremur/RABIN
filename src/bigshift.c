#include <stddef.h>
#include <string.h>

#include "../include/bignum.h"

void bn_lshift1(bignum* r)
{
  if (r->size == 0) return;

  u64 carry = 0;

  for (u64 i = 0; i < r->size; i++) {
    u64 new_carry = r->limbs[i] >> 63;
    r->limbs[i] = (r->limbs[i] << 1) | carry;
    carry = new_carry;
  }

  if (carry) {
    bn_alloc(r, r->size + 1);
    r->limbs[r->size++] = carry;
  }
}

void bn_rshift1(bignum* r)
{
  if (r->size == 0) return;
  u64 carry = 0;
  for (i64 i = r->size - 1; i >= 0; i--) {
    u64 next = (r->limbs[i] & 1) << 63;
    r->limbs[i] = (r->limbs[i] >> 1) | carry;
    carry = next;
  }

  bn_trim(r);
}

void bn_lshift(bignum* r, const bignum* a, int shift)
{
  // aliasing
  if (r == a) {
    bignum tmp;
    bn_init(&tmp);

    bn_lshift(&tmp, a, shift);
    bn_copy(r, &tmp);

    bn_free(&tmp);
    return;
  }

  // early exit
  if (shift == 0) {
    bn_copy(r, a);
    return;
  }

  u64 words = shift / 64;
  u64 bits = shift % 64;

  u64 max_size = a->size + words + 1;
  // if alloc fails just return 0
  if (!bn_alloc(r, max_size)) {
    bn_set_u64(r, 0);
    return;
  }

  // zero the limbs
  memset(r->limbs, 0, max_size * sizeof(u64));

  r->size = max_size;
  r->is_neg = a->is_neg;

  u64 carry = 0;

  // shift left limb by limb
  for (u64 i = 0; i < a->size; i++) {
    u64 src = a->limbs[i];
    u64 dst = i + words;

    // use 128 bit sum, maybe optimize with assembly
    unsigned __int128 sum = (unsigned __int128)src << bits;
    sum += carry;

    // top 64 bits
    r->limbs[dst] = (u64)sum;
    // lower 64 bits go to the carry
    carry = (u64)(sum >> 64);
  }

  // handle remaining carry
  if (carry) {
    r->limbs[a->size + words] = carry;
    r->size = a->size + words + 1;
  } else {
    r->size = a->size + words;
  }

  bn_trim(r);
}

void bn_rshift(bignum* r, const bignum* a, int shift)
{
  // aliasing
  if (r == a) {
    bignum tmp;
    bn_init(&tmp);

    bn_rshift(&tmp, a, shift);
    bn_copy(r, &tmp);

    bn_free(&tmp);
    return;
  }

  // early exit
  if (shift == 0) {
    bn_copy(r, a);
    return;
  }

  u64 words = shift / 64;
  u64 bits = shift % 64;

  // we shift more than the number of limbs
  if (words >= a->size) {
    bn_set_u64(r, 0);
    return;
  }

  u64 new_size = a->size - words;

  if (!bn_alloc(r, new_size)) {
    bn_set_u64(r, 0);
    return;
  }

  u64* new_limbs = r->limbs;

  if (bits == 0) {
    for (u64 i = 0; i < new_size; i++) {
      new_limbs[i] = a->limbs[i + words];
    }
  } else {
    for (u64 i = 0; i < new_size; i++) {
      u64 curr = a->limbs[i + words];
      u64 next = (i + words + 1 < a->size) ? a->limbs[i + words + 1] : 0;

      r->limbs[i] = (curr >> bits) | (next << (64 - bits));
    }
  }

  r->size = new_size;

  if (r->capacity > r->size) {
    memset(r->limbs + r->size, 0, (r->capacity - r->size) * sizeof(u64));
  }
  bn_trim(r);
}

// Shift r left by 1 and add 0 or 1
void bn_lshift1_add(bignum* r, int bit)
{
  u64 carry = (bit != 0);
  for (u64 i = 0; i < r->size; i++) {
    u64 tmp = r->limbs[i] >> 63;
    r->limbs[i] = (r->limbs[i] << 1) | carry;
    carry = tmp;
  }

  if (carry) {
    bn_alloc(r, r->size + 1);
    r->limbs[r->size++] = carry;
  }
}
