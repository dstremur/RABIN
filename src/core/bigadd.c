#include <ctype.h>
#include <immintrin.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>

#include "../../include/bignum.h"

extern uint64_t bn_add_inner(uint64_t* r, const uint64_t* a, uint64_t a_size,
                             const uint64_t* b, uint64_t b_size);

// TODO: sign handling
void bn_add_abs(bignum* r, const bignum* a, const bignum* b)
{
  // aliasing
  if (r == a || r == b) {
    bignum tmp;
    bn_init(&tmp);
    bn_add_abs(&tmp, a, b);
    bn_copy(r, &tmp);
    bn_free(&tmp);
    return;
  }

  // ensure a is always larger
  if (a->size < b->size) {
    const bignum* tmp = a;
    a = b;
    b = tmp;
  }

  if (r->capacity < a->size + 1) {
    bn_alloc(r, a->size + 1);
  }

  u64 carry = bn_add_inner(r->limbs, a->limbs, a->size, b->limbs, b->size);

  r->limbs[a->size] = carry;
  r->size = a->size + carry;
}

void bn_add(bignum* r, const bignum* a, const bignum* b)
{
  if (a->is_neg == b->is_neg) {
    bn_add_abs(r, a, b);
    r->is_neg = a->is_neg;
    return;
  }

  i64 cmp = bn_cmp_abs(a, b);

  if (cmp == 0) {
    bn_set_u64(r, 0);
    r->is_neg = false;
    return;
  }

  if (cmp > 0) {
    bn_sub_abs(r, a, b);
    r->is_neg = a->is_neg;
  } else {
    bn_sub_abs(r, b, a);
    r->is_neg = b->is_neg;
  }
}

void bn_add_u64(bignum* r, const bignum* a, u64 b)
{
  if (b == 0) {
    bn_copy(r, a);
    return;
  }

  if (!a->is_neg) {
    // POSITIVE + POSITIVE
    bn_copy(r, a);
    u64 carry = b;
    for (u64 i = 0; i < r->size; i++) {
      u64 old = r->limbs[i];
      r->limbs[i] += carry;
      if (r->limbs[i] >= old) {
        carry = 0;
        break;
      }  // Optimization: exit early
      carry = 1;
    }
    if (carry) {
      if (r->capacity < r->size + 1) bn_alloc(r, r->size + 1);
      r->limbs[r->size++] = 1;
    }
  } else {
    // NEGATIVE + POSITIVE (e.g., -80 + 16)
    // This is Magnitude Subtraction: |a| - b
    if (a->size == 1 && a->limbs[0] <= b) {
      bn_set_u64(r, b - a->limbs[0]);
      r->is_neg = false;
    } else {
      bn_copy(r, a);
      u64 borrow = b;
      for (u64 i = 0; i < r->size && borrow; i++) {
        u64 old = r->limbs[i];
        r->limbs[i] -= borrow;
        borrow = (old < borrow) ? 1 : 0;
      }
      r->is_neg = true;
      bn_trim(r);
      if (r->size == 0) r->is_neg = false;
    }
  }
}

void bn_add_at_offset(bignum* r, const bignum* a, u64 offset)
{
  if (a->size == 0) return;

  u64 carry = 0;
  // Ensure r has enough limbs to hold a + carry at this offset
  u64 required_size = offset + a->size + 1;
  if (r->size < required_size) {
    bn_alloc(r, required_size);
    r->size = required_size;
  }

  for (u64 i = 0; (i < a->size || carry > 0); i++) {
    u64 r_idx = offset + i;
    u64 av = (i < a->size) ? a->limbs[i] : 0;
    u64 rv = r->limbs[r_idx];

    unsigned __int128 sum = (unsigned __int128)av + rv + carry;
    r->limbs[r_idx] = (u64)sum;
    carry = (u64)(sum >> 64);
  }
  bn_trim(r);
}
