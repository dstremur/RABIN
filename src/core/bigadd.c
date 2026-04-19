#include <ctype.h>
#include <immintrin.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>

#include "../../include/bignum.h"

extern uint64_t bn_add_inner(uint64_t* r, const uint64_t* a, uint64_t a_size,
                             const uint64_t* b, uint64_t b_size);

// TODO: sign handling
void bn_add(bignum* r, const bignum* a, const bignum* b)
{
  // aliasing
  if (r == a || r == b) {
    bignum tmp;
    bn_init(&tmp);
    bn_add(&tmp, a, b);
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
  /*
    u64 a_size = a->size;
    u64 b_size = b->size;

    char carry = 0;
    u64 i = 0;

    for (; i < b_size; i++) {
      carry = _addcarry_u64(carry, a->limbs[i], b->limbs[i],
                            (unsigned long long*)&r->limbs[i]);
    }

    for (; i < a_size; i++) {
      if (carry == 0) {
        if (r != a) {
          memcpy(&r->limbs[i], &a->limbs[i], (a_size - i) * sizeof(u64));
        }
        break;
      }
      carry =
          _addcarry_u64(carry, a->limbs[i], 0, (unsigned long
    long*)&r->limbs[i]);
    }
  */

  u64 carry = bn_add_inner(r->limbs, a->limbs, a->size, b->limbs, b->size);

  r->limbs[a->size] = carry;
  r->size = a->size + carry;
}

void bn_add_u64(bignum* r, const bignum* a, u64 b)
{
  if (b == 0) {
    bn_copy(r, a);
    return;
  }

  bn_copy(r, a);

  u64 carry = b;

  for (u64 i = 0; i < r->size && carry > 0; i++) {
    u64 old = r->limbs[i];
    r->limbs[i] += carry;

    if (r->limbs[i] < old) {
      carry = 1;
    } else {
      carry = 0;
    }
  }

  if (carry) {
    if (bn_alloc(r, r->size + 1)) {
      r->limbs[r->size] = carry;
      r->size++;
    }
  }
}

void bn_add_at_offset(bignum* r, const bignum* a, u64 offset)
{
  if (a->size == 0) return;

  u64 carry = 0;
  u64 i = 0;

  // Use a single loop that continues as long as there is
  // either data in 'a' OR a carry to propagate.
  for (i = 0; (i < a->size || carry > 0); i++) {
    u64 r_idx = offset + i;

    // If we exceed r's capacity, we must stop to avoid a crash
    if (r_idx >= r->size) break;

    u64 av = (i < a->size) ? a->limbs[i] : 0;
    u64 rv = r->limbs[r_idx];

    unsigned __int128 sum = (unsigned __int128)av + rv + carry;

    r->limbs[r_idx] = (u64)sum;
    carry = (u64)(sum >> 64);
  }
}
