#include <ctype.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>

#include "../include/bignum.h"

void bn_mod(bignum* r, const bignum* a, const bignum* b)
{
  if (bn_cmp(a, b) < 0) {
    bn_copy(r, a);
    return;
  }
  if (r == a || r == b) {
    bignum tmp;
    bn_init(&tmp);
    bn_mod(&tmp, a, b);
    bn_copy(r, &tmp);
    bn_free(&tmp);
    return;
  }

  bn_set_u64(r, 0);

  int nbits = bn_bit_length(a);

  for (int i = nbits - 1; i >= 0; i--) {
    bn_lshift1_add(r, bn_get_bit(a, i));

    if (bn_cmp(r, b) >= 0) {
      bn_sub_abs(r, r, b);  // safe subtraction
    }
  }

  bn_trim(r);
}

// Fast modular inverse for a single 64-bit limb (Newton's method)
uint64_t mod_inverse_u64(uint64_t n)
{
  uint64_t inv = 1;
  for (int i = 0; i < 6; i++) {
    inv *= (2 - n * inv);
  }
  return -inv;
}

uint64_t bn_divmod_u64(bignum* q, const bignum* a, uint64_t d)
{
  if (q == a) {
    bignum tmp;
    bn_init(&tmp);
    uint64_t rem = bn_divmod_u64(&tmp, a, d);
    bn_copy(q, &tmp);
    bn_free(&tmp);
    return rem;
  }
  bn_alloc(q, a->size);
  q->size = a->size;

  unsigned __int128 rem = 0;

  for (i64 i = a->size - 1; i >= 0; i--) {
    unsigned __int128 cur = (rem << 64) | a->limbs[i];

    q->limbs[i] = (u64)(cur / d);
    rem = cur % d;
  }

  bn_trim(q);
  return (uint64_t)rem;
}
