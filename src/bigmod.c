#include <ctype.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>

#include "../include/bignum.h"

void bn_mod2(bignum* r, const bignum* a, const bignum* b)
{
  if (bn_cmp(a, b) < 0) {
    bn_copy(r, a);
    return;
  }
  // aliasing
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
    bn_lshift1(r);
    if (bn_get_bit(a, i)) {
      bn_set_bit(r, 0);
    }

    if (bn_cmp(r, b) >= 0) {
      bn_sub(r, r, b);
    }
  }

  bn_trim(r);
}

void bn_mod(bignum* r, const bignum* a, const bignum* b)
{
  if (b->size == 0 || (b->size == 1 && b->limbs[0] == 0)) return;

  if (r == a || r == b) {
    bignum tmp;
    bn_init(&tmp);
    bn_div(&tmp, a, b);
    bn_copy(r, &tmp);
    bn_free(&tmp);
    return;
  }

  u64 n = b->size;
  u64 m = a->size - n;

  bignum u, v;
  bn_init_multi(&u, &v, NULL);

  // Normalize, d = 2^s
  u64 s = __builtin_clzll(b->limbs[n - 1]);

  // multiply by d
  bn_lshift(&v, b, s);
  bn_lshift(&u, a, s);

  // introduce new digit position
  if (u.size == a->size) {
    bn_alloc(&u, u.size + 1);
    u.limbs[u.size] = 0;
    u.size++;
  }

  u64 vn1 = v.limbs[n - 1];
  u64 vn2 = v.limbs[n - 2];

  // main loop
  for (i64 j = m; j >= 0; j--) {
    u64 q_hat = 0;
    u64 r_hat = 0;
    u64 uj_n = u.limbs[j + n];
    u64 uj_n1 = u.limbs[j + n - 1];
    u64 uj_n2 = u.limbs[j + n - 2];

    // check if quotient fits in 64 bit
    if (uj_n == vn1) {
      q_hat = ~0ULL;
      r_hat = uj_n1 + vn1;

      // refine q_hat
      while ((__uint128_t)q_hat * vn2 > (((__uint128_t)r_hat << 64) + uj_n2)) {
        q_hat--;
        r_hat += vn1;
        if (r_hat < vn1) break;
      }
    } else {
      __uint128_t temp = ((__uint128_t)uj_n << 64) + uj_n1;
      q_hat = (u64)(temp / vn1);
      r_hat = (u64)(temp % vn1);
    }

    // refine q_hat
    while ((__uint128_t)q_hat * vn2 > (((__uint128_t)r_hat << 64) + uj_n2)) {
      q_hat--;
      r_hat += vn1;
      if (r_hat < vn1) break;
    }

    // multiply and subtract
    u64 borrow = 0;
    for (u64 i = 0; i < n; i++) {
      __uint128_t p = (__uint128_t)q_hat * v.limbs[i];
      __uint128_t sub = p + borrow;

      if (u.limbs[j + i] < (u64)sub) {
        borrow = (sub >> 64) + 1;
        u.limbs[j + i] -= (u64)sub;
      } else {
        borrow = (sub >> 64);
        u.limbs[j + i] -= (u64)sub;
      }
    }

    bool is_neg = u.limbs[j + n] < borrow;
    u.limbs[j + n] -= borrow;

    // add back
    if (__builtin_expect(is_neg, 0)) {
      u64 carry = 0;
      for (u64 i = 0; i < n; i++) {
        __uint128_t sum = (__uint128_t)u.limbs[j + i] + v.limbs[i] + carry;
        u.limbs[j + i] = (u64)sum;
        carry = sum >> 64;
      }
      u.limbs[j + n] += carry;
    }
  }

  bn_rshift(r, &u, s);
  r->size = n;
  bn_trim(r);

  bn_free(&u);
  bn_free(&v);
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

uint64_t bn_mod_u64(const bignum* a, uint64_t d)
{
  unsigned __int128 rem = 0;

  // Purely mathematical reduction with zero memory allocation
  for (i64 i = a->size - 1; i >= 0; i--) {
    unsigned __int128 cur = (rem << 64) | a->limbs[i];
    rem = cur % d;
  }

  return (uint64_t)rem;
}

uint64_t bn_divmod_u64(bignum* q, const bignum* a, uint64_t d)
{
  // aliasing
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
  return (u64)rem;
}
