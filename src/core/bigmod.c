#include <ctype.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>

#include "../../include/bighelper.h"
#include "../../include/bignum.h"

void bn_mod1(bignum* r, const bignum* a, const bignum* b)
{
  if (bn_cmp_abs(a, b) < 0) {
    bn_copy(r, a);
    return;
  }
  // aliasing
  if (r == a || r == b) {
    bignum tmp;
    bn_init(&tmp);
    bn_mod1(&tmp, a, b);
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
      bn_sub_abs(r, r, b);
    }
  }

  if (r->is_neg) {
    bn_add_abs(r, r, b);
  }

  bn_trim(r);
}

void bn_mod2(bignum* r, const bignum* a, const bignum* b)
{
  if (b->size == 0 || (b->size == 1 && b->limbs[0] == 0)) return;

  if (r == a || r == b) {
    bignum tmp;
    bn_init(&tmp);
    bn_mod(&tmp, a, b);
    bn_copy(r, &tmp);
    bn_free(&tmp);
    return;
  }

  if (bn_cmp(a, b) < 0) {
    bn_copy(r, a);

    if (r->is_neg) {
      r->is_neg = false;
      bn_add_abs(r, r, b);
    }
    return;
  }

  if (b->size == 1) {
    u64 rem = bn_mod_u64(a, b->limbs[0]);

    bn_set_u64(r, rem);

    // euclidian division
    if (a->is_neg && rem != 0) {
      bignum tmp;
      bn_init(&tmp);

      bn_copy(&tmp, b);
      tmp.is_neg = 0;

      bn_sub_abs(r, &tmp, r);

      bn_free(&tmp);
    }

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

    // estimate q_hat
    if (uj_n == vn1) {
      q_hat = ~0ULL;
      r_hat = uj_n1;
    } else {
      __uint128_t numerator = ((__uint128_t)uj_n << 64) | uj_n1;

      q_hat = (u64)(numerator / vn1);
      r_hat = (u64)(numerator % vn1);
    }

    // refine q_hat
    while ((__uint128_t)q_hat * vn2 > (((__uint128_t)r_hat << 64) + uj_n2)) {
      q_hat--;
      r_hat += vn1;
      if (r_hat < vn1) break;
    }

    // 6. Multiply and Subtract: u = u - q_hat * v
    u64 borrow = limbs_submul_1(u.limbs + j, v.limbs, n, q_hat);

    // Final borrow check against the "extra" limb
    u64 top = u.limbs[j + n];
    u.limbs[j + n] = top - borrow;
    unsigned char borrow_sub = (top < borrow);

    // 7. Add Back (if q_hat was 1 too large)
    if (borrow_sub) {
      u64 carry = limbs_add_n(u.limbs + j, v.limbs, n);
      u.limbs[j + n] += carry;
    }
  }

  u.size = n;
  bn_rshift(r, &u, s);
  r->is_neg = false;

  if (a->is_neg && !bn_is_zero(r)) {
    bn_sub_abs(r, b, r);
    r->is_neg = false;
  }

  bn_trim(r);
  bn_free(&u);
  bn_free(&v);
}

uint64_t bn_mod_u64(const bignum* a, uint64_t d)
{
  unsigned __int128 rem = 0;

  // Purely mathematical reduction with zero memory allocation
  for (i64 i = a->size - 1; i >= 0; i--) {
    unsigned __int128 cur = (rem << 64) | a->limbs[i];
    rem = cur % d;
  }

  u64 res = (u64)rem;

  if (a->is_neg && res != 0) {
    return d - res;
  }

  return res;
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

bool bn_mod_inverse(bignum* res, const bignum* a, const bignum* m)
{
  // If a is 0 or modulus is <= 1, no inverse exists
  if (bn_is_zero(a) || bn_is_zero(m) || bn_cmp(m, &BN_ONE) == 0) {
    return false;
  }

  bignum u, v, x1, x2;
  bn_init(&u);
  bn_init(&v);
  bn_init(&x1);
  bn_init(&x2);

  bn_copy(&u, a);
  bn_copy(&v, m);
  bn_set_u64(&x1, 1);
  bn_set_u64(&x2, 0);

  while (!bn_is_zero(&u) && !bn_is_zero(&v)) {
    // Eliminate powers of 2 in u
    while (bn_is_even(&u)) {
      bn_rshift1(&u);  // u = u / 2
      if (bn_is_even(&x1)) {
        bn_rshift1(&x1);
      } else {
        bn_add(&x1, &x1, m);
        bn_rshift1(&x1);  // x1 = (x1 + m) / 2
      }
    }

    // Eliminate powers of 2 in v
    while (bn_is_even(&v)) {
      bn_rshift1(&v);  // v = v / 2
      if (bn_is_even(&x2)) {
        bn_rshift1(&x2);
      } else {
        bn_add(&x2, &x2, m);
        bn_rshift1(&x2);  // x2 = (x2 + m) / 2
      }
    }

    // Step-down subtraction
    if (bn_cmp(&u, &v) >= 0) {
      bn_sub(&u, &u, &v);
      // Simulating signed subtraction under unsigned bignum bounds
      if (bn_cmp(&x1, &x2) < 0) {
        bn_add(&x1, &x1, m);
      }
      bn_sub(&x1, &x1, &x2);
    } else {
      bn_sub(&v, &v, &u);
      if (bn_cmp(&x2, &x1) < 0) {
        bn_add(&x2, &x2, m);
      }
      bn_sub(&x2, &x2, &x1);
    }
  }

  bool success = false;
  // If gcd is 1, the matching variable holds the modular inverse
  if (bn_cmp(&u, &BN_ONE) == 0) {
    bn_copy(res, &x1);
    success = true;
  } else if (bn_cmp(&v, &BN_ONE) == 0) {
    bn_copy(res, &x2);
    success = true;
  }

  bn_free(&u);
  bn_free(&v);
  bn_free(&x1);
  bn_free(&x2);

  return success;
}
