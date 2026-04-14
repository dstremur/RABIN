#include <stddef.h>
#include <string.h>

#include "../include/bignum.h"

u64 estimate_iterations(u64 bits)
{
  if (bits <= 64) return 2;

  u64 log2 = 64 - __builtin_clzll(bits);
  u64 its = (log2 * 2) / 3;

  return its + 1;
}

void bn_newton_div(bignum* q, const bignum* a, const bignum* d)
{
  if (bn_is_zero(d)) return;
  if (bn_cmp(a, d) < 0) {
    bn_set_u64(q, 0);
    return;
  }

  // estimate the reciprocal of d using newton iterations
  // find roots of f(x) = 1/x - d => y = x + x * (1 - x * d)

  bignum x, two_p, tmp, c1, c2, error;
  bn_init_multi(&x, &two_p, &tmp, &c1, &c2, &error, NULL);

  // choose initial estimate
  u64 a_bits = bn_bit_length(a);
  u64 d_bits = bn_bit_length(d);
  // tune if precision is lacking
  u64 P = a_bits + 32;

  bn_set_u64(&x, 1);
  bn_lshift(&x, &x, P - d_bits);

  u64 iters = estimate_iterations(d_bits);

  for (u64 i = 0; i < iters; i++) {
    // knuth: x_n+1 = x_n (1 + (1 - v*x_n)(1+ (1 - v*x_n)))
    // x = x + (x * (2^(P) - d * x)) >> P

    // tmp = d * x
    bn_mul(&tmp, &x, d);

    // two_p = 2^{P}
    bn_set_u64(&two_p, 1);
    bn_lshift(&two_p, &two_p, P);

    // Calculate (2^P - d*x) and update x
    if (bn_cmp(&two_p, &tmp) >= 0) {
      // x is to small or perfect
      bn_sub(&error, &two_p, &tmp);

      // c1 = (x * error) >> P
      bn_mul(&c1, &x, &error);
      bn_rshift(&c1, &c1, P);

      // c2 = (c1 * error) >> P
      bn_mul(&c2, &c1, &error);
      bn_rshift(&c2, &c2, P);

      // x = x + c1 + c2
      bn_add(&x, &x, &c1);
      bn_add(&x, &x, &c2);
    } else {
      // x was slightly too large, subtract correction instead
      bn_sub(&error, &tmp, &two_p);

      // c1 = (x * error) >> P
      bn_mul(&c1, &x, &error);
      bn_rshift(&c1, &c1, P);

      // c2 = (c1 * error) >> P
      bn_mul(&c2, &c1, &error);
      bn_rshift(&c2, &c2, P);

      bn_sub(&x, &x, &c1);
      bn_add(&x, &x, &c2);
    }

    if (bn_is_zero(&error)) break;
  }

  bn_mul(q, a, &x);
  bn_rshift(q, q, P);

  bignum r;
  bn_init(&r);
  bn_mul(&tmp, q, d);
  bn_sub(&r, a, &tmp);

  while (bn_cmp(&r, d) >= 0) {
    bn_add_u64(q, q, 1);
    bn_sub(&r, &r, d);
  }

  bn_free_multi(&x, &tmp, &two_p, &r, &error, &c1, &c2, NULL);
}

static void bn_div_single_limb(bignum* q, const bignum* a, u64 b_limb)
{
  if (b_limb == 0) return;  // Divide by zero error

  bn_alloc(q, a->size);
  q->size = a->size;

  u64 remainder = 0;
  for (i64 i = a->size - 1; i >= 0; i--) {
    __uint128_t temp = ((__uint128_t)remainder << 64) | a->limbs[i];
    q->limbs[i] = (u64)(temp / b_limb);
    remainder = (u64)(temp % b_limb);
  }
  bn_trim(q);
}

void bn_div(bignum* q, const bignum* a, const bignum* b)
{
  if (b->size == 0 || (b->size == 1 && b->limbs[0] == 0)) return;

  if (q == a || q == b) {
    bignum tmp;
    bn_init(&tmp);
    bn_div(&tmp, a, b);
    bn_copy(q, &tmp);
    bn_free(&tmp);
    return;
  }

  if (bn_cmp(a, b) < 0) {
    bn_set_u64(q, 0);
    return;
  }

  if (b->size == 1) {
    bn_div_single_limb(q, a, b->limbs[0]);
    q->is_neg = a->is_neg ^ b->is_neg;
    return;
  }

  u64 n = b->size;
  u64 m = a->size - n;

  bignum u, v;
  bn_init_multi(&u, &v, NULL);

  // Normalize
  u64 s = __builtin_clzll(b->limbs[n - 1]);

  bn_lshift(&v, b, s);
  bn_lshift(&u, a, s);

  if (u.size == a->size) {
    bn_alloc(&u, u.size + 1);
    u.limbs[u.size] = 0;
    u.size++;
  }

  bn_alloc(q, m + 1);
  q->size = m + 1;
  memset(q->limbs, 0, q->size * sizeof(u64));
  q->is_neg = a->is_neg ^ b->is_neg;

  u64 vn1 = v.limbs[n - 1];
  u64 vn2 = v.limbs[n - 2];

  // main loop
  for (i64 j = m; j >= 0; j--) {
    u64 q_hat = 0;
    u64 r_hat = 0;
    u64 uj_n = u.limbs[j + n];
    u64 uj_n1 = u.limbs[j + n - 1];

    if (uj_n == vn1) {
      q_hat = ~0ULL;
      // Knuth: if (uj_n1 + vn1) < base, r_hat = uj_n1 + vn1 and check
      // refinement.
      // Otherwise, r_hat >= base, so we stop refinement immediately.
      if (uj_n1 <= (~0ULL - vn1)) {
        r_hat = uj_n1 + vn1;
        // Refine q_hat using the 3rd limb
        while ((__uint128_t)q_hat * vn2 >
               (((__uint128_t)r_hat << 64) | u.limbs[j + n - 2])) {
          q_hat--;
          r_hat += vn1;
          if (r_hat < vn1)
            break;  // r_hat >= base (overflowed u64), stop refinement
        }
      }
    } else {
      __uint128_t temp = ((__uint128_t)uj_n << 64) | uj_n1;
      q_hat = (u64)(temp / vn1);
      r_hat = (u64)(temp % vn1);
    }

    // Refine q_hat using the 3rd limb
    while ((__uint128_t)q_hat * vn2 >
           (((__uint128_t)r_hat << 64) | u.limbs[j + n - 2])) {
      q_hat--;
      r_hat += vn1;
      if (r_hat < vn1)
        break;  // r_hat >= base (overflowed u64), stop refinement
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

    q->limbs[j] = q_hat;

    // add back

    if (is_neg) {
      q->limbs[j]--;
      u64 carry = 0;
      for (u64 i = 0; i < n; i++) {
        __uint128_t sum = (__uint128_t)u.limbs[j + i] + v.limbs[i] + carry;
        u.limbs[j + i] = (u64)sum;
        carry = sum >> 64;
      }
      u.limbs[j + n] += carry;
    }
  }

  bn_trim(q);
  bn_free(&u);
  bn_free(&v);
}

// algorithm D knuth
void bn_div_long(bignum* q, const bignum* a, const bignum* b)
{
  if (b->size == 0 || (b->size == 1 && b->limbs[0] == 0)) {
    return;
  }
  // aliasing
  if (q == a || q == b) {
    bignum tmp;
    bn_init(&tmp);
    bn_div(&tmp, a, b);
    bn_copy(q, &tmp);
    bn_free(&tmp);
    return;
  }

  bignum r;
  bn_init(&r);

  bn_alloc(q, a->size);
  memset(q->limbs, 0, q->size * sizeof(u64));
  q->size = a->size;
  q->is_neg = a->is_neg ^ b->is_neg;

  i64 nbits = bn_bit_length(a);

  // binary long division
  for (i64 i = nbits - 1; i >= 0; i--) {
    bn_lshift1(&r);
    if (bn_get_bit(a, i)) {
      bn_set_bit(&r, 0);
    }

    if (bn_cmp(&r, b) >= 0) {
      bn_sub(&r, &r, b);
      bn_set_bit(q, i);
    }
  }
  bn_trim(q);
  bn_free(&r);
}
