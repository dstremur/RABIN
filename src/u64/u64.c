#include "../include/u64.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../../include/bigrns.h"

typedef unsigned __int128 u128;

inline u64 mod_add(u64 a, u64 b, u64 p)
{
  u64 res = a + b;

  u64 mask = -(u64)((res >= p) | (res < a));
  return res - (p & mask);
}

u64 mod_sub(u64 a, u64 b, u64 p)
{
  u64 res = a - b;
  // If a < b, a borrow occurred.
  // (a < b) evaluates to 1 or 0. Negating it creates a mask of all 1s or all
  // 0s.
  u64 mask = -(u64)(a < b);
  return res + (p & mask);
}

inline u64 mod_mul(u64 a, u64 b, u64 p)
{
  unsigned __int128 res = (unsigned __int128)a * b;
  return (u64)(res % p);
}

u64 mod_pow(u64 base, u64 exp, u64 p)
{
  u64 res = 1;
  base %= p;
  while (exp > 0) {
    if (exp % 2 == 1) res = mod_mul(res, base, p);
    base = mod_mul(base, base, p);
    exp /= 2;
  }
  return res;
}

// p needs to be prime
u64 mod_inverse_euclid(u64 a, u64 p)
{
  if (a == 0) return 0;  // Should not happen with primes

  __int128 t = 0, newt = 1;
  __int128 r = p, newr = a;

  while (newr != 0) {
    __int128 q = r / newr;

    __int128 tmp_t = newt;
    newt = t - q * newt;
    t = tmp_t;

    __int128 tmp_r = newr;
    newr = r - q * newr;
    r = tmp_r;
  }

  if (t < 0) t += p;

  return (u64)t;
}

// p needs to be prime
u64 mod_inverse(u64 n, u64 p) { return mod_pow(n, p - 2, p); }

u64 compute_mu(u64 q)
{
  u128 dividend = ~((u128)0);
  return (u64)(dividend / q);
}

u64 barrett_reduction(u128 c, u64 q, u64 mu)
{
  unsigned __int128 q_est = (unsigned __int128)((c * mu) >> 128);

  u64 r = (u64)(c - q_est * q);

  while (r >= q) {
    r -= q;
  }

  return r;
}

// q = 2^64 - 2^32 + 1
u64 goldilock_red(u128 c)
{
  u128 X_3 = c >> 96;
  u128 X_2 = (c >> 64) & ((1ULL << 32) - 1);
  u128 X_1 = c & ((1ULL >> 64) - 1);
  u128 C_out = X_1 + (X_2 * (1ULL << 32) - 1) - X_3;
  if (C_out >= 18446744069414584321ULL) {
    C_out -= 18446744069414584321ULL;
  }

  return (u64)C_out;
}
