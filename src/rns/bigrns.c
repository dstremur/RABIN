#include "../../include/bigrns.h"

u64 mod_add(u64 a, u64 b, u64 p)
{
  u64 res = a + b;
  // If res >= p or if an overflow occurred (res < a)
  if (res >= p || res < a) {
    res -= p;
  }
  return res;
}

u64 mod_sub(u64 a, u64 b, u64 p)
{
  if (a < b) {
    return a + p - b;
  }
  return a - b;
}

u64 mod_mul(u64 a, u64 b, u64 p)
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
u64 mod_inverse(u64 n, u64 p) { return mod_pow(n, p - 2, p); }

// estimate how my primes are needed
u64 bigrns_estimate_primes(const bignum* a)
{
  u64 k = bn_bit_length(a);
  u64 res = (k + 1) / 62;

  return res;
}
