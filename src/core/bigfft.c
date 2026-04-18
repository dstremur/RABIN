#include <stdio.h>

#include "../include/bignum.h"

const u64 NTT_PRIME = 998244353;
const u64 NTT_ROOT = 3;  // Primitive root modulo 998244353

// Basic modular exponentiation (x^y mod p)
u64 power(u64 base, u64 exp, u64 mod)
{
  u64 res = 1;
  base %= mod;
  while (exp > 0) {
    if (exp % 2 == 1) res = (res * base) % mod;
    base = (base * base) % mod;
    exp /= 2;
  }
  return res;
}

// Modular inverse
u64 modInverse(u64 n, u64 mod)
{
  return power(n, mod - 2, mod);  // Fermat's Little Theorem
}

// The NTT algorithm. 'a' is an array of size 'n' (must be power of 2).
void ntt(u64* a, int n, bool invert)
{
  // 1. Bit-reversal permutation (puts array in correct order for bottom-up
  // processing)
  for (int i = 1, j = 0; i < n; i++) {
    int bit = n >> 1;
    for (; j & bit; bit >>= 1) j ^= bit;
    j ^= bit;
    if (i < j) {
      u64 temp = a[i];
      a[i] = a[j];
      a[j] = temp;
    }
  }

  // 2. Cooley-Tukey Butterfly operations
  for (int len = 2; len <= n; len <<= 1) {
    // Calculate the root of unity for this level
    u64 wlen = power(NTT_ROOT, (NTT_PRIME - 1) / len, NTT_PRIME);
    if (invert) {
      wlen = modInverse(wlen, NTT_PRIME);
    }

    for (int i = 0; i < n; i += len) {
      u64 w = 1;
      for (int j = 0; j < len / 2; j++) {
        u64 u = a[i + j];
        u64 v = (a[i + j + len / 2] * w) % NTT_PRIME;

        a[i + j] = (u + v < NTT_PRIME ? u + v : u + v - NTT_PRIME);
        a[i + j + len / 2] = (u >= v ? u - v : u + NTT_PRIME - v);

        w = (w * wlen) % NTT_PRIME;
      }
    }
  }

  // 3. If inverse transform, scale by 1/n
  if (invert) {
    u64 n_inv = modInverse(n, NTT_PRIME);
    for (int i = 0; i < n; i++) {
      a[i] = (a[i] * n_inv) % NTT_PRIME;
    }
  }
}

void schoolbook_mul(u64* a, int an, u64* b, int bn, u64* res)
{
  for (int i = 0; i < an; i++) {
    for (int j = 0; j < bn; j++) {
      res[i + j] = (res[i + j] + (a[i] * b[j])) % NTT_PRIME;
    }
  }
}

void fftest()
{
  int n = 8;  // Must be power of 2
  u64 original[8] = {1, 2, 3, 4, 0, 0, 0, 0};
  u64 a[8] = {1, 2, 3, 4, 0, 0, 0, 0};

  printf("--- Test 1: Identity (Forward then Inverse) ---\n");
  ntt(a, n, false);  // Forward
  ntt(a, n, true);   // Inverse

  bool identity_pass = true;
  for (int i = 0; i < n; i++) {
    if (a[i] != original[i]) identity_pass = false;
    printf("Index %d: Expected %lu, Got %lu\n", i, original[i], a[i]);
  }
  printf("Result: %s\n\n", identity_pass ? "PASS" : "FAIL");

  printf("--- Test 2: Polynomial Multiplication ---\n");
  // Reset arrays
  u64 poly_a[8] = {1, 2, 3, 4, 0, 0, 0, 0};
  u64 poly_b[8] = {5, 6, 7, 8, 0, 0, 0, 0};
  u64 expected[8] = {0};

  // 1. Get expected result via schoolbook
  schoolbook_mul(poly_a, 4, poly_b, 4, expected);

  // 2. Get NTT result
  ntt(poly_a, n, false);
  ntt(poly_b, n, false);
  for (int i = 0; i < n; i++) {
    poly_a[i] = (poly_a[i] * poly_b[i]) % NTT_PRIME;
  }
  ntt(poly_a, n, true);

  bool mul_pass = true;
  for (int i = 0; i < n; i++) {
    if (poly_a[i] != expected[i]) mul_pass = false;
    printf("Index %d: Schoolbook %lu, NTT %lu\n", i, expected[i], poly_a[i]);
  }
  printf("Result: %s\n", mul_pass ? "PASS" : "FAIL");
}
