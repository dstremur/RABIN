#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../include/u64.h"
#include "../../include/bigrns.h"


void barrett_init(barrett_ctx* ctx, u64 p)
{

  ctx->p = p;

  // max = 2^128 - 1
  unsigned __int128 max = ~((unsigned __int128)0);

  ctx->mu = max / p;
}

u64 barrett_reduce(unsigned __int128 a, const barrett_ctx* ctx) {
  // 1. Estimate quotient: q = (a * mu) >> 128
  unsigned __int128 q = mul128_high(a, ctx->mu);

  // 2. Calculate remainder: r = a - q * p
  // We can safely cast to u64 here because we know r < 2 * p, 
  // which fits easily within 64 bits. Modulo 2^64 arithmetic handles this perfectly.
  u64 r = (u64)a - (u64)(q * ctx->p);

  // 3. Correction step
  // Barrett quotient approximation can be slightly too small (usually by 0 or 1)
  while (r >= ctx->p) {
    r -= ctx->p;
  }

  return r;
}

u64 mod_mul_barrett(u64 a, u64 b, const barrett_ctx* ctx) {
  unsigned __int128 res = (unsigned __int128)a * b;
  return barrett_reduce(res, ctx);
}

u64 mod_add(u64 a, u64 b, u64 p)
{
  u64 res = a + b;
  if (res >= p || res < a) {
    res -= p;
  }
  return res;
}

u64 mod_sub(u64 a, u64 b, u64 p)
{
  if (a < b) {
    return (p - b) + a;
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

u64 matrix_u64_det(matrix_u64* M)
{
  if (M->r_size != M->c_size) {
    printf("Not square\n");
  }
  u64 p = M->modulus;
  u64 det = 1;
  u64 n = M->r_size;

  u64* mat = malloc(sizeof(u64) * n * n);
  memcpy(mat, M->data, sizeof(u64) * n * n);

  for (u64 i = 0; i < n; i++) {
    // find pivot
    u64 pivot = i;
    while (pivot < n && mat[pivot * n + i] == 0) {
      pivot++;
    }

    if (pivot == n) {
      free(mat);
      return 0;
    }

    // swap rows
    if (pivot != i) {
      for (u64 j = i; j < n; j++) {
        u64 tmp = mat[i * n + j];
        mat[i * n + j] = mat[pivot * n + j];
        mat[pivot * n + j] = tmp;
      }
      // swapping multiplies det by -1 or p-1 mod p
      det = mod_sub(0, det, p);
    }

    // multiply det by pivot
    u64 pivot_val = mat[i * n + i];
    det = mod_mul(det, pivot_val, p);

    // eliminate below pivot
    u64 inv = mod_inverse_euclid(pivot_val, p);
    for (u64 j = i + 1; j < n; j++) {
      u64 factor = mod_mul(mat[j * n + i], inv, p);
      for (u64 k = i; k < n; k++) {
        u64 sub = mod_mul(factor, mat[i * n + k], p);
        mat[j * n + k] = mod_sub(mat[j * n + k], sub, p);
      }
    }
  }
  free(mat);
  return det;
}

u64 matrix_u64_det_optimized(u64* data, u64 n, const barrett_ctx* ctx)
{
  u64 det = 1;
  u64 p = ctx->p;
  // Work on a copy to avoid destroying the original reduced matrix
  u64* mat = malloc(sizeof(u64) * n * n);
  memcpy(mat, data, sizeof(u64) * n * n);

  for (u64 i = 0; i < n; i++) {
    u64 pivot = i;
    while (pivot < n && mat[pivot * n + i] == 0) pivot++;

    if (pivot == n) {
      free(mat);
      return 0;
    }

    if (pivot != i) {
      for (u64 j = i; j < n; j++) {
        u64 tmp = mat[i * n + j];
        mat[i * n + j] = mat[pivot * n + j];
        mat[pivot * n + j] = tmp;
      }
      det = mod_sub(0, det, p);
    }

    u64 pivot_val = mat[i * n + i];
    det = mod_mul_barrett(det, pivot_val, ctx);

    u64 inv = mod_inverse_euclid(pivot_val, p);
    for (u64 j = i + 1; j < n; j++) {
      u64 factor = mod_mul_barrett(mat[j * n + i], inv, ctx);
      for (u64 k = i + 1; k < n; k++) {  // Start k from i+1
        u64 sub = mod_mul_barrett(factor, mat[i * n + k], ctx);
        mat[j * n + k] = mod_sub(mat[j * n + k], sub, p);
      }
    }
  }
  free(mat);
  return det;
}