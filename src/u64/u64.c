#include "../include/u64.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../../include/bigrns.h"

void mont_init(mont_ctx* ctx, u64 p)
{
  ctx->p = p;

  ctx->p_inv = mod_inverse_u64(p);

  // Calculate 2^128 mod p
  unsigned __int128 r2 = ((unsigned __int128)1 << 64) % p;
  r2 = (r2 * r2) % p;
  ctx->r2_mod_p = (u64)r2;
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
  u64 res = a - b;
    return res + (p & (u64)((i64)res >> 63));
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

u64 mont_redc(unsigned __int128 T, mont_ctx* ctx)
{
  u64 p = ctx->p; 

  // T mod R
  u64 T_modR = (u64)T; 

  u64 m = T_modR * ctx->p_inv;

  unsigned __int128 t_wide = (T + (unsigned __int128)m * p);
  u64 t = (u64)(t_wide >> 64);
        
  return (t >= p) ? (t - p) : t;
}

u64 mont_mul(u64 a, u64 b, mont_ctx* ctx)
{
  unsigned __int128 R = (unsigned __int128)a * b; 

  return mont_redc(R, ctx);
}

u64 mont_in(u64 a, mont_ctx* ctx)
{
  return mont_mul(a, ctx->r2_mod_p, ctx);
}

u64 mont_out(u64 a_hat, mont_ctx* ctx)
{
  return mont_mul(a_hat, 1, ctx);
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

u64 mont_inverse(u64 a_mont, const mont_ctx* ctx)
{
  u64 a = mont_out(a_mont, ctx);

  u64 inv = mod_inverse_euclid(a, ctx->p);

  return mont_in(inv, ctx);
}
#define TILE_SIZE 64

// works in place, need to pass copy
u64 matrix_u64_det_optimized(u64* mat, u64 n, const mont_ctx* ctx)
{
  u64 det = mont_in(1, ctx);
  u64 p = ctx->p;

  // convert matrix to montgomery form
  for (u64 i = 0; i < n * n; i++) {
    mat[i] = mont_in(mat[i], ctx);
  }

  for (u64 i = 0; i < n; i++) {
    // find pivot
    u64 pivot = i;
    while (pivot < n && mat[pivot * n + i] == 0) pivot++;

    if (pivot == n) {
      return 0;
    }

    if (pivot != i) {
      u64* row_i = mat + i * n;
      u64* row_p = mat + pivot * n;
      for (u64 j = i; j < n; j++) {
        u64 tmp = row_i[j];
        row_i[j] = row_p[j];
        row_p[j] = tmp;
      }
      det = mod_sub(0, det, p);
    }

    u64 pivot_val = mat[i * n + i];
    det = mont_mul(det, pivot_val, ctx);

    u64 pivot_real = mont_out(pivot_val, ctx);
    u64 inv_real = mod_inverse_euclid(pivot_real, p);
    u64 inv = mont_in(inv_real, ctx);

    u64* row_i = mat + i * n;

    for (u64 j = i + 1; j < n; j++) {
      u64 factor = mont_mul(mat[j * n + i], inv, ctx);
      if (factor == 0) continue;
      u64* row_j = mat + j * n;

      // use tiling 
      u64 k = i + 1;
      for (; k <= (n >= TILE_SIZE ? n - TILE_SIZE: 0); k += TILE_SIZE){
        for (u64 tk = 0; tk < TILE_SIZE; tk++){
          u64 idx = k + tk;
          u64 prod = mont_mul(factor, row_i[idx], ctx);
          row_j[idx] = mod_sub(row_j[idx], prod, p);
        }
      }

      // process remaining elements
            for (; k < n; k++) {
                u64 prod = mont_mul(factor, row_i[k], ctx);
                row_j[k] = mod_sub(row_j[k], prod, p);
            }
    }
  }

  return mont_out(det, ctx);
}
