#include "../../include/u64.h"
#define TILE_SIZE 64

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
      for (; k <= (n >= TILE_SIZE ? n - TILE_SIZE : 0); k += TILE_SIZE) {
        for (u64 tk = 0; tk < TILE_SIZE; tk++) {
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