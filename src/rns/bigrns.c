#include "../../include/bigrns.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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

// estimate how my primes are needed
u64 rns_estimate_primes(const bignum* a)
{
  u64 k = bn_bit_length(a);
  u64 res = (k + 62) / 62;

  return res;
}

void rns_context_init(ctx_rns* ctx, const u64* primes, u64 count)
{
  ctx->primes = malloc(sizeof(u64) * count);
  ctx->crt_weights = malloc(sizeof(bignum) * count);
  memcpy(ctx->primes, primes, sizeof(u64) * count);
  ctx->count = count;

  bignum tmp;
  bn_init_multi(&tmp, NULL);
  bn_init(&ctx->prod);
  bn_set_u64(&ctx->prod, 1);

  // calculate produt
  for (u64 i = 0; i < count; i++) {
    bn_set_u64(&tmp, primes[i]);
    bn_mul(&ctx->prod, &ctx->prod, &tmp);
  }

  // calculate crt weights
  for (u64 i = 0; i < count; i++) {
    bn_init(&ctx->crt_weights[i]);

    bignum M_div_tmp, inv_bn;
    bn_init_multi(&M_div_tmp, &inv_bn, NULL);

    bn_set_u64(&tmp, primes[i]);
    bn_div(&M_div_tmp, &ctx->prod, &tmp);

    // find inverse
    u64 m_mod_p = bn_mod_u64(&M_div_tmp, primes[i]);
    u64 inv = mod_inverse_euclid(m_mod_p, primes[i]);

    bn_set_u64(&inv_bn, inv);

    bn_mul(&ctx->crt_weights[i], &M_div_tmp, &inv_bn);

    bn_free_multi(&M_div_tmp, &inv_bn, NULL);
  }

  bn_free_multi(&tmp, NULL);
}

void rns_context_free(ctx_rns* ctx)
{
  if (!ctx) return;

  // 1. Free the individual bignum weights
  if (ctx->crt_weights) {
    for (u64 i = 0; i < ctx->count; i++) {
      bn_free(&ctx->crt_weights[i]);
    }
    free(ctx->crt_weights);
    ctx->crt_weights = NULL;
  }

  // 2. Free the big product
  bn_free(&ctx->prod);

  // 3. Free the primes array
  if (ctx->primes) {
    free(ctx->primes);
    ctx->primes = NULL;
  }

  ctx->count = 0;
}

void bignum_to_rns(rns_num* r, const bignum* a, ctx_rns* ctx)
{
  if (!r->residues) {
    r->residues = malloc(sizeof(u64) * ctx->count);
  }
  r->size = ctx->count;

  for (u64 i = 0; i < r->size; i++) {
    r->residues[i] = bn_mod_u64(a, ctx->primes[i]);
  }
}

void rns_add(rns_num* r, const rns_num* a, const rns_num* b, const ctx_rns* ctx)
{
  for (u64 i = 0; i < r->size; i++) {
    r->residues[i] = mod_add(a->residues[i], b->residues[i], ctx->primes[i]);
  }
}

void rns_to_bignum1(bignum* a, const rns_num* r, ctx_rns* ctx)
{
  bignum sum, tmp;
  bn_init_multi(&sum, &tmp, NULL);

  bn_set_u64(a, 0);

  for (u64 i = 0; i < r->size; i++) {
    bn_set_u64(&tmp, r->residues[i]);
    bn_mul(&sum, &ctx->crt_weights[i], &tmp);

    bn_add(a, a, &sum);
  }

  bn_mod(a, a, &ctx->prod);

  bn_free_multi(&sum, &tmp, NULL);
}

void rns_to_bignum(bignum* a, const rns_num* r, ctx_rns* ctx)
{
  u64 n = r->size;
  u64* mixed_radix = malloc(sizeof(u64) * n);

  for (u64 i = 0; i < n; i++) {
    u64 temp = r->residues[i];
    for (u64 j = 0; j < i; j++) {
      u64 inv = mod_inverse_euclid(ctx->primes[j], ctx->primes[i]);
      u64 diff = mod_sub(temp, mixed_radix[j], ctx->primes[i]);
      temp = mod_mul(diff, inv, ctx->primes[i]);
    }
    mixed_radix[i] = temp;
  }

  bn_set_u64(a, 0);
  bignum term, p_prod, p_val;
  bn_init_multi(&term, &p_prod, &p_val, NULL);
  bn_set_u64(&p_prod, 1);

  for (u64 i = 0; i < n; i++) {
    bn_set_u64(&term, mixed_radix[i]);
    bn_mul(&term, &term, &p_prod);
    bn_add(a, a, &term);

    bn_set_u64(&p_val, ctx->primes[i]);
    bn_mul(&p_prod, &p_prod, &p_val);
  }

  bn_free_multi(&term, &p_prod, &p_val, NULL);
  free(mixed_radix);
}

// Estimates primes for a matrix determinant (The most robust way)
u64 rns_estimate_determinant(bigmatrix* A)
{
  bignum res;
  bn_init(&res);

  // hadamard bound
  bigmatrix_hadamard(&res, A);

  // double to get to range [-M, M]
  bn_lshift1(&res);

  return rns_estimate_primes(&res);
}

void bigmatrix_reduce(matrix_u64* R, const bigmatrix* A, u64 m)
{
  R->data = malloc(sizeof(u64) * A->c_size * A->r_size);
  R->c_size = A->c_size;
  R->r_size = A->r_size;
  R->modulus = m;

  for (u64 i = 0; i < A->r_size; i++) {
    for (u64 j = 0; j < A->c_size; j++) {
      R->data[i * R->c_size + j] = bn_mod_u64(GET(A, i, j), m);
    }
  }
}

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

u64 matrix_u64_det_optimized(u64* data, u64 n, u64 p)
{
  u64 det = 1;
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
    det = mod_mul(det, pivot_val, p);

    u64 inv = mod_inverse_euclid(pivot_val, p);
    for (u64 j = i + 1; j < n; j++) {
      u64 factor = mod_mul(mat[j * n + i], inv, p);
      for (u64 k = i + 1; k < n; k++) {  // Start k from i+1
        u64 sub = mod_mul(factor, mat[i * n + k], p);
        mat[j * n + k] = mod_sub(mat[j * n + k], sub, p);
      }
    }
  }
  free(mat);
  return det;
}

void bigmatrix_det_rns(bignum* det, const bigmatrix* A, const ctx_rns* ctx)
{
  // create array of n matrices mod p_i

  u64 n = A->r_size;
  u64* residues = malloc(sizeof(u64) * ctx->count);

#pragma omp parallel for schedule(static)
  for (u64 i = 0; i < ctx->count; i++) {
    u64 p = ctx->primes[i];

    u64* reduced_data = malloc(sizeof(u64) * n * n);

#pragma omp simd
    for (u64 r = 0; r < n; r++) {
#pragma omp simd
      for (u64 c = 0; c < n; c++) {
        reduced_data[r * n + c] = bn_mod_u64(GET(A, r, c), p);
      }
    }
    residues[i] = matrix_u64_det_optimized(reduced_data, n, p);
    free(reduced_data);
  }

  rns_num r;
  r.residues = residues;
  r.size = ctx->count;

  // reconstruct x in [0, M - 1]
  rns_to_bignum(det, &r, ctx);

  // if d > M / 2 det is negative

  bignum M_half;
  bn_init_multi(&M_half, NULL);
  bn_copy(&M_half, &ctx->prod);
  bn_rshift1(&M_half);

  if (bn_cmp(det, &M_half) > 0) {
    bn_sub(det, det, &ctx->prod);
  }

  bn_free(&M_half);
  free(residues);
}
