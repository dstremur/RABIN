/*
 * rmat_kernel.c
 *
 * Integer kernel computations for rmat_t: LLL-based kernel, kernel of the image
 * and the rational inverse helper.
 *
 * Copyright (C) 2026 Diego Strebel
 * SPDX-License-Identifier: GPL-3.0-or-later
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#include <math.h>
#include <string.h>

#include "../../include/rmat.h"
#include "../../include/rq.h"
#include "../../include/rzmath.h"
#include "rmat_internal.h"
/* ------------- Kernel and image of a matrix (Cohen 2.7) ------------- */

/*
 * Extract columns [c0, c0 + len) of src into a fresh matrix dst.
 */
static rabin_err_t rmat_extract_cols(rmat_t* dst, const rmat_t* src, u64 c0,
                                     u64 len)
{
  if (dst == NULL || src == NULL) return RABIN_ERR_NULL_PTR;
  if (c0 + len > src->cols) return RABIN_ERR_INVALID_ARG;

  rabin_err_t err = rmat_init(dst, src->rows, len);
  if (err != RABIN_SUCCESS) return err;
  for (u64 x = 0; x < src->rows && err == RABIN_SUCCESS; x++) {
    for (u64 y = 0; y < len; y++) {
      err = rz_copy(RMAT_GET(dst, x, y), RMAT_GET(src, x, c0 + y));
    }
  }
  return err;
}

/*
 * Kernel-only LLL (Cohen 2.7.2): integer Gram-Schmidt with flags iota[]
 * marking the vectors that are nonzero in the Gram-Schmidt sense (iota = 1
 * iff B_k != 0). When u == 0 the convention d_k <- d_{k-1} keeps every d_k
 * nonzero, so no division ever hits zero and dependent vectors are simply
 * pushed to the front by SWAPK until iota[0..r-1] == 0 exactly delimits the
 * kernel: the first r columns of H then form a Z-basis of ker(A).
 */

/*
 * Exact integer Gram-Schmidt for one column (1-based row i of lam): fills
 * lambda[i][j] (j < i), d_i and iota[i-1] from the vectors in V and the GS
 * data of the rows below i. Mirrors the step-2 loop of rmat_lll_run.
 */
static rabin_err_t lll_gsrow(rmat_t* V, rmat_t* lam, rz_t* dv, bool* iota,
                             u64 i, lll_scr_t* s)
{
  rabin_err_t err;
  for (u64 j = 1; j <= i; j++) {
    if (!iota[j - 1] && j < i) {
      if ((err = rz_set_u64(RMAT_GET(lam, i - 1, j - 1), 0)) != RABIN_SUCCESS)
        return err;
      continue;
    }
    if ((err = lll_dot(&s->u, V, i - 1, j - 1, &s->v)) != RABIN_SUCCESS)
      return err;
    for (u64 x = 1; x < j; x++) {
      if (!iota[x - 1]) continue;
      if ((err = rz_mul(&s->w, &dv[x], &s->u)) != RABIN_SUCCESS) return err;
      if ((err = rz_mul(&s->v, RMAT_GET(lam, i - 1, x - 1),
                        RMAT_GET(lam, j - 1, x - 1))) != RABIN_SUCCESS)
        return err;
      if ((err = rz_sub(&s->w, &s->w, &s->v)) != RABIN_SUCCESS) return err;
      if ((err = rz_div_exact(&s->u, &s->w, &dv[x - 1])) != RABIN_SUCCESS)
        return err;
    }
    if (j < i) {
      if ((err = rz_copy(RMAT_GET(lam, i - 1, j - 1), &s->u)) != RABIN_SUCCESS)
        return err;
    } else {
      if (!rz_is_zero(&s->u)) {
        if ((err = rz_copy(&dv[i], &s->u)) != RABIN_SUCCESS) return err;
        iota[i - 1] = true;
      } else {
        if ((err = rz_copy(&dv[i], &dv[i - 1])) != RABIN_SUCCESS) return err;
        iota[i - 1] = false;
      }
    }
  }
  return RABIN_SUCCESS;
}

/*
 * SWAPK (Cohen 2.7.2): exchange the vectors at (1-based) positions k-1 and k.
 * On entry iota[k-2] is true and iota[k-1] is false (the vector at position
 * k is Gram-Schmidt-zero). Instead of updating lambda/d incrementally, the
 * GS data of rows k-1..kmax is recomputed from the (already swapped) vectors
 * with lll_gsrow(); this is exact by construction and keeps iota in sync.
 */
static rabin_err_t lll_swapk(rmat_t* V, rmat_t* H, rmat_t* lam, rz_t* dv,
                             bool* iota, u64 k, u64 kmax, lll_scr_t* s)
{
  rabin_err_t err;

  // swap the columns k-1 and k of V and H
  for (u64 i = 0; i < V->rows; i++) {
    if ((err = rz_copy(&s->u, RMAT_GET(V, i, k - 1))) != RABIN_SUCCESS)
      return err;
    if ((err = rz_copy(RMAT_GET(V, i, k - 1), RMAT_GET(V, i, k - 2))) !=
        RABIN_SUCCESS)
      return err;
    if ((err = rz_copy(RMAT_GET(V, i, k - 2), &s->u)) != RABIN_SUCCESS)
      return err;
  }
  for (u64 i = 0; i < H->rows; i++) {
    if ((err = rz_copy(&s->u, RMAT_GET(H, i, k - 1))) != RABIN_SUCCESS)
      return err;
    if ((err = rz_copy(RMAT_GET(H, i, k - 1), RMAT_GET(H, i, k - 2))) !=
        RABIN_SUCCESS)
      return err;
    if ((err = rz_copy(RMAT_GET(H, i, k - 2), &s->u)) != RABIN_SUCCESS)
      return err;
  }

  // recompute the GS data of all rows at or below the swap (ascending, so
  // each row sees the freshly updated data of the rows above it)
  for (u64 i = k - 1; i <= kmax; i++) {
    if ((err = lll_gsrow(V, lam, dv, iota, i, s)) != RABIN_SUCCESS) return err;
  }
  return RABIN_SUCCESS;
}

static rabin_err_t rmat_kernel_lll_run(rmat_t* V, rmat_t* H, rmat_t* lam,
                                       rz_t* dv, bool* iota, u64 n)
{
  lll_scr_t s;
  if (lll_scr_init(&s) != RABIN_SUCCESS) return RABIN_ERR_OUT_OF_MEMORY;

  rabin_err_t err = rz_set_u64(&dv[0], 1);  // d_0 = 1
  if (err != RABIN_SUCCESS) {
    lll_scr_clear(&s);
    return err;
  }

  u64 k = 1;
  u64 kmax = 0;

step2:
  // [Initialize Gram-Schmidt data for vector k]
  if (k > kmax) {
    kmax = k;
    TRY(lll_gsrow(V, lam, dv, iota, k, &s));
  }

step3:
  // [Test]
  if (k >= 2) {
    if (iota[k - 2]) {
      TRY(lll_red1(V, H, lam, dv, k - 1, k - 2, &s));
      if (!iota[k - 1]) {
        TRY(lll_swapk(V, H, lam, dv, iota, k, kmax, &s));
        if (k > 2) k--;  // k <- max(2, k-1)
        goto step3;
      }
    }
    // descending: RED(k,l) updates lambda[k][i] for i < l
    for (u64 l = k - 2; l >= 1; l--) {
      if (!iota[l - 1]) continue;
      TRY(lll_red1(V, H, lam, dv, k - 1, l - 1, &s));
    }
  }
  k++;
  if (k <= n) goto step2;

  err = RABIN_SUCCESS;
cleanup:
  lll_scr_clear(&s);
  return err;
}

/*
 * Nearest integer of the rational x (ties go up): with den > 0 the invariant
 * of rq_t,
 *
 *   round(num/den) = floor((2*num + den) / (2*den)),
 *
 * evaluated with the exact Euclidean quotient. t1, t2 are scratch.
 */
static rabin_err_t rz_round_nearest(rz_t* q, const rq_t* x, rz_t* t1, rz_t* t2)
{
  rabin_err_t err = rz_lshift(t1, &x->num, 1);  // 2*num
  if (err == RABIN_SUCCESS) err = rz_add(t1, t1, &x->den);
  if (err == RABIN_SUCCESS) err = rz_lshift(t2, &x->den, 1);  // 2*den > 0
  if (err == RABIN_SUCCESS) err = rz_div_euclid(q, t1, t2);
  return err;
}

/*
 * Invert the r x r rational matrix G by Gauss-Jordan elimination, returning
 * the inverse as a freshly allocated row-major array of r*r rationals
 * (freed by the caller with rq_clear + free). Fails with
 * RABIN_ERR_INVALID_ARG when G is singular.
 */
static rabin_err_t rq_mat_inverse(rq_t** Dout, const rmat_t* G, u64 r)
{
  if (Dout == NULL || G == NULL) return RABIN_ERR_NULL_PTR;
  *Dout = NULL;
  if (G->rows != r || G->cols != r) return RABIN_ERR_INVALID_ARG;
  if (r == 0) return RABIN_SUCCESS;

  rabin_err_t err = RABIN_SUCCESS;
  rq_t* A = calloc(r * 2 * r, sizeof(rq_t));
  if (A == NULL) return RABIN_ERR_OUT_OF_MEMORY;
  for (u64 x = 0; x < r * 2 * r; x++) {
    rq_init(&A[x]);
  }

  rq_t pv, f, t;
  rq_init_multi(&pv, &f, &t, NULL);

  rq_t* D = NULL;

  for (u64 x = 0; x < r && err == RABIN_SUCCESS; x++) {
    for (u64 y = 0; y < 2 * r; y++) {
      if (y < r) {
        err = rq_set(&A[x * 2 * r + y], RMAT_GET(G, x, y));
      } else {
        err = rq_set_u64(&A[x * 2 * r + y], y - r == x ? 1 : 0);
      }
    }
  }

  for (u64 c = 0; c < r && err == RABIN_SUCCESS; c++) {
    u64 p = c;
    while (p < r && rz_is_zero(&A[p * 2 * r + c].num)) p++;
    if (p == r) {
      err = RABIN_ERR_INVALID_ARG;  // singular
      break;
    }
    if (p != c) {
      for (u64 y = 0; y < 2 * r; y++) {
        err = rq_copy(&t, &A[p * 2 * r + y]);
        if (err == RABIN_SUCCESS)
          err = rq_copy(&A[p * 2 * r + y], &A[c * 2 * r + y]);
        if (err == RABIN_SUCCESS) err = rq_copy(&A[c * 2 * r + y], &t);
        if (err != RABIN_SUCCESS) break;
      }
    }
    if (err != RABIN_SUCCESS) break;

    // normalize row c so that the pivot becomes 1
    err = rq_copy(&pv, &A[c * 2 * r + c]);
    for (u64 y = c + 1; err == RABIN_SUCCESS && y < 2 * r; y++) {
      err = rq_div(&A[c * 2 * r + y], &A[c * 2 * r + y], &pv);
    }
    if (err == RABIN_SUCCESS) err = rq_set_u64(&A[c * 2 * r + c], 1);

    // clear the other rows
    for (u64 x = 0; x < r && err == RABIN_SUCCESS; x++) {
      if (x == c) continue;
      if (rz_is_zero(&A[x * 2 * r + c].num)) continue;
      err = rq_copy(&f, &A[x * 2 * r + c]);
      if (err == RABIN_SUCCESS) err = rq_set_u64(&A[x * 2 * r + c], 0);
      for (u64 y = c + 1; err == RABIN_SUCCESS && y < 2 * r; y++) {
        if (rz_is_zero(&A[c * 2 * r + y].num)) continue;
        err = rq_mul(&t, &f, &A[c * 2 * r + y]);
        if (err == RABIN_SUCCESS)
          err = rq_sub(&A[x * 2 * r + y], &A[x * 2 * r + y], &t);
      }
    }
  }
  if (err != RABIN_SUCCESS) goto cleanup;

  D = malloc(r * r * sizeof(rq_t));
  if (D == NULL) {
    err = RABIN_ERR_OUT_OF_MEMORY;
    goto cleanup;
  }
  for (u64 x = 0; x < r * r; x++) {
    rq_init(&D[x]);
  }
  for (u64 x = 0; x < r * r; x++) {
    if (rq_copy(&D[x], &A[(x / r) * 2 * r + r + (x % r)]) != RABIN_SUCCESS) {
      err = RABIN_ERR_OUT_OF_MEMORY;
      goto cleanup;
    }
  }
  *Dout = D;
  D = NULL;  // ownership transferred to the caller

cleanup:
  if (D != NULL) {
    for (u64 x = 0; x < r * r; x++) {
      rq_clear(&D[x]);
    }
    free(D);
  }
  rq_clear_multi(&pv, &f, &t, NULL);
  for (u64 x = 0; x < r * 2 * r; x++) {
    rq_clear(&A[x]);
  }
  free(A);
  return err;
}

/*
 * Steps 3-5 of Cohen 2.7.1: size-reduce every column of H with index
 * >= r against the (fixed) kernel basis formed by the first r columns.
 * H is updated in place.
 */
static rabin_err_t rmat_reduce_against_kernel(rmat_t* H, u64 n, u64 r)
{
  if (r == 0) return RABIN_SUCCESS;

  rabin_err_t err = RABIN_SUCCESS;
  rmat_t G = {0};
  rz_t* mr = NULL;
  rq_t* vq = NULL;
  rq_t* D = NULL;
  u64 mr_done = 0;
  rq_t acc, tmpq;
  rz_t t1, t2;
  bool temps = false;

  if ((err = rmat_init(&G, r, r)) != RABIN_SUCCESS) return err;
  mr = calloc(r, sizeof(rz_t));
  vq = calloc(r, sizeof(rq_t));
  if (mr == NULL || vq == NULL) {
    err = RABIN_ERR_OUT_OF_MEMORY;
    goto cleanup;
  }
  for (u64 j = 0; j < r; j++) {
    if (rz_init(&mr[j]) != RABIN_SUCCESS) goto cleanup;
    mr_done++;
  }
  for (u64 j = 0; j < r; j++) {
    if (rq_init(&vq[j]) != RABIN_SUCCESS) goto cleanup;
  }
  TRY(rq_init_multi(&acc, &tmpq, NULL));
  TRY(rz_init_multi(&t1, &t2, NULL));
  temps = true;

  // e = Gram matrix of the kernel basis columns
  for (u64 j = 0; j < r && err == RABIN_SUCCESS; j++) {
    for (u64 k = 0; k < r; k++) {
      err = lll_dot(RMAT_GET(&G, j, k), H, j, k, &t1);
    }
  }
  if (err == RABIN_SUCCESS) err = rq_mat_inverse(&D, &G, r);
  if (err != RABIN_SUCCESS) goto cleanup;

  for (u64 i = r; i < n && err == RABIN_SUCCESS; i++) {
    // v_j = H_i . H_j
    for (u64 j = 0; j < r; j++) {
      err = lll_dot(&t1, H, i, j, &t2);
      if (err != RABIN_SUCCESS) break;
      err = rq_set(&vq[j], &t1);
    }
    if (err != RABIN_SUCCESS) break;

    // m_j = nearest integer of x_j = sum_k D[j][k] * v_k
    for (u64 j = 0; j < r; j++) {
      err = rq_set_u64(&acc, 0);
      for (u64 k = 0; err == RABIN_SUCCESS && k < r; k++) {
        err = rq_mul(&tmpq, &D[j * r + k], &vq[k]);
        if (err == RABIN_SUCCESS) err = rq_add(&acc, &acc, &tmpq);
      }
      if (err != RABIN_SUCCESS) break;
      err = rz_round_nearest(&mr[j], &acc, &t1, &t2);
    }
    if (err != RABIN_SUCCESS) break;

    // H_i <- H_i - sum_j m_j * H_j
    for (u64 j = 0; j < r; j++) {
      if (rz_is_zero(&mr[j])) continue;
      for (u64 x = 0; x < n; x++) {
        err = rz_mul(&t1, &mr[j], RMAT_GET(H, x, j));
        if (err == RABIN_SUCCESS)
          err = rz_sub(RMAT_GET(H, x, i), RMAT_GET(H, x, i), &t1);
        if (err != RABIN_SUCCESS) break;
      }
      if (err != RABIN_SUCCESS) break;
    }
  }

cleanup:
  if (temps) {
    rq_clear_multi(&acc, &tmpq, NULL);
    rz_clear_multi(&t1, &t2, NULL);
  }
  if (D != NULL) {
    for (u64 x = 0; x < r * r; x++) rq_clear(&D[x]);
    free(D);
  }
  for (u64 j = 0; j < r && vq != NULL; j++) rq_clear(&vq[j]);
  free(vq);
  for (u64 j = 0; j < mr_done; j++) rz_clear(&mr[j]);
  free(mr);
  rmat_clear(&G);
  return err;
}

rabin_err_t rmat_kernel(const rmat_t* A, rmat_t* K)
{
  if (A == NULL || K == NULL) return RABIN_ERR_NULL_PTR;

  const u64 n = A->cols;
  if (n == 0) return rmat_init(K, 0, 0);

  rmat_t W = {0}, U = {0}, Uk = {0};
  rmat_t Kraw = {0};

  // no coordinates: every vector of Z^n is in the kernel
  if (A->rows == 0) {
    rabin_err_t err = rmat_init(K, n, n);
    if (err != RABIN_SUCCESS) return err;
    return rmat_id(K, n);
  }

  rabin_err_t err = hnf_with_transform(&W, &U, A);
  if (err != RABIN_SUCCESS) goto cleanup;

  const u64 r = n - W.cols;
  if (r == 0) {
    err = rmat_init(K, n, 0);
    goto cleanup;
  }

  err = rmat_extract_cols(&Kraw, &U, 0, r);
  if (err != RABIN_SUCCESS) goto cleanup;
  if (r > 1) {
    err = rmat_lll(&Kraw, r, 0.75, &Uk);
    if (err != RABIN_SUCCESS) goto cleanup;
  }
  rmat_swap(K, &Kraw);

cleanup:
  rmat_clear(&Kraw);
  rmat_clear(&Uk);
  rmat_clear(&U);
  rmat_clear(&W);
  return err;
}

rabin_err_t rmat_kernel_lll(const rmat_t* A, rmat_t* K)
{
  if (A == NULL || K == NULL) return RABIN_ERR_NULL_PTR;

  const u64 n = A->cols;
  if (n == 0) return rmat_init(K, 0, 0);

  rmat_t V = {0}, H = {0}, lam = {0}, Kraw = {0}, Uk = {0};
  rz_t* dv = NULL;
  bool* iota = NULL;

  // no coordinates: every vector of Z^n is in the kernel
  if (A->rows == 0) {
    rabin_err_t err = rmat_init(K, n, n);
    if (err != RABIN_SUCCESS) return err;
    return rmat_id(K, n);
  }

  rabin_err_t err = rmat_init(&V, A->rows, n);
  if (err != RABIN_SUCCESS) return err;
  err = rmat_copy(&V, A);
  if (err != RABIN_SUCCESS) goto cleanup;
  TRY(rmat_init(&H, n, n));
  TRY(rmat_id(&H, n));
  TRY(rmat_init(&lam, n, n));

  u64 dv_done = 0;
  dv = malloc((n + 1) * sizeof(rz_t));
  iota = calloc(n, sizeof(bool));
  if (dv == NULL || iota == NULL) {
    err = RABIN_ERR_OUT_OF_MEMORY;
    goto cleanup;
  }
  for (u64 x = 0; x <= n; x++) {
    if (rz_init(&dv[x]) != RABIN_SUCCESS) goto cleanup;
    dv_done++;
  }

  err = rmat_kernel_lll_run(&V, &H, &lam, dv, iota, n);
  if (err != RABIN_SUCCESS) goto cleanup;

  // r+1 is the first index with a nonzero flag: the first r columns of H
  // form a Z-basis of the kernel (they are independent)
  u64 r = 0;
  while (r < n && !iota[r]) r++;

  if (r == 0) {
    err = rmat_init(K, n, 0);
    goto cleanup;
  }
  err = rmat_extract_cols(&Kraw, &H, 0, r);
  if (err != RABIN_SUCCESS) goto cleanup;
  if (r > 1) {
    err = rmat_lll(&Kraw, r, 0.75, &Uk);
    if (err != RABIN_SUCCESS) goto cleanup;
  }
  rmat_swap(K, &Kraw);

cleanup:
  if (dv != NULL) {
    for (u64 x = 0; x < dv_done; x++) rz_clear(&dv[x]);
    free(dv);
  }
  free(iota);
  rmat_clear(&Kraw);
  rmat_clear(&Uk);
  rmat_clear(&lam);
  rmat_clear(&H);
  rmat_clear(&V);
  return err;
}

rabin_err_t rmat_kernel_image(const rmat_t* A, rmat_t* H, u64* p)
{
  if (A == NULL || H == NULL || p == NULL) return RABIN_ERR_NULL_PTR;
  *p = 0;

  const u64 n = A->cols;
  if (n == 0) return rmat_init(H, 0, 0);

  rmat_t W = {0}, U = {0}, Q = {0}, Uk = {0}, Kraw = {0}, Usub = {0};
  rmat_t Tail = {0}, Htmp = {0};
  rz_t det, one;
  bool det_live = false;

  // no coordinates: the whole Z^n is the kernel, the image is trivial
  if (A->rows == 0) {
    rabin_err_t err = rmat_init(H, n, n);
    if (err != RABIN_SUCCESS) return err;
    return rmat_id(H, n);
  }

  rabin_err_t err = hnf_with_transform(&W, &U, A);
  if (err != RABIN_SUCCESS) goto cleanup;
  const u64 rank = W.cols;
  const u64 r = n - rank;

  // LLL-reduce the HNF image basis W (in place), keeping the transform Q
  if (rank >= 1 && (err = rmat_lll(&W, rank, 0.75, &Q)) != RABIN_SUCCESS)
    goto cleanup;

  // LLL-reduce the kernel basis (columns 1..r of U, independent)
  if (r > 0) {
    err = rmat_extract_cols(&Kraw, &U, 0, r);
    if (err != RABIN_SUCCESS) goto cleanup;
    if (r > 1) err = rmat_lll(&Kraw, r, 0.75, &Uk);
    if (err != RABIN_SUCCESS) goto cleanup;
  }

  // image part: columns r..n-1 of U times Q
  if (rank >= 1) {
    err = rmat_extract_cols(&Usub, &U, r, rank);
    if (err == RABIN_SUCCESS) err = rmat_init(&Tail, n, rank);
    if (err == RABIN_SUCCESS) err = rmat_mul(&Tail, &Usub, &Q);
    if (err != RABIN_SUCCESS) goto cleanup;
  }

  // assemble H = [ K | U_tail * Q ] and commit only when every check passed
  TRY(rmat_init(&Htmp, n, n));
  for (u64 y = 0; y < n && err == RABIN_SUCCESS; y++) {
    for (u64 x = 0; x < r; x++) {
      err = rz_copy(RMAT_GET(&Htmp, y, x), RMAT_GET(&Kraw, y, x));
    }
    for (u64 x = 0; x < rank; x++) {
      if (err == RABIN_SUCCESS)
        err = rz_copy(RMAT_GET(&Htmp, y, r + x), RMAT_GET(&Tail, y, x));
    }
  }
  if (err != RABIN_SUCCESS) goto cleanup;

  // H must be unimodular
  TRY(rz_init_multi(&det, &one, NULL));
  det_live = true;
  TRY(rz_set_u64(&one, 1));
  TRY(rmat_det(&det, &Htmp));
  if (rz_cmp_abs(&det, &one) != 0) {
    err = RABIN_ERR_INVALID_ARG;
    goto cleanup;
  }

  // steps 3-5: keep the image columns small against the kernel basis
  TRY(rmat_reduce_against_kernel(&Htmp, n, r));

  rmat_swap(H, &Htmp);
  *p = rank;

cleanup:
  if (det_live) {
    rz_clear(&det);
    rz_clear(&one);
  }
  rmat_clear(&Htmp);
  rmat_clear(&Tail);
  rmat_clear(&Usub);
  rmat_clear(&Kraw);
  rmat_clear(&Uk);
  rmat_clear(&Q);
  rmat_clear(&U);
  rmat_clear(&W);
  return err;
}
