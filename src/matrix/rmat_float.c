/*
 * rmat_float.c
 *
 * Floating-point entry detection: approximate linear dependence and algebraic
 * relation search via LLL.
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

#include "../../include/rmat.h"
#include "../../include/rzmath.h"
/* ----------------- Dependence relations (Cohen 2.7.4) -------------- */

/* round scale*z to an integer, rejecting values outside the int64 range */
static rabin_err_t rz_round_scaled(rz_t* q, double z, double scale)
{
  if (!isfinite(z) || !(scale > 0.0) || !isfinite(scale) ||
      fabs(z) > 4.0e18 / scale)
    return RABIN_ERR_INVALID_ARG;
  return rz_set_i64(q, llround(scale * z));
}

/* best-effort double conversion of an integer, used for error residuals */
static double rz_as_double(const rz_t* z)
{
  double d = 0.0;
  for (size_t i = z->size; i > 0; i--) {
    d = d * 18446744073709551616.0 + (double)z->limbs[i - 1];
  }
  return z->is_neg ? -d : d;
}

rabin_err_t rmat_linear_dependence(const double* re, const double* im, u64 n,
                                   double scale, rmat_t* REL, double* err)
{
  if (re == NULL || REL == NULL) return RABIN_ERR_NULL_PTR;
  if (n < 2) return RABIN_ERR_INVALID_ARG;
  if (!(scale > 0.0) || !isfinite(scale)) return RABIN_ERR_INVALID_ARG;

  const bool complex_in = (im != NULL);
  for (u64 i = 0; i < n; i++) {
    if (!isfinite(re[i]) || (complex_in && !isfinite(im[i])))
      return RABIN_ERR_INVALID_ARG;
  }

  const u64 m = n + (complex_in ? 2 : 1);
  rmat_t M = {0};
  rabin_err_t e = rmat_init(&M, m, n);
  if (e != RABIN_SUCCESS) return e;

  // symmetric construction: rows 1..n hold the identity, the last rows the
  // scaled real (and imaginary) parts
  for (u64 i = 0; i < n && e == RABIN_SUCCESS; i++) {
    for (u64 j = 0; j < n; j++) {
      e = rz_set_u64(RMAT_GET(&M, i, j), i == j ? 1 : 0);
    }
  }
  for (u64 j = 0; j < n && e == RABIN_SUCCESS; j++) {
    e = rz_round_scaled(RMAT_GET(&M, n, j), re[j], scale);
  }
  if (e == RABIN_SUCCESS && complex_in) {
    for (u64 j = 0; j < n; j++) {
      e = rz_round_scaled(RMAT_GET(&M, n + 1, j), im[j], scale);
      if (e != RABIN_SUCCESS) break;
    }
  }

  // rmat_lll reduces M in place and returns the transformation matrix
  rmat_t Q = {0};
  if (e == RABIN_SUCCESS) e = rmat_lll(&M, n, 0.75, &Q);
  if (e == RABIN_SUCCESS) {
    e = rmat_init(REL, m, n);
    if (e == RABIN_SUCCESS) e = rmat_copy(REL, &M);
  }
  rmat_clear(&Q);
  if (e != RABIN_SUCCESS) {
    rmat_clear(&M);
    return e;
  }

  if (err != NULL) {
    // residual of every candidate relation: |sum_j REL[j][c] * z_j|
    for (u64 c = 0; c < n; c++) {
      long double sr = 0.0L, si = 0.0L;
      for (u64 j = 0; j < n; j++) {
        const double cj = rz_as_double(RMAT_GET(REL, j, c));
        sr += (long double)cj * (long double)re[j];
        if (complex_in) si += (long double)cj * (long double)im[j];
      }
      err[c] = (double)sqrtl(sr * sr + si * si);
    }
  }

  rmat_clear(&M);
  return RABIN_SUCCESS;
}

rabin_err_t rmat_algebraic_relation(double re, double im, u64 degree,
                                    double scale, rmat_t* poly, double* err)
{
  if (poly == NULL) return RABIN_ERR_NULL_PTR;
  if (degree < 1 || degree > 64) return RABIN_ERR_INVALID_ARG;
  if (!isfinite(re) || !isfinite(im) || !(scale > 0.0) || !isfinite(scale))
    return RABIN_ERR_INVALID_ARG;

  const u64 n = degree + 1;
  const bool complex_in = (im != 0.0);

  double* zre = malloc(n * sizeof(double));
  double* zim = complex_in ? malloc(n * sizeof(double)) : NULL;
  if (zre == NULL || (complex_in && zim == NULL)) {
    free(zre);
    free(zim);
    return RABIN_ERR_OUT_OF_MEMORY;
  }

  // z_i = a^i (1-based i, powers of the number)
  zre[0] = 1.0;
  if (complex_in) zim[0] = 0.0;
  for (u64 i = 1; i < n; i++) {
    if (complex_in) {
      const double pr = zre[i - 1], pi = zim[i - 1];
      zre[i] = pr * re - pi * im;
      zim[i] = pr * im + pi * re;
    } else {
      zre[i] = zre[i - 1] * re;
    }
    if (!isfinite(zre[i]) || (complex_in && !isfinite(zim[i]))) {
      free(zre);
      free(zim);
      return RABIN_ERR_INVALID_ARG;
    }
  }

  const rabin_err_t e = rmat_linear_dependence(zre, zim, n, scale, poly, err);
  free(zre);
  free(zim);
  return e;
}
