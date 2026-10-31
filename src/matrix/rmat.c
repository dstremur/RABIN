/*
 * rmat.c
 *
 * Matrix arithmetic over bignums: lifecycle, copying, element/row/column
 * access, structural predicates and printing.
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

#include "../../include/rmat.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../core/strparse_internal.h"

#define RMAT_PRINT_MAX_WIDTH 16
rabin_err_t rmat_init(rmat_t* M, u64 r, u64 c)
{
  if (M == NULL) return RABIN_ERR_NULL_PTR;

  // overflow-safe row * cols
  if (r != 0 && c > UINT64_MAX / r) return RABIN_ERR_OVERFLOW;

  // idempotent init: release any previous allocation
  if (M->data != NULL) {
    rmat_clear(M);
  }

  if (r * c == 0) {
    M->data = NULL;
    M->rows = r;
    M->cols = c;
    return RABIN_SUCCESS;
  }

  M->data = malloc(r * c * sizeof(rz_t));
  if (M->data == NULL) return RABIN_ERR_OUT_OF_MEMORY;

  M->rows = r;
  M->cols = c;

  for (u64 i = 0; i < r * c; i++) {
    rz_init(&M->data[i]);
  }
  return RABIN_SUCCESS;
}

rabin_err_t rmat_clear(rmat_t* A)
{
  if (A == NULL) return RABIN_ERR_NULL_PTR;
  if (!A->data) return RABIN_SUCCESS;

  for (u64 i = 0; i < A->rows * A->cols; i++) {
    rz_clear(&A->data[i]);
  }

  free(A->data);

  A->data = NULL;
  A->rows = 0;
  A->cols = 0;
  return RABIN_SUCCESS;
}

rabin_err_t rmat_print_python(const rmat_t* A)
{
  if (A == NULL) return RABIN_ERR_NULL_PTR;

  printf("[");  // Start outer list
  for (u64 i = 0; i < A->rows; i++) {
    printf("[");  // Start row
    for (u64 j = 0; j < A->cols; j++) {
      rz_print(RMAT_GET(A, i, j));
      if (j < A->cols - 1) printf(", ");
    }
    printf("]");  // End row
    if (i < A->rows - 1) printf(",\n ");
  }
  printf("]\n");  // End outer list
  return RABIN_SUCCESS;
}

rabin_err_t rmat_copy(rmat_t* R, const rmat_t* A)
{
  if (R == NULL || A == NULL) return RABIN_ERR_NULL_PTR;

  // check if sizes match
  if (A->cols != R->cols || A->rows != R->rows) return RABIN_ERR_MATRIX_DIM;

  for (u64 i = 0; i < A->cols * A->rows; i++) {
    rabin_err_t err = rz_copy(&R->data[i], &A->data[i]);
    if (err != RABIN_SUCCESS) return err;
  }
  return RABIN_SUCCESS;
}

rabin_err_t rmat_get(rz_t* R, const rmat_t* A, u64 r, u64 c)
{
  if (R == NULL || A == NULL) return RABIN_ERR_NULL_PTR;
  if (r >= A->rows || c >= A->cols) return RABIN_ERR_INVALID_ARG;

  return rz_copy(R, RMAT_GET(A, r, c));
}

rabin_err_t rmat_set(rmat_t* A, const rz_t* a, u64 r, u64 c)
{
  if (A == NULL || a == NULL) return RABIN_ERR_NULL_PTR;
  if (r >= A->rows || c >= A->cols) return RABIN_ERR_INVALID_ARG;

  return rz_copy(RMAT_GET(A, r, c), a);
}

rabin_err_t rmat_get_col(rvec_t* c, const rmat_t* A, u64 col)
{
  if (c == NULL || A == NULL) return RABIN_ERR_NULL_PTR;
  if (col >= A->cols) return RABIN_ERR_INVALID_ARG;
  if (c->size != A->rows) return RABIN_ERR_MATRIX_DIM;

  for (u64 i = 0; i < A->rows; i++) {
    rabin_err_t err = rvec_set(c, RMAT_GET(A, i, col), i);
    if (err != RABIN_SUCCESS) return err;
  }
  return RABIN_SUCCESS;
}

rabin_err_t rmat_get_row(rvec_t* r, const rmat_t* A, u64 row)
{
  if (r == NULL || A == NULL) return RABIN_ERR_NULL_PTR;
  if (row >= A->rows) return RABIN_ERR_INVALID_ARG;
  if (r->size != A->cols) return RABIN_ERR_MATRIX_DIM;

  for (u64 i = 0; i < A->cols; i++) {
    rabin_err_t err = rvec_set(r, RMAT_GET(A, row, i), i);
    if (err != RABIN_SUCCESS) return err;
  }
  return RABIN_SUCCESS;
}

bool rmat_is_square(const rmat_t* A)
{
  if (A == NULL) return false;
  return A->rows == A->cols;
}

bool rmat_is_zero(const rmat_t* A)
{
  if (A == NULL) return false;
  // single sequential pass over the flat row-major storage; rz_is_zero is
  // an O(1) size/limb inspection, so no temporary bignums are needed
  u64 total = A->rows * A->cols;
  for (u64 i = 0; i < total; i++) {
    if (!rz_is_zero(&A->data[i])) {
      return false;
    }
  }
  return true;
}

bool rmat_is_identity(const rmat_t* A)
{
  if (A == NULL) return false;
  // dimension short-circuit: a non-square matrix cannot be the identity,
  // bail out before touching any entry
  if (!rmat_is_square(A)) {
    return false;
  }

  u64 n = A->rows;
  for (u64 i = 0; i < n; i++) {
    rz_t* row = &A->data[i * A->cols];
    for (u64 j = 0; j < n; j++) {
      if (i == j) {
        // diagonal entry must be exactly 1 (rz_is_one rejects -1)
        if (!rz_is_one(&row[j])) {
          return false;
        }
      } else {
        // every off-diagonal entry must be 0
        if (!rz_is_zero(&row[j])) {
          return false;
        }
      }
    }
  }
  return true;
}

bool rmat_is_diagonal(const rmat_t* A)
{
  if (A == NULL) return false;
  // only the off-diagonal sub-regions can violate: in row i those are the
  // columns left of the diagonal, j < i (capped at cols), and the
  // columns right of it, j > i (empty once i >= cols). The diagonal
  // entry (i, i) itself is never inspected
  for (u64 i = 0; i < A->rows; i++) {
    rz_t* row = &A->data[i * A->cols];

    // left of the diagonal: j in [0, min(i, cols))
    u64 left = RZ_MIN(i, A->cols);
    for (u64 j = 0; j < left; j++) {
      if (!rz_is_zero(&row[j])) {
        return false;
      }
    }

    // right of the diagonal: j in (i, cols)
    for (u64 j = i + 1; j < A->cols; j++) {
      if (!rz_is_zero(&row[j])) {
        return false;
      }
    }
  }
  return true;
}

bool rmat_is_upper_triangular(const rmat_t* A)
{
  if (A == NULL) return false;
  // only the strict lower triangle (row > col) can violate; row 0 has no
  // entries below the diagonal, so the scan starts at i = 1
  for (u64 i = 1; i < A->rows; i++) {
    rz_t* row = &A->data[i * A->cols];

    // below the diagonal in row i: j in [0, min(i-1, cols-1)], i.e.
    // j < min(i, cols) — the cap handles tall (c < r) matrices
    u64 limit = RZ_MIN(i, A->cols);
    for (u64 j = 0; j < limit; j++) {
      if (!rz_is_zero(&row[j])) {
        return false;
      }
    }
  }
  return true;
}

bool rmat_is_lower_triangular(const rmat_t* A)
{
  if (A == NULL) return false;
  // only the strict upper triangle (col > row) can violate; in row i the
  // first candidate is column i+1, and the inner loop is empty once
  // i >= cols, so rows below the diagonal region cost nothing
  for (u64 i = 0; i < A->rows; i++) {
    rz_t* row = &A->data[i * A->cols];
    for (u64 j = i + 1; j < A->cols; j++) {
      if (!rz_is_zero(&row[j])) {
        return false;
      }
    }
  }
  return true;
}

bool rmat_is_symmetric(const rmat_t* A)
{
  if (A == NULL) return false;
  // dimension short-circuit: a non-square matrix cannot be symmetric,
  // bail out before touching any entry
  if (!rmat_is_square(A)) {
    return false;
  }

  u64 n = A->rows;
  // walk only the strict upper triangle: each mirror pair (i, j)/(j, i)
  // with i < j is compared exactly once; row i is swept sequentially,
  // and RMAT_GET(A, j, i) is a constant offset into row j's contiguous block
  for (u64 i = 0; i < n; i++) {
    rz_t* row_i = &A->data[i * A->cols];
    for (u64 j = i + 1; j < n; j++) {
      if (rz_cmp(&row_i[j], RMAT_GET(A, j, i)) != 0) {
        return false;
      }
    }
  }
  return true;
}

// upper bound on the printed width (digits plus sign) of a single entry,
// computed without a decimal conversion: ceil(bits * log10(2)) using the
// rational over-approximation 30103/100000 > log10(2), so the result is
// always at least the true digit count
static size_t rmat_entry_width(const rz_t* x)
{
  if (rz_is_zero(x)) return 1;

  size_t bits = (size_t)rz_bit_length(x);
  size_t w = (bits * 30103 + 99999) / 100000;
  if (x->is_neg) w++;
  return w;
}

static rabin_err_t rmat_print_impl(const rmat_t* A, bool full, bool tail)
{
  if (A == NULL) return RABIN_ERR_NULL_PTR;
  if (A->rows == 0 || A->cols == 0) return RABIN_SUCCESS;

  // per-column width: widest entry of the column, capped in the truncated
  // modes. A single outlier (e.g. the last invariant factor of a Smith form
  // sitting in a column of zeros) must not inflate the column, so in the
  // truncated modes outliers are excluded from the column width; they are
  // rendered at the capped width and simply overflow the column to the
  // right
  size_t cap = RMAT_PRINT_MAX_WIDTH;
  size_t* colw = calloc(A->cols, sizeof(size_t));
  if (colw == NULL) return RABIN_ERR_OUT_OF_MEMORY;
  for (u64 i = 0; i < A->rows; i++) {
    for (u64 j = 0; j < A->cols; j++) {
      size_t xw = rmat_entry_width(RMAT_GET(A, i, j));
      if (!full && xw > cap) continue;
      if (xw > colw[j]) colw[j] = xw;
    }
  }

  rabin_err_t err = RABIN_SUCCESS;
  for (u64 i = 0; i < A->rows && err == RABIN_SUCCESS; i++) {
    for (u64 j = 0; j < A->cols && err == RABIN_SUCCESS; j++) {
      char* s = rz_to_string(RMAT_GET(A, i, j));
      if (s == NULL) {
        err = RABIN_ERR_OUT_OF_MEMORY;
        break;
      }
      size_t len = strlen(s);
      size_t w = colw[j];

      if (full || len <= w) {
        printf("%*s", (int)w, s);
      } else if (len <= cap) {
        // fits within the cap: print it whole, overflowing by its own width
        printf("%s", s);
      } else if (tail) {
        size_t sign = (s[0] == '-') ? 1 : 0;
        if (sign) fputc('-', stdout);
        size_t keep = cap - 1 - sign;
        printf("\u2026%.*s", (int)keep, s + len - keep);
      } else {
        printf("%.*s\u2026", (int)(cap - 1), s);
      }

      free(s);
      printf(" ");
    }
    printf("\n");
  }

  free(colw);
  return err;
}

rabin_err_t rmat_print(const rmat_t* A)
{
  return rmat_print_impl(A, false, false);
}

rabin_err_t rmat_print_full(const rmat_t* A)
{
  return rmat_print_impl(A, true, false);
}

rabin_err_t rmat_print_tail(const rmat_t* A)
{
  return rmat_print_impl(A, false, true);
}

size_t rmat_to_string_len(const rmat_t* A)
{
  if (A == NULL) return 0;

  size_t len = 3; /* outer '[' + ']' + NUL */
  for (u64 r = 0; r < A->rows; r++) {
    if (r + 1 < A->rows) len++;
    len += 2; /* row '[' + ']' */
    for (u64 c = 0; c < A->cols; c++) {
      if (c + 1 < A->cols) len++;
      char* ns = rz_to_string(RMAT_GET(A, r, c));
      if (ns == NULL) return 0;
      len += strlen(ns);
      free(ns);
    }
  }
  return len;
}

rabin_err_t rmat_to_string(const rmat_t* A, char* out, size_t cap)
{
  if (A == NULL || out == NULL) return RABIN_ERR_NULL_PTR;

  size_t need = rmat_to_string_len(A);
  if (need == 0) return RABIN_ERR_OUT_OF_MEMORY;
  if (need > cap) return RABIN_ERR_INVALID_ARG;

  size_t w = 0;
  out[w++] = '[';
  for (u64 r = 0; r < A->rows; r++) {
    if (r > 0) out[w++] = ',';
    out[w++] = '[';
    for (u64 c = 0; c < A->cols; c++) {
      if (c > 0) out[w++] = ',';
      char* ns = rz_to_string(RMAT_GET(A, r, c));
      if (ns == NULL) return RABIN_ERR_OUT_OF_MEMORY;
      const size_t l = strlen(ns);
      memcpy(out + w, ns, l);
      w += l;
      free(ns);
    }
    out[w++] = ']';
  }
  out[w++] = ']';
  out[w] = '\0';
  return RABIN_SUCCESS;
}

rabin_err_t rmat_init_val(rmat_t* A, const char* str)
{
  if (A == NULL || str == NULL) return RABIN_ERR_NULL_PTR;

  strpx_lexer lx;
  strpx_init(&lx, str);
  if (!strpx_eat(&lx, '[')) return RABIN_ERR_INVALID_ARG;

  const u64 total = A->rows * A->cols;
  rz_t* tmp = NULL;
  u64 ninit = 0, rows_got = 0;
  rabin_err_t err = RABIN_SUCCESS;

  if (total > 0) {
    tmp = malloc(total * sizeof(rz_t));
    if (tmp == NULL) return RABIN_ERR_OUT_OF_MEMORY;
    for (u64 i = 0; i < total; i++) {
      if (rz_init(&tmp[i]) != RABIN_SUCCESS) {
        err = RABIN_ERR_OUT_OF_MEMORY;
        goto cleanup;
      }
      ninit++;
    }
  }

  bool closed = strpx_eat(&lx, ']'); /* empty row list */
  while (!closed) {
    if (rows_got == A->rows) {
      err = RABIN_ERR_MATRIX_DIM; /* more rows than fit */
      goto cleanup;
    }
    if (!strpx_eat(&lx, '[')) {
      err = RABIN_ERR_INVALID_ARG;
      goto cleanup;
    }

    u64 cols_got = 0;
    if (!strpx_eat(&lx, ']')) {
      for (;;) {
        if (cols_got == A->cols) {
          err = RABIN_ERR_MATRIX_DIM; /* more columns than fit */
          goto cleanup;
        }
        if (strpx_int(&lx, &tmp[rows_got * A->cols + cols_got]) !=
            RABIN_SUCCESS) {
          err = RABIN_ERR_INVALID_ARG;
          goto cleanup;
        }
        cols_got++;
        if (strpx_eat(&lx, ',')) continue;
        break;
      }
      if (!strpx_eat(&lx, ']')) {
        err = RABIN_ERR_INVALID_ARG;
        goto cleanup;
      }
    }
    if (cols_got != A->cols) {
      err = RABIN_ERR_MATRIX_DIM;
      goto cleanup;
    }

    rows_got++;
    if (strpx_eat(&lx, ',')) continue;
    if (!strpx_eat(&lx, ']')) {
      err = RABIN_ERR_INVALID_ARG;
      goto cleanup;
    }
    closed = true;
  }

  if (rows_got != A->rows) {
    err = RABIN_ERR_MATRIX_DIM;
    goto cleanup;
  }
  if (!strpx_done(&lx)) {
    err = RABIN_ERR_INVALID_ARG;
    goto cleanup;
  }

  for (u64 i = 0; i < total; i++) {
    TRY(rz_copy(&A->data[i], &tmp[i]));
  }

cleanup:
  if (tmp != NULL) {
    for (u64 i = 0; i < ninit; i++) rz_clear(&tmp[i]);
    free(tmp);
  }
  return err;
}
