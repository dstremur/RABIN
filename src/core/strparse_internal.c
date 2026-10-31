/*
 * strparse_internal.c
 *
 * Shared lexer primitives for the rvec, rmat and rpol string grammars.
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

#include "strparse_internal.h"

#include <stdlib.h>
#include <string.h>

static bool strpx_is_ws(char c)
{
  return c == ' ' || c == '\t' || c == '\n' || c == '\r';
}

static bool strpx_is_digit(char c) { return c >= '0' && c <= '9'; }

void strpx_init(strpx_lexer* lx, const char* str)
{
  lx->s = str;
  lx->i = 0;
  lx->n = strlen(str);
}

void strpx_skip_ws(strpx_lexer* lx)
{
  while (lx->i < lx->n && strpx_is_ws(lx->s[lx->i])) lx->i++;
}

bool strpx_done(strpx_lexer* lx)
{
  strpx_skip_ws(lx);
  return lx->i == lx->n;
}

bool strpx_eat(strpx_lexer* lx, char c)
{
  strpx_skip_ws(lx);
  if (lx->i < lx->n && lx->s[lx->i] == c) {
    lx->i++;
    return true;
  }
  return false;
}

rabin_err_t strpx_int(strpx_lexer* lx, rz_t* out)
{
  strpx_skip_ws(lx);

  bool neg = false;
  if (lx->i < lx->n && (lx->s[lx->i] == '+' || lx->s[lx->i] == '-')) {
    neg = (lx->s[lx->i] == '-');
    lx->i++;
  }

  const size_t ds = lx->i;
  while (lx->i < lx->n && strpx_is_digit(lx->s[lx->i])) lx->i++;
  if (lx->i == ds) return RABIN_ERR_INVALID_ARG;

  char* tmp = strndup(lx->s + ds, lx->i - ds);
  if (tmp == NULL) return RABIN_ERR_OUT_OF_MEMORY;

  rabin_err_t err = rz_init_val(out, tmp);
  free(tmp);
  if (err != RABIN_SUCCESS) return err;

  if (neg) return rz_neg(out, out);
  return RABIN_SUCCESS;
}
