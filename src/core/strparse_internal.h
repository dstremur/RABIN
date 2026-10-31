/*
 * strparse_internal.h
 *
 * Shared private lexer for the bracketed decimal-list grammars used by
 * rvec, rmat and rpol parsing. Not a public header.
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

#ifndef STRPARSE_INTERNAL_H
#define STRPARSE_INTERNAL_H

#include <stdbool.h>
#include <stddef.h>

#include "../../include/rabin_errors.h"
#include "../../include/rz.h"

/*
 * Token grammar shared by the container parsers:
 *   ws   := ' ' | '\t' | '\n' | '\r'   (allowed between tokens only)
 *   int  := ['+' | '-'] digit+         (signed decimal, no leading forms)
 * Elements must be separated by ',' with no trailing comma or empty slot.
 */
typedef struct
{
  const char* s;
  size_t i;
  size_t n;
} strpx_lexer;

void strpx_init(strpx_lexer* lx, const char* str);

/* Advance past whitespace. */
void strpx_skip_ws(strpx_lexer* lx);

/* True when every remaining character is whitespace. */
bool strpx_done(strpx_lexer* lx);

/* Consume c when it is the next non-whitespace token; false otherwise. */
bool strpx_eat(strpx_lexer* lx, char c);

/*
 * Parse one int token into out. RABIN_ERR_INVALID_ARG on a missing or
 * malformed token, RABIN_ERR_OUT_OF_MEMORY on allocation failure.
 */
rabin_err_t strpx_int(strpx_lexer* lx, rz_t* out);

#endif /* STRPARSE_INTERNAL_H */
