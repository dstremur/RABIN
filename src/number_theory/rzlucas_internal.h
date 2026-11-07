/*
 * rzlucas_internal.h
 *
 * Private cross-TU glue for the Lucas sequence module. Not a public header.
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

#ifndef RZLUCAS_INTERNAL_H
#define RZLUCAS_INTERNAL_H

#include "../../include/rabin_errors.h"
#include "../../include/rz.h"

/*
 * Strong Lucas test driven by an existing Montgomery context; the test
 * never leaves the Montgomery domain. Requires n + 1 = d * 2^s with
 * d >= 1 odd, and ctx->n equal to n. Sets *res to the test verdict.
 */
rabin_err_t rz_strong_lucas_mont_ctx(bool* res, rz_mont_ctx* ctx,
                                     const rz_t* n, const rz_t* P,
                                     const rz_t* Q, const rz_t* d, u64 s);

#endif
