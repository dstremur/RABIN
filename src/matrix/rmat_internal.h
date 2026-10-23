/*
 * rmat_internal.h
 *
 * Private cross-TU glue for the matrix module. Not a public header.
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

#ifndef RMAT_INTERNAL_H
#define RMAT_INTERNAL_H

#include "../../include/rmat.h"

/* scratch used by the RED/SWAP steps */
typedef struct {
  rz_t u, v, w;  // general temporaries
  rz_t one, two;
} lll_scr_t;

rabin_err_t lll_scr_init(lll_scr_t* s);
void        lll_scr_clear(lll_scr_t* s);
rabin_err_t lll_dot(rz_t* r, const rmat_t* V, u64 a, u64 b, rz_t* t);
rabin_err_t lll_red1(rmat_t* V, rmat_t* H, rmat_t* lam, const rz_t* dv,
                     u64 k1, u64 k2, lll_scr_t* s);


/* shared between rmat_hnf.c and the LLL/kernel TUs */
rabin_err_t hnf_euclid_step(rmat_t* A_work, rmat_t* U, rz_t* buf, rz_t* C,
                           u64 i, u64 j, u64 k, rz_t* u, rz_t* v, rz_t* d,
                           rz_t* q_k, rz_t* q_j, rz_t* temp, rz_t* temp2);
rabin_err_t hnf_with_transform(rmat_t* W, rmat_t* U, const rmat_t* A);

#endif  // RMAT_INTERNAL_H
