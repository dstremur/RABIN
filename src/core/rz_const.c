/*
 * rz_const.c
 *
 * Global rz_t constants.
 *
 * This file defines the shared constant bignums RZ_ZERO, RZ_ONE, and
 * RZ_TWO, along with the functions to initialize and free them. The
 * constants must be initialized (rz_init_constants) before use and
 * freed (rz_clear_constants) on shutdown.
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

#include "../../include/rabin.h"

// Shared constant bignums. Initialized by rz_init_constants(), freed by
// rz_clear_constants(). Treat as read-only after initialization.
rz_t RZ_ZERO;
rz_t RZ_ONE;
rz_t RZ_TWO;

rabin_err_t rz_init_constants()
{
  rabin_err_t err = rz_init(&RZ_ZERO);
  if (err == RABIN_SUCCESS) err = rz_set_u64(&RZ_ZERO, 0);
  if (err == RABIN_SUCCESS) err = rz_init(&RZ_ONE);
  if (err == RABIN_SUCCESS) err = rz_set_u64(&RZ_ONE, 1);
  if (err == RABIN_SUCCESS) err = rz_init(&RZ_TWO);
  if (err == RABIN_SUCCESS) err = rz_set_u64(&RZ_TWO, 2);

  if (err != RABIN_SUCCESS) {
    rz_clear_multi(&RZ_ZERO, &RZ_ONE, &RZ_TWO, NULL);
  }
  return err;
}

rabin_err_t rz_clear_constants()
{
  return rz_clear_multi(&RZ_ZERO, &RZ_ONE, &RZ_TWO, NULL);
}
