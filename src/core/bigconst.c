/*
 * bigconst.c
 *
 * Global bignum constants.
 *
 * This file defines the shared constant bignums BN_ZERO, BN_ONE, and
 * BN_TWO, along with the functions to initialize and free them. The
 * constants must be initialized (bn_init_constants) before use and
 * freed (bn_free_constants) on shutdown.
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

#include "../../include/bignum.h"

// Shared constant bignums. Initialized by bn_init_constants(), freed by
// bn_free_constants(). Treat as read-only after initialization.
bignum BN_ZERO;
bignum BN_ONE;
bignum BN_TWO;

void bn_init_constants()
{
  bn_init(&BN_ZERO);
  bn_set_u64(&BN_ZERO, 0);

  bn_init(&BN_ONE);
  bn_set_u64(&BN_ONE, 1);

  bn_init(&BN_TWO);
  bn_set_u64(&BN_TWO, 2);
}

void bn_free_constants() { bn_free_multi(&BN_ZERO, &BN_ONE, &BN_TWO, NULL); }
