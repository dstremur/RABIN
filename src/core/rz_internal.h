#ifndef RZ_INTERNAL_H
#define RZ_INTERNAL_H

/*
 * Internal helpers for the rz (integer) module.
 *
 * These are NOT part of the public RABIN API. Public mathematical
 * functions never allocate the outer rz_t struct (Caller Allocates);
 * they only grow the limb buffer of an already-initialized output
 * through rz_alloc().
 */

#include "../../include/rz.h"

/**
 * @brief Grow the limb buffer of an initialized rz_t in place.
 *
 * Grows exponentially (doubling capacity) unless a larger capacity is
 * requested. Newly allocated limbs are zeroed. The logical size is not
 * changed.
 *
 * @param[in,out] r rz_t whose limb buffer grows.
 * @param[in]     capacity Minimum number of limbs to hold.
 *
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR if r is NULL,
 *         or RABIN_ERR_OUT_OF_MEMORY if the realloc fails.
 */
rabin_err_t rz_alloc(rz_t* r, u64 capacity);

#endif
