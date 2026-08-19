/*===========================================================================
 *  bighelper.c
 *
 *  Raw-limb multiplication kernels and multi-limb add, moved here from
 *  bigmul.c so that all low-level limb helpers live in one place.
 *  See bighelper.h for the full list of helpers.
 *===========================================================================*/

#include "../../include/bighelper.h"

#include <string.h>

/*---------------------------------------------------------------------------
 *  add
 *-------------------------------------------------------------------------*/

u64 limbs_add_raw(u64* r, const u64* a, u64 a_len, const u64* b, u64 b_len)
{
  u64 carry;
  u64 max_len;

  // bn_add_inner requires the first operand to be the longer one
  if (a_len >= b_len) {
    carry = bn_add_inner(r, a, a_len, b, b_len);
    max_len = a_len;
  } else {
    carry = bn_add_inner(r, b, b_len, a, a_len);
    max_len = b_len;
  }

  // if carry present
  if (carry) {
    r[max_len] = carry;
    return max_len + 1;
  }
  return max_len;
}

/*---------------------------------------------------------------------------
 *  multiplication
 *-------------------------------------------------------------------------*/

// schoolbook multiplication directly on the limbs
void limbs_mul_school(u64* r, const u64* a, u64 a_size, const u64* b,
                      u64 b_size)
{
  memset(r, 0, (a_size + b_size) * sizeof(u64));
  for (u64 i = 0; i < a_size; i++) {
    if (a[i] == 0) continue;
    bn_mul_add_inner(&r[i], b, a[i], b_size);
  }
}

void limbs_mul_karatsuba(u64* r, const u64* a, const u64* b, u64 n,
                         u64* scratch)
{
  // if below limit, switch to school multiplication
  if (n < KARATSUBA_LIMIT) {
    limbs_mul_school(r, a, n, b, n);
    return;
  }

  // split number, upper can be larger by 1
  u64 m = n / 2;
  u64 high_len = n - m;

  // pointer arithmetic a = [low...][...high]
  const u64* a0 = a;
  const u64* a1 = a + m;
  const u64* b0 = b;
  const u64* b1 = b + m;

  // partition scratchpad memory, scratch = [s_a][s_b][z1][next_scratch]
  // s_a = a0 + a1, s_b = b0 + b1
  // z1 = (a0 + a1)(b0 + b1)
  u64* s_a = scratch;
  u64* s_b = scratch + (high_len + 1);
  u64* z1 = s_b + (high_len + 1);
  u64* next_scratch = z1 + (2 * high_len + 2);

  // store output buffer
  u64* z0 = r;
  u64* z2 = r + 2 * m;

  memset(r, 0, 2 * n * sizeof(u64));

  // recursive calls z0 = a0*b0, z1 = a1*b1
  limbs_mul_karatsuba(z0, a0, b0, m, next_scratch);
  limbs_mul_karatsuba(z2, a1, b1, high_len, next_scratch);

  // computing sums s_a = a0 + a1
  u64 a_len = limbs_add_raw(s_a, a0, m, a1, high_len);
  u64 b_len = limbs_add_raw(s_b, b0, m, b1, high_len);

  // normalize lengths, shorter is zero-padded
  u64 max_len = MAX(a_len, b_len);
  if (max_len > a_len) {
    memset(s_a + a_len, 0, (max_len - a_len) * sizeof(u64));
  }
  if (max_len > b_len) {
    memset(s_b + b_len, 0, (max_len - b_len) * sizeof(u64));
  }

  // recursive middle multiplication
  // z1 = s_a * s_b
  limbs_mul_karatsuba(z1, s_a, s_b, max_len, next_scratch);

  // intermediate subtraction step: z1 = z1 - z0 - z2
  bn_sub_inner(z1, z1, 2 * max_len, z0, 2 * m);
  bn_sub_inner(z1, z1, 2 * max_len, z2, 2 * high_len);

  // add middle term z1 shifted by m limbs = z1 * B^m
  bn_add_inner(r + m, r + m, 2 * n - m, z1, 2 * max_len);
}

void limbs_sqr_karatsuba(u64* r, const u64* a, u64 a_len, u64* scratch)
{
  if (a_len < KARATSUBA_LIMIT) {
    // Fine as a base case, though a dedicated sqr_school would be faster.
    limbs_mul_school(r, a, a_len, a, a_len);
    return;
  }

  u64 m = a_len / 2;
  u64 high_len = a_len - m;

  const u64* a0 = a;
  const u64* a1 = a + m;

  // scratch = [s_a][z1][next_scratch]
  // s_a = a0 + a1
  // z1  = (a0 + a1)^2
  u64* s_a = scratch;              // size: high_len + 1
  u64* z1 = s_a + (high_len + 1);  // size: 2*high_len + 2
  u64* next_scratch = z1 + (2 * high_len + 2);

  // output layout
  u64* z0 = r;
  u64* z2 = r + 2 * m;

  memset(r, 0, 2 * a_len * sizeof(u64));

  // z0 = a0^2
  limbs_sqr_karatsuba(z0, a0, m, next_scratch);

  // z2 = a1^2
  limbs_sqr_karatsuba(z2, a1, high_len, next_scratch);

  // s_a = a0 + a1
  u64 sa_len = limbs_add_raw(s_a, a0, m, a1, high_len);

  // z1 = (a0 + a1)^2
  limbs_sqr_karatsuba(z1, s_a, sa_len, next_scratch);

  // z1 = z1 - z0 - z2
  bn_sub_inner(z1, z1, 2 * sa_len, z0, 2 * m);
  bn_sub_inner(z1, z1, 2 * sa_len, z2, 2 * high_len);

  // add middle term shifted by m limbs
  bn_add_inner(r + m, r + m, 2 * a_len - m, z1, 2 * sa_len);
}
