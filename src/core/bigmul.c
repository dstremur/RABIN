#include <stddef.h>
#include <stdio.h>
#include <string.h>

#include "../../include/bignum.h"
#include "../../include/u64.h"
#define KARATSUBA_LIMIT 128

extern uint64_t bn_add_inner(uint64_t* r, const uint64_t* a, uint64_t a_size,
                             const uint64_t* b, uint64_t b_size);

extern uint64_t bn_sub_inner(uint64_t* r, const uint64_t* a, uint64_t a_size,
                             const uint64_t* b, uint64_t b_size);

extern void bn_mul_add_inner(u64* r, const u64* b, u64 a_limb, u64 len);

// schoolbook multiplication directly on the limbs
void limbs_mul_school(u64* r, const u64* a, const u64 a_size, const u64* b,
                      const u64 b_size)
{
  memset(r, 0, (a_size + b_size) * sizeof(u64));
  for (u64 i = 0; i < a_size; i++) {
    if (a[i] == 0) continue;
    bn_mul_add_inner(&r[i], b, a[i], b_size);
  }
}

// helper function to add directly on the limbs
static u64 limbs_add_raw(u64* r, const u64* a, const u64 a_len, const u64* b,
                         const u64 b_len)
{
  u64 carry;
  u64 max_len;
  // choose max_len
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
  for (u64 i = a_len; i < max_len; i++) s_a[i] = 0;
  for (u64 i = b_len; i < max_len; i++) s_b[i] = 0;

  // recursive middle multiplication
  // z1 = s_a * s_b
  limbs_mul_karatsuba(z1, s_a, s_b, max_len, next_scratch);

  // intermediate subtraction step: z1 = z1 - z0 - z2
  bn_sub_inner(z1, z1, 2 * max_len, z0, 2 * m);
  bn_sub_inner(z1, z1, 2 * max_len, z2, 2 * high_len);

  // add middle term z1 shifted by m limbs = z1 * B^m
  bn_add_inner(r + m, r + m, 2 * n - m, z1, 2 * max_len);
}

void bn_mul(bignum* r, const bignum* a, const bignum* b)
{
  if (a->size == 0 || b->size == 0) {
    r->size = 0;
    if (r->capacity > 0) r->limbs[0] = 0;
    return;
  }

  if (r == a || r == b) {
    bignum tmp;
    bn_init(&tmp);
    bn_mul(&tmp, a, b);
    bn_copy(r, &tmp);
    bn_free(&tmp);
    return;
  }

  u64 max_len = MAX(a->size, b->size);
  bn_alloc(r, 2 * max_len);

  // use karasuba if both inputs are larger then the limit
  if (a->size >= KARATSUBA_LIMIT && b->size >= KARATSUBA_LIMIT) {
    // Normalize asymmetric arrays with zero-padding up front
    uint64_t* pad_a = calloc(max_len, sizeof(uint64_t));
    uint64_t* pad_b = calloc(max_len, sizeof(uint64_t));
    memcpy(pad_a, a->limbs, a->size * sizeof(uint64_t));
    memcpy(pad_b, b->limbs, b->size * sizeof(uint64_t));

    // allocate scratchpad
    uint64_t scratch_size = 8 * max_len;
    uint64_t* scratch = malloc(scratch_size * sizeof(uint64_t));

    limbs_mul_karatsuba(r->limbs, pad_a, pad_b, max_len, scratch);

    free(scratch);
    free(pad_a);
    free(pad_b);
  } else {
    limbs_mul_school(r->limbs, a->limbs, a->size, b->limbs, b->size);
  }

  r->size = a->size + b->size;
  r->is_neg = a->is_neg ^ b->is_neg;
  bn_trim(r);
}

void bn_mul_school(bignum* r, const bignum* a, const bignum* b)
{
  u64 max = a->size + b->size;
  bn_alloc(r, max);

  memset(r->limbs, 0, max * sizeof(u64));
  r->size = max;
  r->is_neg = a->is_neg ^ b->is_neg;

  for (u64 i = 0; i < a->size; i++) {
    if (a->limbs[i] == 0) continue;

    bn_mul_add_inner(&r->limbs[i], b->limbs, a->limbs[i], b->size);
  }

  bn_trim(r);
}

void bn_mul_karatsuba(bignum* r, const bignum* a, const bignum* b)
{
  if (r == a || r == b) {
    bignum tmp;
    bn_init(&tmp);
    bn_mul_karatsuba(&tmp, a, b);
    bn_copy(r, &tmp);
    bn_free(&tmp);
    return;
  }

  if (a->size < KARATSUBA_LIMIT || b->size < KARATSUBA_LIMIT) {
    bn_mul(r, a, b);
    return;
  }

  u64 max = MAX(a->size, b->size);

  u64 m = max / 2;

  bignum a0, a1, b0, b1;

  bn_init_multi(&a0, &a1, &b0, &b1, NULL);

  if (a->size > m) {
    bn_alloc(&a0, m);
    memcpy(a0.limbs, a->limbs, m * sizeof(u64));
    a0.size = m;
    bn_trim(&a0);

    bn_alloc(&a1, a->size - m);
    memcpy(a1.limbs, a->limbs + m, (a->size - m) * sizeof(u64));
    a1.size = a->size - m;
    bn_trim(&a1);
  } else {
    bn_copy(&a0, a);
  }

  a0.is_neg = 0;
  a1.is_neg = 0;

  if (b->size > m) {
    bn_alloc(&b0, m);
    memcpy(b0.limbs, b->limbs, m * sizeof(u64));
    b0.size = m;
    bn_trim(&b0);

    bn_alloc(&b1, b->size - m);
    memcpy(b1.limbs, b->limbs + m, (b->size - m) * sizeof(u64));
    b1.size = b->size - m;
    bn_trim(&b1);
  } else {
    bn_copy(&b0, b);
  }

  b0.is_neg = 0;
  b1.is_neg = 0;

  bignum z0, z1, z2, s_a, s_b, tmp;

  bn_init_multi(&z0, &z1, &z2, &s_a, &s_b, &tmp, NULL);

  bn_mul(&z0, &a0, &b0);

  bn_mul(&z2, &a1, &b1);

  bn_add_abs(&s_a, &a0, &a1);
  bn_add_abs(&s_b, &b0, &b1);
  bn_mul(&z1, &s_a, &s_b);

  bn_sub_abs(&tmp, &z1, &z0);
  bn_sub_abs(&z1, &tmp, &z2);

  u64 max_len = a->size + b->size;
  bn_alloc(r, max_len);
  memset(r->limbs, 0, max_len * sizeof(u64));
  r->size = max_len;
  r->is_neg = a->is_neg ^ b->is_neg;

  bn_add_at_offset(r, &z0, 0);
  bn_add_at_offset(r, &z1, m);
  bn_add_at_offset(r, &z2, 2 * m);

  bn_free_multi(&a0, &a1, &b0, &b1, &z0, &z1, &z2, &s_a, &s_b, &tmp, NULL);

  bn_trim(r);
}
