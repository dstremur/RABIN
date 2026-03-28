#ifndef BIGNUM_H
#define BIGNUM_H

#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>

#include "stdint.h"

typedef struct bignum {
  uint64_t* limbs;  // 64-bit limbs (Base 2^64)
  size_t size;
  size_t capacity;
  bool is_neg;
} bignum;

typedef struct {
  bignum n;
  uint64_t n_inv;
  bignum r_square;
  bignum one_mont;
} bn_mont_ctx;

typedef uint64_t u64;
typedef int64_t i64;

#define MAX(a, b) ((a) > (b) ? (a) : (b));

#define MAX_LIMBS 512

// bignum.c
void bn_init(bignum* r);
int bn_is_even(const bignum* a);
bool bn_alloc(bignum* r, u64 capacity);
void bn_init_val(bignum* n, const char* str);
void bn_print(bignum* n);
void bn_free(bignum* r);
void bn_trim(bignum* r);
int bn_get_bit(const bignum* a, int i);
void bn_copy(bignum* dest, const bignum* src);
void bn_set_u64(bignum* r, uint64_t val);
void bn_set_bit(bignum* a, int i);
int bn_bitlen(const bignum* a);
bool bn_is_zero(const bignum* a);
int bn_cmp(const bignum* a, const bignum* b);
int bn_bit_length(const bignum* a);

// bigadd.c
void bn_add(bignum* r, const bignum* a, const bignum* b);

// bigsub.c
void bn_sub(bignum* r, const bignum* a, const bignum* b);
void bn_sub_abs(bignum* r, const bignum* a, const bignum* b);

// bigmul.c
void bn_mul(bignum* r, const bignum* a, const bignum* b);

// bigdiv.c
void bn_div(bignum* q, const bignum* a, const bignum* b);

// bigmod.c
void bn_mod(bignum* r, const bignum* a, const bignum* n);
uint64_t bn_divmod_u64(bignum* q, const bignum* a, uint64_t d);
uint64_t mod_inverse_u64(uint64_t n);

// bigexp.c
void bn_pow(bignum* r, bignum* a, bignum* b);

// bigshift.c
void bn_lshift1(bignum* r);
void bn_lshift1_add(bignum* r, int bit);
void bn_lshift(bignum* r, const bignum* a, int shift);
void bn_rshift(bignum* r, const bignum* a, int shift);
void bn_rshift1(bignum* r);
void bn_rshift(bignum* r, const bignum* a, int shift);

#endif
