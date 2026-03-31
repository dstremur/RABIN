#ifndef BIGNUM_H
#define BIGNUM_H

#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>

typedef struct bignum {
  uint64_t* limbs;  // 64-bit limbs (Base 2^64)
  size_t size;
  size_t capacity;
  bool is_neg;
} bignum;

typedef struct {
    bignum n;          /* modulus (odd, > 1)               */
    uint64_t n_inv;    /* -n^{-1} mod 2^64                  */
    bignum r_square;   /* R^2 mod N   (R = 2^{64 * n.limbs})*/
    bignum one_mont;   /* R   mod N   (representation of 1) */
} bn_mont_ctx;

typedef uint64_t u64;
typedef int64_t i64;

#define MAX(a, b) ((a) > (b) ? (a) : (b));

#define MAX_LIMBS 512

// bignum.c
void bn_init(bignum* r);
void bn_init_multi(bignum* first, ...);
int bn_is_even(const bignum* a);
bool bn_alloc(bignum* r, u64 capacity);
void bn_init_val(bignum* n, const char* str);
void bn_print(bignum* n);
void bn_free(bignum* r);
void bn_trim(bignum* r);
int bn_get_bit(const bignum* a, int i);
bool bn_is_eq_i64(const bignum* n, i64 a);
void bn_copy(bignum* dest, const bignum* src);
void bn_set_u64(bignum* r, uint64_t val);
void bn_set_i64(bignum* r, int64_t val); 
void bn_set_bit(bignum* a, int i);
int bn_bitlen(const bignum* a);
bool bn_is_zero(const bignum* a);
int bn_cmp(const bignum* a, const bignum* b);
int bn_bit_length(const bignum* a);
bool bn_gen_random(bignum* r, int bits);
bool bn_gen_prime(bignum* p, int bits);
u64 bn_cnt_trailing_zeros(const bignum* a);
// bigadd.c
void bn_add(bignum* r, const bignum* a, const bignum* b);
void bn_add_u64(bignum* r, const bignum* a, u64 b); 
void bn_add_at_offset(bignum* r, const bignum* a, u64 offset);
// bigsub.c
void bn_sub(bignum* r, const bignum* a, const bignum* b);
void bn_sub_abs(bignum* r, const bignum* a, const bignum* b);

// bigmul.c
void bn_mul(bignum* r, bignum* a, bignum* b);
void bn_mul_karatsuba(bignum* r, bignum* a, bignum* b);
void bn_mul_school(bignum* r, bignum* a, bignum* b); 

// bigdiv.c
void bn_div(bignum* q, const bignum* a, const bignum* b);

// bigmod.c
void bn_mod(bignum* r, const bignum* a, const bignum* n);
uint64_t bn_divmod_u64(bignum* q, const bignum* a, uint64_t d);
uint64_t mod_inverse_u64(uint64_t n);

// bigexp.c
void bn_pow(bignum* r, bignum* a, bignum* b);
void bn_mod_exp(bignum* r, bignum* a, bignum* b, bignum* m); 

// bigshift.c
void bn_lshift1(bignum* r);
void bn_lshift1_add(bignum* r, int bit);
void bn_lshift(bignum* r, const bignum* a, int shift);
void bn_rshift(bignum* r, const bignum* a, int shift);
void bn_rshift1(bignum* r);
void bn_rshift(bignum* r, const bignum* a, int shift);
void bn_lshift(bignum* r, const bignum* a, int shift);
//bigrabin.c 
bool bn_rabin(bignum* n, bignum* a); 

// bigmath.c 
i64 bn_jacobi(bignum* a, bignum* m);








void bn_mont_ctx_init(bn_mont_ctx *ctx, const bignum *N);

void bn_mont_ctx_free(bn_mont_ctx *ctx);

void bn_mont_redc(bn_mont_ctx *ctx, bignum *t);

void bn_mont_mul(const bn_mont_ctx *ctx,
                 const bignum *x,
                 const bignum *y,
                 bignum *r);

void bn_mont_in(const bn_mont_ctx *ctx,
                const bignum *x,
                bignum *r);

void bn_mont_out(const bn_mont_ctx *ctx,
                 const bignum *x,
                 bignum *r);

void bn_mod_exp_mont(bignum* r, const bignum* a, const bignum* b, const bignum* m, bn_mont_ctx* ctx);

#endif
