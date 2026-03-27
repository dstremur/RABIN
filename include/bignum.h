#ifndef BIGNUM_H
#define BIGNUM_H

#include "stdint.h"
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>

typedef struct {
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

void bn_init_val(bignum* n, const char* str);
void bn_print(bignum* n); 
void bn_rshift(bignum* r, const bignum* a, int shift);
int bn_is_even(const bignum* a);
void bn_init(bignum* r);
void bn_alloc(bignum* r, size_t capacity);
void bn_free(bignum* r);
void bn_copy(bignum* dest, const bignum* src);

void bn_div(bignum* q, const bignum* a, const bignum* b);

void bn_add(bignum* r, const bignum* a, const bignum* b);
void bn_mul(bignum* r, const bignum* a, const bignum* b);

void bn_set_u64(bignum* r, uint64_t val);
void bn_print_hex(const char* label, const bignum* a);
int bn_bit_length(const bignum* a);
void bn_lshift1(bignum* r);
void bn_lshift1_add(bignum* r, int bit);

int bn_cmp(const bignum* a, const bignum* b);
void bn_copy(bignum* r, const bignum* a);

void bn_sub(bignum* r, const bignum* a, const bignum* b);
void bn_sub_abs(bignum* r, const bignum* a, const bignum* b);
uint64_t bn_divmod_u64(bignum* q, const bignum* a, uint64_t d);
int bn_bitlen(const bignum* a);
void bn_lshift(bignum* r, const bignum* a, int shift);
void bn_rshift(bignum* r, const bignum* a, int shift);
void bn_mod(bignum* r, const bignum* a, const bignum* n);
uint64_t bn_mod_u64(const bignum* a, uint64_t d);

void bn_mod_mul(bignum* r, const bignum* a, const bignum* b, const bignum* n);

void bn_mod_exp(bignum* r, const bignum* base, const bignum* exp,
                const bignum* n);

int bn_get_bit(const bignum* a, int i);
uint32_t bn_get_bits(const bignum* a, int i, int n);

// Update this signature
bool bn_millerRabin2(const bignum* n, const bignum* base, const bignum* d, int s, const bn_mont_ctx* ctx);
bool bn_millerRabin(const bignum* n, uint64_t base);


static inline bool bn_is_zero(const bignum* a) {
  return (a->size == 0 || (a->size == 1 && a->limbs[0] == 0));
}

void bn_set_bit(bignum* a, int i);
void bn_mont_init(bn_mont_ctx* ctx, const bignum* n);
void bn_mont_free(bn_mont_ctx* ctx);
void bn_redc_raw(uint64_t* res, uint64_t* T, const bn_mont_ctx* ctx);

// Utility
void bn_trim(bignum* r);

// --- Function Prototypes ---

// Core Montgomery Operations
uint64_t mod_inverse_u64(uint64_t n);
void bn_redc(bignum* r, uint64_t* T, const bn_mont_ctx* ctx);

// Exponentiation
void bn_mod_exp_mont(bignum* r, const bignum* base, const bignum* exp,
                     const bn_mont_ctx* ctx);
void bn_mod_exp_mont_window(bignum* r, const bignum* base, const bignum* exp,
                            const bn_mont_ctx* ctx);

// Raw Limb Operations (No struct overhead)
void bn_mul_raw(uint64_t* res, const uint64_t* a, size_t an, const uint64_t* b,
                size_t bn);
void bn_redc_raw(
    uint64_t* r, uint64_t* T,
    const bn_mont_ctx* ctx);  // Added: Called in C code but missing prototype

// Specialized / Hardware-Accelerated Functions
void lucas_double_mont(uint64_t* V, uint64_t* Qk, const bn_mont_ctx* ctx);
void bn_add512_unrolled(uint64_t* dest, const uint64_t* src);
void bn_redc_2048_fast(uint64_t* T, const bn_mont_ctx* ctx);

int bn_jacobi(bignum* a_in, const bignum* n);

/**
 * Performs the Strong Lucas Primality Test on n.
 * Returns true if n is a probable prime, false if composite.
 */
bool bn_strongLucasTest(const bignum* n);

/**
 * Modular subtraction: r = (a - b) % n.
 * Handles cases where a < b by adding n.
 */
void bn_sub_mod(bignum* r, const bignum* a, const bignum* b, const bignum* n);

bool bn_millerRabin_optimized(const bignum* n, const bignum* base, const bn_mont_ctx* ctx);

#endif

