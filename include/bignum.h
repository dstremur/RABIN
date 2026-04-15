#ifndef BIGNUM_H
#define BIGNUM_H

#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>

/* Ideas
 Catalan pseudoprime

Hensel lifting

berlenkamp algo

LLL 

Discrete log problem

suntherlands algorithm

smith normal form

 */

typedef uint64_t u64;
typedef int64_t i64;

typedef struct bignum {
  uint64_t* limbs;  // 64-bit limbs (Base 2^64)
  size_t size;
  size_t capacity;
  bool is_neg;
} bignum;

typedef struct {
  bignum n;        /* modulus (odd, > 1)               */
  uint64_t n_inv;  /* -n^{-1} mod 2^64                  */
  bignum r_square; /* R^2 mod N   (R = 2^{64 * n.limbs})*/
  bignum one_mont; /* R   mod N   (representation of 1) */
} bn_mont_ctx;

typedef struct {
	u64 n; // transform size (power of 2)
	bignum q; // modulus
	bn_mont_ctx ctx;
	bignum* psi_pow;
	bignum* inv_psi_pow;
	bignum n_inv; 
} bn_ntt;

typedef struct {
	u64* limbs;
	u64 size;
} bn52;

void nntest(); 

#define MAX(a, b) ((a) > (b) ? (a) : (b));
#define MIN(a, b) ((a) < (b) ? (a) : (b)); 

#define MAX_LIMBS 512


void bn_mul_512(bignum* r, bignum* a, bignum* b);

void fftest();
void bn_provable_prime(bignum* p, u64 k);
// bignum.c
void bn_init(bignum* r);
void bn_init_multi(bignum* first, ...);
int bn_is_even(const bignum* a);
bool bn_alloc(bignum* r, u64 capacity);
void bn_init_val(bignum* n, const char* str);
char* bn_to_string(const bignum* n);
void bn_print(const bignum* n);
void bn_println(const bignum* n);
void bn_free(bignum* r);
void bn_free_multi(bignum* r, ...);
void bn_trim(bignum* r);
void bn_swap(bignum* a, bignum* b); 
int bn_get_bit(const bignum* a, int i);
bool bn_is_eq_i64(const bignum* n, i64 a);
void bn_copy(bignum* dest, const bignum* src);
void bn_set_u64(bignum* r, uint64_t val);
void bn_set_i64(bignum* r, int64_t val);
void bn_set_bit(bignum* a, int i);
void bn_clear_bit(bignum* a, int i);
int bn_bitlen(const bignum* a);
bool bn_is_zero(const bignum* a);
int bn_cmp(const bignum* a, const bignum* b);
int bn_bit_length(const bignum* a);

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
void bn_mul(bignum* r, const bignum* a, const bignum* b);
void bn_mul_raw(bignum* r, const bignum* a, const bignum* b);
void bn_mul_karatsuba(bignum* r, const bignum* a, const bignum* b);
void bn_mul_school(bignum* r, const bignum* a, const bignum* b);

// bigdiv.c
void bn_div(bignum* q, const bignum* a, const bignum* b);
void bn_newton_div(bignum* q, const bignum* a, const bignum* d);
void bn_div_knuth(bignum* q, const bignum* a, const bignum* b);
// bigmod.c
void bn_mod(bignum* r, const bignum* a, const bignum* n);
uint64_t bn_divmod_u64(bignum* q, const bignum* a, uint64_t d);
uint64_t mod_inverse_u64(uint64_t n);
uint64_t bn_mod_u64(const bignum* a, uint64_t d);

// bigexp.c
void bn_pow(bignum* r, const bignum* a, const bignum* b);
void bn_mod_exp(bignum* r, const bignum* a, const bignum* b, const bignum* m);
void bn_mont_exp(bignum* r_bar, const bignum* a_bar, const bignum* d, bn_mont_ctx* ctx);

void bn_mod_exp_mont(bignum* r, const bignum* a, const bignum* b, const bignum* m, bn_mont_ctx* ctx);

// bigshift.c
void bn_lshift1(bignum* r);
void bn_lshift1_add(bignum* r, int bit);
void bn_lshift(bignum* r, const bignum* a, int shift);
void bn_rshift(bignum* r, const bignum* a, int shift);
void bn_rshift1(bignum* r);
void bn_rshift(bignum* r, const bignum* a, int shift);
void bn_lshift(bignum* r, const bignum* a, int shift);
// bigrabin.c
bool bn_rabin(const bignum* n, const bignum* a);
bool bn_rabin_mont(const bignum* n, const bignum* a);

// bigrand.c 
bool bn_gen_random(bignum* r, u64 bits);
bool bn_gen_random_with_fd(bignum* r, u64 bits, int fd); 
void bn_gen_random_range(bignum* r, const bignum* low, const bignum* high);

// bigpseudo.c 
bool bn_gen_strps(bignum* p, u64 k); 

// bigmath.c
i64 bn_jacobi(const bignum* a, const bignum* m);
void tonelli_shanks(bignum* r, const bignum* n, const bignum* p);
void bn_gcd(bignum* d, const bignum* a, const bignum* b);
// biglucas.c
void bn_lucas(bignum* u, bignum* v, const bignum* p, const bignum* q, const bignum* n);
void bn_lucas_mod(bignum* u, bignum* v, const bignum* p, const bignum* q, const bignum* n,
                  const bignum* m);
void bn_lucas_solve_mod(bignum* u, bignum* v, const bignum* p, const bignum* q, bignum* qn, const bignum* n, const bignum* m);

// bigfactor.ctx
bool bn_pollard_rho(bignum* f, const bignum* n);

// bigprime.c
bool bn_bpsw(const bignum* n);

void bn_mont_ctx_init(bn_mont_ctx* ctx, const bignum* n);

void bn_mont_ctx_free(bn_mont_ctx* ctx);

void bn_mont_redc(bignum* r, bignum* t, bn_mont_ctx* ctx);

void bn_mont_mul(bignum* result, bignum* A_bar, bignum* B_bar,
                 bn_mont_ctx* ctx);
void bn_mont_mul_raw(bignum* result, bignum* A_bar, bignum* B_bar,
                     bn_mont_ctx* ctx);

void bn_mont_in(bignum* A_bar, const bignum* A, bn_mont_ctx* ctx);

void bn_mont_out(bignum* A, const bignum* A_bar, bn_mont_ctx* ctx);



void avxtest();

#endif
