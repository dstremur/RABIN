#include "../include/bignum.h"

#include <ctype.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>

#define BASE_10_19 10000000000000000000ULL

// Clear oject
void bn_init(bignum* r) {
  r->limbs = NULL;
  r->size = 0;
  r->capacity = 0;
  r->is_neg = false;
}

int bn_is_even(const bignum* n) {
	if (n->size == 0 || n->limbs == NULL) {
        return 1; 
    }

	return (n->limbs[0] & 1) == 0; 

}

void bn_init_val(bignum* n, const char* str) {
  bn_init(n);

  n->is_neg = (str[0] == '-');
  const char* s = n->is_neg ? str + 1 : str;

  // start with zero
  n->limbs = calloc(1, sizeof(u64));
  n->size = 1;
  n->capacity = 1;

  for (size_t i = 0; s[i]; i++) {
    int digit = s[i] - '0';

    // n = n * 10
    u64 carry = 0;
    for (u64 j = 0; j < n->size; j++) {
      unsigned __int128 prod = (unsigned __int128)n->limbs[j] * 10 + carry;

      n->limbs[j] = (u64)prod;
      carry = (u64)(prod >> 64);
    }

    if (carry) {
      n->limbs = realloc(n->limbs, (n->size + 1) * sizeof(u64));
      n->limbs[n->size++] = carry;
    }

    // n = n + digit
    unsigned __int128 sum = (unsigned __int128)n->limbs[0] + digit;

    n->limbs[0] = (u64)sum;
    carry = (u64)(sum >> 64);

    u64 k = 1;
    while (carry && k < n->size) {
      unsigned __int128 s = (unsigned __int128)n->limbs[k] + carry;

      n->limbs[k] = (u64)s;
      carry = (u64)(s >> 64);
      k++;
    }

    if (carry) {
      n->limbs = realloc(n->limbs, (n->size + 1) * sizeof(u64));
      n->limbs[n->size++] = carry;
    }
  }

  bn_trim(n);
}

uint64_t bn_mod_u64(const bignum* a, uint64_t d) {}

// Fast modular inverse for a single 64-bit limb (Newton's method)
uint64_t mod_inverse_u64(uint64_t n) {
  uint64_t inv = n;
  for (int i = 0; i < 5; i++) {
    inv *= (2 - n * inv);
  }
  return -inv;
}



// print in base 10
void bn_print(bignum* n) {
  if (n->size == 0 || (n->size == 1 && n->limbs[0] == 0)) {
    printf("0\n");
    return;
  }

  if (n->is_neg) printf("-");

  bignum tmp;
  bn_init(&tmp);

  // copy n → tmp
  tmp.size = n->size;
  tmp.capacity = n->size;
  tmp.limbs = malloc(tmp.size * sizeof(u64));
  memcpy(tmp.limbs, n->limbs, tmp.size * sizeof(u64));

  const uint64_t BASE10 = 10000000000000000000ULL;  // 10^19

  uint64_t* parts = NULL;
  size_t parts_count = 0;

  while (!(tmp.size == 1 && tmp.limbs[0] == 0)) {
    bignum q;
    bn_init(&q);

    uint64_t rem = bn_divmod_u64(&q, &tmp, BASE10);

    parts = realloc(parts, (parts_count + 1) * sizeof(uint64_t));
    parts[parts_count++] = rem;

    bn_free(&tmp);
    tmp = q;
  }

  // print most significant chunk normally
  printf("%llu", parts[parts_count - 1]);

  // remaining chunks padded with leading zeros
  for (i64 i = parts_count - 2; i >= 0; i--) {
    printf("%019llu", parts[i]);
  }

  printf("\n");

  free(parts);
  bn_free(&tmp);
}

void bn_alloc(bignum* r, size_t capacity) {}

void bn_free(bignum* r) {
  if (r->limbs) {
    free(r->limbs);
    r->limbs = NULL;
  }
  r->size = 0;
  r->capacity = 0;
}

void bn_trim(bignum* r) {
  while (r->size > 1 && r->limbs[r->size - 1] == 0) {
    r->size--;
  }
}

// Returns the value of the i-th bit (0 or 1)
int bn_get_bit(const bignum* a, int i) {
  int limb = i / 64;
  int offset = i % 64;

  if (limb >= a->size) return 0;

  return (a->limbs[limb] >> offset) & 1;
}

// Extracts 'n' bits starting from index 'i' (for windowed exp)
uint32_t bn_get_bits(const bignum* a, int i, int n) {}



// Compare two bignums: returns 1 if a > b, -1 if a < b, 0 if a == b
int bn_cmp(const bignum* a, const bignum* b) {
  if (a->size != b->size) return (a->size < b->size) ? -1 : 1;

  for (i64 i = a->size - 1; i >= 0; i--) {
    if (a->limbs[i] != b->limbs[i]) return (a->limbs[i] < b->limbs[i]) ? -1 : 1;
  }
  return 0;
}
void bn_copy(bignum* r, const bignum* a) {
	    if (r == a) return;

    bn_free(r);

    r->size = a->size;
    r->capacity = a->size;
    r->is_neg = a->is_neg;

    r->limbs = malloc(r->size * sizeof(u64));
    memcpy(r->limbs, a->limbs, r->size * sizeof(u64));
}



// Set to 64 bit unsigned integer
void bn_set_u64(bignum* n, uint64_t val) {
  bn_free(n);  // free any existing limbs
  n->limbs = calloc(1, sizeof(u64));
  n->limbs[0] = val;
  n->size = 1;
  n->capacity = 1;
  n->is_neg = false;
}
void bn_set_bit(bignum* a, int i) {
  int limb = i / 64;
  int offset = i % 64;

  if (limb >= a->size) {
    a->limbs = realloc(a->limbs, (limb + 1) * sizeof(u64));
    for (u64 j = a->size; j <= limb; j++) {
      a->limbs[j] = 0;
    }
    a->size = limb + 1;
  }

  a->limbs[limb] |= ((u64)1 << offset);
}

// Bit length of a bignum (useful for division)
int bn_bit_length(const bignum* a) {
  if (a->size == 0) return 0;

  u64 top = a->limbs[a->size - 1];
  int bits = 0;

  while (top) {
    top >>= 1;
    bits++;
  }

  return (a->size - 1) * 64 + bits;
}


// Left shift: r = a << shift
void bn_lshift(bignum* r, const bignum* a, int shift) {}
// Modulo: r = a % n (using binary long division)

// r = (a * b) % n
void bn_mod_mul(bignum* r, const bignum* a, const bignum* b, const bignum* n) {}

// r = (base ^ exp) % n
void bn_mod_exp(bignum* r, const bignum* base, const bignum* exp,
                const bignum* n) {}

void bn_mont_init(bn_mont_ctx* ctx, const bignum* n) {}

void bn_mont_free(bn_mont_ctx* ctx) {}

void bn_redc(bignum* r, uint64_t* T, const bn_mont_ctx* ctx) {}

// r = (base ^ exp) % n using Montgomery Reduction
void bn_mod_exp_mont(bignum* r, const bignum* base, const bignum* exp,
                     const bn_mont_ctx* ctx) {}

int bn_jacobi(bignum* a_in, const bignum* n_in) {}

// Helper for modular subtraction: (a - b) % n
void bn_sub_mod(bignum* r, const bignum* a, const bignum* b, const bignum* n) {}

