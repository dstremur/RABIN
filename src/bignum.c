#include "../include/bignum.h"

#include <ctype.h>
#include <immintrin.h>
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
void bn_rshift(bignum* r, const bignum* a, int shift) {}

// Fast modular inverse for a single 64-bit limb (Newton's method)
uint64_t mod_inverse_u64(uint64_t n) {
  uint64_t inv = n;
  for (int i = 0; i < 5; i++) {
    inv *= (2 - n * inv);
  }
  return -inv;
}

uint64_t bn_divmod_u64(bignum* q, const bignum* a, uint64_t d) {
  q->limbs = calloc(a->size, sizeof(u64));
  q->size = a->size;

  __int128 rem = 0;

  for (i64 i = a->size - 1; i >= 0; i--) {
    __int128 cur = (rem << 64) | a->limbs[i];

    q->limbs[i] = (u64)(cur / d);
    rem = cur % d;
  }

  bn_trim(q);
  return (uint64_t)rem;
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
#define MAX(a, b) ((a) > (b) ? (a) : (b))

void bn_add(bignum* r, const bignum* a, const bignum* b) {
  u64 max = MAX(a->size, b->size);

  r->limbs = calloc(max + 1, sizeof(u64));

  u64 carry = 0;

  for (u64 i = 0; i < max; i++) {
    u64 av = (i < a->size) ? a->limbs[i] : 0;
    u64 bv = (i < b->size) ? b->limbs[i] : 0;

    unsigned __int128 sum = (unsigned __int128)av + bv + carry;

    r->limbs[i] = (u64)sum;    // low 64 bits
    carry = (u64)(sum >> 64);  // high 64 bits
  }

  r->limbs[max] = carry;
  r->size = max + (carry ? 1 : 0);
}
void bn_mul(bignum* r, const bignum* a, const bignum* b) {
  u64 max = a->size + b->size;
  u64* temp = calloc(max, sizeof(u64));

  for (u64 i = 0; i < a->size; i++) {
    u64 carry = 0;

    for (u64 j = 0; j < b->size; j++) {
      u64 idx = i + j;

      unsigned __int128 prod =
          (unsigned __int128)a->limbs[i] * b->limbs[j] + temp[idx] + carry;

      temp[idx] = (u64)prod;
      carry = (u64)(prod >> 64);
    }

    temp[i + b->size] += carry;
  }

  free(r->limbs);
  r->limbs = temp;
  r->size = max;
  bn_trim(r);
}

// Compare two bignums: returns 1 if a > b, -1 if a < b, 0 if a == b
int bn_cmp(const bignum* a, const bignum* b) {
  if (a->size != b->size) return (a->size < b->size) ? -1 : 1;

  for (i64 i = a->size - 1; i >= 0; i--) {
    if (a->limbs[i] != b->limbs[i]) return (a->limbs[i] < b->limbs[i]) ? -1 : 1;
  }
  return 0;
}
void bn_copy(bignum* r, const bignum* a) {}

void bn_sub_abs(bignum* r, const bignum* a, const bignum* b) {
  r->limbs = calloc(a->size, sizeof(u64));

  u64 borrow = 0;

  for (u64 i = 0; i < a->size; i++) {
    u64 av = a->limbs[i];
    u64 bv = (i < b->size) ? b->limbs[i] : 0;

    unsigned __int128 diff = (unsigned __int128)av - bv - borrow;

    r->limbs[i] = (u64)diff;
    borrow = (diff >> 127) & 1;  // detect underflow
  }

  r->size = a->size;
  bn_trim(r);
}
// Subtract: r = a - b
void bn_sub(bignum* r, const bignum* a, const bignum* b) {
  // a - (-b) = a + b
  if (!a->is_neg && b->is_neg) {
    bn_add(r, a, b);
    r->is_neg = false;
    return;
  }

  // (-a) - b = -(a + b)
  if (a->is_neg && !b->is_neg) {
    bn_add(r, a, b);
    r->is_neg = true;
    return;
  }

  i64 cmp = bn_cmp(a, b);

  if (cmp == 0) {
    bn_init(r);
    return;
  }

  if (cmp > 0) {
    // |a| > |b|
    bn_sub_abs(r, a, b);
    r->is_neg = a->is_neg;
  } else {
    // |b| > |a|
    bn_sub_abs(r, b, a);
    r->is_neg = !a->is_neg;
  }
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
void bn_lshift1(bignum* r) {
  u64 carry = 0;

  for (u64 i = 0; i < r->size; i++) {
    u64 new_carry = r->limbs[i] >> 63;
    r->limbs[i] = (r->limbs[i] << 1) | carry;
    carry = new_carry;
  }

  if (carry) {
    r->limbs = realloc(r->limbs, (r->size + 1) * sizeof(u64));
    r->limbs[r->size++] = carry;
  }
}

// Shift r left by 1 and add 0 or 1
void bn_lshift1_add(bignum* r, int bit) {
  u64 carry = bit;
  for (u64 i = 0; i < r->size; i++) {
    __uint128_t tmp = ((__uint128_t)r->limbs[i] << 1) | carry;
    r->limbs[i] = (u64)tmp;
    carry = tmp >> 64;
  }

  if (carry) {
    r->limbs = realloc(r->limbs, (r->size + 1) * sizeof(u64));
    r->limbs[r->size++] = carry;
  }
}

// Left shift: r = a << shift
void bn_lshift(bignum* r, const bignum* a, int shift) {}
// Modulo: r = a % n (using binary long division)
void bn_mod(bignum* r, const bignum* a, const bignum* b) {
  bn_init(r);
  r->limbs = calloc(1, sizeof(u64));
  r->size = 1;

  int nbits = bn_bit_length(a);

  for (int i = nbits - 1; i >= 0; i--) {
    bn_lshift1_add(r, bn_get_bit(a, i));

    if (bn_cmp(r, b) >= 0) {
      bignum tmp;
      bn_init(&tmp);
      bn_sub_abs(&tmp, r, b);  // safe subtraction
      bn_free(r);
      *r = tmp;
    }
  }

  bn_trim(r);
}
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

void bn_add512_unrolled(uint64_t* dest, const uint64_t* src) {
  __asm__ __volatile__(
      "clc\n\t"  // Clear Carry Flag (CF) to start the chain

      // Limb 0
      "movq 0(%0), %%r8\n\t"   // Load dest[0]
      "adcxq 0(%1), %%r8\n\t"  // Add src[0] + CF
      "movq %%r8, 0(%0)\n\t"   // Store result

      // Limb 1
      "movq 8(%0), %%r9\n\t"   // Load dest[1]
      "adcxq 8(%1), %%r9\n\t"  // Add src[1] + CF
      "movq %%r9, 8(%0)\n\t"   // Store result

      // Limb 2
      "movq 16(%0), %%r10\n\t"
      "adcxq 16(%1), %%r10\n\t"
      "movq %%r10, 16(%0)\n\t"

      // Limb 3
      "movq 24(%0), %%r11\n\t"
      "adcxq 24(%1), %%r11\n\t"
      "movq %%r11, 24(%0)\n\t"

      // Limb 4
      "movq 32(%0), %%r12\n\t"
      "adcxq 32(%1), %%r12\n\t"
      "movq %%r12, 32(%0)\n\t"

      // Limb 5
      "movq 40(%0), %%r13\n\t"
      "adcxq 40(%1), %%r13\n\t"
      "movq %%r13, 40(%0)\n\t"

      // Limb 6
      "movq 48(%0), %%r14\n\t"
      "adcxq 48(%1), %%r14\n\t"
      "movq %%r14, 48(%0)\n\t"

      // Limb 7
      "movq 56(%0), %%r15\n\t"
      "adcxq 56(%1), %%r15\n\t"
      "movq %%r15, 56(%0)\n\t"

      :
      : "r"(dest), "r"(src)
      : "r8", "r9", "r10", "r11", "r12", "r13", "r14", "r15", "memory", "cc");
}

int bn_jacobi(bignum* a_in, const bignum* n_in) {}

// Helper for modular subtraction: (a - b) % n
void bn_sub_mod(bignum* r, const bignum* a, const bignum* b, const bignum* n) {}

