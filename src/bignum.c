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

bool bn_is_zero(const bignum* a) {
  return (a->size == 0 || (a->size == 1 && a->limbs[0] == 0));
}

bool bn_alloc(bignum* r, u64 capacity) {
  if (capacity <= r->capacity) return true;

  // Grow memory exponentially
  u64 new_cap = r->capacity * 2;
  if (new_cap < capacity) new_cap = capacity;

  u64* new_limbs = realloc(r->limbs, new_cap * sizeof(u64));
  if (!new_limbs) return false;

  memset(new_limbs + r->capacity, 0, (new_cap - r->capacity) * sizeof(u64));

  r->limbs = new_limbs;
  r->capacity = new_cap;
  return true;
}

void bn_init_val(bignum* n, const char* str) {
  bn_init(n);
  if (!str) return;

  n->is_neg = (str[0] == '-');
  const char* s = n->is_neg ? str + 1 : str;

  // Start with a value of 0
  bn_set_u64(n, 0);

  for (size_t i = 0; s[i]; i++) {
    if (!isdigit(s[i])) continue;
    int digit = s[i] - '0';

    // Correct way to do n = n * 10 + digit:
    // Perform a carry-propagation multiplication across ALL limbs
    u64 carry = digit;
    for (u64 j = 0; j < n->size; j++) {
      unsigned __int128 prod = (unsigned __int128)n->limbs[j] * 10 + carry;
      n->limbs[j] = (u64)prod;
      carry = (u64)(prod >> 64);
    }

    // If there's still a carry after the last limb, grow the bignum
    while (carry) {
      bn_alloc(n, n->size + 1);
      n->limbs[n->size++] = carry % 0xFFFFFFFFFFFFFFFFULL;  // Simplified
      // Actually, with base 10, the carry will never exceed a single u64
      n->limbs[n->size - 1] = carry;
      carry = 0;
    }
  }
  bn_trim(n);
}

void bn_init_val2(bignum* n, const char* str) {
  bn_init(n);

  n->is_neg = (str[0] == '-');
  const char* s = n->is_neg ? str + 1 : str;

  // start with zero
  /*  n->limbs = calloc(1, sizeof(u64));
    n->size = 1;
    n->capacity = 1;
  */
  bn_alloc(n, 1);
  n->size = 1;

  for (size_t i = 0; s[i]; i++) {
    if (!isdigit(s[i])) continue;
    int digit = s[i] - '0';

    // n = n * 10
    u64 carry = 0;
    for (u64 j = 0; j < n->size; j++) {
      unsigned __int128 prod = (unsigned __int128)n->limbs[j] * 10 + carry;

      n->limbs[j] = (u64)prod;
      carry = (u64)(prod >> 64);
    }

    if (carry) {
      bn_alloc(n, n->size + 1);
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

void bn_print(bignum* n) {
  if (n->size == 0 || (n->size == 1 && n->limbs[0] == 0)) {
    printf("0\n");
    return;
  }

  if (n->is_neg) printf("-");

  bignum tmp;
  bn_init(&tmp);
  bn_copy(&tmp, n);

  const uint64_t BASE10 = 10000000000000000000ULL;  // 10^19
  uint64_t* parts = NULL;
  size_t parts_count = 0;

  // Extract 19-digit chunks
  while (!(tmp.size == 1 && tmp.limbs[0] == 0)) {
    bignum q;
    bn_init(&q);
    // Ensure bn_divmod_u64 is correctly updating 'q' and returning 'rem'
    uint64_t rem = bn_divmod_u64(&q, &tmp, BASE10);

    parts = realloc(parts, (parts_count + 1) * sizeof(uint64_t));
    parts[parts_count++] = rem;

    bn_free(&tmp);
    tmp = q;
  }

  // 1. Print the most significant chunk (no leading zeros)
  printf("%llu", parts[parts_count - 1]);

  // 2. Print all other chunks (MUST have 19 digits, pad with zeros)
  for (int64_t i = (int64_t)parts_count - 2; i >= 0; i--) {
    printf("%019llu", parts[i]);
  }

  printf("\n");

  free(parts);
  bn_free(&tmp);
}

// print in base 10
void bn_print2(bignum* n) {
  if (n->size == 0 || (n->size == 1 && n->limbs[0] == 0)) {
    printf("0\n");
    return;
  }

  if (n->is_neg) printf("-");

  bignum tmp;
  bn_init(&tmp);
  bn_copy(&tmp, n);

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

  if (a->size == 0) {
    r->size = 0;
    return;
  }

  r->size = a->size;
  r->is_neg = a->is_neg;

  bn_alloc(r, a->size);
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
    bn_alloc(a, limb + 1);
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
