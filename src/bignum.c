#include "../include/bignum.h"

#include <ctype.h>
#include <fcntl.h>
#include <stdarg.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

#define BASE_10_19 10000000000000000000ULL

// Clear oject
void bn_init(bignum* r)
{
  r->limbs = NULL;
  r->size = 0;
  r->capacity = 0;
  r->is_neg = false;
}

// initializes multiple bignums
// IMPORTANT: terminate with NULL
void bn_init_multi(bignum* r, ...)
{
  if (r == NULL) return;

  bn_init(r);

  va_list arg;
  va_start(arg, r);

  bignum* next;
  while ((next = va_arg(arg, bignum*)) != NULL) {
    bn_init(next);
  }

  va_end(arg);
}

// returns 1 if a bignum is even, 0 else
int bn_is_even(const bignum* n)
{
  if (n->size == 0 || n->limbs == NULL) {
    return 1;
  }

  return (n->limbs[0] & 1) == 0;
}

// returns true if a is even, false else
bool bn_is_zero(const bignum* a)
{
  return (a->size == 0 || (a->size == 1 && a->limbs[0] == 0));
}

// returns true if n is equal to the 64 bit signed integer a
bool bn_is_eq_i64(const bignum* n, i64 a)
{
  bignum test;
  bn_init(&test);
  bn_set_i64(&test, a);
  if (bn_cmp(n, &test) == 0) {
    return true;
  } else {
    return false;
  }
}

// allocates more memory for r, at least the specified capacity
bool bn_alloc(bignum* r, u64 capacity)
{
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

// reads a string input in base 10 and sets n to that value
void bn_init_val(bignum* n, const char* str)
{
  bn_init(n);
  if (!str) return;

  n->is_neg = (str[0] == '-');
  const char* s = n->is_neg ? str + 1 : str;

  // Start with a value of 0
  bn_set_u64(n, 0);

  for (size_t i = 0; s[i]; i++) {
    if (!isdigit(s[i])) continue;

    int digit = s[i] - '0';

    // go digits by digit propagating carry
    u64 carry = digit;
    for (u64 j = 0; j < n->size; j++) {
      unsigned __int128 prod = (unsigned __int128)n->limbs[j] * 10 + carry;
      n->limbs[j] = (u64)prod;
      carry = (u64)(prod >> 64);
    }

    // If there's still a carry after the last limb, grow the bignum
    while (carry) {
      bn_alloc(n, n->size + 1);
      n->limbs[n->size++] = carry % 0xFFFFFFFFFFFFFFFFULL;
      n->limbs[n->size - 1] = carry;
      carry = 0;
    }
  }
  bn_trim(n);
}

// returns the n as a string
char* bn_to_string(bignum* n)
{
  if (n->size == 0 || (n->size == 1 && n->limbs[0] == 0)) {
    return strdup("0");
  }

  bignum tmp;
  bn_init(&tmp);
  bn_copy(&tmp, n);

  uint64_t* parts = NULL;
  size_t parts_count = 0;

  // Extract 19-digit chunks
  while (!(tmp.size == 1 && tmp.limbs[0] == 0)) {
    bignum q;
    bn_init(&q);
    // Ensure bn_divmod_u64 is correctly updating 'q' and returning 'rem'
    uint64_t rem = bn_divmod_u64(&q, &tmp, BASE_10_19);

    parts = realloc(parts, (parts_count + 1) * sizeof(uint64_t));
    parts[parts_count++] = rem;

    bn_free(&tmp);
    tmp = q;
  }

  u64 buffer_size = (parts_count * 19) + 2;
  char* result = (char*)malloc(buffer_size);
  if (!result) return NULL;

  char* ptr = result;

  if (n->is_neg) {
    *ptr++ = '-';
  }

  // 1. Print the most significant chunk (no leading zeros)
  ptr += sprintf(ptr, "%llu", (unsigned long long)parts[parts_count - 1]);

  // 2. Print all other chunks (MUST have 19 digits, pad with zeros)
  for (int64_t i = (int64_t)parts_count - 2; i >= 0; i--) {
    ptr += sprintf(ptr, "%019llu", (unsigned long long)parts[i]);
  }

  free(parts);
  bn_free(&tmp);

  return result;
}

// prints a bignum to the console in base 10
void bn_print(bignum* n)
{
  if (n->size == 0 || (n->size == 1 && n->limbs[0] == 0)) {
    printf("0");
    return;
  }

  if (n->is_neg) printf("-");

  bignum tmp;
  bn_init(&tmp);
  bn_copy(&tmp, n);

  uint64_t* parts = NULL;
  size_t parts_count = 0;

  // Extract 19-digit chunks
  while (!(tmp.size == 1 && tmp.limbs[0] == 0)) {
    bignum q;
    bn_init(&q);
    // Ensure bn_divmod_u64 is correctly updating 'q' and returning 'rem'
    uint64_t rem = bn_divmod_u64(&q, &tmp, BASE_10_19);

    parts = realloc(parts, (parts_count + 1) * sizeof(uint64_t));
    parts[parts_count++] = rem;

    bn_free(&tmp);
    tmp = q;
  }

  // 1. Print the most significant chunk (no leading zeros)
  printf("%llu", (unsigned long long)parts[parts_count - 1]);

  // 2. Print all other chunks (MUST have 19 digits, pad with zeros)
  for (int64_t i = (int64_t)parts_count - 2; i >= 0; i--) {
    printf("%019llu", (unsigned long long)parts[i]);
  }

  free(parts);
  bn_free(&tmp);
}

// same as bn_print with a newline
void bn_println(bignum* n)
{
  if (n->size == 0 || (n->size == 1 && n->limbs[0] == 0)) {
    printf("0\n");
    return;
  }
  bn_print(n);

  printf("\n");
}

// frees a bignum
void bn_free(bignum* r)
{
  if (r->limbs) {
    free(r->limbs);
    r->limbs = NULL;
  }
  r->size = 0;
  r->capacity = 0;
}

// frees multiple bignums
// IMPORTANT: terminate with NULL
void bn_free_multi(bignum* r, ...)
{
  if (r == NULL) return;

  va_list arg;
  va_start(arg, r);

  bignum* next;
  while ((next = va_arg(arg, bignum*)) != NULL) {
    bn_free(next);
  }

  va_end(arg);

  bn_free(r);
}

// removes empty limbs from a bignum
void bn_trim(bignum* r)
{
  while (r->size > 1 && r->limbs[r->size - 1] == 0) {
    r->size--;
  }
}

// Returns the value of the i-th bit (0 or 1)
int bn_get_bit(const bignum* a, int i)
{
  u64 limb = i / 64;
  u64 offset = i % 64;

  if (limb >= a->size) return 0;

  return (a->limbs[limb] >> offset) & 1;
}

// Compare two bignums: returns 1 if a > b, -1 if a < b, 0 if a == b
int bn_cmp(const bignum* a, const bignum* b)
{
  if (!a->is_neg && b->is_neg) return 1;

  if (a->is_neg && !b->is_neg) return -1;

  int cmp = 0;

  if (a->size != b->size)
    cmp = (a->size < b->size) ? -1 : 1;
  else {
    for (i64 i = a->size - 1; i >= 0; i--) {
      if (a->limbs[i] != b->limbs[i]) {
        cmp = (a->limbs[i] < b->limbs[i]) ? -1 : 1;
        break;
      }
    }
  }
  if (a->is_neg && b->is_neg) {
    cmp = -cmp;
  }

  return cmp;
}

// deep copies one a into r
void bn_copy(bignum* r, const bignum* a)
{
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
void bn_set_u64(bignum* n, uint64_t val)
{
  bn_free(n);  // free any existing limbs
  bn_alloc(n, 1);
  n->limbs[0] = val;
  n->size = 1;
  n->capacity = 1;
  n->is_neg = false;
}

// Set to a 64 bit signed integer
void bn_set_i64(bignum* n, int64_t val)
{
  if (val >= 0) {
    bn_set_u64(n, (u64)val);
  } else {
    u64 abs = (u64) - (val - 1) + 1;
    bn_set_u64(n, abs);
    n->is_neg = true;
  }
}

// sets the i-th bit of a
void bn_set_bit(bignum* a, int i)
{
  u64 limb = i / 64;
  u64 offset = i % 64;

  if (limb >= a->size) {
    bn_alloc(a, limb + 1);
    a->size = limb + 1;
  }

  a->limbs[limb] |= ((u64)1 << offset);
}

// Bit length of a bignum
int bn_bit_length(const bignum* a)
{
  if (a->size == 0) return 0;

  u64 top = a->limbs[a->size - 1];
  int bits = 0;

  while (top) {
    top >>= 1;
    bits++;
  }

  return (a->size - 1) * 64 + bits;
}

static inline u64 count_trailing_zeros_u64(u64 val)
{
  if (val == 0) return 64;
  return (u64)__builtin_ctzll(val);
}

static inline u64 count_leading_zeros_u64(u64 val)
{
  if (val == 0) return 64;
  return (u64)__builtin_clzll(val);
}

// returns the number of trailing zeros of a
u64 bn_cnt_trailing_zeros(const bignum* a)
{
  if (bn_is_zero(a)) return 0;

  u64 zeros = 0;
  u64 i = 0;

  while (i < a->size && a->limbs[i] == 0) {
    zeros += 64;
    i++;
  }

  if (i == a->size) {
    return 0;
  }

  zeros += count_trailing_zeros_u64(a->limbs[i]);

  return zeros;
}

// returns the number of leading zeros of a
u64 bn_cnt_leading_zeros(const bignum* a)
{
  if (bn_is_zero(a)) return 0;

  u64 zeros = 0;
  u64 i = 0;

  while (i < a->size && a->limbs[i] == 0) {
    zeros += 64;
    i++;
  }

  if (i == a->size) {
    return 0;
  }

  zeros += count_leading_zeros_u64(a->limbs[i]);

  return zeros;
}

// generates a random odd number using a hardware random number generator (on
// Linux)
bool bn_gen_random(bignum* r, int bits)
{
  int limbs_needed = (bits + 63) / 64;

  if (!bn_alloc(r, limbs_needed)) return false;
  r->size = limbs_needed;

  int fd = open("/dev/urandom", O_RDONLY);
  if (fd < 0) return false;

  // Read random bytes directly into the limb memory
  if (read(fd, r->limbs, limbs_needed * sizeof(u64)) !=
      (ssize_t)(limbs_needed * sizeof(u64))) {
    close(fd);
    return false;
  }
  close(fd);

  // Mask the top limb to fit the exact bit length
  int top_bits = bits % 64;
  if (top_bits != 0) {
    u64 mask = ((u64)1 << top_bits) - 1;
    r->limbs[r->size - 1] &= mask;
  }

  // Ensure it's exactly 'bits' long by setting the MSB
  r->limbs[r->size - 1] |= ((u64)1 << ((bits - 1) % 64));

  // Ensure it's odd by setting the LSB
  r->limbs[0] |= 1;

  return true;
}

bool bn_gen_random_with_fd(bignum* r, int bits, int fd)
{
  int limbs_needed = (bits + 63) / 64;

  if (!bn_alloc(r, limbs_needed)) return false;
  r->size = limbs_needed;

  // Read random bytes directly using the open file descriptor
  if (read(fd, r->limbs, limbs_needed * sizeof(u64)) !=
      (ssize_t)(limbs_needed * sizeof(u64))) {
    return false;
  }

  // Mask the top limb to fit the exact bit length
  int top_bits = bits % 64;
  if (top_bits != 0) {
    u64 mask = ((u64)1 << top_bits) - 1;
    r->limbs[r->size - 1] &= mask;
  }

  // Ensure it's exactly 'bits' long by setting the MSB
  r->limbs[r->size - 1] |= ((u64)1 << ((bits - 1) % 64));

  // Ensure it's odd by setting the LSB
  r->limbs[0] |= 1;

  return true;
}

// generates a random prime p with bits length using a variety of primality
// tests
bool bn_gen_prime(bignum* p, int bits)
{
  bignum a;
  bn_init(&a);

  int fd = open("/dev/urandom", O_RDONLY);
  if (fd < 0) {
    bn_free(&a);
    return false;
  }

  // Small primes to check for quick trial division
  u64 small_primes[] = {
      2,   3,   5,   7,   11,  13,  17,  19,  23,  29,  31,  37,  41,  43,
      47,  53,  59,  61,  67,  71,  73,  79,  83,  89,  97,  101, 103, 107,
      109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181,
      191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263,
      269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349,
      353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431, 433,
      439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503, 509, 521,
      523, 541, 547, 557, 563, 569, 571, 577, 587, 593, 599, 601, 607, 613,
      617, 619, 631, 641, 643, 647, 653, 659, 661, 673, 677, 683, 691, 701,
      709, 719, 727, 733, 739, 743, 751, 757, 761, 769, 773, 787, 797, 809,
      811, 821, 823, 827, 829, 839, 853, 857, 859, 863, 877, 881, 883, 887,
      907, 911, 919, 929, 937, 941, 947, 953, 967, 971, 977, 983, 991, 997};

  int num_primes = sizeof(small_primes) / sizeof(small_primes[0]);
  while (true) {
    bool composite = false;

    if (!bn_gen_random_with_fd(p, bits, fd)) {
      bn_free(&a);
      return false;
    }

    uint64_t multi_prime = 3ULL * 5 * 7 * 11 * 13 * 17;
    uint64_t rem = bn_mod_u64(p, multi_prime);

    if (rem % 3 == 0 || rem % 5 == 0 || rem % 7 == 0 || rem % 11 == 0 ||
        rem % 13 == 0 || rem % 17 == 0)
      continue;

    for (int i = 0; i < num_primes; i++) {
      if (bn_mod_u64(p, small_primes[i]) == 0) {
        composite = true;
        break;
      }
    }

    // 2. BPSW deterministic :)
    if (!composite) {
      if (bn_bpsw(p)) {
        close(fd);
        return true;
      }
    }
  }
}
