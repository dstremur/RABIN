#include <ctype.h>
#include <fcntl.h>
#include <stdarg.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

#include "../include/bignum.h"

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
  int fd = open("/dev/urandom", O_RDONLY);
  if (fd < 0) {
    close(fd);
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
      return false;
    }

    bn_set_bit(p, bits - 1);
    bn_set_bit(p, 0);

    // uint64_t multi_prime = 3ULL * 5 * 7 * 11 * 13 * 17;
    // uint64_t rem = bn_mod_u64(p, multi_prime);

    // if (rem % 3 == 0 || rem % 5 == 0 || rem % 7 == 0 || rem % 11 == 0 ||
    //     rem % 13 == 0 || rem % 17 == 0)
    //   continue;

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

bool bn_is_perfect_square(bignum* n)
{
  if (bn_is_zero(n)) return true;

  if (n->is_neg) return false;

  // Newton-Raphson method
  bignum x, y, tmp, rem;
  bn_init_multi(&x, &y, &tmp, &rem, NULL);

  // Initial guess: x = 2^(bits/2)
  int bits = bn_bit_length(n);
  bn_set_u64(&x, 1);
  bn_lshift(&x, &x, (bits + 1) / 2);

  // Newton-Raphson iteration
  while (true) {
    // y = (x + n/x) / 2
    bn_newton_div(&tmp, n, &x);
    bn_add(&y, &x, &tmp);
    bn_rshift1(&y);

    if (bn_cmp(&y, &x) >= 0) {
      break;
    }

    bn_copy(&x, &y);
  }

  // Check if x * x == n
  bn_mul(&tmp, &x, &x);
  bool is_square = (bn_cmp(&tmp, n) == 0);

  bn_free(&x);
  bn_free(&y);
  bn_free(&tmp);
  bn_free(&rem);

  return is_square;
}

bool bn_stronglucas(bignum* n, bignum* P, bignum* Q)
{
  bignum p, q, d, u, v, qn, tmp, n_plus_1;
  bn_init_multi(&p, &q, &d, &u, &v, &qn, &tmp, &n_plus_1, NULL);

  bool prime = false;

  int s = 0;

  bn_copy(&n_plus_1, n);
  bn_add_u64(&n_plus_1, &n_plus_1, 1);
  bn_copy(&d, &n_plus_1);

  while (bn_is_even(&d)) {
    bn_rshift1(&d);
    s++;
  }

  bn_lucas_solve_mod(&u, &v, P, Q, &qn, &d, n);
  // 1. condition
  if (bn_is_zero(&u)) {
    prime = true;
    goto cleanup;
  }
  // 2. condition
  if (bn_is_zero(&v)) {
    prime = true;
    goto cleanup;
  }
  // check 2. condition for all d * 2^r
  for (int r = 1; r < s; r++) {
    bn_mul(&tmp, &v, &v);
    bn_mod(&tmp, &tmp, n);

    bn_add(&tmp, &tmp, n);
    bn_sub(&tmp, &tmp, &qn);
    bn_add(&tmp, &tmp, n);
    bn_sub(&tmp, &tmp, &qn);
    bn_mod(&v, &tmp, n);

    bn_mul(&qn, &qn, &qn);
    bn_mod(&qn, &qn, n);

    if (bn_is_eq_i64(&v, 0)) {
      prime = true;
      goto cleanup;
    }
  }

cleanup:
  bn_free(&p);
  bn_free(&q);
  bn_free(&d);
  bn_free(&u);
  bn_free(&v);
  bn_free(&qn);
  bn_free(&tmp);
  bn_free(&n_plus_1);
  return prime;
}

// check n for primality using a Baillie-PSW test
// Assumes sufficient trial division was done previously
bool bn_bpsw(const bignum* n)
{
  // 1. run miller rabin base 2
  bignum two;
  bn_init(&two);
  bn_set_u64(&two, 2);
  if (!bn_rabin_mont(n, &two)) {
    bn_free(&two);
    return false;
  }

  bignum D, magnitude;
  bn_init(&D);
  bn_init(&magnitude);
  bn_set_i64(&magnitude, 5);

  bool negative = false;

  int rounds = 0;
  // finds a D using Selfridges method A*
  while (1) {
    // check if n is perfect square after 5 rounds
    if (rounds == 5 && bn_is_perfect_square(n)) {
      bn_free(&two);
      bn_free(&D);
      bn_free(&magnitude);
      return false;
    }

    bn_copy(&D, &magnitude);
    D.is_neg = negative;

    i64 jacobi = bn_jacobi(&D, n);
    if (jacobi == -1) {
      break;
    }

    if (jacobi == 0) {
      bn_free(&D);
      bn_free(&magnitude);
      bn_free(&two);
      return false;
    }

    bn_add_u64(&magnitude, &magnitude, 2);
    negative = !negative;

    rounds++;
  }

  bignum P, Q;
  bn_init_multi(&P, &Q, NULL);
  bn_set_u64(&P, 1);
  bn_set_u64(&Q, 1);

  if (negative) {
    bn_copy(&Q, &magnitude);
    bn_add_u64(&Q, &Q, 1);
    bn_divmod_u64(&Q, &Q, 4);
  } else {
    bn_copy(&Q, &magnitude);
    bn_sub(&Q, &Q, &P);
    bn_divmod_u64(&Q, &Q, 4);
    bn_sub(&Q, n, &Q);
  }

  if (bn_is_eq_i64(&D, 5)) {
    bn_set_u64(&P, 5);
    bn_set_u64(&Q, 5);
  }

  bool res = bn_stronglucas(n, &P, &Q);

  bn_free(&P);
  bn_free(&Q);
  bn_free(&D);
  bn_free(&magnitude);
  bn_free(&two);
  return res;
}
