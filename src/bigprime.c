#include <fcntl.h>
#include <math.h>
#include <stdarg.h>
#include <stddef.h>
#include <stdlib.h>
#include <unistd.h>

#include "../include/bignum.h"
#include "../include/primes.h"

// generates a random prime p with bits length using a variety of primality
// tests
bool bn_gen_prime(bignum* p, int bits)
{
  int fd = open("/dev/urandom", O_RDONLY);
  if (fd < 0) {
    close(fd);
    return false;
  }

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

    for (int i = 0; i < 1500; i++) {
      if (bn_mod_u64(p, primes[i]) == 0) {
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

bool bn_is_perfect_square(const bignum* n)
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

bool bn_stronglucas(const bignum* n, bignum* P, bignum* Q)
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

// Implementation of Maurers algorithm
void bn_provable_prime(bignum* p, u64 k)
{
  // ---- Small base case ----
  if (k <= 20) {
    do {
      bn_gen_random(p, k);
      bn_set_bit(p, k - 1);
    } while (!bn_bpsw(p));
    return;
  }

  // ---- Step 1: choose r ----
  double r;
  const u64 m = 20;

  if (k > 2 * m) {
    do {
      double s = (double)rand() / RAND_MAX;
      r = pow(2.0, s - 1.0);
    } while ((u64)(r * k) <= m);
  } else {
    r = 0.5;
  }

  // trial division bound
  u64 B = k * k;

  // ---- Step 2: recursively generate q ----
  bignum q;
  bn_init(&q);
  bn_provable_prime(&q, (u64)(r * k));

  // ---- Precompute ----
  bignum I, R, n;
  bignum a, b, d;
  bignum tmp, two, one;
  bignum n_minus_one, n_minus_two;

  bn_init_multi(&I, &R, &n, &a, &b, &d, &tmp, &two, &one, &n_minus_one,
                &n_minus_two, NULL);

  bn_set_u64(&two, 2);
  bn_set_u64(&one, 1);

  // I = 2^(k-1) / (2q)
  bignum pow2, two_q;
  bn_init_multi(&pow2, &two_q, NULL);

  bn_set_u64(&pow2, 1);
  bn_lshift(&pow2, &pow2, k - 1);

  bn_div(&I, &pow2, &q);

  bn_free(&pow2);
  bn_free(&two_q);

  // Precompute bounds: [I+1, 2I]
  bignum low, high;
  bn_init_multi(&low, &high, NULL);

  bn_add_u64(&low, &I, 1);
  bn_copy(&high, &I);
  bn_lshift1(&high);

  // Small primes for trial division

  int attempts = 0;

  while (1) {
    attempts++;
    if (attempts > 100000) {
      // fallback safety (should basically never happen)
      bn_gen_prime(p, k);
      break;
    }

    // ---- Step 3: pick R ----
    bn_gen_random_range(&R, &low, &high);

    // ---- Step 4: n = 2Rq + 1 ----
    bn_mul(&n, &R, &q);
    bn_mul(&n, &n, &two);
    bn_add_u64(&n, &n, 1);

    // Ensure n is k bits
    if (bn_bit_length(&n) != (int)k) continue;

    // ---- Step 5: trial division ----
    int composite = 0;
    u64 i = 0;
    for (; primes[i] < B; i++) {
      if (bn_mod_u64(&n, primes[i]) == 0) {
        composite = 1;
        break;
      }
    }
    if (composite) continue;

    // ---- Step 6: enforce q > n^(1/3) ----
    if (3 * bn_bit_length(&q) <= bn_bit_length(&n)) {
      continue;
    }

    // ---- Step 7: prepare values ----
    bn_copy(&n_minus_one, &n);
    bn_sub(&n_minus_one, &n_minus_one, &one);

    bn_copy(&n_minus_two, &n);
    bn_sub(&n_minus_two, &n_minus_two, &two);

    // ---- Step 8: pick random a ----
    bn_gen_random_range(&a, &two, &n_minus_two);

    // ---- Step 9: check a^(n-1) mod n == 1 ----
    bn_mod_exp(&b, &a, &n_minus_one, &n);
    if (!bn_is_eq_i64(&b, 1)) continue;

    // ---- Step 10: Pocklington test ----
    // exponent = (n-1)/q = 2R
    bn_mul(&tmp, &R, &two);
    bn_mod_exp(&b, &a, &tmp, &n);

    bn_sub(&b, &b, &one);
    bn_gcd(&d, &b, &n);

    if (bn_is_eq_i64(&d, 1)) {
      bn_copy(p, &n);
      break;
    }
  }

  // ---- Cleanup ----
  bn_free(&q);
  bn_free_multi(&I, &R, &n, &a, &b, &d, &tmp, &two, &one, &n_minus_one,
                &n_minus_two, &low, &high, NULL);
}
