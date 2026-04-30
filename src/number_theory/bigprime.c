#include <fcntl.h>
#include <math.h>
#include <stdarg.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "../../include/bignum.h"
#include "../../include/primes.h"

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

    uint64_t multi_prime = 3ULL * 5 * 7 * 11 * 13 * 17;
    uint64_t rem = bn_mod_u64(p, multi_prime);

    if (rem % 3 == 0 || rem % 5 == 0 || rem % 7 == 0 || rem % 11 == 0 ||
        rem % 13 == 0 || rem % 17 == 0)
      continue;

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

bool trialdiv(bignum* n, u64 g)
{
  u64 i = 0;

  while (i < 100000 && primes[i] < g) {
    if (bn_mod_u64(n, primes[i]) == 0) {
      return false;
    }
    i++;
  }

  return true;
}

// simple version for r = 1
bool checkLemma1(bignum* n, bignum* n_min1, bignum* a, bignum* q)
{
  bignum tmp, exp, X, gcd;
  bn_init_multi(&tmp, &exp, &X, &gcd, NULL);

  bool result = false;

  bn_copy(&exp, n_min1);
  bn_div(&exp, &exp, q);
  bn_copy(&tmp, a);

  bn_mod_exp(&X, &tmp, &exp, n);

  if (bn_is_eq_i64(&X, 1)) {
    goto cleanup;
  }

  bn_set_u64(&tmp, 1);
  bn_sub(&tmp, &X, &tmp);

  bn_gcd(&gcd, &tmp, n);
  if (!bn_is_eq_i64(&gcd, 1)) {
    goto cleanup;
  }

  bn_mod_exp(&X, &X, q, n);

  result = bn_is_eq_i64(&X, 1);

cleanup:
  bn_free_multi(&tmp, &exp, &X, &gcd, NULL);
  return result;
}

double gen_rel_size()
{
  // Generate uniform random variable u in [0, 1]
  double u = (double)rand() / RAND_MAX;

  return pow(2.0, u - 1.0);
}

void bn_provable_prime(bignum* p, u64 k)
{
  // Keep trying from the absolute top until it succeeds.
  // Every time it fails, it will have cleanly unwound the stack and freed all
  // memory.
  while (!bn_provable_prime_inner(p, k)) {
    printf("Restarting prime generation from scratch...\n");
  }
}

// Implementation of Maurers simpler algorithm
bool bn_provable_prime_inner(bignum* p, u64 k)
{
  //  base case k <= 20
  //  use Baillie-PSW instead of trial factoring
  if (k <= 20) {
    do {
      // generates a random odd k-bit integer
      bn_gen_random(p, k);
    } while (!bn_bpsw(p));

    return true;
  }

  // constants
  const double c_opt = 0.1;
  u64 margin = k / 6;

  bignum a, n, q, I, R, twoI, n_min1, two, two_q, tmp;
  bn_init_multi(&a, &n, &q, &I, &R, &twoI, &n_min1, &two, &two_q, &tmp, NULL);
  i64 i, g;
  bool success;
  bn_set_u64(&two, 2);

  // trial division bound
  g = (u64)(c_opt * k * k + 1);

restart:

  double rel_size;
  do {
    rel_size = gen_rel_size();
  } while ((k * rel_size >= (k - margin)));

  printf("new size %llu \n", (u64)(rel_size * k));
  // recursive call
  if (!bn_provable_prime_inner(&q, (u64)(rel_size * k))) {
    bn_free_multi(&a, &n, &q, &I, &R, &twoI, &n_min1, &two, &two_q, &tmp, NULL);
    return false;
  }

  bn_copy(&two_q, &q);
  bn_lshift1(&two_q);

  // I = 2^(k-1) / 2q
  bn_set_u64(&I, 1);
  bn_lshift(&I, &I, k - 1);
  bn_div(&I, &I, &two_q);

  // twoI = 2^k / 2q
  bn_set_u64(&twoI, 1);
  bn_lshift(&twoI, &twoI, k);
  bn_div(&twoI, &twoI, &two_q);

  success = false;
  u64 attempts = 0;
  while (!success) {
    attempts++;
    if (attempts > 1000) {
      bn_free_multi(&a, &n, &q, &I, &R, &twoI, &n_min1, &two, &two_q, &tmp,
                    NULL);
      return false;
    }

    bn_gen_random_range(&R, &I, &twoI);

    // n = 2 * rand(I, 2I) * q + 1
    bn_mul(&n, &R, &q);
    bn_lshift1(&n);
    bn_copy(&n_min1, &n);
    bn_add_u64(&n, &n, 1);

    if (trialdiv(&n, g)) {
      for (int j = 0; j < 200; j++) {
        bn_gen_random_range(&a, &two, &n_min1);

        if (checkLemma1(&n, &n_min1, &a, &q)) {
          success = true;
          break;
        }
      }

      if (success) {
        printf("   [FOUND] %llu-bit prime\n", k);
        break;
      }
    }
  }

  bn_copy(p, &n);

  bn_free_multi(&a, &n, &q, &I, &R, &twoI, &n_min1, &two, &two_q, &tmp, NULL);
  return true;
}
