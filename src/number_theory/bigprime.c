/*
 * bigprime.c
 *
 * Primality testing and prime generation.
 *
 * This file implements the Baillie-PSW primality test (Miller-Rabin
 * base 2 plus a strong Lucas test with Selfridge's parameter search),
 * random prime generation, provable prime generation via Maurer's
 * algorithm, and generators of Proth primes / RNS prime tables.
 *
 * Copyright (C) 2026 Diego Strebel
 * SPDX-License-Identifier: GPL-3.0-or-later
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#include <fcntl.h>
#include <math.h>
#include <stdarg.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "../../include/bignum.h"
#include "../../include/primes.h"

/*
 * Generate a random prime of the given bit length into p.
 *
 * Let k = bits.
 *
 * Repeatedly draws random k-bit candidates from /dev/urandom with the
 * top and bottom bits set, rejects candidates divisible by 3, 5, 7,
 * 11, 13, or 17, trial-divides against the first 1500 primes, and
 * accepts the first candidate that passes BPSW.
 *
 * Returns true on success, false if /dev/urandom cannot be opened or
 * reading from it fails.
 *
 * Complexity:
 *   Time: O(k log k) expected - about k / ln(k) candidates, each
 *         costing trial division plus one BPSW test
 *   Auxiliary memory: O(k/64) limbs
 *   Output memory: O(k/64) limbs
 */
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
/**
 * @brief Generates a safe prime number p of a specified bit length.
 *
 * A prime p is safe if p = 2q + 1, where q is also a prime (a Sophie Germain
 * prime). This function repeatedly generates a random prime q of length (bits -
 * 1) until the calculated p passes the BPSW primality test.
 *
 * @param[out] p    Pointer to the bignum structure where the generated safe
 * prime will be stored.
 * @param[in]  bits The desired total bit length of the safe prime p.
 *
 * @return True  If the safe prime was successfully generated.
 * @return False If prime generation failed (e.g., maximum iterations reached or
 * invalid bit size).
 */
bool bn_gen_safe_prime(bignum* p, int bits)
{
  bignum q;
  bn_init(&q);
  do {
    bn_gen_prime(&q, bits - 1);

    // p = 2q + 1
    bn_lshift(p, &q, 1);
    bn_add_u64(p, p, 1);

  } while (!bn_bpsw(p));

  return true;
}

/*
 * Test whether n is a perfect square.
 *
 * Let n_l = n->size, measured in 64-bit limbs.
 *
 * Computes the integer square root with Newton-Raphson iteration
 * (x <- (x + n/x) / 2, seeded at 2^(bits/2)) and checks whether
 * x^2 == n. Zero is considered a square; negative numbers are not.
 *
 * Complexity:
 *   Time: O(n_l^2 log n_l) - O(log n_l) iterations, each dominated by
 *         a division
 *   Auxiliary memory: O(n_l) limbs for temporaries
 *   Output memory: O(1)
 */
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

/*
 * Strong Lucas test of n with Lucas parameters (P, Q).
 *
 * Let n_l = n->size, measured in 64-bit limbs.
 *
 * Writes n + 1 = d * 2^s with d odd, computes (U_d, V_d, Q^d) mod n,
 * and n passes if:
 *
 *   - U_d == 0 mod n, or
 *   - V_d == 0 mod n, or
 *   - V_{d * 2^r} == 0 mod n for some 1 <= r < s
 *
 * (V is iterated with the doubling identity V_2k = V_k^2 - 2 Q^k.)
 *
 * Returns true if n passes the test (probably prime), false otherwise.
 *
 * Complexity:
 *   Time: O(n_l^3) for the Lucas sequence computation (Montgomery
 *         context init), then O(s * n_l^2) for the V iteration
 *   Auxiliary memory: O(n_l) limbs for temporaries
 *   Output memory: O(1)
 */
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

/*
 * Baillie-PSW primality test.
 *
 * Let n_l = n->size, measured in 64-bit limbs.
 *
 * n passes if all of the following hold:
 *
 *   1. n is not divisible by any prime below 10000 (n itself is
 *      accepted if it is one of them),
 *   2. n passes Miller-Rabin to base 2,
 *   3. n passes the strong Lucas test with parameters (P, Q) =
 *      (1, (1 - D)/4) or (1, (1 + D)/4), where D is the first value
 *      in the sequence 5, -7, 9, -11, ... (Selfridge's method A*)
 *      with Jacobi symbol (D / n) = -1.
 *
 * If no such D is found within 15 rounds, n is rejected (this also
 * catches perfect squares). No known Baillie-PSW pseudoprime exists.
 *
 * Returns true if n is probably prime, false if n is even or a
 * witness for compositeness was found.
 *
 * Complexity:
 *   Time: O(n_l^3) - dominated by the Miller-Rabin and strong Lucas
 *         tests (each paying a Montgomery context initialization)
 *   Auxiliary memory: O(n_l) limbs for temporaries
 *   Output memory: O(1)
 */
bool bn_bpsw(const bignum* n)
{
  if (bn_is_even(n)) return false;

  // Trial division by small primes to quickly reject composites (mirrors
  // GMP's fast path). ~88% of random odd numbers have a factor < 10000.
  for (u64 i = 0; primes[i] < 10000; i++) {
    u64 p = primes[i];
    if (n->size == 1 && n->limbs[0] == p) return true;
    if (bn_mod_u64(n, p) == 0) return false;
  }

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
    if (rounds == 15 && bn_is_perfect_square(n)) {
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
  }

  bool res = bn_stronglucas(n, &P, &Q);

  bn_free(&P);
  bn_free(&Q);
  bn_free(&D);
  bn_free(&magnitude);
  bn_free(&two);
  return res;
}

/*
 * Lemma 1 check of Maurer's algorithm for the case r = 1.
 *
 * Let n_l = n->size, measured in 64-bit limbs.
 *
 * Given n = 2 R q + 1 with q prime, this checks the primality
 * criterion for a random base a:
 *
 *   X = a^((n-1)/q) mod n
 *
 * n is certified if X != 1, gcd(X - 1, n) == 1, and X^q == 1 mod n.
 *
 * Returns true if the check passes (n is prime), false otherwise.
 *
 * Complexity:
 *   Time: O(n_l^2) - two modular exponentiations
 *   Auxiliary memory: O(n_l) limbs for temporaries
 *   Output memory: O(1)
 */
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

/*
 * Draw a relative size for Maurer's algorithm: 2^u with u uniform in
 * [0, 1), i.e. a value in [1/2, 1) biased toward 1.
 *
 * Complexity:
 *   Time: O(1)
 *   Auxiliary memory: O(1)
 *   Output memory: O(1)
 */
double gen_rel_size()
{
  // Generate uniform random variable u in [0, 1]
  double u = (double)rand() / RAND_MAX;

  return pow(2.0, u - 1.0);
}

/*
 * Generate a provable prime of k bits into p (Maurer's algorithm).
 *
 * Repeatedly invokes bn_provable_prime_inner() until it succeeds; each
 * failed attempt unwinds its recursion and frees all memory, so the
 * loop is safe to retry indefinitely.
 *
 * Complexity:
 *   Time: O(k^3) expected per successful attempt (recursive provable
 *         prime generation plus trial division and Lemma 1 checks)
 *   Auxiliary memory: O(k^2 / 64) limbs on the recursion stack
 *   Output memory: O(k/64) limbs
 */
void bn_provable_prime(bignum* p, u64 k)
{
  // Keep trying from the absolute top until it succeeds.
  // Every time it fails, it will have cleanly unwound the stack and freed all
  // memory.
  while (!bn_provable_prime_inner(p, k)) {
    printf("Restarting prime generation from scratch...\n");
  }
}

/*
 * One attempt at Maurer's simpler algorithm for a provable k-bit prime.
 *
 * Let k = bit length of the target prime.
 *
 * Strategy:
 *
 *   - base case k <= 20: accept a random k-bit candidate that passes
 *     BPSW (treated as proven at this size),
 *   - otherwise: pick a relative size rel_size in [1/2, 1) with
 *     rel_size * k < k - k/6, recursively generate a provable prime q
 *     of that size, then search for n = 2 R q + 1 with R in
 *     [2^(k-1)/(2q), 2^k/(2q)) that survives trial division up to
 *     0.1 k^2 + 1 and passes the Lemma 1 check for some random base a
 *     (up to 200 bases).
 *
 * Returns true and stores the prime in p on success. Returns false
 * (having freed all temporaries) if the recursion fails or the search
 * exceeds 1000 candidates.
 *
 * Complexity:
 *   Time: O(k^3) expected - dominated by the recursive call and the
 *         Lemma 1 modular exponentiations
 *   Auxiliary memory: O(k^2 / 64) limbs on the recursion stack
 *   Output memory: O(k/64) limbs
 */
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

/*
 * Generate `count` Proth primes of the form p = c * 2^k + 1 and print
 * them as a C array initializer.
 *
 * Let k = exponent of the power of two.
 *
 * Scans odd c starting from the given c (rounded up to odd), testing
 * each candidate c * 2^k + 1 with BPSW, until `count` primes are
 * found or the 64-bit search space is exhausted. The output is a
 * `static const uint64_t RNS_PRIMES[]` table suitable for pasting
 * into primes.h.
 *
 * Complexity:
 *   Time: O(count * k^3) expected - one BPSW test per candidate
 *   Auxiliary memory: O(k/64) limbs
 *   Output memory: O(count) printed entries
 */
void bn_gen_proth_primes(u64 count, u64 k, u64 c)
{
  bignum p_bn, c_bn, two_k;
  bn_init_multi(&p_bn, &c_bn, &two_k, NULL);

  // 1. Calculate 2^k ONCE outside the loop
  bignum k_bn;
  bn_init(&k_bn);
  bn_set_u64(&k_bn, k);
  bn_pow(&two_k, &BN_TWO, &k_bn);
  bn_free(&k_bn);

  u64 curr_c = (c & 1) ? c : c + 1;
  u64 found = 0;

  printf("/* Generated %llu Proth Primes with k=%llu */\n", count, k);
  printf("static const uint64_t RNS_PRIMES[] = {\n");

  while (found < count) {
    // 2. p = curr_c * (precomputed 2^k)
    bn_set_u64(&c_bn, curr_c);
    bn_mul(&p_bn, &two_k, &c_bn);

    // 3. p = p + 1
    bn_add(&p_bn, &p_bn, &BN_ONE);

    if (bn_bpsw(&p_bn)) {
      // Using limbs[0] works IF your limbs are 64-bit.
      u64 prime = p_bn.limbs[0];
      printf("    %lluULL, // c=%llu\n", prime, curr_c);
      found++;
    }

    curr_c += 2;

    // 64-bit boundary check
    if (curr_c > (0xFFFFFFFFFFFFFFFF >> k)) {
      fprintf(stderr, "\nError: Search space exhausted.\n");
      break;
    }
  }

  printf("};\n");
  bn_free_multi(&p_bn, &c_bn, &two_k, NULL);
}

/*
 * Generate `count` primes just below 2^64 and print them as a C array
 * initializer.
 *
 * Scans downward from 2^64 - 47 in steps of 2 (so all candidates are
 * of the form 2^64 - 1 - 2i, i.e. close to the top of the 64-bit
 * range), testing each candidate with BPSW until `count` primes are
 * found. The output is a `static const u64 RNS_PRIMES[]` table
 * suitable for pasting into primes.h.
 *
 * Complexity:
 *   Time: O(count * 1024^3) expected - one BPSW test per candidate
 *   Auxiliary memory: O(1) limbs (single-limb candidates)
 *   Output memory: O(count) printed entries
 */
void gen_rns_primes(u64 count)
{
  u64 found = 0;

  u64 cand = 0xFFFFFFFFFFFFFFBULL;

  bignum tmp;
  bn_init(&tmp);

  bn_set_u64(&tmp, cand);

  printf("static const u64 RNS_PRIMES[] = {\n");
  while (found < count) {
    if (bn_bpsw(&tmp)) {
      bn_print(&tmp);
      printf(", \n");
      found++;
    }
    bn_sub(&tmp, &tmp, &BN_TWO);
  }

  printf("}; \n");
}
