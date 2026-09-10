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

#include "../../include/bigprime.h"

#include <fcntl.h>
#include <math.h>
#include <stdarg.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "../../include/bigcert.h"
#include "../../include/bigcore.h"
#include "../../include/bigfactor.h"
#include "../../include/biglucas.h"
#include "../../include/bigmath.h"
#include "../../include/bigprime.h"
#include "../../include/bigrabin.h"
#include "../../include/bigrand.h"
#include "../../include/primes.h"
#include "../../include/u64.h"

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
      if (p->limbs[0] == primes[i]) {
        close(fd);
        return true;
      }

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
    if (rounds == 5 && bn_is_square(NULL, n)) {
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
  i64 g;
  bool success;
  bn_set_u64(&two, 2);

  // trial division bound
  g = (u64)(c_opt * k * k + 1);

  double rel_size;
  do {
    rel_size = gen_rel_size();
  } while ((k * rel_size >= (k - margin)));

  printf("new size %lu \n", (u64)(rel_size * k));
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
        printf("   [FOUND] %lu-bit prime\n", k);
        break;
      }
    }
  }

  bn_copy(p, &n);

  bn_free_multi(&a, &n, &q, &I, &R, &twoI, &n_min1, &two, &two_q, &tmp, NULL);
  return true;
}

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

  printf("/* Generated %lu Proth Primes with k=%lu */\n", count, k);
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
      printf("    %luULL, // c=%lu\n", prime, curr_c);
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

u64 optimal_table_len(u64 n, u64 B)
{
  double k = 0.4;
  u64 m = n / B;
  return (u64)(k * n * m) / log(n * n * m);
}

void gen_provable_primes_arithmetic(bignum* p, u64 n,
                                    pocklington_cert** cert_out)
{
  const u64 num_bases = 20;

  if (n < 64) {
    bn_gen_prime(p, n);

    if (cert_out) {
      *cert_out = malloc(sizeof(pocklington_cert));
      bn_init(&(*cert_out)->N);
      bn_copy(&(*cert_out)->N, p);
      (*cert_out)->size = 0;
      (*cert_out)->capacity = 0;
      (*cert_out)->data = NULL;
    }

    return;
  }

  // F is a proven prime roughly of size 2^(n/2)
  bignum F;
  bn_init(&F);
  bool found_prime = false;

  pocklington_cert* F_cert = NULL;

  gen_provable_primes_arithmetic(&F, (n / 2) + 2, &F_cert);

  u64 s = optimal_table_len(n, 64);
  // printf("Table length: %lu \n", s);

  bignum t, A, B, exp, temp;
  bn_init_multi(&t, &A, &B, &exp, &temp, NULL);

  bignum a, N0, N, N_minus_1, test_val, gcd_val, pock_pow, X;
  bn_init_multi(&a, &N0, &N, &N_minus_1, &test_val, &gcd_val, &pock_pow, &X,
                NULL);

  bignum term, I, alpha;
  bn_init_multi(&term, &I, &alpha, NULL);

  while (!found_prime) {
    // Step 3: draw random number t in (2^(n-2) / F, 2^(n-1) / F - sn)

    bn_set_u64(&exp, n - 2);
    // A = 2^(n - 2) / F
    bn_pow(&A, &BN_TWO, &exp);
    bn_div(&A, &A, &F);

    // B = 2^(n - 1) / F - sn
    bn_mul(&B, &A, &BN_TWO);
    bn_set_u64(&temp, s * n);
    bn_sub(&B, &B, &temp);

    bn_gen_random_range(&t, &A, &B);

    // Step 4: Find a prime in the arithmetic progression
    // P = {N | N = N_0 + ia; N_0 = ta + 1; a = 2F; i <= i <= s}

    // a = 2F
    bn_mul(&a, &F, &BN_TWO);

    // N0 = t * a + 1
    bn_mul(&N0, &t, &a);
    bn_add_u64(&N0, &N0, 1);

    // Part 1: Trial division by primes < T_bound

    // tab[i] = 1 iff N'_p + ia'_p = 0 mod p forall p < T_bound
    bool tab[s + 1];
    memset(tab, 0, sizeof(tab));

    u64 T_bound = MAX(10 * n, 1000);

    for (u64 i = 0; i < 70000; i++) {
      u64 p = primes[i];
      if (p >= T_bound) {
        // printf("i: %lu \n", i);
        break;
      }

      u64 a_prime = bn_mod_u64(&a, p);
      u64 N0_prime = bn_mod_u64(&N0, p);

      // Mark solutions to N0_prime + i * a_prime == 0 (mod p)
      for (u64 j = 0; j <= s; j++) {
        if (tab[j]) continue;
        if ((N0_prime + j * a_prime) % p == 0) {
          tab[j] = true;
        }
      }
    }

    // Part II & III: Compositeness Test and Primality Proof
    for (u64 i = 0; i <= s; i++) {
      if (found_prime || tab[i]) continue;  // Skip sieved candidates

      // N = N0 + i * a
      bn_set_u64(&I, i);
      bn_mul(&term, &a, &I);
      bn_add(&N, &N0, &term);
      bn_sub(&N_minus_1, &N, &BN_ONE);

      // Initialize the ctx once
      bn_mont_ctx ctx;
      bn_mont_ctx_init(&ctx, &N);

      // Part II: Rabin-Miller test with base 2
      if (!bn_rabin_mont_ctx(&N, &BN_TWO, &ctx)) {
        bn_mont_ctx_free(&ctx);
        continue;
      }

      // Part III: Primality proof using Pocklington lemma
      // Since F is a prime from Step 2, q = F. We seek a base alpha_q.
      bn_div(&pock_pow, &N_minus_1, &F);

      for (u64 b = 0; b < num_bases; b++) {
        bn_set_u64(&alpha, primes[b]);

        // X = alpha^((N- 1)/ F) mod N
        bn_mod_exp_mont(&X, &alpha, &pock_pow, &N, &ctx);

        // alpha^(N-1) = X^F == 1 mod N (Little Fermat)
        bn_mod_exp_mont(&test_val, &X, &F, &N, &ctx);
        if (bn_cmp(&test_val, &BN_ONE) != 0) {
          break;
        }

        // gcd(alpha^((N-1)/F) - 1, N) = gcd(X - 1, N) == 1
        bn_sub(&test_val, &X, &BN_ONE);
        bn_gcd(&gcd_val, &test_val, &N);

        if (bn_cmp(&gcd_val, &BN_ONE) == 0) {
          found_prime = true;
          break;  // Base found, N is verified prime
        }
      }
      bn_mont_ctx_free(&ctx);

      if (found_prime) {
        bn_copy(p, &N);

        if (cert_out) {
          // Build the certificate for current N
          *cert_out = malloc(sizeof(pocklington_cert));
          bn_init(&(*cert_out)->N);
          bn_copy(&(*cert_out)->N, &N);

          // Allocate space for the single factor F we are proving against
          (*cert_out)->size = 1;
          (*cert_out)->capacity = 1;
          (*cert_out)->data = malloc(sizeof(pocklington_cert_elem));

          bn_init(&(*cert_out)->data[0].q);
          bn_copy(&(*cert_out)->data[0].q, &F);

          bn_init(&(*cert_out)->data[0].alpha_q);
          bn_copy(&(*cert_out)->data[0].alpha_q, &alpha);

          // Link the recursive proof for F
          (*cert_out)->data[0].q_cert = F_cert;
        } else {
          pocklington_cert_free(F_cert);
        }

        break;
      }
    }
  }

  bn_free_multi(&term, &I, &alpha, NULL);
  bn_free_multi(&t, &A, &B, &exp, &temp, NULL);
  bn_free_multi(&a, &N0, &N, &N_minus_1, &test_val, &gcd_val, &pock_pow, &X,
                NULL);

  bn_free(&F);
}
