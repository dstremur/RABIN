/*
 * rzprime.c
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

#include "../../include/rzprime.h"

#include <fcntl.h>
#include <limits.h>
#include <math.h>
#include <stdarg.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "../../include/rns_primes.h"
#include "../../include/rz.h"
#include "../../include/rzcert.h"
#include "../../include/rzfactor.h"
#include "../../include/rzlucas.h"
#include "../../include/rzmath.h"
#include "../../include/rzprime.h"
#include "../../include/rzrabin.h"
#include "../../include/rzrand.h"
#include "../../include/u64.h"

rabin_err_t rz_gen_prime(rz_t* p, int bits)
{
  if (p == NULL) return RABIN_ERR_NULL_PTR;
  if (bits < 2) return RABIN_ERR_INVALID_ARG;

  int fd = open("/dev/urandom", O_RDONLY);
  if (fd < 0) {
    return RABIN_ERR_INVALID_ARG;
  }

  rabin_err_t err = RABIN_SUCCESS;
  while (true) {
    bool composite = false;

    err = rz_gen_random_with_fd(p, (u64)bits, fd);
    if (err != RABIN_SUCCESS) break;

    rz_set_bit(p, bits - 1);
    rz_set_bit(p, 0);

    uint64_t multi_prime = 3ULL * 5 * 7 * 11 * 13 * 17;
    uint64_t rem = rz_mod_u64(p, multi_prime);

    if (rem % 3 == 0 || rem % 5 == 0 || rem % 7 == 0 || rem % 11 == 0 ||
        rem % 13 == 0 || rem % 17 == 0)
      continue;

    for (int i = 0; i < 1500; i++) {
      if (p->limbs[0] == SMALL_PRIMES[i]) {
        goto done;
      }

      if (rz_mod_u64(p, SMALL_PRIMES[i]) == 0) {
        composite = true;
        break;
      }
    }

    // 2. BPSW deterministic :)
    if (!composite) {
      if (rz_bpsw(p)) {
        goto done;
      }
    }
  }

done:
  close(fd);
  return err;
}

rabin_err_t rz_gen_safe_prime(rz_t* p, int bits)
{
  if (p == NULL) return RABIN_ERR_NULL_PTR;
  if (bits < 3) return RABIN_ERR_INVALID_ARG;

  rz_t q;
  rz_init(&q);
  rabin_err_t err = RABIN_SUCCESS;
  do {
    if ((err = rz_gen_prime(&q, bits - 1)) != RABIN_SUCCESS) break;

    // p = 2q + 1
    if ((err = rz_lshift(p, &q, 1)) != RABIN_SUCCESS) break;
    if ((err = rz_add_u64(p, p, 1)) != RABIN_SUCCESS) break;

  } while (!rz_bpsw(p));

  rz_clear(&q);
  return err;
}

bool rz_stronglucas(const rz_t* n, const rz_t* P, const rz_t* Q)
{
  if (n == NULL || P == NULL || Q == NULL) return false;

  rz_t p, q, d, u, v, qn, tmp, n_plus_1;
  rz_init_multi(&p, &q, &d, &u, &v, &qn, &tmp, &n_plus_1, NULL);

  bool prime = false;

  int s = 0;

  if (rz_copy(&n_plus_1, n) != RABIN_SUCCESS) goto cleanup;
  if (rz_add_u64(&n_plus_1, &n_plus_1, 1) != RABIN_SUCCESS) goto cleanup;
  if (rz_copy(&d, &n_plus_1) != RABIN_SUCCESS) goto cleanup;

  while (rz_is_even(&d)) {
    if (rz_rshift1(&d) != RABIN_SUCCESS) goto cleanup;
    s++;
  }

  if (rz_lucas_solve_mod(&u, &v, P, Q, &qn, &d, n) != RABIN_SUCCESS)
    goto cleanup;
  // 1. condition
  if (rz_is_zero(&u)) {
    prime = true;
    goto cleanup;
  }
  // 2. condition
  if (rz_is_zero(&v)) {
    prime = true;
    goto cleanup;
  }
  // check 2. condition for all d * 2^r
  for (int r = 1; r < s; r++) {
    if (rz_mul(&tmp, &v, &v) != RABIN_SUCCESS) goto cleanup;
    if (rz_mod(&tmp, &tmp, n) != RABIN_SUCCESS) goto cleanup;

    if (rz_add(&tmp, &tmp, n) != RABIN_SUCCESS) goto cleanup;
    if (rz_sub(&tmp, &tmp, &qn) != RABIN_SUCCESS) goto cleanup;
    if (rz_add(&tmp, &tmp, n) != RABIN_SUCCESS) goto cleanup;
    if (rz_sub(&tmp, &tmp, &qn) != RABIN_SUCCESS) goto cleanup;
    if (rz_mod(&v, &tmp, n) != RABIN_SUCCESS) goto cleanup;

    if (rz_mul(&qn, &qn, &qn) != RABIN_SUCCESS) goto cleanup;
    if (rz_mod(&qn, &qn, n) != RABIN_SUCCESS) goto cleanup;

    if (rz_is_eq_i64(&v, 0)) {
      prime = true;
      goto cleanup;
    }
  }

cleanup:
  rz_clear(&p);
  rz_clear(&q);
  rz_clear(&d);
  rz_clear(&u);
  rz_clear(&v);
  rz_clear(&qn);
  rz_clear(&tmp);
  rz_clear(&n_plus_1);
  return prime;
}

bool rz_bpsw(const rz_t* n)
{
  if (n == NULL) return false;
  if (rz_is_even(n)) return false;

  // Trial division by small primes to quickly reject composites (mirrors
  // GMP's fast path). ~88% of random odd numbers have a factor < 10000.
  for (u64 i = 0; SMALL_PRIMES[i] < 10000; i++) {
    u64 p = SMALL_PRIMES[i];
    if (n->size == 1 && n->limbs[0] == p) return true;
    if (rz_mod_u64(n, p) == 0) return false;
  }

  // 1. run miller rabin base 2
  rz_t two;
  rz_init(&two);
  if (rz_set_u64(&two, 2) != RABIN_SUCCESS) {
    rz_clear(&two);
    return false;
  }
  if (!rz_rabin_mont(n, &two)) {
    rz_clear(&two);
    return false;
  }

  rz_t D, magnitude;
  rz_init(&D);
  rz_init(&magnitude);
  if (rz_set_i64(&magnitude, 5) != RABIN_SUCCESS) goto reject;

  bool negative = false;

  int rounds = 0;
  // finds a D using Selfridges method A*
  while (1) {
    // check if n is perfect square after 5 rounds
    if (rounds == 5 && rz_is_square(NULL, n)) {
      goto reject;
    }

    if (rz_copy(&D, &magnitude) != RABIN_SUCCESS) goto reject;
    D.is_neg = negative;

    // rz_jacobi() returns 0 for both "no solution" and failure; either way
    // this candidate is rejected
    i64 jacobi = rz_jacobi(&D, n);
    if (jacobi == -1) {
      break;
    }

    if (jacobi == 0) {
      goto reject;
    }

    if (rz_add_u64(&magnitude, &magnitude, 2) != RABIN_SUCCESS) goto reject;
    negative = !negative;

    rounds++;
  }

  // zero-initialized so the reject path can clear them even when the jump
  // lands before rz_init_multi() ran
  rz_t P = {0}, Q = {0};
  rz_init_multi(&P, &Q, NULL);
  if (rz_set_u64(&P, 1) != RABIN_SUCCESS) goto reject;
  if (rz_set_u64(&Q, 1) != RABIN_SUCCESS) goto reject;

  if (negative) {
    if (rz_copy(&Q, &magnitude) != RABIN_SUCCESS) goto reject;
    if (rz_add_u64(&Q, &Q, 1) != RABIN_SUCCESS) goto reject;
    rz_divmod_u64(&Q, &Q, 4);
  } else {
    if (rz_copy(&Q, &magnitude) != RABIN_SUCCESS) goto reject;
    if (rz_sub(&Q, &Q, &P) != RABIN_SUCCESS) goto reject;
    rz_divmod_u64(&Q, &Q, 4);
  }

  bool res = rz_stronglucas(n, &P, &Q);

  rz_clear(&P);
  rz_clear(&Q);
  rz_clear(&D);
  rz_clear(&magnitude);
  rz_clear(&two);
  return res;
reject:
  rz_clear(&P);
  rz_clear(&Q);
  rz_clear(&D);
  rz_clear(&magnitude);
  rz_clear(&two);
  return false;
}

bool rz_check_lemma1(const rz_t* n, const rz_t* n_min1, const rz_t* a,
                     const rz_t* q)
{
  if (n == NULL || n_min1 == NULL || a == NULL || q == NULL) return false;

  rz_t tmp, exp, X, gcd;

  rz_init_multi(&tmp, &exp, &X, &gcd, NULL);

  bool result = false;

  if (rz_copy(&exp, n_min1) != RABIN_SUCCESS) goto cleanup;
  if (rz_div(&exp, &exp, q) != RABIN_SUCCESS) goto cleanup;
  if (rz_copy(&tmp, a) != RABIN_SUCCESS) goto cleanup;

  if (rz_mod_exp(&X, &tmp, &exp, n) != RABIN_SUCCESS) goto cleanup;

  if (rz_is_eq_i64(&X, 1)) {
    goto cleanup;
  }

  if (rz_set_u64(&tmp, 1) != RABIN_SUCCESS) goto cleanup;
  if (rz_sub(&tmp, &X, &tmp) != RABIN_SUCCESS) goto cleanup;

  if (rz_gcd(&gcd, &tmp, n) != RABIN_SUCCESS) goto cleanup;
  if (!rz_is_eq_i64(&gcd, 1)) {
    goto cleanup;
  }

  if (rz_mod_exp(&X, &X, q, n) != RABIN_SUCCESS) goto cleanup;

  result = rz_is_eq_i64(&X, 1);

cleanup:
  rz_clear_multi(&tmp, &exp, &X, &gcd, NULL);
  return result;
}

double rz_gen_rel_size()
{
  // Generate uniform random variable u in [0, 1]
  double u = (double)rand() / RAND_MAX;

  return pow(2.0, u - 1.0);
}

rabin_err_t rz_provable_prime(rz_t* p, u64 k)
{
  if (p == NULL) return RABIN_ERR_NULL_PTR;
  if (k < 2) return RABIN_ERR_INVALID_ARG;

  // Keep trying from the absolute top until it succeeds.
  // Every rejected attempt has cleanly unwound the stack and freed all
  // memory, so the retry is safe.
  rabin_err_t err;
  while ((err = rz_provable_prime_inner(p, k)) == RABIN_ERR_INVALID_ARG) {
    // candidate rejected; retry from scratch
  }
  return err;
}

rabin_err_t rz_provable_prime_inner(rz_t* p, u64 k)
{
  if (p == NULL) return RABIN_ERR_NULL_PTR;

  //  base case k <= 20
  //  use Baillie-PSW instead of trial factoring
  if (k <= 20) {
    rabin_err_t err = RABIN_SUCCESS;
    do {
      // generates a random odd k-bit integer
      if ((err = rz_gen_random(p, k)) != RABIN_SUCCESS) break;
    } while (!rz_bpsw(p));
    return err;
  }

  // constants
  const double c_opt = 0.1;
  u64 margin = k / 6;

  rz_t a, n, q, I, R, twoI, n_min1, two, two_q, tmp;
  rz_init_multi(&a, &n, &q, &I, &R, &twoI, &n_min1, &two, &two_q, &tmp, NULL);

  rabin_err_t err = rz_set_u64(&two, 2);
  if (err != RABIN_SUCCESS) goto out;

  // trial division bound
  u64 g = (u64)(c_opt * k * k + 1);

  double rel_size;
  do {
    rel_size = rz_gen_rel_size();
  } while ((k * rel_size >= (k - margin)));

  // recursive call
  // RABIN_ERR_INVALID_ARG means the recursive candidate was rejected, which
  // rejects this whole attempt
  err = rz_provable_prime_inner(&q, (u64)(rel_size * k));
  if (err != RABIN_SUCCESS) goto out;

  if ((err = rz_copy(&two_q, &q)) != RABIN_SUCCESS) goto out;
  if ((err = rz_lshift1(&two_q)) != RABIN_SUCCESS) goto out;

  // I = 2^(k-1) / 2q
  if ((err = rz_set_u64(&I, 1)) != RABIN_SUCCESS) goto out;
  if ((err = rz_lshift(&I, &I, k - 1)) != RABIN_SUCCESS) goto out;
  if ((err = rz_div(&I, &I, &two_q)) != RABIN_SUCCESS) goto out;

  // twoI = 2^k / 2q
  if ((err = rz_set_u64(&twoI, 1)) != RABIN_SUCCESS) goto out;
  if ((err = rz_lshift(&twoI, &twoI, k)) != RABIN_SUCCESS) goto out;
  if ((err = rz_div(&twoI, &twoI, &two_q)) != RABIN_SUCCESS) goto out;

  bool success = false;
  u64 attempts = 0;
  while (!success) {
    attempts++;
    if (attempts > 1000) {
      err = RABIN_ERR_INVALID_ARG;
      goto out;
    }

    if ((err = rz_gen_random_range(&R, &I, &twoI)) != RABIN_SUCCESS) goto out;

    // n = 2 * rand(I, 2I) * q + 1
    if ((err = rz_mul(&n, &R, &q)) != RABIN_SUCCESS) goto out;
    if ((err = rz_lshift1(&n)) != RABIN_SUCCESS) goto out;
    if ((err = rz_copy(&n_min1, &n)) != RABIN_SUCCESS) goto out;
    if ((err = rz_add_u64(&n, &n, 1)) != RABIN_SUCCESS) goto out;

    if (rz_trialdiv(&n, g)) {
      for (int j = 0; j < 200; j++) {
        if ((err = rz_gen_random_range(&a, &two, &n_min1)) != RABIN_SUCCESS)
          goto out;

        if (rz_check_lemma1(&n, &n_min1, &a, &q)) {
          success = true;
          break;
        }
      }
    }
  }

  if ((err = rz_copy(p, &n)) != RABIN_SUCCESS) goto out;

  err = RABIN_SUCCESS;
out:
  rz_clear_multi(&a, &n, &q, &I, &R, &twoI, &n_min1, &two, &two_q, &tmp, NULL);
  return err;
}

rabin_err_t rz_gen_proth_primes(u64 count, u64 k, u64 c)
{
  if (count == 0) return RABIN_ERR_INVALID_ARG;
  if (k >= 64) return RABIN_ERR_INVALID_ARG;

  rz_t p_bn, c_bn, two_k;
  rz_init_multi(&p_bn, &c_bn, &two_k, NULL);

  // 1. Calculate 2^k ONCE outside the loop
  rz_t k_bn;
  rz_init(&k_bn);
  rabin_err_t err = rz_set_u64(&k_bn, k);
  if (err == RABIN_SUCCESS) {
    err = rz_pow(&two_k, &RZ_TWO, &k_bn);
  }
  rz_clear(&k_bn);
  if (err != RABIN_SUCCESS) goto out;

  u64 curr_c = (c & 1) ? c : c + 1;
  u64 found = 0;

  printf("/* Generated %lu Proth Primes with k=%lu */\n", count, k);
  printf("static const uint64_t RNS_PRIMES[] = {\n");

  while (found < count) {
    // 2. p = curr_c * (precomputed 2^k)
    if ((err = rz_set_u64(&c_bn, curr_c)) != RABIN_SUCCESS) break;
    if ((err = rz_mul(&p_bn, &two_k, &c_bn)) != RABIN_SUCCESS) break;

    // 3. p = p + 1
    if ((err = rz_add(&p_bn, &p_bn, &RZ_ONE)) != RABIN_SUCCESS) break;

    if (rz_bpsw(&p_bn)) {
      // Using limbs[0] works IF your limbs are 64-bit.
      u64 prime = p_bn.limbs[0];
      printf("    %luULL, // c=%lu\n", prime, curr_c);
      found++;
    }

    curr_c += 2;

    // 64-bit boundary check
    if (curr_c > (UINT64_MAX >> k)) {
      printf("};\n");
      err = RABIN_ERR_OVERFLOW;
      goto out;
    }
  }

  printf("};\n");
  err = RABIN_SUCCESS;
out:
  rz_clear_multi(&p_bn, &c_bn, &two_k, NULL);
  return err;
}

rabin_err_t rns_gen_primes(u64 count)
{
  if (count == 0) return RABIN_ERR_INVALID_ARG;

  u64 found = 0;

  u64 cand = 0xFFFFFFFFFFFFFFBULL;

  rz_t tmp;
  rz_init(&tmp);
  rabin_err_t err = rz_set_u64(&tmp, cand);
  if (err != RABIN_SUCCESS) goto out;

  printf("static const u64 RNS_PRIMES[] = {\n");
  while (found < count) {
    if (tmp.size == 0 || tmp.limbs[0] < 2) {
      // search space exhausted (underflow guard)
      err = RABIN_ERR_OVERFLOW;
      break;
    }
    if (rz_bpsw(&tmp)) {
      rz_print(&tmp);
      printf(", \n");
      found++;
    }
    if ((err = rz_sub(&tmp, &tmp, &RZ_TWO)) != RABIN_SUCCESS) break;
  }

  printf("}; \n");
out:
  rz_clear(&tmp);
  return err;
}

u64 rz_optimal_table_len(u64 n, u64 B)
{
  double k = 0.4;
  u64 m = n / B;
  return (u64)(k * n * m) / log(n * n * m);
}

rabin_err_t rz_gen_provable_arithmetic(rz_t* p, u64 n, rzcert_t** cert_out)
{
  const u64 num_bases = 20;
  if (p == NULL) return RABIN_ERR_NULL_PTR;
  if (cert_out != NULL) *cert_out = NULL;

  if (n < 64) {
    if (n > (u64)INT_MAX) return RABIN_ERR_INVALID_ARG;
    rabin_err_t err = rz_gen_prime(p, (int)n);
    if (err != RABIN_SUCCESS) return err;

    if (cert_out != NULL) {
      *cert_out = malloc(sizeof(rzcert_t));
      if (*cert_out == NULL) return RABIN_ERR_OUT_OF_MEMORY;
      **cert_out = (rzcert_t){0};
      rz_init(&(*cert_out)->N);
      if ((err = rz_copy(&(*cert_out)->N, p)) != RABIN_SUCCESS) {
        rz_clear(&(*cert_out)->N);
        free(*cert_out);
        *cert_out = NULL;
        return err;
      }
      (*cert_out)->size = 0;
      (*cert_out)->capacity = 0;
      (*cert_out)->data = NULL;
    }

    return RABIN_SUCCESS;
  }

  // F is a proven prime roughly of size 2^(n/2)
  rz_t F;
  rz_init(&F);

  rzcert_t* F_cert = NULL;

  rabin_err_t err = rz_gen_provable_arithmetic(&F, (n / 2) + 2, &F_cert);
  if (err != RABIN_SUCCESS) {
    rz_clear(&F);
    return err;
  }

  u64 s = rz_optimal_table_len(n, 64);
  // the sieve table is a VLA; cap its size so the stack stays sane
  if (s == 0 || s > 4096 * 1024) {
    rzcert_clear(F_cert);
    rz_clear(&F);
    return RABIN_ERR_INVALID_ARG;
  }
  if (s > UINT64_MAX / n) {
    rzcert_clear(F_cert);
    rz_clear(&F);
    return RABIN_ERR_OVERFLOW;
  }

  rz_t t, A, B, exp, temp;
  rz_init_multi(&t, &A, &B, &exp, &temp, NULL);

  rz_t a, N0, N, N_minus_1, test_val, gcd_val, pock_pow, X;
  rz_init_multi(&a, &N0, &N, &N_minus_1, &test_val, &gcd_val, &pock_pow, &X,
                NULL);

  rz_t term, I, alpha;
  rz_init_multi(&term, &I, &alpha, NULL);

  bool found_prime = false;
  while (!found_prime) {
    // Step 3: draw random number t in (2^(n-2) / F, 2^(n-1) / F - sn)

    if ((err = rz_set_u64(&exp, n - 2)) != RABIN_SUCCESS) goto out;
    // A = 2^(n - 2) / F
    if ((err = rz_pow(&A, &RZ_TWO, &exp)) != RABIN_SUCCESS) goto out;
    if ((err = rz_div(&A, &A, &F)) != RABIN_SUCCESS) goto out;

    // B = 2^(n - 1) / F - sn
    if ((err = rz_mul(&B, &A, &RZ_TWO)) != RABIN_SUCCESS) goto out;
    if ((err = rz_set_u64(&temp, s * n)) != RABIN_SUCCESS) goto out;
    if ((err = rz_sub(&B, &B, &temp)) != RABIN_SUCCESS) goto out;

    if ((err = rz_gen_random_range(&t, &A, &B)) != RABIN_SUCCESS) goto out;

    // Step 4: Find a prime in the arithmetic progression
    // P = {N | N = N_0 + ia; N_0 = ta + 1; a = 2F; i <= i <= s}

    // a = 2F
    if ((err = rz_mul(&a, &F, &RZ_TWO)) != RABIN_SUCCESS) goto out;

    // N0 = t * a + 1
    if ((err = rz_mul(&N0, &t, &a)) != RABIN_SUCCESS) goto out;
    if ((err = rz_add_u64(&N0, &N0, 1)) != RABIN_SUCCESS) goto out;

    // Part 1: Trial division by primes < T_bound

    // tab[i] = 1 iff N'_p + ia'_p = 0 mod p forall p < T_bound
    bool tab[s + 1];
    memset(tab, 0, sizeof(tab));

    u64 T_bound = RZ_MAX(10 * n, 1000);

    for (u64 i = 0; i < 70000; i++) {
      u64 p_small = SMALL_PRIMES[i];
      if (p_small >= T_bound) {
        break;
      }

      u64 a_prime = rz_mod_u64(&a, p_small);
      u64 N0_prime = rz_mod_u64(&N0, p_small);

      // Mark solutions to N0_prime + i * a_prime == 0 (mod p)
      for (u64 j = 0; j <= s; j++) {
        if (tab[j]) continue;
        if ((N0_prime + j * a_prime) % p_small == 0) {
          tab[j] = true;
        }
      }
    }

    // Part II & III: Compositeness Test and Primality Proof
    for (u64 i = 0; i <= s; i++) {
      if (found_prime || tab[i]) continue;  // Skip sieved candidates

      // N = N0 + i * a
      if ((err = rz_set_u64(&I, i)) != RABIN_SUCCESS) goto out;
      if ((err = rz_mul(&term, &a, &I)) != RABIN_SUCCESS) goto out;
      if ((err = rz_add(&N, &N0, &term)) != RABIN_SUCCESS) goto out;
      if ((err = rz_sub(&N_minus_1, &N, &RZ_ONE)) != RABIN_SUCCESS) goto out;

      // Initialize the ctx once
      rz_mont_ctx ctx;
      if ((err = rz_mont_ctx_init(&ctx, &N)) != RABIN_SUCCESS) {
        if (err == RABIN_ERR_OUT_OF_MEMORY) goto out;
        // candidate rejected (e.g. even modulus); try the next one
        continue;
      }

      // Part II: Rabin-Miller test with base 2
      if (!rz_rabin_mont_ctx(&N, &RZ_TWO, &ctx)) {
        rz_mont_ctx_clear(&ctx);
        continue;
      }

      // Part III: Primality proof using Pocklington lemma
      // Since F is a prime from Step 2, q = F. We seek a base alpha_q.
      if ((err = rz_div(&pock_pow, &N_minus_1, &F)) != RABIN_SUCCESS) {
        rz_mont_ctx_clear(&ctx);
        goto out;
      }

      for (u64 b = 0; b < num_bases; b++) {
        if ((err = rz_set_u64(&alpha, SMALL_PRIMES[b])) != RABIN_SUCCESS) break;

        // X = alpha^((N- 1)/ F) mod N
        if ((err = rz_mod_exp_mont(&X, &alpha, &pock_pow, &N, &ctx)) !=
            RABIN_SUCCESS)
          break;

        // alpha^(N-1) = X^F == 1 mod N (Little Fermat)
        if ((err = rz_mod_exp_mont(&test_val, &X, &F, &N, &ctx)) !=
            RABIN_SUCCESS)
          break;
        if (rz_cmp(&test_val, &RZ_ONE) != 0) {
          break;
        }

        // gcd(alpha^((N-1)/F) - 1, N) = gcd(X - 1, N) == 1
        if ((err = rz_sub(&test_val, &X, &RZ_ONE)) != RABIN_SUCCESS) break;
        if ((err = rz_gcd(&gcd_val, &test_val, &N)) != RABIN_SUCCESS) break;

        if (rz_cmp(&gcd_val, &RZ_ONE) == 0) {
          found_prime = true;
          break;  // Base found, N is verified prime
        }
      }
      rz_mont_ctx_clear(&ctx);
      if (err != RABIN_SUCCESS) goto out;

      if (found_prime) {
        if (cert_out != NULL) {
          // Build the certificate for current N before touching p
          *cert_out = malloc(sizeof(rzcert_t));
          if (*cert_out == NULL) {
            err = RABIN_ERR_OUT_OF_MEMORY;
            goto out;
          }
          **cert_out = (rzcert_t){0};
          rz_init(&(*cert_out)->N);
          if ((err = rz_copy(&(*cert_out)->N, &N)) != RABIN_SUCCESS)
            goto fail_cert;

          // Allocate space for the single factor F we are proving against
          (*cert_out)->size = 1;
          (*cert_out)->capacity = 1;
          (*cert_out)->data = malloc(sizeof(rzcert_elem_t));
          if ((*cert_out)->data == NULL) {
            err = RABIN_ERR_OUT_OF_MEMORY;
            goto fail_cert;
          }
          rz_init(&(*cert_out)->data[0].q);
          rz_init(&(*cert_out)->data[0].alpha_q);
          if ((err = rz_copy(&(*cert_out)->data[0].q, &F)) != RABIN_SUCCESS)
            goto fail_cert;

          if ((err = rz_copy(&(*cert_out)->data[0].alpha_q, &alpha)) !=
              RABIN_SUCCESS)
            goto fail_cert;

          // Link the recursive proof for F (ownership transfers to the
          // certificate)
          (*cert_out)->data[0].q_cert = F_cert;
          F_cert = NULL;
        } else {
          rzcert_clear(F_cert);
          F_cert = NULL;
        }

        if ((err = rz_copy(p, &N)) != RABIN_SUCCESS) goto out;

        err = RABIN_SUCCESS;
        break;
      }
    }
  }

  err = RABIN_SUCCESS;
out:
  if (F_cert != NULL) {
    rzcert_clear(F_cert);
  }
  rz_clear_multi(&term, &I, &alpha, NULL);
  rz_clear_multi(&t, &A, &B, &exp, &temp, NULL);
  rz_clear_multi(&a, &N0, &N, &N_minus_1, &test_val, &gcd_val, &pock_pow, &X,
                 NULL);

  rz_clear(&F);
  return err;
fail_cert:
  // partially built certificate; q_cert has not been linked yet
  if ((*cert_out)->data != NULL) {
    if ((*cert_out)->size > 0) {
      rz_clear(&(*cert_out)->data[0].q);
      rz_clear(&(*cert_out)->data[0].alpha_q);
    }
    free((*cert_out)->data);
  }
  rz_clear(&(*cert_out)->N);
  free(*cert_out);
  *cert_out = NULL;
  goto out;
}
