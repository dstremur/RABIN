/*
 * main.c
 *
 * Interactive demo driver for the rz_t library.
 *
 * Reads two bignums from stdin and exercises the library: basic
 * arithmetic (add, sub, mul, fast mul, div, mod, shifts), square
 * root and natural log, polynomial self-tests, factorization,
 * Pollard rho / p-1, BPSW primality, Jacobi symbol, Tonelli-Shanks
 * square roots, and 2048-bit prime generation. Also prints a debug
 * dump of a rz_t NTT context.
 *
 * This is a manual playground, not part of the library API.
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
#include <assert.h>
#include <fcntl.h>
#include <stdint.h>
#include <stdio.h>
#include <unistd.h>

#include "../include/rabin.h"
#include "../include/rmat.h"
#include "../include/rntt.h"
#include "ctype.h"
#include "string.h"

typedef struct {
  const char* n_hex;
  const char* base_hex;
  bool expected;
  const char* description;
} test_case;

/**
 * @brief Print a debug dump of a rz_t NTT context.
 *
 * Shows the transform length, modulus, n inverse, selected omega /
 * omega_inv / psi powers (with expected values annotated), and the
 * first 8 bit-reversal indices. Intended for manual verification of
 * context initialization.
 *
 * Complexity:
 *   Time: O(n) where n is the transform length (printing the tables)
 *   Auxiliary memory: O(1)
 *   Output memory: O(n) characters written
 *
 * @param[in] ctx NTT context to print.
 */
void rntt_ctx_print_debug(const rntt_ctx_t* ctx)
{
  if (!ctx) {
    printf("Context is NULL\n");
    return;
  }

  printf("=========================================\n");
  printf("         NTT CONTEXT DEBUG LOG           \n");
  printf("=========================================\n");
  printf("Transform Length (n) : %lu\n", ctx->n);
  printf("Modulus (q)          : ");
  rz_println(&ctx->q);
  printf("N Inverse (n_inv)    : ");
  rz_println(&ctx->n_inv);

  printf("\n--- OMEGA POWERS (w^i mod q) ---\n");
  printf("omega[0]       (Expected: 1)   : ");
  rz_println(&ctx->omega_powers[0]);
  printf("omega[1]       (Base Root)     : ");
  rz_println(&ctx->omega_powers[1]);
  printf("omega[n/2]     (Expected: q-1) : ");
  rz_println(&ctx->omega_powers[ctx->n / 2]);
  printf("omega[n-1]                     : ");
  rz_println(&ctx->omega_powers[ctx->n - 1]);

  printf("\n--- OMEGA INVERSE POWERS (w^-i mod q) ---\n");
  printf("omega_inv[0]   (Expected: 1)   : ");
  rz_println(&ctx->omega_inv_powers[0]);
  printf("omega_inv[1]   (Base Inv Root) : ");
  rz_println(&ctx->omega_inv_powers[1]);
  printf("omega_inv[n/2] (Expected: q-1) : ");
  rz_println(&ctx->omega_inv_powers[ctx->n / 2]);

  if (ctx->psi_powers != NULL) {
    printf("\n--- PSI POWERS (psi^i mod q) ---\n");
    printf("psi[0]         (Expected: 1)   : ");
    rz_println(&ctx->psi_powers[0]);
    printf("psi[1]         (Base Psi Root) : ");
    rz_println(&ctx->psi_powers[1]);
    printf("psi[n/2]       (Square rt of -1): ");
    rz_println(&ctx->psi_powers[ctx->n / 2]);
  }

  printf("\n--- BIT REVERSAL INDICES (First 8 Samples) ---\n");
  u64 max_samples = (ctx->n < 8) ? ctx->n : 8;
  for (u64 i = 0; i < max_samples; i++) {
    printf("  Index [%lu] ---> Bit-Reversed Index [%lu]\n", i,
           ctx->bit_rev_indices[i]);
  }
  printf("=========================================\n");
}

/**
 * @brief Demo entry point.
 *
 * Reads two bignums a and b from stdin, then runs a sequence of
 * library calls printing each result. Also initializes a Goldilocks
 * NTT context (size 2^8) for the debug dump, factorizes a, runs
 * Pollard rho and p-1 on a, tests a with BPSW, computes the Jacobi
 * symbol (a, b), a Tonelli-Shanks square root of a mod b, and
 * generates a 2048-bit prime.
 *
 * Complexity:
 *   Time: dominated by the 2048-bit prime generation and the
 *         factorization attempts on a
 *   Auxiliary memory: O(size of a and b) limbs
 *   Output memory: O(1)
 *
 * @return 0 On success, 1 if input reading fails.
 */
int main()
{
  rz_init_constants();

  // rns_gen_primes(500);
  rz_t a, b, c, psi, omega;
  rz_init_multi(&a, &b, &c, &psi, &omega, NULL);
  char buf1[1024];
  char buf2[1024];

  printf("Enter first big number (a): ");
  if (!fgets(buf1, sizeof(buf1), stdin)) return 1;
  buf1[strcspn(buf1, "\n")] = 0;  // Remove newline

  printf("Enter second big number (b, p): ");
  if (!fgets(buf2, sizeof(buf2), stdin)) return 1;
  buf2[strcspn(buf2, "\n")] = 0;  // Remove newline

  rz_init_val(&a, buf1);
  rz_init_val(&b, buf2);
  printf("Inputs: \n");
  printf("a: ");
  rz_println(&a);
  printf("b: ");
  rz_println(&b);
  printf("ADD: ");
  rz_add(&c, &a, &b);
  rz_println(&c);
  rz_sub(&c, &a, &b);
  printf("SUB: ");
  rz_println(&c);
  rz_isqrt(&c, &a);
  printf("ISQRT (a): ");
  rz_println(&c);
  rz_ln(&c, &a);
  printf("LN (a): ");
  rz_println(&c);
  rz_mul(&c, &a, &b);
  printf("MUL: ");
  rz_println(&c);
  rz_mul_fast(&c, &a, &b);
  printf("MUL NTT: ");
  rz_println(&c);
  rz_div(&c, &a, &b);
  printf("DIV: ");
  rz_println(&c);
  rz_div_euclid(&c, &a, &b);
  printf("DIV euclid: ");
  rz_println(&c);
  rz_div_exact(&c, &a, &b);
  printf("DIV exact: ");
  rz_println(&c);
  rz_mod(&c, &a, &b);
  printf("MOD: ");
  rz_println(&c);
  printf("RSHIFT: ");
  rz_rshift(&c, &a, 64);
  rz_println(&c);
  printf("LSHIFT: ");
  rz_lshift(&c, &a, 32);
  rz_println(&c);
  printf("Pow: ");
  // rz_pow(&c, &a, &b);
  rz_println(&c);
  rz_mul_i64(&c, &a, 3);
  printf("3 * a: ");
  rz_println(&c);

  rz_t x, y, z;
  rz_init_multi(&x, &y, &z);

  rz_gcd_extended(&x, &y, &z, &a, &b);

  printf("u: ");
  rz_println(&x);
  printf("v: ");
  rz_println(&y);
  printf("d: ");
  rz_println(&z);

  if (rz_is_square(&c, &a)) {
    printf("a is square \n");
  } else {
    printf("a is not square \n");
  }

  if (rz_is_prime_power(&c, &a)) {
    printf("a is prime power \n");
  } else {
    printf("a is not a prime power\n");
  }

  i64 j = rz_jacobi(&a, &b);

  printf("JACOBI: %lli\n", (long long)j);

  j = rz_kronecker(&a, &b);

  printf("KRONECKER: %lli\n", (long long)j);

  rz_cornacchia(&x, &y, &a, &b);
  printf("Diophantine eq x^2 + by^2 = a solutions\n");
  printf("x: ");
  rz_println(&x);
  printf("y: ");
  rz_println(&y);

  rz_clear_multi(&x, &y, &z);

  rz_gen_safe_prime(&c, 512);

  printf("Safe prime: ");
  rz_println(&c);

  rpol_test();

  printf("factorize a: ");
  rvec_t v = {0};
  rvec_init_dynamic(&v);

  rz_factorize(&v, &a);

  rvec_println(&v);
  rvec_clear(&v);

  rntt_ctx_t ctx = {0};

  // rntt_ctx_init(&ctx, 512, &n1, &omega, &psi);
  rntt_ctx_init_golden(&ctx, 8);
  // rntt_ctx_init_simple(&ctx, 8, 4294967295);
  rntt_ctx_print_debug(&ctx);

  rntt_ctx_clear(&ctx);

  if (rz_pollard_rho(&c, &a)) {
    printf("Factor found: ");
    rz_println(&c);
  } else {
    printf("No non-trivial factor found (Number might be prime).\n");
  }

  if (rz_pollard_p_minus_one(&c, &a)) {
    printf("Factor found p - 1: ");
    rz_println(&c);
  } else {
    printf("p - 1 fail \n");
  }

  rz_t n, d;
  rz_init_multi(&n, &d, NULL);
  /*
  u64 n_1 = 10;
    for (u64 i = 0; i < n_1; i++) {
      rz_set_u64(&n, i);
      rz_lucas(&c, &d, &a, &b, &n);
      printf("Lucas\n");
      rz_println(&c);
      rz_println(&d);
    }
  */
  if (rz_bpsw(&a)) {
    printf("a is Prime \n");
  } else {
    printf("a is not prime to base b \n");
  }

  rz_t r;
  rz_init(&r);

  rz_tonelli_shanks(&r, &a, &b);

  rz_print(&r);
  printf("^2 = ");
  rz_print(&a);
  printf(" mod ");
  rz_println(&b);
  char* s = rz_to_string(&b);
  printf("Result: %s\n", s);
  free(s);

  rz_t my_prime;
  rz_init(&my_prime);
  printf("Generating a 2048-bit prime... (this may take a minute)\n");

  if (rz_gen_prime(&my_prime, 2048) == RABIN_SUCCESS) {
    printf("\nSuccess! Your 2048-bit prime is:\n");
    printf(
        "----------------------------------------------------------------\n");
    rz_print(&my_prime);
    printf(
        "----------------------------------------------------------------\n");
  } else {
    fprintf(stderr, "Failed to generate prime or read from /dev/urandom\n");
  }

  //  printf("Generating 512-bit provable prime...\n");
  //// rz_provable_prime(&my_prime, 512);

  // printf("Result: ");
  // rz_println(&my_prime);

  // avxtest();

  //
  //
  rpol_test();
  //
  rpol_t p = {0};
  rpol_init(&p);

  rz_clear(&my_prime);
  rz_clear_multi(&a, &b, &c, &n, &d, &r, &omega, &psi, NULL);
  rz_clear_constants();
  return 0;
}
