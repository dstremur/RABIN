/*
 * test_all.c
 *
 * Property-based smoke test over the core rz_t API: arithmetic, bit logic,
 * powers, primes, factoring, and vector operations, all gated on
 * self-consistency assertions and known answers.
 *
 * Copyright (C) 2026 Diego Strebel
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <stdio.h>

#include "../include/rabin.h"
#include "../include/rvec.h"
#include "../include/rzlogic.h"
#include "test_util.h"

static const char* A_STR =
    "214240923082520938209829358290858352903593285435345345346346346456456456"
    "456456464256324636346346346";
static const char* B_STR = "2344";

TEST(arith_consistency)
{
  rz_t a, b, sum, diff, quo, rem, tmp, tmp2;
  rz_init_multi(&a, &b, &sum, &diff, &quo, &rem, &tmp, &tmp2, NULL);
  TEST_ASSERT(rz_init_val(&a, A_STR) == RABIN_SUCCESS, "init_val a");
  TEST_ASSERT(rz_init_val(&b, B_STR) == RABIN_SUCCESS, "init_val b");

  TEST_ASSERT(rz_add(&sum, &a, &b) == RABIN_SUCCESS, "add");
  TEST_ASSERT(rz_add(&tmp, &b, &a) == RABIN_SUCCESS, "add commuted");
  TEST_ASSERT(rz_cmp(&sum, &tmp) == 0, "addition is commutative");
  TEST_ASSERT(rz_sub(&tmp, &sum, &b) == RABIN_SUCCESS, "sub");
  TEST_ASSERT(rz_cmp(&tmp, &a) == 0, "sub undoes add");

  TEST_ASSERT(rz_mul(&tmp, &a, &b) == RABIN_SUCCESS, "mul");
  TEST_ASSERT(rz_mul(&tmp2, &b, &a) == RABIN_SUCCESS, "mul commuted");
  TEST_ASSERT(rz_cmp(&tmp, &tmp2) == 0, "multiplication is commutative");

  TEST_ASSERT(rz_div(&quo, &a, &b) == RABIN_SUCCESS, "div");
  TEST_ASSERT(rz_mod(&rem, &a, &b) == RABIN_SUCCESS, "mod");
  TEST_ASSERT(rz_mul(&tmp, &quo, &b) == RABIN_SUCCESS, "div * b");
  TEST_ASSERT(rz_add(&tmp, &tmp, &rem) == RABIN_SUCCESS, "div*b + rem");
  TEST_ASSERT(rz_cmp(&tmp, &a) == 0, "a == (a/b)*b + a%b");
  TEST_ASSERT(rz_cmp(&rem, &b) < 0 && !rem.is_neg, "0 <= a%b < |b|");

  rz_clear_multi(&a, &b, &sum, &diff, &quo, &rem, &tmp, &tmp2, NULL);
}

TEST(bitlogic_identities)
{
  rz_t a, b, and_, or_, xor_, tmp;
  rz_init_multi(&a, &b, &and_, &or_, &xor_, &tmp, NULL);
  TEST_ASSERT(rz_init_val(&a, A_STR) == RABIN_SUCCESS, "init_val a");
  TEST_ASSERT(rz_init_val(&b, B_STR) == RABIN_SUCCESS, "init_val b");

  TEST_ASSERT(rz_and(&and_, &a, &b) == RABIN_SUCCESS, "and");
  TEST_ASSERT(rz_or(&or_, &a, &b) == RABIN_SUCCESS, "or");
  TEST_ASSERT(rz_xor(&xor_, &a, &b) == RABIN_SUCCESS, "xor");

  /* For non-negative operands: (a|b) == (a&b) + (a^b). */
  TEST_ASSERT(rz_add(&tmp, &and_, &xor_) == RABIN_SUCCESS, "and + xor");
  TEST_ASSERT(rz_cmp(&tmp, &or_) == 0, "(a|b) == (a&b) + (a^b)");

  TEST_ASSERT(rz_xor(&tmp, &a, &a) == RABIN_SUCCESS, "a xor a");
  TEST_ASSERT(rz_is_zero(&tmp), "a ^ a == 0");
  TEST_ASSERT(rz_and(&tmp, &a, &a) == RABIN_SUCCESS, "a and a");
  TEST_ASSERT(rz_cmp(&tmp, &a) == 0, "a & a == a");

  rz_clear_multi(&a, &b, &and_, &or_, &xor_, &tmp, NULL);
}

TEST(pow_isqrt_and_modexp)
{
  rz_t base, e, r, a, m, one;
  rz_init_multi(&base, &e, &r, &a, &m, &one, NULL);
  TEST_ASSERT(rz_init_val(&a, A_STR) == RABIN_SUCCESS, "init_val a");
  rz_set_u64(&base, 2);
  rz_set_u64(&e, 10);
  TEST_ASSERT(rz_set_u64(&one, 1) == RABIN_SUCCESS, "set one");

  TEST_ASSERT(rz_pow(&r, &base, &e) == RABIN_SUCCESS, "2^10");
  rz_set_u64(&m, 1024);
  TEST_ASSERT(rz_cmp(&r, &m) == 0, "2^10 == 1024");

  rz_set_u64(&e, 0);
  TEST_ASSERT(rz_pow(&r, &base, &e) == RABIN_SUCCESS, "2^0");
  TEST_ASSERT(rz_cmp(&r, &one) == 0, "2^0 == 1");

  TEST_ASSERT(rz_isqrt(&r, &a) == RABIN_SUCCESS, "isqrt(a)");
  /* r^2 <= a < (r+1)^2 */
  rz_t sq, sqp1, r1;
  rz_init_multi(&sq, &sqp1, &r1, NULL);
  TEST_ASSERT(rz_mul(&sq, &r, &r) == RABIN_SUCCESS, "r*r");
  TEST_ASSERT(rz_cmp(&sq, &a) <= 0, "isqrt^2 <= a");
  TEST_ASSERT(rz_add(&r1, &r, &one) == RABIN_SUCCESS, "r+1");
  TEST_ASSERT(rz_mul(&sqp1, &r1, &r1) == RABIN_SUCCESS, "(r+1)^2");
  TEST_ASSERT(rz_cmp(&sqp1, &a) > 0, "(isqrt+1)^2 > a");
  rz_clear_multi(&sq, &sqp1, &r1, NULL);

  /* a^1 mod m == a mod m, a^0 mod m == 1 (m > 1). */
  rz_set_u64(&m, 1000000007ULL);
  rz_set_u64(&e, 1);
  TEST_ASSERT(rz_mod_exp(&r, &a, &e, &m) == RABIN_SUCCESS, "a^1 mod m");
  TEST_ASSERT(rz_mod(&sq, &a, &m) == RABIN_SUCCESS, "a mod m");
  TEST_ASSERT(rz_cmp(&r, &sq) == 0, "a^1 mod m == a mod m");
  rz_set_u64(&e, 0);
  TEST_ASSERT(rz_mod_exp(&r, &a, &e, &m) == RABIN_SUCCESS, "a^0 mod m");
  TEST_ASSERT(rz_cmp(&r, &one) == 0, "a^0 mod m == 1");

  rz_clear_multi(&base, &e, &r, &a, &m, &one, NULL);
}

TEST(prime_and_provable)
{
  rz_t p, q;
  rz_init_multi(&p, &q, NULL);

  TEST_ASSERT(rz_gen_prime(&p, 1024) == RABIN_SUCCESS, "rz_gen_prime 1024");
  TEST_ASSERT(rz_bpsw(&p), "generated 1024-bit prime passes BPSW");

  TEST_ASSERT(rz_provable_prime(&q, 256) == RABIN_SUCCESS,
              "rz_provable_prime 256");
  TEST_ASSERT(rz_bpsw(&q), "provable 256-bit prime passes BPSW");

  rz_clear_multi(&p, &q, NULL);
}

TEST(factoring_jacobi_tonelli)
{
  rz_t n, f, r, p, four, mod;
  rz_init_multi(&n, &f, &r, &p, &four, &mod, NULL);

  /* 1000000007 * 1000000009 = 1000000016000000063 */
  TEST_ASSERT(rz_init_val(&n, "1000000016000000063") == RABIN_SUCCESS,
              "init_val n");
  bool found = rz_pollard_rho(&f, &n);
  TEST_ASSERT(found, "pollard_rho finds a factor of 1e18 semiprime");
  if (found) {
    TEST_ASSERT(rz_mod(&mod, &n, &f) == RABIN_SUCCESS, "n mod f");
    TEST_ASSERT(rz_is_zero(&mod), "factor divides n");
    TEST_ASSERT(rz_cmp(&f, &RZ_ONE) > 0, "factor > 1");
    TEST_ASSERT(rz_cmp(&f, &n) < 0, "factor < n");
  }

  i64 j = rz_jacobi(&f, &n);
  TEST_ASSERT(j >= -1 && j <= 1, "jacobi in {-1,0,1}");
  TEST_ASSERT(rz_jacobi(&n, &n) == 0, "jacobi(n,n) == 0");

  TEST_ASSERT(rz_gen_prime(&p, 256) == RABIN_SUCCESS, "gen 256-bit prime");
  rz_set_u64(&four, 4);
  TEST_ASSERT(rz_tonelli_shanks(&r, &four, &p) == RABIN_SUCCESS,
              "tonelli_shanks(4, p) succeeds");
  TEST_ASSERT(rz_mod(&mod, &r, &p) == RABIN_SUCCESS, "reduce r");
  TEST_ASSERT(rz_mul(&f, &mod, &mod) == RABIN_SUCCESS, "r^2");
  TEST_ASSERT(rz_mod(&mod, &f, &p) == RABIN_SUCCESS, "r^2 mod p");
  TEST_ASSERT(rz_cmp(&mod, &four) == 0, "tonelli root squares to n mod p");

  rz_clear_multi(&n, &f, &r, &p, &four, &mod, NULL);
}

TEST(vector_operations)
{
  rz_t a, b, dot, norm, expected, term;
  rz_init_multi(&a, &b, &dot, &norm, &expected, &term, NULL);
  TEST_ASSERT(rz_init_val(&a, A_STR) == RABIN_SUCCESS, "init_val a");
  TEST_ASSERT(rz_init_val(&b, B_STR) == RABIN_SUCCESS, "init_val b");

  rvec_t v = {0};
  TEST_ASSERT(rvec_init(&v, 5) == RABIN_SUCCESS, "rvec_init");
  TEST_ASSERT(rvec_set(&v, &a, 0) == RABIN_SUCCESS, "set 0");
  TEST_ASSERT(rvec_set(&v, &b, 1) == RABIN_SUCCESS, "set 1");
  TEST_ASSERT(rvec_set(&v, &a, 2) == RABIN_SUCCESS, "set 2");
  TEST_ASSERT(rvec_set(&v, &b, 3) == RABIN_SUCCESS, "set 3");
  TEST_ASSERT(rvec_set(&v, &b, 4) == RABIN_SUCCESS, "set 4");

  /* expected = 2*a^2 + 3*b^2 */
  TEST_ASSERT(rz_mul(&expected, &a, &a) == RABIN_SUCCESS, "a^2");
  TEST_ASSERT(rz_mul(&expected, &expected, &RZ_TWO) == RABIN_SUCCESS, "2a^2");
  TEST_ASSERT(rz_mul(&term, &b, &b) == RABIN_SUCCESS, "b^2");
  rz_t three;
  rz_init(&three);
  rz_set_u64(&three, 3);
  TEST_ASSERT(rz_mul(&term, &term, &three) == RABIN_SUCCESS, "3b^2");
  TEST_ASSERT(rz_add(&expected, &expected, &term) == RABIN_SUCCESS, "sum");

  TEST_ASSERT(rvec_dot(&dot, &v, &v) == RABIN_SUCCESS, "rvec_dot");
  TEST_ASSERT(rz_cmp(&dot, &expected) == 0, "dot(v,v) == 2a^2 + 3b^2");
  TEST_ASSERT(rvec_norm(&norm, &v) == RABIN_SUCCESS, "rvec_norm");
  TEST_ASSERT(rz_isqrt(&term, &expected) == RABIN_SUCCESS, "isqrt(expected)");
  TEST_ASSERT(rz_cmp(&norm, &term) == 0, "norm(v) == isqrt(dot(v,v))");

  rvec_clear(&v);
  rz_clear_multi(&a, &b, &dot, &norm, &expected, &term, &three, NULL);
}

int main(void)
{
  rz_init_constants();

  RUN_TEST(arith_consistency);
  RUN_TEST(bitlogic_identities);
  RUN_TEST(pow_isqrt_and_modexp);
  RUN_TEST(prime_and_provable);
  RUN_TEST(factoring_jacobi_tonelli);
  RUN_TEST(vector_operations);

  rz_clear_constants();
  return TEST_REPORT();
}
