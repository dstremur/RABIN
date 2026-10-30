/*
 * test_lucas.c
 *
 * Differential tests for the Lucas sequence solvers: an independent
 * straight-line recurrence over __int128 (written from the definition
 * U_k = P*U_{k-1} - Q*U_{k-2}, V_k = P*V_{k-1} - Q*V_{k-2}) is compared
 * against rz_lucas_solve() (exact) and rz_lucas_solve_mod() (Montgomery
 * modular). Also covers the discriminant identity, error-domain guards
 * and oracle-verified rz_stronglucas() vectors.
 *
 * Copyright (C) 2026 Diego Strebel
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <stdio.h>

#include "../include/rabin.h"
#include "test_util.h"

typedef __int128 i128;

/* Independent oracle: straight iteration of the defining recurrence.
 * All (P, Q, k) below are sized so |U_k|, |V_k|, |Q^k| fit in i128. */
static void lucas_naive(long P, long Q, int k, i128* U, i128* V, i128* Qn)
{
  i128 um2 = 0, um1 = 1;       /* U_0, U_1 */
  i128 vm2 = 2, vm1 = (i128)P; /* V_0, V_1 */

  if (k == 0) {
    *U = 0;
    *V = 2;
    *Qn = 1;
    return;
  }
  for (int i = 1; i < k; i++) {
    i128 u = (i128)P * um1 - (i128)Q * um2;
    i128 v = (i128)P * vm1 - (i128)Q * vm2;
    um2 = um1;
    um1 = u;
    vm2 = vm1;
    vm1 = v;
  }
  i128 qn = 1;
  for (int i = 0; i < k; i++) qn *= (i128)Q;
  *U = um1;
  *V = vm1;
  *Qn = qn;
}

static void rz_from_i128(rz_t* out, i128 v)
{
  int neg = v < 0;
  unsigned __int128 m = neg ? (unsigned __int128)(-v) : (unsigned __int128)v;
  u64 lo = (u64)m;
  u64 hi = (u64)(m >> 64);

  if (hi) {
    rz_set_u64(out, hi);
    rz_lshift(out, out, 64);
    rz_add_u64(out, out, lo);
  } else {
    rz_set_u64(out, lo);
  }
  out->is_neg = neg && ((lo | hi) != 0);
}

typedef struct {
  long P;
  long Q;
  int max_k;
} pair_t;

static const pair_t PAIRS[] = {
    {1, -1, 30}, {2, 1, 25}, {1, 2, 25},  {5, -2, 28},
    {-3, 1, 20}, {3, 1, 20}, {0, -1, 20}, {-1, -2, 22},
};
#define NUM_PAIRS (sizeof(PAIRS) / sizeof(PAIRS[0]))

TEST(exact_vs_naive)
{
  rz_t p, q, k, u, v, qn, want;
  rz_init_multi(&p, &q, &k, &u, &v, &qn, &want, NULL);

  for (unsigned i = 0; i < NUM_PAIRS; i++) {
    rz_set_i64(&p, PAIRS[i].P);
    rz_set_i64(&q, PAIRS[i].Q);

    for (int kk = 0; kk <= PAIRS[i].max_k; kk++) {
      rz_set_u64(&k, (u64)kk);
      if (rz_lucas_solve(&u, &v, &p, &q, &qn, &k) != RABIN_SUCCESS) {
        TEST_ASSERT(false, "rz_lucas_solve failed");
        continue;
      }
      i128 U, V, Qn;
      lucas_naive(PAIRS[i].P, PAIRS[i].Q, kk, &U, &V, &Qn);

      rz_from_i128(&want, U);
      if (rz_cmp(&u, &want) != 0)
        TEST_ASSERT(false, "U_k mismatch vs naive recurrence");
      rz_from_i128(&want, V);
      if (rz_cmp(&v, &want) != 0)
        TEST_ASSERT(false, "V_k mismatch vs naive recurrence");
      rz_from_i128(&want, Qn);
      if (rz_cmp(&qn, &want) != 0)
        TEST_ASSERT(false, "Q^k mismatch vs naive recurrence");
    }
  }

  rz_clear_multi(&p, &q, &k, &u, &v, &qn, &want, NULL);
}

TEST(mod_vs_exact)
{
  /* regression: rz_mont_in dropped the sign of negative values, so the
   * modular solver disagreed with (exact mod m) whenever Q < 0 */
  static const u64 mods[] = {9ULL, 2047ULL, 10403ULL, 1000003ULL,
                             1099511627775ULL}; /* 2^40 - 1 */
  rz_t p, q, k, u, v, qn, m, r, exact, eu, ev, eqn;
  rz_init_multi(&p, &q, &k, &u, &v, &qn, &m, &r, &exact, &eu, &ev, &eqn, NULL);

  for (unsigned i = 0; i < NUM_PAIRS; i++) {
    rz_set_i64(&p, PAIRS[i].P);
    rz_set_i64(&q, PAIRS[i].Q);

    /* sample a few indices per pair to keep the suite fast */
    for (int kk = 0; kk <= PAIRS[i].max_k; kk += 7) {
      rz_set_u64(&k, (u64)kk);
      for (unsigned j = 0; j < 5; j++) {
        rz_set_u64(&m, mods[j]);
        if (rz_lucas_solve(&eu, &ev, &p, &q, &eqn, &k) != RABIN_SUCCESS) {
          TEST_ASSERT(false, "exact solve failed");
          continue;
        }
        if (rz_lucas_solve_mod(&u, &v, &p, &q, &qn, &k, &m) != RABIN_SUCCESS) {
          TEST_ASSERT(false, "mod solve failed");
          continue;
        }

        /* expected: floor(exact mod m) in [0, m) */
        rz_mod(&r, &eu, &m);
        if (r.is_neg) rz_add(&r, &r, &m);
        if (rz_cmp(&u, &r) != 0) TEST_ASSERT(false, "U mod m mismatch");

        rz_mod(&r, &ev, &m);
        if (r.is_neg) rz_add(&r, &r, &m);
        if (rz_cmp(&v, &r) != 0) TEST_ASSERT(false, "V mod m mismatch");

        rz_mod(&r, &eqn, &m);
        if (r.is_neg) rz_add(&r, &r, &m);
        if (rz_cmp(&qn, &r) != 0) TEST_ASSERT(false, "Q^k mod m mismatch");
      }
    }
  }

  rz_clear_multi(&p, &q, &k, &u, &v, &qn, &m, &r, &exact, &eu, &ev, &eqn, NULL);
}

TEST(discriminant_identity)
{
  /* V_k^2 - D * U_k^2 = 4 Q^k with D = P^2 - 4Q, both D signs */
  struct {
    long P, Q, k;
  } vec[] = {{1, -1, 25}, {1, 2, 20}, {5, -2, 14}, {-3, 1, 15}};

  rz_t p, q, k, u, v, qn, t1, t2, d, lhs, rhs, want;
  rz_init_multi(&p, &q, &k, &u, &v, &qn, &t1, &t2, &d, &lhs, &rhs, &want, NULL);

  for (unsigned i = 0; i < 4; i++) {
    rz_set_i64(&p, vec[i].P);
    rz_set_i64(&q, vec[i].Q);
    rz_set_u64(&k, (u64)vec[i].k);

    if (rz_lucas_solve(&u, &v, &p, &q, &qn, &k) != RABIN_SUCCESS) {
      TEST_ASSERT(false, "solve failed in identity test");
      continue;
    }

    /* D = P^2 - 4Q */
    rz_mul(&d, &p, &p);
    rz_mul(&t2, &q, &q); /* placeholder overwritten below */
    rz_set_u64(&t1, 4);
    rz_mul(&t2, &t1, &q);
    rz_sub(&d, &d, &t2);

    /* lhs = V^2 - D*U^2 */
    rz_mul(&t1, &v, &v);
    rz_mul(&t2, &u, &u);
    rz_mul(&t2, &t2, &d);
    rz_sub(&lhs, &t1, &t2);

    /* rhs = 4 * Q^k */
    i128 U, V, Qn;
    lucas_naive(vec[i].P, vec[i].Q, vec[i].k, &U, &V, &Qn);
    rz_from_i128(&rhs, (i128)4 * Qn);

    if (rz_cmp(&lhs, &rhs) != 0) TEST_ASSERT(false, "identity violated");
    /* also cross-check exact U against the naive oracle here */
    rz_from_i128(&want, U);
    if (rz_cmp(&u, &want) != 0) TEST_ASSERT(false, "U identity cross-check");
  }

  rz_clear_multi(&p, &q, &k, &u, &v, &qn, &t1, &t2, &d, &lhs, &rhs, &want,
                 NULL);
}

TEST(domain_guards)
{
  rz_t p, q, k, m, u, v, qn;
  rz_init_multi(&p, &q, &k, &m, &u, &v, &qn, NULL);
  rz_set_u64(&p, 1);
  rz_set_u64(&q, 1);
  rz_set_u64(&m, 10007);

  rz_set_i64(&k, -5);
  TEST_ASSERT(rz_lucas_solve(&u, &v, &p, &q, &qn, &k) == RABIN_ERR_INVALID_ARG,
              "exact solver must reject negative index");
  TEST_ASSERT(
      rz_lucas_solve_mod(&u, &v, &p, &q, &qn, &k, &m) == RABIN_ERR_INVALID_ARG,
      "mod solver must reject negative index");

  rz_set_u64(&k, 5);
  rz_set_u64(&m, 1);
  TEST_ASSERT(rz_lucas_solve_mod(&u, &v, &p, &q, &qn, &k, &m) != RABIN_SUCCESS,
              "mod solver must reject m = 1");
  rz_set_u64(&m, 2);
  TEST_ASSERT(rz_lucas_solve_mod(&u, &v, &p, &q, &qn, &k, &m) != RABIN_SUCCESS,
              "mod solver must reject even m");
  rz_set_i64(&m, -10007);
  TEST_ASSERT(rz_lucas_solve_mod(&u, &v, &p, &q, &qn, &k, &m) != RABIN_SUCCESS,
              "mod solver must reject negative m");

  rz_clear_multi(&p, &q, &k, &m, &u, &v, &qn, NULL);
}

TEST(strong_lucas_vectors)
{
  /* Oracle-verified: primes pass strong Lucas under parameters with
   * (D|n) = -1; the six composites are strong MR base-2 pseudoprimes that
   * only the Lucas stage catches (independent u128 matrix-powering oracle). */
  struct {
    u64 n;
    long q;
    int want;
  } vec[] = {
      {7ULL, -1, 1},
      {13ULL, -1, 1},
      {101ULL, 2, 1},
      {10007ULL, -1, 1},
      {104729ULL, 4, 1},
      {2305843009213693951ULL, -4, 1},
      {19ULL, 2, 1},
      {2047ULL, -1, 0},
      {1373653ULL, -1, 0},
      {3215031751ULL, 3, 0},
      {2152302898747ULL, -1, 0},
      {3474749660383ULL, -1, 0},
      {3825123056546413051ULL, 2, 0},
  };

  rz_t n, P, Q;
  rz_init_multi(&n, &P, &Q, NULL);
  rz_set_u64(&P, 1);

  for (unsigned i = 0; i < sizeof(vec) / sizeof(vec[0]); i++) {
    rz_set_u64(&n, vec[i].n);
    rz_set_i64(&Q, vec[i].q);
    if (rz_stronglucas(&n, &P, &Q) != (bool)vec[i].want)
      TEST_ASSERT(false, "strong Lucas disagrees with oracle vector");
  }

  /* domain guards: 0, 1, 2, even, negative are outside the test's domain */
  rz_set_u64(&Q, 1);
  rz_set_u64(&n, 0);
  TEST_ASSERT(rz_stronglucas(&n, &P, &Q) == false, "stronglucas(0)");
  rz_set_u64(&n, 1);
  TEST_ASSERT(rz_stronglucas(&n, &P, &Q) == false, "stronglucas(1)");
  rz_set_u64(&n, 2);
  TEST_ASSERT(rz_stronglucas(&n, &P, &Q) == false, "stronglucas(2)");
  rz_set_u64(&n, 4);
  TEST_ASSERT(rz_stronglucas(&n, &P, &Q) == false, "stronglucas(4)");
  rz_set_i64(&n, -9);
  TEST_ASSERT(rz_stronglucas(&n, &P, &Q) == false, "stronglucas(-9)");

  rz_clear_multi(&n, &P, &Q, NULL);
}

int main(void)
{
  RUN_TEST(exact_vs_naive);
  RUN_TEST(mod_vs_exact);
  RUN_TEST(discriminant_identity);
  RUN_TEST(domain_guards);
  RUN_TEST(strong_lucas_vectors);
  return TEST_REPORT();
}
