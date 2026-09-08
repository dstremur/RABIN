/*
 * test_div_gcd.c
 *
 * Cross-fuzz division and GCD against GMP.
 *
 * Covers:
 *   - bn_divmod / bn_div / bn_mod (truncated semantics, all sign combos)
 *     with random sizes, the qn == 1 (equal-size) fast path, single-limb
 *     divisors, power-of-two divisors, exact-quotient +/- 1 correction
 *     cases, and the kernel boundary sizes (Karatsuba / NTT / Newton
 *     dispatch thresholds in limbs).
 *   - bn_gcd (Euclid chain with 2-adic extraction) against mpz_gcd,
 *     including sign handling (result must be non-negative) and the
 *     Fibonacci worst case.
 *
 * Copyright (C) 2026 Diego Strebel
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "../include/bignum.h"

static int failures = 0;

// export a bignum to an mpz (exact, no decimal round-trip)
static void bn_to_mpz(mpz_t out, const bignum* in)
{
  u64 n = 0;
  if (in->size > 0) n = limbs_norm(in->limbs, in->size);
  if (n == 0) {
    mpz_set_ui(out, 0);
    return;
  }
  mpz_import(out, (size_t)n, -1, 8, 0, 0, in->limbs);
  if (in->is_neg) mpz_neg(out, out);
}

static void report_mismatch(const char* op, const char* detail,
                            const bignum* res, const mpz_t gres)
{
  char* rs = bn_to_string(res);
  char* gs = mpz_get_str(NULL, 10, gres);
  fprintf(stderr, "[FAIL] %s (%s)\n  ours: %s\n  gmp : %s\n", op, detail, rs,
          gs);
  free(rs);
  free(gs);
  failures++;
}

static int check_bignum_mpz(const char* op, const char* detail,
                            const bignum* res, const mpz_t gres)
{
  mpz_t t;
  mpz_init(t);
  bn_to_mpz(t, res);
  int ok = mpz_cmp(t, gres) == 0;
  if (!ok) report_mismatch(op, detail, res, gres);
  mpz_clear(t);
  return ok;
}

// a, b given as mpz; build bignums, run op, verify against gmp_q/gmp_r
static void check_divmod_case(mpz_t a, mpz_t b, int tag)
{
  char* sa = mpz_get_str(NULL, 10, a);
  char* sb = mpz_get_str(NULL, 10, b);

  bignum x, y, q, r;
  mpz_t gq, gr;
  mpz_inits(gq, gr, NULL);

  bn_init_multi(&x, &y, &q, &r, NULL);
  bn_init_val(&x, sa);
  bn_init_val(&y, sb);

  bn_divmod(&q, &r, &x, &y);
  mpz_tdiv_qr(gq, gr, a, b);

  char detail[64];
  snprintf(detail, sizeof(detail), "case %d: a=%ld bits, b=%ld bits", tag,
           (long)mpz_sizeinbase(a, 2), (long)mpz_sizeinbase(b, 2));

  int ok = check_bignum_mpz("q", detail, &q, gq) &&
           check_bignum_mpz("r", detail, &r, gr);

  // invariant: q*b + r == a
  if (ok) {
    bignum t, u;
    mpz_t g;
    mpz_init(g);
    bn_init_multi(&t, &u, NULL);
    bn_mul(&t, &q, &y);
    bn_add(&u, &t, &r);
    mpz_mul(g, gq, b);
    mpz_add(g, g, gr);
    if (!check_bignum_mpz("q*b+r", detail, &u, g)) ok = 0;
    mpz_clear(g);
    bn_free(&t);
    bn_free(&u);
  }

  bn_free_multi(&x, &y, &q, &r, NULL);
  mpz_clears(gq, gr, NULL);
  free(sa);
  free(sb);
}

static void check_gcd_case(mpz_t a, mpz_t b, int tag)
{
  char* sa = mpz_get_str(NULL, 10, a);
  char* sb = mpz_get_str(NULL, 10, b);

  bignum x, y, d;
  mpz_t g;
  mpz_init(g);

  bn_init_multi(&x, &y, &d, NULL);
  bn_init_val(&x, sa);
  bn_init_val(&y, sb);

  bn_gcd(&d, &x, &y);
  mpz_gcd(g, a, b);

  // gcd must be non-negative
  if (d.is_neg && !bn_is_zero(&d)) {
    fprintf(stderr, "[FAIL] gcd (%d): negative result\n", tag);
    failures++;
  }

  char detail[64];
  snprintf(detail, sizeof(detail), "case %d: a=%ld bits, b=%ld bits", tag,
           (long)mpz_sizeinbase(a, 2), (long)mpz_sizeinbase(b, 2));
  check_bignum_mpz("gcd", detail, &d, g);

  bn_free_multi(&x, &y, &d, NULL);
  mpz_clear(g);
  free(sa);
  free(sb);
}

int main()
{
  gmp_randstate_t rs;
  gmp_randinit_default(rs);
  gmp_randseed_ui(rs, (unsigned long)time(NULL));
  srand(time(NULL));

  int tag = 0;

  /*-----------------------------------------------------------------------
   * 1. random div/mod, mixed sizes and signs (fast, small-to-medium)
   *---------------------------------------------------------------------*/
  for (int i = 0; i < 20000; i++) {
    mpz_t a, b;
    mpz_inits(a, b, NULL);
    u64 ba = 1 + (u64)(rand() % 2048);
    u64 bb = 1 + (u64)(rand() % 1024);
    if (bb > ba && (rand() & 1)) bb = 1 + (u64)(rand() % ba);
    mpz_urandomb(a, rs, (size_t)ba);
    mpz_urandomb(b, rs, (size_t)bb);
    if (mpz_cmp_ui(b, 0) == 0) mpz_set_ui(b, 1);
    if (rand() & 1) mpz_neg(a, a);
    if (rand() & 1) mpz_neg(b, b);
    check_divmod_case(a, b, tag++);
    mpz_clears(a, b, NULL);
  }
  printf("PASS 1/4: 20000 random div/mod cases\n");

  /*-----------------------------------------------------------------------
   * 2. structured division cases: kernel boundaries and corrections
   *---------------------------------------------------------------------*/
  {
    // limb-boundary divisor sizes (Karatsuba 128, NTT min 256, Newton 512)
    int m_limbs[] = {1,   2,   63,  64,  65,  127, 128, 129,
                     255, 256, 257, 511, 512, 513, 514};
    // quotient sizes: qn = 1 (fast path) and the Newton/other boundaries
    int qn_list[] = {1, 15, 16, 17, 509, 511, 512, 513, 514};

    for (int mi = 0; mi < (int)(sizeof(m_limbs) / sizeof(m_limbs[0])); mi++) {
      for (int qi = 0; qi < (int)(sizeof(qn_list) / sizeof(qn_list[0])); qi++) {
        int m = m_limbs[mi];
        int qn = qn_list[qi];
        int an = m + qn - 1;
        // keep the total size bounded (a few 10k-bit cases are enough)
        if (an > 1500) continue;

        mpz_t a, b;
        mpz_inits(a, b, NULL);
        mpz_urandomb(a, rs, (size_t)(64L * an));
        mpz_setbit(a, (mp_bitcnt_t)(64L * an - 1));
        mpz_urandomb(b, rs, (size_t)(64L * m));
        mpz_setbit(b, (mp_bitcnt_t)(64L * m - 1));
        if (mpz_cmp(a, b) < 0) mpz_swap(a, b);
        check_divmod_case(a, b, tag++);
        mpz_clears(a, b, NULL);
      }
    }

    // a = b * 2^t +/- 1: forces quotient +/- 1 corrections
    for (int t = 0; t <= 65; t++) {
      for (int bits = 32; bits <= 4096; bits *= 4) {
        mpz_t a, b;
        mpz_inits(a, b, NULL);
        mpz_urandomb(b, rs, (size_t)bits);
        mpz_setbit(b, (mp_bitcnt_t)(bits - 1));
        mpz_mul_ui(a, b, 1);
        mpz_mul_2exp(a, a, (size_t)t);
        int plus = rand() & 1;
        if (plus)
          mpz_add_ui(a, a, 1);
        else
          mpz_sub_ui(a, a, 1);
        check_divmod_case(a, b, tag++);
        mpz_clears(a, b, NULL);
      }
    }

    // powers of two as divisor, all sizes
    for (int k = 0; k <= 4096; k += k ? k / 64 + 1 : 1) {
      mpz_t a, b;
      mpz_inits(a, b, NULL);
      mpz_urandomb(a, rs, (size_t)(k * 2 + 128));
      mpz_setbit(b, (mp_bitcnt_t)k);
      check_divmod_case(a, b, tag++);
      mpz_clears(a, b, NULL);
    }

    // degenerate pairs
    {
      const char* pairs[][2] = {
          {"1", "1"},
          {"1", "2"},
          {"2", "1"},
          {"0", "7"},
          {"7", "0"},
          {"0", "0"},
          {"18446744073709551615", "1"},
          {"18446744073709551615", "18446744073709551615"}, /* 2^64-1 */
          {"18446744073709551616",
           "18446744073709551615"}, /* 2^64 / (2^64-1) */
          {"-5", "3"},
          {"5", "-3"},
          {"-5", "-3"},
      };
      for (int pi = 0; pi < (int)(sizeof(pairs) / sizeof(pairs[0])); pi++) {
        mpz_t a, b;
        mpz_inits(a, b, NULL);
        mpz_set_str(a, pairs[pi][0], 10);
        mpz_set_str(b, pairs[pi][1], 10);
        if (mpz_sgn(b) != 0) {
          check_divmod_case(a, b, tag++);
        }
        // gcd is defined for 0 too
        check_gcd_case(a, b, tag++);
        mpz_clears(a, b, NULL);
      }
    }
  }
  printf("PASS 2/4: structured div/mod cases (boundaries, powers of 2)\n");

  /*-----------------------------------------------------------------------
   * 3. gcd fuzz: random + adversarial
   *---------------------------------------------------------------------*/
  for (int i = 0; i < 2000; i++) {
    mpz_t a, b;
    mpz_inits(a, b, NULL);
    u64 ba = 1 + (u64)(rand() % 4096);
    u64 bb = 1 + (u64)(rand() % 4096);
    mpz_urandomb(a, rs, (size_t)ba);
    mpz_urandomb(b, rs, (size_t)bb);
    if (rand() & 1) mpz_neg(a, a);
    if (rand() & 1) mpz_neg(b, b);
    check_gcd_case(a, b, tag++);
    mpz_clears(a, b, NULL);
  }

  // Fibonacci worst case for the Euclid chain (~3000 bits)
  {
    mpz_t fa, fb;
    mpz_inits(fa, fb, NULL);
    mpz_set_ui(fa, 1);
    mpz_set_ui(fb, 1);
    while (mpz_sizeinbase(fb, 2) < 3000) {
      mpz_t t;
      mpz_init(t);
      mpz_add(t, fa, fb);
      mpz_swap(fa, fb);
      mpz_set(fb, t);
      mpz_clear(t);
    }
    check_gcd_case(fa, fb, tag++);
    mpz_clears(fa, fb, NULL);
  }

  // gcd with powers of two (2-adic extraction)
  for (int k = 1; k <= 4096; k *= 2) {
    mpz_t a, b;
    mpz_inits(a, b, NULL);
    mpz_urandomb(a, rs, (size_t)(k + 128));
    mpz_setbit(b, (mp_bitcnt_t)k);
    // make a share the power-of-two factor (half the time keep it odd)
    if (rand() & 1) mpz_mul_2exp(a, a, (size_t)k);
    check_gcd_case(a, b, tag++);
    mpz_clears(a, b, NULL);
  }
  printf("PASS 3/4: gcd fuzz + adversarial\n");

  /*-----------------------------------------------------------------------
   * 4. large sizes (Newton division + big gcd chains)
   *---------------------------------------------------------------------*/
  {
    struct {
      int abits, bbits;
    } big[] = {
        {65536, 32768},  {65539, 1024},  {131072, 65536},
        {131072, 2048},  {65539, 65539},  // equal size -> qn == 1 path at scale
        {100000, 50000},
    };
    for (int i = 0; i < (int)(sizeof(big) / sizeof(big[0])); i++) {
      mpz_t a, b;
      mpz_inits(a, b, NULL);
      mpz_urandomb(a, rs, (size_t)big[i].abits);
      mpz_setbit(a, (mp_bitcnt_t)(big[i].abits - 1));
      mpz_urandomb(b, rs, (size_t)big[i].bbits);
      mpz_setbit(b, (mp_bitcnt_t)(big[i].bbits - 1));
      if (mpz_cmp(a, b) < 0) mpz_swap(a, b);
      check_divmod_case(a, b, tag++);
      check_gcd_case(a, b, tag++);
      mpz_clears(a, b, NULL);
    }
  }
  printf("PASS 4/4: large div/mod + gcd\n");

  gmp_randclear(rs);

  if (failures == 0) {
    printf("ALL TESTS PASSED (%d cases)\n", tag);
    return 0;
  }
  fprintf(stderr, "%d FAILURES\n", failures);
  return 1;
}
