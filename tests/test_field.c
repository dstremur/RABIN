// tests/test_field.c
//
// Correctness tests for the Z_m field and Z_m[x]/(q) ring layers in
// bigfield.c. Field operations are cross-checked against GMP; the
// polynomial division and ring operations are checked against algebraic
// identities.

#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "../include/bigfield.h"
#include "../include/bignum.h"

static int failures = 0;
static int passes = 0;

// secp256k1 prime (odd -> exercises the Montgomery path)
static const char* P256 =
    "11579208923731619542357098500868790785326998466564056403945758400790883467"
    "1663";

// ---- conversion helpers ---------------------------------------------------

static void mpz_to_bn(bignum* out, const mpz_t m)
{
  char* s = mpz_get_str(NULL, 10, m);
  bn_init(out);
  bn_init_val(out, s);
  free(s);
}

static void bn_to_mpz(mpz_t out, const bignum* bn)
{
  char* s = bn_to_string(bn);
  mpz_set_str(out, s, 10);
  free(s);
}

static void check_eq(const char* name, const bignum* got, const mpz_t want)
{
  mpz_t got_mpz;
  mpz_init(got_mpz);
  bn_to_mpz(got_mpz, got);
  if (mpz_cmp(got_mpz, want) != 0) {
    char* gs = mpz_get_str(NULL, 10, got_mpz);
    char* ws = mpz_get_str(NULL, 10, want);
    fprintf(stderr, "[FAIL] %s: got %s, want %s\n", name, gs, ws);
    free(gs);
    free(ws);
    failures++;
  } else {
    passes++;
  }
  mpz_clear(got_mpz);
}

// Fill p with a random ring element of degree < 2 (coeffs in [0, m)).
static void fill_ring_elem(bigpoly* p, field_ctx* ctx, gmp_randstate_t state,
                           mpz_t M)
{
  bigpoly_alloc(p, 1);
  for (u64 i = 0; i <= 1; i++) {
    mpz_t c;
    mpz_init(c);
    mpz_urandomb(c, state, 8);
    mpz_mod(c, c, M);
    bignum coeff;
    mpz_to_bn(&coeff, c);
    field_in(&p->coeff[i], &coeff, ctx);
    bn_free(&coeff);
    mpz_clear(c);
  }
  p->deg = 1;
  bigpoly_trim(p);
}

// ---- GF(p) field ops vs GMP (Montgomery path) ----------------------------

static void test_field_mont(void)
{
  printf("--- GF(p) field ops (Montgomery path, p = secp256k1) ---\n");

  bignum m;
  bn_init(&m);
  bn_init_val(&m, P256);

  field_ctx ctx;
  field_ctx_init(&ctx, &m);
  if (!ctx.mont) {
    fprintf(stderr, "[FAIL] expected Montgomery path for odd prime\n");
    failures++;
  }

  gmp_randstate_t state;
  gmp_randinit_default(state);
  gmp_randseed_ui(state, 0xC0FFEE);

  mpz_t A, B, E, M, res;
  mpz_inits(A, B, E, M, res, NULL);
  bn_to_mpz(M, &m);

  for (int iter = 0; iter < 50; iter++) {
    mpz_urandomb(A, state, 256);
    mpz_urandomb(B, state, 256);
    mpz_urandomb(E, state, 256);
    mpz_mod(A, A, M);
    mpz_mod(B, B, M);

    bignum a, b, e;
    mpz_to_bn(&a, A);
    mpz_to_bn(&b, B);
    mpz_to_bn(&e, E);

    bignum abar, bbar;
    bn_init_multi(&abar, &bbar, NULL);
    field_in(&abar, &a, &ctx);
    field_in(&bbar, &b, &ctx);

    // add
    bignum r, rp;
    bn_init_multi(&r, &rp, NULL);
    field_add(&r, &abar, &bbar, &ctx);
    field_out(&rp, &r, &ctx);
    mpz_add(res, A, B);
    mpz_mod(res, res, M);
    check_eq("field_add", &rp, res);

    // sub
    field_sub(&r, &abar, &bbar, &ctx);
    field_out(&rp, &r, &ctx);
    mpz_sub(res, A, B);
    mpz_mod(res, res, M);
    check_eq("field_sub", &rp, res);

    // mul
    field_mul(&r, &abar, &bbar, &ctx);
    field_out(&rp, &r, &ctx);
    mpz_mul(res, A, B);
    mpz_mod(res, res, M);
    check_eq("field_mul", &rp, res);

    // pow
    field_pow(&r, &abar, &e, &ctx);
    field_out(&rp, &r, &ctx);
    mpz_powm(res, A, E, M);
    check_eq("field_pow", &rp, res);

    // inv (A nonzero)
    if (mpz_cmp_ui(A, 0) != 0) {
      bool ok = field_inv(&r, &abar, &ctx);
      if (!ok) {
        fprintf(stderr, "[FAIL] field_inv failed for nonzero A\n");
        failures++;
      } else {
        field_out(&rp, &r, &ctx);
        mpz_invert(res, A, M);
        check_eq("field_inv", &rp, res);
      }
    }

    // div (B nonzero)
    if (mpz_cmp_ui(B, 0) != 0) {
      bool ok = field_div(&r, &abar, &bbar, &ctx);
      if (!ok) {
        fprintf(stderr, "[FAIL] field_div failed for nonzero B\n");
        failures++;
      } else {
        field_out(&rp, &r, &ctx);
        mpz_invert(res, B, M);
        mpz_mul(res, A, res);
        mpz_mod(res, res, M);
        check_eq("field_div", &rp, res);
      }
    }

    bn_free_multi(&r, &rp, NULL);
    bn_free_multi(&abar, &bbar, NULL);
    bn_free_multi(&a, &b, &e, NULL);
  }

  mpz_clears(A, B, E, M, res, NULL);
  gmp_randclear(state);
  field_ctx_free(&ctx);
  bn_free(&m);
}

// ---- Z_m ops vs GMP (plain path, even modulus) ---------------------------

static void test_field_even(void)
{
  printf("--- Z_m field ops (plain path, m = 2^128) ---\n");

  bignum m;
  bn_init(&m);
  bn_set_u64(&m, 1);
  bn_lshift(&m, &m, 128);  // m = 2^128

  field_ctx ctx;
  field_ctx_init(&ctx, &m);
  if (ctx.mont) {
    fprintf(stderr, "[FAIL] expected plain path for even modulus\n");
    failures++;
  }

  gmp_randstate_t state;
  gmp_randinit_default(state);
  gmp_randseed_ui(state, 0xDEAD);

  mpz_t A, B, E, M, res;
  mpz_inits(A, B, E, M, res, NULL);
  bn_to_mpz(M, &m);

  for (int iter = 0; iter < 50; iter++) {
    mpz_urandomb(A, state, 128);
    mpz_urandomb(B, state, 128);
    mpz_urandomb(E, state, 64);
    mpz_mod(A, A, M);
    mpz_mod(B, B, M);

    bignum a, b, e;
    mpz_to_bn(&a, A);
    mpz_to_bn(&b, B);
    mpz_to_bn(&e, E);

    bignum abar, bbar;
    bn_init_multi(&abar, &bbar, NULL);
    field_in(&abar, &a, &ctx);
    field_in(&bbar, &b, &ctx);

    bignum r, rp;
    bn_init_multi(&r, &rp, NULL);

    // add
    field_add(&r, &abar, &bbar, &ctx);
    field_out(&rp, &r, &ctx);
    mpz_add(res, A, B);
    mpz_mod(res, res, M);
    check_eq("even_add", &rp, res);

    // sub
    field_sub(&r, &abar, &bbar, &ctx);
    field_out(&rp, &r, &ctx);
    mpz_sub(res, A, B);
    mpz_mod(res, res, M);
    check_eq("even_sub", &rp, res);

    // mul
    field_mul(&r, &abar, &bbar, &ctx);
    field_out(&rp, &r, &ctx);
    mpz_mul(res, A, B);
    mpz_mod(res, res, M);
    check_eq("even_mul", &rp, res);

    // pow
    field_pow(&r, &abar, &e, &ctx);
    field_out(&rp, &r, &ctx);
    mpz_powm(res, A, E, M);
    check_eq("even_pow", &rp, res);

    // inv of a unit: force A odd (the units of Z_{2^128} are the odds)
    bignum a2, a2bar;
    bn_init_multi(&a2, &a2bar, NULL);
    mpz_t Aunit;
    mpz_init(Aunit);
    mpz_set(Aunit, A);
    if (mpz_even_p(Aunit)) mpz_add_ui(Aunit, Aunit, 1);
    mpz_mod(Aunit, Aunit, M);
    mpz_to_bn(&a2, Aunit);
    field_in(&a2bar, &a2, &ctx);
    bool ok = field_inv(&r, &a2bar, &ctx);
    if (!ok) {
      fprintf(stderr, "[FAIL] even field_inv failed for a unit\n");
      failures++;
    } else {
      field_out(&rp, &r, &ctx);
      mpz_invert(res, Aunit, M);
      check_eq("even_inv (unit)", &rp, res);
    }

    // inv of a non-unit must fail: force A even (gcd(A, 2^128) > 1)
    mpz_t Anon;
    mpz_init(Anon);
    mpz_set(Anon, A);
    if (mpz_odd_p(Anon)) mpz_add_ui(Anon, Anon, 1);  // make it even
    if (mpz_cmp_ui(Anon, 0) == 0) mpz_add_ui(Anon, Anon, 2);
    mpz_mod(Anon, Anon, M);
    if (mpz_cmp_ui(Anon, 0) != 0) {
      bignum an, anbar;
      bn_init_multi(&an, &anbar, NULL);
      mpz_to_bn(&an, Anon);
      field_in(&anbar, &an, &ctx);
      bignum dummy;
      bn_init(&dummy);
      bool ok2 = field_inv(&dummy, &anbar, &ctx);
      if (ok2) {
        fprintf(stderr, "[FAIL] even field_inv succeeded for a non-unit\n");
        failures++;
      } else {
        passes++;
      }
      bn_free(&dummy);
      bn_free_multi(&an, &anbar, NULL);
    }

    mpz_clear(Aunit);
    mpz_clear(Anon);
    bn_free_multi(&a2, &a2bar, NULL);

    bn_free_multi(&r, &rp, NULL);
    bn_free_multi(&abar, &bbar, NULL);
    bn_free_multi(&a, &b, &e, NULL);
  }

  mpz_clears(A, B, E, M, res, NULL);
  gmp_randclear(state);
  field_ctx_free(&ctx);
  bn_free(&m);
}

// ---- polynomial division --------------------------------------------------

static void test_poly_divmod(void)
{
  printf("--- Polynomial division in GF(p)[x] ---\n");

  bignum m;
  bn_init(&m);
  bn_init_val(&m, P256);
  field_ctx ctx;
  field_ctx_init(&ctx, &m);

  gmp_randstate_t state;
  gmp_randinit_default(state);
  gmp_randseed_ui(state, 0xBEEF);

  mpz_t M, c;
  mpz_inits(M, c, NULL);
  bn_to_mpz(M, &m);

  for (int iter = 0; iter < 50; iter++) {
    u64 da = 8 + (u64)(iter % 5);
    u64 db = 2 + (u64)(iter % 3);

    bigpoly a, b, q, r;
    bigpoly_init(&a);
    bigpoly_init(&b);
    bigpoly_init(&q);
    bigpoly_init(&r);

    // build a (random coeffs)
    bigpoly_alloc(&a, da);
    for (u64 i = 0; i <= da; i++) {
      mpz_urandomb(c, state, 256);
      mpz_mod(c, c, M);
      bignum coeff;
      mpz_to_bn(&coeff, c);
      field_in(&a.coeff[i], &coeff, &ctx);
      bn_free(&coeff);
    }
    a.deg = da;
    bigpoly_trim(&a);

    // build b (monic)
    bigpoly_alloc(&b, db);
    for (u64 i = 0; i < db; i++) {
      mpz_urandomb(c, state, 256);
      mpz_mod(c, c, M);
      bignum coeff;
      mpz_to_bn(&coeff, c);
      field_in(&b.coeff[i], &coeff, &ctx);
      bn_free(&coeff);
    }
    field_set_u64(&b.coeff[db], 1, &ctx);
    b.deg = db;

    bool ok = bigpoly_divmod(&q, &r, &a, &b, &ctx);
    if (!ok) {
      fprintf(stderr, "[FAIL] bigpoly_divmod returned false (iter %d)\n", iter);
      failures++;
    } else {
      if (r.deg >= b.deg) {
        fprintf(stderr, "[FAIL] deg(r) >= deg(b) (iter %d)\n", iter);
        failures++;
      } else {
        passes++;
      }

      // a - r must be divisible by b: (a - r) mod b == 0
      poly_ring ring;
      poly_ring_init(&ring, &b, &ctx);
      bigpoly t;
      bigpoly_init(&t);
      poly_ring_sub(&t, &a, &r, &ring);
      if (!poly_ring_is_zero(&t)) {
        fprintf(stderr, "[FAIL] a - r not divisible by b (iter %d)\n", iter);
        failures++;
      } else {
        passes++;
      }
      bigpoly_free(&t);
      poly_ring_free(&ring);
    }

    bigpoly_free(&a);
    bigpoly_free(&b);
    bigpoly_free(&q);
    bigpoly_free(&r);
  }

  mpz_clears(M, c, NULL);
  gmp_randclear(state);
  field_ctx_free(&ctx);
  bn_free(&m);
}

// ---- ring identities ------------------------------------------------------

static void test_ring(void)
{
  printf("--- Ring Z_17[x]/(x^2+3) identities ---\n");

  bignum m;
  bn_init(&m);
  bn_set_u64(&m, 17);
  field_ctx ctx;
  field_ctx_init(&ctx, &m);

  // q = x^2 + 3 (irreducible over GF(17), so the ring is the field GF(289))
  bigpoly q;
  bigpoly_init(&q);
  bigpoly_alloc(&q, 2);
  field_set_u64(&q.coeff[0], 3, &ctx);
  field_set_u64(&q.coeff[1], 0, &ctx);
  field_set_u64(&q.coeff[2], 1, &ctx);
  q.deg = 2;

  poly_ring ring;
  poly_ring_init(&ring, &q, &ctx);

  gmp_randstate_t state;
  gmp_randinit_default(state);
  gmp_randseed_ui(state, 0x1234);
  mpz_t M, c;
  mpz_inits(M, c, NULL);
  bn_to_mpz(M, &m);

  bigpoly one;
  bigpoly_init(&one);
  bigpoly_alloc(&one, 0);
  field_set_u64(&one.coeff[0], 1, &ctx);
  one.deg = 0;

  bignum e0, e1, e2, e12;
  bn_init_multi(&e0, &e1, &e2, &e12, NULL);
  bn_set_u64(&e0, 0);
  bn_set_u64(&e1, 1);
  bn_set_u64(&e2, 5);
  bn_add(&e12, &e1, &e2);  // e12 = 6

  for (int iter = 0; iter < 100; iter++) {
    bigpoly a, b, d;
    bigpoly_init(&a);
    bigpoly_init(&b);
    bigpoly_init(&d);
    fill_ring_elem(&a, &ctx, state, M);
    fill_ring_elem(&b, &ctx, state, M);
    fill_ring_elem(&d, &ctx, state, M);

    bigpoly ab, ad, bd, lhs, rhs, t1, t2, ainv, ae2;
    bigpoly_init(&ab);
    bigpoly_init(&ad);
    bigpoly_init(&bd);
    bigpoly_init(&lhs);
    bigpoly_init(&rhs);
    bigpoly_init(&t1);
    bigpoly_init(&t2);
    bigpoly_init(&ainv);
    bigpoly_init(&ae2);

    // a * 1 == a
    poly_ring_mul(&t1, &a, &one, &ring);
    if (!poly_ring_equal(&t1, &a)) {
      fprintf(stderr, "[FAIL] a*1 != a (iter %d)\n", iter);
      failures++;
    } else {
      passes++;
    }

    // a * a_inv == 1 (a nonzero; ring is a field)
    if (!poly_ring_is_zero(&a)) {
      bool ok = poly_ring_inv(&ainv, &a, &ring);
      if (!ok) {
        fprintf(stderr, "[FAIL] poly_ring_inv failed for nonzero a (iter %d)\n",
                iter);
        failures++;
      } else {
        poly_ring_mul(&t2, &a, &ainv, &ring);
        if (!poly_ring_equal(&t2, &one)) {
          fprintf(stderr, "[FAIL] a*a_inv != 1 (iter %d)\n", iter);
          failures++;
        } else {
          passes++;
        }
      }
    }

    // (a + b) * d == a*d + b*d
    poly_ring_add(&ab, &a, &b, &ring);
    poly_ring_mul(&lhs, &ab, &d, &ring);
    poly_ring_mul(&ad, &a, &d, &ring);
    poly_ring_mul(&bd, &b, &d, &ring);
    poly_ring_add(&rhs, &ad, &bd, &ring);
    if (!poly_ring_equal(&lhs, &rhs)) {
      fprintf(stderr, "[FAIL] (a+b)*d != a*d+b*d (iter %d)\n", iter);
      failures++;
    } else {
      passes++;
    }

    // a * (b + d) == a*b + a*d
    poly_ring_add(&ab, &b, &d, &ring);
    poly_ring_mul(&lhs, &a, &ab, &ring);
    poly_ring_mul(&bd, &b, &a, &ring);
    poly_ring_add(&rhs, &ad, &bd, &ring);
    if (!poly_ring_equal(&lhs, &rhs)) {
      fprintf(stderr, "[FAIL] a*(b+d) != a*b+a*d (iter %d)\n", iter);
      failures++;
    } else {
      passes++;
    }

    // a * (b * d) == (a * b) * d
    poly_ring_mul(&ab, &b, &d, &ring);
    poly_ring_mul(&lhs, &a, &ab, &ring);
    poly_ring_mul(&bd, &a, &b, &ring);
    poly_ring_mul(&rhs, &bd, &d, &ring);
    if (!poly_ring_equal(&lhs, &rhs)) {
      fprintf(stderr, "[FAIL] a*(b*d) != (a*b)*d (iter %d)\n", iter);
      failures++;
    } else {
      passes++;
    }

    // a^0 == 1, a^1 == a
    poly_ring_pow(&t1, &a, &e0, &ring);
    if (!poly_ring_equal(&t1, &one)) {
      fprintf(stderr, "[FAIL] a^0 != 1 (iter %d)\n", iter);
      failures++;
    } else {
      passes++;
    }
    poly_ring_pow(&t2, &a, &e1, &ring);
    if (!poly_ring_equal(&t2, &a)) {
      fprintf(stderr, "[FAIL] a^1 != a (iter %d)\n", iter);
      failures++;
    } else {
      passes++;
    }

    // a^6 == a^1 * a^5
    poly_ring_pow(&t1, &a, &e12, &ring);
    poly_ring_pow(&t2, &a, &e1, &ring);
    poly_ring_pow(&ae2, &a, &e2, &ring);
    poly_ring_mul(&t2, &t2, &ae2, &ring);
    if (!poly_ring_equal(&t1, &t2)) {
      fprintf(stderr, "[FAIL] a^6 != a^1*a^5 (iter %d)\n", iter);
      failures++;
    } else {
      passes++;
    }

    bigpoly_free(&ab);
    bigpoly_free(&ad);
    bigpoly_free(&bd);
    bigpoly_free(&lhs);
    bigpoly_free(&rhs);
    bigpoly_free(&t1);
    bigpoly_free(&t2);
    bigpoly_free(&ainv);
    bigpoly_free(&ae2);
    bigpoly_free(&a);
    bigpoly_free(&b);
    bigpoly_free(&d);
  }

  bn_free_multi(&e0, &e1, &e2, &e12, NULL);
  bigpoly_free(&one);
  mpz_clears(M, c, NULL);
  gmp_randclear(state);
  poly_ring_free(&ring);
  bigpoly_free(&q);
  field_ctx_free(&ctx);
  bn_free(&m);
}

int main(void)
{
  bn_init_constants();

  test_field_mont();
  test_field_even();
  test_poly_divmod();
  test_ring();

  printf("\n=== test_field summary: %d passed, %d failed ===\n", passes,
         failures);

  bn_free_constants();
  return failures ? 1 : 0;
}
