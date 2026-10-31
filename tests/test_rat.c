/*
 * test_rat.c
 *
 * Assert-based checks for rq_t arithmetic invariants (normalized results,
 * positive denominators), the strict decimal parse/print API, randomized
 * algebraic properties, rmatq_t add/sub with known answers, and an
 * optional cross-check against GMP's mpq type (auto-skipped without GMP).
 *
 * Copyright (C) 2026 Diego Strebel
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <errno.h>
#include <stdlib.h>
#include <string.h>

#include "../include/rmat.h"
#include "../include/rmatq.h"
#include "../include/rq.h"
#include "test_util.h"

#if defined(__has_include)
#if __has_include(<gmp.h>)
#include <gmp.h>
#define TEST_RAT_HAVE_GMP 1
#endif
#endif
#ifndef TEST_RAT_HAVE_GMP
#define TEST_RAT_HAVE_GMP 0
#endif

static void rq_set_pair(rq_t* r, i64 num, i64 den)
{
  rz_set_i64(&r->num, num);
  rz_set_i64(&r->den, den);
}

static void rq_set_num_den_str(rq_t* r, i64 num, const char* den)
{
  rz_set_i64(&r->num, num);
  rz_init_val(&r->den, den);
}

/* Compare two rationals independent of their current reduction level. */
static bool rq_equal(const rq_t* a, const rq_t* b)
{
  rq_t x, y;
  rq_init_multi(&x, &y, NULL);
  rq_copy(&x, a);
  rq_copy(&y, b);
  rq_normalize(&x);
  rq_normalize(&y);
  bool eq = (rq_cmp(&x, &y) == 0);
  rq_clear_multi(&x, &y, NULL);
  return eq;
}

/* A normalized rational must have a positive denominator. */
static bool rq_den_positive(const rq_t* r)
{
  return !r->den.is_neg && !rz_is_zero(&r->den);
}

/*
 * Full canonical-form audit: den > 0, gcd(num, den) == 1, and zero stored
 * as 0/1. Canonical form is unique, so rq_cmp() doubles as exact equality
 * for values that pass this check.
 */
static void rq_require_normalized(const rq_t* r, const char* what)
{
  if (!rq_den_positive(r)) {
    TEST_ASSERT(false, what);
    printf("  (%s: non-positive denominator)\n", what);
    return;
  }

  rz_t g;
  if (rz_init(&g) != RABIN_SUCCESS) {
    TEST_ASSERT(false, what);
    return;
  }
  bool ok = (rz_gcd(&g, &r->num, &r->den) == RABIN_SUCCESS) && rz_is_one(&g);
  rz_clear(&g);
  if (!ok) {
    TEST_ASSERT(false, what);
    printf("  (%s: gcd(num, den) != 1)\n", what);
  }

  if (rz_is_zero(&r->num) && !rz_is_one(&r->den)) {
    TEST_ASSERT(false, what);
    printf("  (%s: zero not stored as 0/1)\n", what);
  }
}

/*
 * Exact comparison of a small rational against num/den given as i128,
 * cross-multiplying in __int128. Only valid for values whose decimal
 * forms fit in long long (all inputs here are 30 bits or less).
 */
static bool rq_eq_i128(const rq_t* r, __int128 num, __int128 den)
{
  char* ns = rz_to_string(&r->num);
  char* ds = rz_to_string(&r->den);
  if (ns == NULL || ds == NULL) {
    free(ns);
    free(ds);
    return false;
  }

  errno = 0;
  long long rn = strtoll(ns, NULL, 10);
  long long rd = strtoll(ds, NULL, 10);
  free(ns);
  free(ds);
  if (errno != 0) return false;

  return (__int128)rn * den == num * (__int128)rd;
}

/* Deterministic random 30-bit magnitude with random sign. */
static i64 rand_signed_30(void)
{
  i64 v = (i64)(rand() & ((1 << 30) - 1));
  v += ((i64)(rand() & ((1 << 30) - 1))) << 15; /* wider spread, still < 2^31 */
  v &= (1 << 30) - 1;
  return (rand() & 1) ? -v : v;
}

static i64 rand_pos_30(void)
{
  i64 v = (i64)(rand() & ((1 << 30) - 1)) | 1;
  return v > 0 ? v : 1;
}

TEST(rq_parse_strict)
{
  static const struct {
    const char* str;
    i64 num;
    i64 den;
  } valid[] = {
      {"3/4", 3, 4},   {"-3/4", -3, 4}, {"+3/4", 3, 4},
      {" 2/6 ", 1, 3}, /* whitespace only at the ends, reduces */
      {"0", 0, 1},     {"0/5", 0, 1},   {"17", 17, 1},
      {"-0/7", 0, 1},  {"-1", -1, 1},
  };

  rq_t r;
  rq_init(&r);
  rz_t n, d;
  rz_init_multi(&n, &d, NULL);

  for (size_t i = 0; i < sizeof(valid) / sizeof(valid[0]); i++) {
    rabin_err_t err = rq_init_val(&r, valid[i].str);
    TEST_ASSERT(err == RABIN_SUCCESS, "valid parse");
    if (err != RABIN_SUCCESS) {
      printf("  (\"%s\" -> %d)\n", valid[i].str, err);
      continue;
    }
    rq_require_normalized(&r, "parsed value normalized");

    rq_get_num(&n, &r);
    rq_get_den(&d, &r);
    bool ok = rz_is_eq_i64(&n, valid[i].num) && rz_is_eq_i64(&d, valid[i].den);
    TEST_ASSERT(ok, "parse value matches");
    if (!ok) printf("  (\"%s\" -> wrong num/den)\n", valid[i].str);
  }

  /* Large numerator with denominator: exactness past 64 bits. */
  TEST_ASSERT(rq_init_val(&r, "36893488147419103232/3") == RABIN_SUCCESS,
              "big parse");
  rq_require_normalized(&r, "big parse normalized");
  rq_get_num(&n, &r);
  rq_get_den(&d, &r);
  TEST_ASSERT(rz_is_eq_i64(&d, 3), "big parse den");
  {
    rz_t three, expect;
    rz_init_multi(&three, &expect, NULL);
    rz_set_u64(&three, 3);
    rz_init_val(&expect, "36893488147419103232");
    TEST_ASSERT(rz_cmp(&n, &expect) == 0, "big parse num exact");
    rz_clear_multi(&three, &expect, NULL);
  }

  static const char* bad[] = {"",    "abc",  "3/",  "/4",   "3/4/5",
                              "3.5", "1 2",  "+",   "-",    "3/x",
                              "++3", "3/-4", "1e3", "0x10", " "};

  for (size_t i = 0; i < sizeof(bad) / sizeof(bad[0]); i++) {
    rabin_err_t err = rq_init_val(&r, bad[i]);
    TEST_ASSERT(err == RABIN_ERR_INVALID_ARG, "malformed rejected");
    if (err != RABIN_ERR_INVALID_ARG) printf("  (\"%s\" -> %d)\n", bad[i], err);
    TEST_ASSERT(rz_is_zero(&r.num) && rz_is_zero(&r.den),
                "failed parse leaves 0/0");
  }

  static const char* zero_den[] = {"3/0", "0/0", "-7/0"};
  for (size_t i = 0; i < sizeof(zero_den) / sizeof(zero_den[0]); i++) {
    rabin_err_t err = rq_init_val(&r, zero_den[i]);
    TEST_ASSERT(err == RABIN_ERR_DIV_BY_ZERO, "zero den rejected");
    if (err != RABIN_ERR_DIV_BY_ZERO)
      printf("  (\"%s\" -> %d)\n", zero_den[i], err);
  }

  TEST_ASSERT(rq_init_val(NULL, "1/2") == RABIN_ERR_NULL_PTR, "NULL r");
  TEST_ASSERT(rq_init_val(&r, NULL) == RABIN_ERR_NULL_PTR, "NULL str");

  rq_clear(&r);
  rz_clear_multi(&n, &d, NULL);
}

TEST(rq_to_string_format)
{
  static const struct {
    const char* src; /* NULL: pair built directly below */
    i64 num;
    i64 den;
    const char* expect;
  } cases[] = {
      {"3/4", 0, 0, "3/4"},   {"-3/4", 0, 0, "-3/4"},
      {NULL, 3, -4, "-3/4"}, /* unnormalized input: printed after norm. */
      {"6/8", 0, 0, "3/4"},   {"7/1", 0, 0, "7"}, /* den == 1 omits "/1" */
      {"-7", 0, 0, "-7"},     {"0/9", 0, 0, "0"},
      {"22/7", 0, 0, "22/7"},
  };

  rq_t r;
  rq_init(&r);

  for (size_t i = 0; i < sizeof(cases) / sizeof(cases[0]); i++) {
    if (cases[i].src != NULL)
      TEST_ASSERT(rq_init_val(&r, cases[i].src) == RABIN_SUCCESS, "parse case");
    else {
      rq_set_pair(&r, cases[i].num, cases[i].den);
      TEST_ASSERT(rq_normalize(&r) == RABIN_SUCCESS, "normalize case");
    }

    size_t need = rq_to_string_len(&r);
    TEST_ASSERT(need == strlen(cases[i].expect) + 1, "len exact");
    if (need != strlen(cases[i].expect) + 1)
      printf("  (%s: len %zu, expected %zu)\n", cases[i].src, need,
             strlen(cases[i].expect) + 1);

    char buf[64];
    memset(buf, '#', sizeof buf);
    TEST_ASSERT(rq_to_string(&r, buf, need) == RABIN_SUCCESS, "fit exact cap");
    TEST_ASSERT(strcmp(buf, cases[i].expect) == 0, "string exact");
    if (strcmp(buf, cases[i].expect) != 0)
      printf("  (%s -> \"%s\", expected \"%s\")\n", cases[i].src, buf,
             cases[i].expect);

    /* One byte too small: error and out untouched. */
    memset(buf, '#', sizeof buf);
    TEST_ASSERT(rq_to_string(&r, buf, need - 1) == RABIN_ERR_INVALID_ARG,
                "too-small cap rejected");
    TEST_ASSERT(buf[0] == '#', "too-small cap leaves out untouched");
  }

  /* Round-trip on a big value. */
  TEST_ASSERT(rq_init_val(&r, "36893488147419103232/3") == RABIN_SUCCESS,
              "big parse");
  char big[80];
  size_t need = rq_to_string_len(&r);
  TEST_ASSERT(need <= sizeof big, "big len fits buffer");
  TEST_ASSERT(rq_to_string(&r, big, need) == RABIN_SUCCESS, "big print");
  TEST_ASSERT(strcmp(big, "36893488147419103232/3") == 0, "big round trip");

  /* Round-trip: parse -> print -> parse is identity. */
  rq_t q;
  rq_init(&q);
  TEST_ASSERT(rq_init_val(&q, big) == RABIN_SUCCESS, "reparse printed");
  TEST_ASSERT(rq_cmp(&r, &q) == 0, "round trip exact");

  TEST_ASSERT(rq_to_string_len(NULL) == 0, "NULL len");
  TEST_ASSERT(rq_to_string(NULL, big, sizeof big) == RABIN_ERR_NULL_PTR,
              "NULL print r");
  TEST_ASSERT(rq_to_string(&r, NULL, 10) == RABIN_ERR_NULL_PTR,
              "NULL print out");

  rq_clear_multi(&r, &q, NULL);
}

TEST(rq_inv_and_helpers)
{
  rq_t r, q, expect;
  rq_init_multi(&r, &q, &expect, NULL);

  TEST_ASSERT(rq_init_val(&r, "3/4") == RABIN_SUCCESS, "parse 3/4");
  TEST_ASSERT(rq_inv(&q, &r) == RABIN_SUCCESS, "rq_inv");
  rq_require_normalized(&q, "inv normalized");
  TEST_ASSERT(rq_init_val(&expect, "4/3") == RABIN_SUCCESS, "parse 4/3");
  TEST_ASSERT(rq_cmp(&q, &expect) == 0, "inv(3/4) == 4/3");

  TEST_ASSERT(rq_init_val(&r, "-5/7") == RABIN_SUCCESS, "parse -5/7");
  TEST_ASSERT(rq_inv(&q, &r) == RABIN_SUCCESS, "rq_inv negative");
  rq_require_normalized(&q, "inv negative normalized");
  TEST_ASSERT(rq_init_val(&expect, "-7/5") == RABIN_SUCCESS, "parse -7/5");
  TEST_ASSERT(rq_cmp(&q, &expect) == 0, "inv(-5/7) == -7/5");

  /* Aliasing: r = 1/r in place, twice, is the identity. */
  TEST_ASSERT(rq_init_val(&r, "2/9") == RABIN_SUCCESS, "parse 2/9");
  TEST_ASSERT(rq_copy(&expect, &r) == RABIN_SUCCESS, "save expect");
  TEST_ASSERT(rq_inv(&r, &r) == RABIN_SUCCESS, "aliased inv");
  TEST_ASSERT(rq_inv(&r, &r) == RABIN_SUCCESS, "aliased inv again");
  TEST_ASSERT(rq_cmp(&r, &expect) == 0, "inv(inv(a)) == a (aliased)");

  TEST_ASSERT(rq_set_i64(&r, 0) == RABIN_SUCCESS, "set_i64 0");
  TEST_ASSERT(rq_inv(&q, &r) == RABIN_ERR_DIV_BY_ZERO, "inv(0) rejected");

  TEST_ASSERT(rq_set_i64(&r, -42) == RABIN_SUCCESS, "set_i64 -42");
  rq_require_normalized(&r, "set_i64 normalized");
  TEST_ASSERT(rq_cmp_i64(&r, -42) == 0, "set_i64 exact");

  /* cmp_i64 on both signs. */
  TEST_ASSERT(rq_init_val(&r, "22/7") == RABIN_SUCCESS, "parse 22/7");
  TEST_ASSERT(rq_cmp_i64(&r, 2) > 0, "22/7 > 2");
  TEST_ASSERT(rq_cmp_i64(&r, 3) > 0, "22/7 > 3");
  TEST_ASSERT(rq_cmp_i64(&r, 4) < 0, "22/7 < 4");
  TEST_ASSERT(rq_init_val(&r, "-22/7") == RABIN_SUCCESS, "parse -22/7");
  TEST_ASSERT(rq_cmp_i64(&r, -3) < 0, "-22/7 < -3");
  TEST_ASSERT(rq_cmp_i64(&r, -4) > 0, "-22/7 > -4");
  TEST_ASSERT(rq_set_i64(&r, 5) == RABIN_SUCCESS, "set_i64 5 exact");
  TEST_ASSERT(rq_cmp_i64(&r, 5) == 0, "integer rational cmp exact");

  /* get_num/get_den deep-copy: mutating the source leaves the copy alone. */
  TEST_ASSERT(rq_init_val(&r, "22/7") == RABIN_SUCCESS, "parse 22/7 again");
  rz_t n, d;
  rz_init_multi(&n, &d, NULL);
  TEST_ASSERT(rq_get_num(&n, &r) == RABIN_SUCCESS, "get_num");
  TEST_ASSERT(rq_get_den(&d, &r) == RABIN_SUCCESS, "get_den");
  rz_set_u64(&r.num, 999);
  TEST_ASSERT(rz_is_eq_i64(&n, 22), "get_num is a copy");
  TEST_ASSERT(rz_is_eq_i64(&d, 7), "get_den is a copy");
  rz_clear_multi(&n, &d, NULL);

  rq_clear_multi(&r, &q, &expect, NULL);
}

TEST(rat_random_properties)
{
  srand(42);

  rq_t a, b, c, l, rr;
  rq_init_multi(&a, &b, &c, &l, &rr, NULL);

  const int rounds = 150;
  for (int i = 0; i < rounds; i++) {
    i64 an = rand_signed_30(), ad = rand_pos_30();
    i64 bn = rand_signed_30(), bd = rand_pos_30();
    i64 cn = rand_signed_30(), cd = rand_pos_30();

    rq_set_pair(&a, an, ad);
    rq_set_pair(&b, bn, bd);
    rq_set_pair(&c, cn, cd);
    TEST_ASSERT(rq_normalize(&a) == RABIN_SUCCESS, "normalize a");
    TEST_ASSERT(rq_normalize(&b) == RABIN_SUCCESS, "normalize b");
    TEST_ASSERT(rq_normalize(&c) == RABIN_SUCCESS, "normalize c");

    /* add: commutative, and matches the exact cross-product oracle. */
    TEST_ASSERT(rq_add(&l, &a, &b) == RABIN_SUCCESS, "add ab");
    rq_require_normalized(&l, "a+b normalized");
    TEST_ASSERT(rq_add(&rr, &b, &a) == RABIN_SUCCESS, "add ba");
    TEST_ASSERT(rq_cmp(&l, &rr) == 0, "add commutative");
    TEST_ASSERT(rq_eq_i128(&l, (__int128)an * bd + (__int128)bn * ad,
                           (__int128)ad * bd),
                "a+b matches exact oracle");

    /* sub: a - a == 0 with canonical 0/1; oracle for a - b. */
    TEST_ASSERT(rq_sub(&l, &a, &a) == RABIN_SUCCESS, "sub aa");
    rq_require_normalized(&l, "a-a normalized");
    TEST_ASSERT(rq_is_zero(&l), "a - a == 0");
    TEST_ASSERT(rq_sub(&l, &a, &b) == RABIN_SUCCESS, "sub ab");
    TEST_ASSERT(rq_eq_i128(&l, (__int128)an * bd - (__int128)bn * ad,
                           (__int128)ad * bd),
                "a-b matches exact oracle");

    /* mul: commutative with oracle. */
    TEST_ASSERT(rq_mul(&l, &a, &b) == RABIN_SUCCESS, "mul ab");
    rq_require_normalized(&l, "a*b normalized");
    TEST_ASSERT(rq_mul(&rr, &b, &a) == RABIN_SUCCESS, "mul ba");
    TEST_ASSERT(rq_cmp(&l, &rr) == 0, "mul commutative");
    TEST_ASSERT(rq_eq_i128(&l, (__int128)an * bn, (__int128)ad * bd),
                "a*b matches exact oracle");

    /* sqr == mul. */
    TEST_ASSERT(rq_sqr(&l, &a) == RABIN_SUCCESS, "sqr a");
    rq_require_normalized(&l, "a*a sqr normalized");
    TEST_ASSERT(rq_mul(&rr, &a, &a) == RABIN_SUCCESS, "mul aa");
    TEST_ASSERT(rq_cmp(&l, &rr) == 0, "sqr(a) == a*a");

    /* associativity on normalized canonical forms (exact compare). */
    TEST_ASSERT(rq_add(&l, &a, &b) == RABIN_SUCCESS, "assoc add lhs mul");
    TEST_ASSERT(rq_add(&l, &l, &c) == RABIN_SUCCESS, "assoc add lhs");
    rq_require_normalized(&l, "(a+b)+c normalized");
    TEST_ASSERT(rq_add(&rr, &b, &c) == RABIN_SUCCESS, "assoc add rhs inner");
    TEST_ASSERT(rq_add(&rr, &a, &rr) == RABIN_SUCCESS, "assoc add rhs");
    TEST_ASSERT(rq_cmp(&l, &rr) == 0, "add associative");

    /* a * inv(a) == 1 exactly, and inv(inv(a)) == a. */
    if (!rz_is_zero(&a.num)) {
      TEST_ASSERT(rq_inv(&l, &a) == RABIN_SUCCESS, "inv a");
      rq_require_normalized(&l, "inv(a) normalized");
      TEST_ASSERT(rq_mul(&rr, &a, &l) == RABIN_SUCCESS, "a * inv(a)");
      rq_require_normalized(&rr, "a*inv(a) normalized");
      TEST_ASSERT(rq_cmp_i64(&rr, 1) == 0, "a * inv(a) == 1");

      TEST_ASSERT(rq_inv(&l, &l) == RABIN_SUCCESS, "inv inv a");
      TEST_ASSERT(rq_cmp(&l, &a) == 0, "inv(inv(a)) == a");
    }

    /* division oracle: (a/b) * b == a. */
    if (!rz_is_zero(&b.num)) {
      TEST_ASSERT(rq_div(&l, &a, &b) == RABIN_SUCCESS, "div ab");
      rq_require_normalized(&l, "a/b normalized");
      TEST_ASSERT(rq_eq_i128(&l, (__int128)an * bd, (__int128)ad * bn),
                  "a/b matches exact oracle");
      TEST_ASSERT(rq_mul(&rr, &l, &b) == RABIN_SUCCESS, "re-mul div");
      TEST_ASSERT(rq_cmp(&rr, &a) == 0, "(a/b)*b == a");
    }
  }

  rq_clear_multi(&a, &b, &c, &l, &rr, NULL);
}

TEST(rmatq_negative_den_regression)
{
  rmatq_t A = {0};
  rq_t a, got;
  rq_init_multi(&a, &got, NULL);

  TEST_ASSERT(rmatq_init(&A, 2, 2) == RABIN_SUCCESS, "rmatq_init");

  /* rmatq_set with a negative-denominator input must not poison the
   * shared denominator; rmatq_get must return a normalized value.
   * 3/-4 == -3/4: the sign belongs in the numerator. */
  rq_set_pair(&a, 3, -4);
  TEST_ASSERT(rmatq_set(&A, &a, 0, 0) == RABIN_SUCCESS, "set 3/-4");
  TEST_ASSERT(!A.den.is_neg && !rz_is_zero(&A.den), "shared den stays > 0");
  TEST_ASSERT(rmatq_get(&got, &A, 0, 0) == RABIN_SUCCESS, "get back");
  rq_require_normalized(&got, "rmatq_get entry normalized");
  rq_set_pair(&a, -3, 4);
  TEST_ASSERT(rq_normalize(&a) == RABIN_SUCCESS, "normalize expect");
  TEST_ASSERT(rq_cmp(&got, &a) == 0, "3/-4 stored as -3/4");

  /* Hand-built hostile state: shared den < 0 encoding 1/2 and -3/5.
   * Normalizing must fold the sign into the entries and keep every
   * represented value identical. */
  rmatq_clear(&A);
  TEST_ASSERT(rmatq_init(&A, 2, 2) == RABIN_SUCCESS, "rmatq_init again");
  rz_t e;
  rz_init(&e);
  rz_set_i64(&A.den, -10);
  rz_set_i64(&e, -5);
  rmat_set(&A.M, &e, 0, 0); /* -5/-10 = 1/2 */
  rz_set_i64(&e, 6);
  rmat_set(&A.M, &e, 1, 1); /*  6/-10 = -3/5 */
  TEST_ASSERT(rmatq_normalize(&A) == RABIN_SUCCESS, "normalize negative den");
  TEST_ASSERT(!A.den.is_neg && !rz_is_zero(&A.den),
              "den positive after normalize");
  TEST_ASSERT(rmatq_get(&got, &A, 0, 0) == RABIN_SUCCESS, "get (0,0)");
  rq_set_pair(&a, 1, 2);
  TEST_ASSERT(rq_normalize(&a) == RABIN_SUCCESS, "normalize 1/2");
  TEST_ASSERT(rq_cmp(&got, &a) == 0, "value (0,0) preserved");
  TEST_ASSERT(rmatq_get(&got, &A, 1, 1) == RABIN_SUCCESS, "get (1,1)");
  rq_set_pair(&a, -3, 5);
  TEST_ASSERT(rq_normalize(&a) == RABIN_SUCCESS, "normalize -3/5");
  TEST_ASSERT(rq_cmp(&got, &a) == 0, "value (1,1) preserved");
  rz_clear(&e);

  rmatq_clear(&A);
  rq_clear_multi(&a, &got, NULL);
}

TEST(rq_gmp_crosscheck)
{
#if TEST_RAT_HAVE_GMP
  srand(1337);

  rq_t a, b, r_lib;
  rq_init_multi(&a, &b, &r_lib, NULL);
  mpq_t ga, gb, gr;
  mpq_init(ga);
  mpq_init(gb);
  mpq_init(gr);

  char buf[128];
  const int rounds = 60;
  for (int i = 0; i < rounds; i++) {
    i64 an = rand_signed_30(), ad = rand_pos_30();
    i64 bn = rand_signed_30(), bd = rand_pos_30();

    snprintf(buf, sizeof buf, "%lld/%lld", (long long)an, (long long)ad);
    TEST_ASSERT(rq_init_val(&a, buf) == RABIN_SUCCESS, "gmp case parse a");
    mpq_set_str(ga, buf, 10);
    mpq_canonicalize(ga);

    snprintf(buf, sizeof buf, "%lld/%lld", (long long)bn, (long long)bd);
    TEST_ASSERT(rq_init_val(&b, buf) == RABIN_SUCCESS, "gmp case parse b");
    mpq_set_str(gb, buf, 10);
    mpq_canonicalize(gb);

    /* parse cross-check: string -> both libs -> compare via GMP */
    TEST_ASSERT(rq_to_string(&a, buf, sizeof buf) == RABIN_SUCCESS, "print a");
    mpq_set_str(gr, buf, 10);
    mpq_canonicalize(gr);
    TEST_ASSERT(mpq_cmp(gr, ga) == 0, "parse+print agrees with mpq");

    struct {
      const char* op;
      rabin_err_t (*lib)(rq_t*, const rq_t*, const rq_t*);
      void (*gmp)(mpq_ptr, mpq_srcptr, mpq_srcptr);
    } ops[] = {
        {"add", rq_add, mpq_add},
        {"sub", rq_sub, mpq_sub},
        {"mul", rq_mul, mpq_mul},
        {"div", rq_div, mpq_div},
    };

    for (size_t k = 0; k < sizeof(ops) / sizeof(ops[0]); k++) {
      if (k == 3 && mpq_sgn(gb) == 0) continue; /* skip div by zero */

      TEST_ASSERT(ops[k].lib(&r_lib, &a, &b) == RABIN_SUCCESS, ops[k].op);
      rq_require_normalized(&r_lib, "gmp result normalized");
      ops[k].gmp(gr, ga, gb);

      TEST_ASSERT(rq_to_string(&r_lib, buf, sizeof buf) == RABIN_SUCCESS,
                  "print result");
      mpq_t check;
      mpq_init(check);
      mpq_set_str(check, buf, 10);
      mpq_canonicalize(check);
      TEST_ASSERT(mpq_cmp(check, gr) == 0, ops[k].op);
      if (mpq_cmp(check, gr) != 0)
        printf("  (%s: %s vs %s)\n", ops[k].op, buf, mpq_get_str(NULL, 10, gr));
      mpq_clear(check);
    }
  }

  mpq_clear(ga);
  mpq_clear(gb);
  mpq_clear(gr);
  rq_clear_multi(&a, &b, &r_lib, NULL);
#else
  printf("GMP not available: rq_gmp_crosscheck skipped\n");
#endif
}

TEST(rat_arithmetic)
{
  rq_t r, q, p, expect;
  rq_init_multi(&r, &q, &p, &expect, NULL);

  rq_set_pair(&r, 3, 8905493208);
  rq_set_pair(&q, 3, 23534563456);

  rq_print(&r);
  TEST_ASSERT(rq_mul(&p, &r, &q) == RABIN_SUCCESS, "rq_mul");
  rq_normalize(&p);
  TEST_ASSERT(rq_den_positive(&p), "product has positive denominator");
  rq_set_num_den_str(&expect, 9, "209586895010653006848");
  TEST_ASSERT(rq_equal(&p, &expect), "3/8905493208 * 3/23534563456 correct");

  rq_print(&p);
  rq_print(&r);

  rq_set_pair(&r, -1, 4);
  rq_set_pair(&q, 1, 4);

  TEST_ASSERT(rq_div(&p, &r, &q) == RABIN_SUCCESS, "rq_div");
  rq_normalize(&p);
  rq_set_pair(&expect, -1, 1);
  TEST_ASSERT(rq_equal(&p, &expect), "(-1/4) / (1/4) == -1");

  TEST_ASSERT(rq_add(&p, &r, &q) == RABIN_SUCCESS, "rq_add");
  rq_normalize(&p);
  TEST_ASSERT(rq_is_zero(&p), "(-1/4) + (1/4) == 0");

  /* Roundtrip: (x*y)/y == x with non-trivial rationals. */
  rq_set_pair(&r, 5, 6);
  rq_set_pair(&q, -7, 15);
  TEST_ASSERT(rq_mul(&p, &r, &q) == RABIN_SUCCESS, "roundtrip mul");
  TEST_ASSERT(rq_div(&p, &p, &q) == RABIN_SUCCESS, "roundtrip div");
  TEST_ASSERT(rq_equal(&p, &r), "(x*y)/y == x");

  /* Negation twice is identity. */
  TEST_ASSERT(rq_neg(&p) == RABIN_SUCCESS, "rq_neg");
  TEST_ASSERT(rq_neg(&p) == RABIN_SUCCESS, "rq_neg again");
  TEST_ASSERT(rq_equal(&p, &r), "-(-x) == x");

  rq_clear_multi(&r, &q, &p, &expect, NULL);
}

TEST(rat_matrix_addsub)
{
  rmatq_t A = {0}, B = {0}, R_add = {0}, R_sub = {0};
  rq_t entry, expect;
  rq_init_multi(&entry, &expect, NULL);

  TEST_ASSERT(rmatq_init(&A, 5, 5) == RABIN_SUCCESS, "rmatq_init A");
  TEST_ASSERT(rmatq_init(&B, 5, 5) == RABIN_SUCCESS, "rmatq_init B");
  TEST_ASSERT(rmatq_init(&R_add, 5, 5) == RABIN_SUCCESS, "rmatq_init R_add");
  TEST_ASSERT(rmatq_init(&R_sub, 5, 5) == RABIN_SUCCESS, "rmatq_init R_sub");

  rz_set_u64(&A.den, 2);
  rz_set_u64(&B.den, 3);

  rz_t val;
  rz_init(&val);
  for (u64 i = 0; i < 5; i++) {
    for (u64 j = 0; j < 5; j++) {
      rz_set_u64(&val, i + j + 1);
      rmat_set(&A.M, &val, i, j);
      rz_set_u64(&val, i * j + 2);
      rmat_set(&B.M, &val, i, j);
    }
  }
  rz_clear(&val);

  TEST_ASSERT(rmatq_add(&R_add, &A, &B) == RABIN_SUCCESS, "rmatq_add");
  TEST_ASSERT(rmatq_normalize(&R_add) == RABIN_SUCCESS, "normalize R_add");
  TEST_ASSERT(rmatq_sub(&R_sub, &A, &B) == RABIN_SUCCESS, "rmatq_sub");
  TEST_ASSERT(rmatq_normalize(&R_sub) == RABIN_SUCCESS, "normalize R_sub");

  printf("Matrix A + B:\n");
  rmatq_print(&R_add);
  printf("Matrix A - B:\n");
  rmatq_print(&R_sub);

  for (u64 i = 0; i < 5; i++) {
    for (u64 j = 0; j < 5; j++) {
      /* A[i][j] + B[i][j] = (3(i+j+1) + 2(ij+2)) / 6 */
      TEST_ASSERT(rmatq_get(&entry, &R_add, i, j) == RABIN_SUCCESS, "get add");
      TEST_ASSERT(rq_den_positive(&entry), "add entry den > 0");
      rq_set_pair(&expect, (i64)(3 * (i + j + 1) + 2 * (i * j + 2)), 6);
      TEST_ASSERT(rq_equal(&entry, &expect), "add entry matches known value");

      /* A[i][j] - B[i][j] = (3(i+j+1) - 2(ij+2)) / 6 */
      TEST_ASSERT(rmatq_get(&entry, &R_sub, i, j) == RABIN_SUCCESS, "get sub");
      rq_set_pair(&expect, (i64)(3 * (i + j + 1) - 2 * (i * j + 2)), 6);
      TEST_ASSERT(rq_equal(&entry, &expect), "sub entry matches known value");
    }
  }

  rq_clear_multi(&entry, &expect, NULL);
  rmatq_clear(&A);
  rmatq_clear(&B);
  rmatq_clear(&R_add);
  rmatq_clear(&R_sub);
}

int main(void)
{
  rz_init_constants();

  RUN_TEST(rq_parse_strict);
  RUN_TEST(rq_to_string_format);
  RUN_TEST(rq_inv_and_helpers);
  RUN_TEST(rat_random_properties);
  RUN_TEST(rmatq_negative_den_regression);
  RUN_TEST(rq_gmp_crosscheck);
  RUN_TEST(rat_arithmetic);
  RUN_TEST(rat_matrix_addsub);

  rz_clear_constants();
  return TEST_REPORT();
}
