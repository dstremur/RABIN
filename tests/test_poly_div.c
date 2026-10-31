/*
 * test_poly_div.c
 *
 * Correctness suite for the Z[x] division and gcd stack:
 * rpol_pseudo_div, rpol_div_exact, rpol_content, rpol_primitive_part,
 * rpol_gcd, rpol_eval and rpol_derive.
 *
 * Randomized gcd results are cross-checked against two independent
 * oracles: a naive pseudo-remainder sequence written on small __int128
 * coefficient arrays (this file only), and FLINT's fmpz_poly_gcd when
 * libflint is available.
 *
 * Copyright (C) 2026 Diego Strebel
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <limits.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "../include/rabin.h"
#include "../include/rpol.h"

#if __has_include(<flint/fmpz_poly.h>)
#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>
#define HAVE_FLINT 1
#endif

static int failures = 0;

#define CHECK(cond, msg)                                       \
  do {                                                         \
    if (!(cond)) {                                             \
      printf("FAILED: %s (%s:%d)\n", msg, __FILE__, __LINE__); \
      failures++;                                              \
    }                                                          \
  } while (0)

/* ------------------------------------------------------------------ */
/* small helpers                                                       */
/* ------------------------------------------------------------------ */

static void mk(rpol_t* p, const i64* c, u64 deg)
{
  CHECK(rpol_set_i64(p, c, deg) == RABIN_SUCCESS, "rpol_set_i64");
}

/* true if p is the zero polynomial, untrimmed leading zeros included */
static bool poly_is_zero(const rpol_t* p)
{
  if (p->size == 0) return true;
  for (u64 i = 0; i <= p->deg; i++) {
    if (!rz_is_zero(&p->coeff[i])) return false;
  }
  return true;
}

static bool poly_lc_positive(const rpol_t* p)
{
  if (poly_is_zero(p)) return true;

  rz_t zero;
  rz_init(&zero);
  bool ok = rz_cmp(&p->coeff[p->deg], &zero) > 0;
  rz_clear(&zero);
  return ok;
}

/* random coefficient in [-bound, bound] */
static i64 rnd_coef(i64 bound)
{
  return (i64)(rand() % (int)(2 * bound + 1)) - bound;
}

/* random nonzero polynomial of degree <= max_deg, coefficients in [-b, b] */
static void rnd_poly(rpol_t* p, u64 max_deg, i64 b)
{
  i64 c[16];
  do {
    u64 deg = (u64)(rand() % (int)(max_deg + 1));
    for (u64 i = 0; i <= deg; i++) c[i] = rnd_coef(b);
    CHECK(rpol_set_i64(p, c, deg) == RABIN_SUCCESS, "rnd_poly set");
  } while (poly_is_zero(p));
}

/* multiply every coefficient of src by s, store in dst (dst != src) */
static void poly_scale_by(rpol_t* dst, const rpol_t* src, i64 s)
{
  rz_t scale;
  rpol_t scale_p = {0};
  rz_init(&scale);
  rpol_init(&scale_p);

  CHECK(rz_set_i64(&scale, s) == RABIN_SUCCESS, "scale set");
  CHECK(rpol_set(&scale_p, &scale, 0) == RABIN_SUCCESS, "scale poly");
  CHECK(rpol_mul(dst, &scale_p, src) == RABIN_SUCCESS, "scale mul");

  rpol_clear(&scale_p);
  rz_clear(&scale);
}

/* fresh trimmed copy of src into an uninitialized dst */
static void trim_copy(rpol_t* dst, const rpol_t* src)
{
  rpol_init(dst);
  CHECK(rpol_set(dst, src->coeff, src->deg) == RABIN_SUCCESS, "trim copy set");
  CHECK(rpol_trim(dst) == RABIN_SUCCESS, "trim copy trim");
}

/* true if p is the constant polynomial +-1 */
static bool poly_is_unit(const rpol_t* p)
{
  return p->deg == 0 &&
         (rz_is_eq_i64(&p->coeff[0], 1) || rz_is_eq_i64(&p->coeff[0], -1));
}

/* ------------------------------------------------------------------ */
/* pseudo-division                                                     */
/* ------------------------------------------------------------------ */

/* verify the defining identity delta * a = q * b + r */
static void check_pseudo_identity(const rpol_t* a, const rpol_t* b,
                                  const rpol_t* q, const rpol_t* r)
{
  rpol_t a_t = {0}, b_t = {0}, q_t = {0}, r_t = {0};
  rpol_t scale_p = {0}, lhs = {0}, qb = {0}, rhs = {0};
  rz_t lc, expo, scale;
  trim_copy(&a_t, a);
  trim_copy(&b_t, b);
  trim_copy(&q_t, q);
  trim_copy(&r_t, r);
  rpol_init(&scale_p);
  rpol_init(&lhs);
  rpol_init(&qb);
  rpol_init(&rhs);
  rz_init_multi(&lc, &expo, &scale, NULL);

  bool a_zero = poly_is_zero(&a_t);
  if (!a_zero && a_t.deg >= b_t.deg) {
    CHECK(rz_copy(&lc, &b_t.coeff[b_t.deg]) == RABIN_SUCCESS, "lc copy");
    CHECK(rz_set_u64(&expo, a_t.deg - b_t.deg + 1) == RABIN_SUCCESS,
          "expo set");
    CHECK(rz_pow(&scale, &lc, &expo) == RABIN_SUCCESS, "scale pow");
  } else {
    CHECK(rz_set_u64(&scale, 1) == RABIN_SUCCESS, "scale one");
  }

  CHECK(rpol_set(&scale_p, &scale, 0) == RABIN_SUCCESS, "scale poly set");
  CHECK(rpol_mul(&lhs, &scale_p, &a_t) == RABIN_SUCCESS, "lhs mul");
  CHECK(rpol_mul(&qb, &q_t, &b_t) == RABIN_SUCCESS, "qb mul");
  CHECK(rpol_add(&rhs, &qb, &r_t) == RABIN_SUCCESS, "rhs add");

  CHECK(rpol_equal(&lhs, &rhs), "pseudo identity delta*a == q*b + r");
  if (!rpol_equal(&lhs, &rhs)) {
    printf("  a: ");
    rpol_print(&a_t);
    printf("  b: ");
    rpol_print(&b_t);
    printf("  q: ");
    rpol_print(&q_t);
    printf("  r: ");
    rpol_print(&r_t);
    printf("  delta*a: ");
    rpol_print(&lhs);
    printf("  q*b+r: ");
    rpol_print(&rhs);
  }

  rpol_clear(&a_t);
  rpol_clear(&b_t);
  rpol_clear(&q_t);
  rpol_clear(&r_t);
  rpol_clear(&scale_p);
  rpol_clear(&lhs);
  rpol_clear(&qb);
  rpol_clear(&rhs);
  rz_clear_multi(&lc, &expo, &scale, NULL);
}

static void test_pseudo_div_hand(void)
{
  rpol_t a = {0}, b = {0}, q = {0}, r = {0}, eq = {0};
  rpol_init(&a);
  rpol_init(&b);
  rpol_init(&q);
  rpol_init(&r);
  rpol_init(&eq);

  /* a = 3x^2 + 5x + 2, b = x + 1 -> q = 3x + 2, r = 0 */
  i64 av[] = {2, 5, 3};
  i64 bv[] = {1, 1};
  mk(&a, av, 2);
  mk(&b, bv, 1);
  CHECK(rpol_pseudo_div(&q, &r, &a, &b) == RABIN_SUCCESS, "pseudo 1");
  i64 qv[] = {2, 3};
  mk(&eq, qv, 1);
  CHECK(rpol_equal(&q, &eq), "pseudo q == 3x + 2");
  CHECK(poly_is_zero(&r), "pseudo r == 0");
  check_pseudo_identity(&a, &b, &q, &r);

  /* same a, b = 2x + 1 -> scale 2^2 = 4, q = 6x + 7, r = 1 */
  i64 b2[] = {1, 2};
  mk(&b, b2, 1);
  CHECK(rpol_pseudo_div(&q, &r, &a, &b) == RABIN_SUCCESS, "pseudo 2");
  i64 q2[] = {7, 6};
  mk(&eq, q2, 1);
  CHECK(rpol_equal(&q, &eq), "pseudo q == 6x + 7");
  i64 r2[] = {1};
  mk(&eq, r2, 0);
  CHECK(rpol_equal(&r, &eq), "pseudo r == 1");
  check_pseudo_identity(&a, &b, &q, &r);

  /* negative leading coeff: b = -2x + 1 -> q = -6x - 13, r = 21 */
  i64 b3[] = {1, -2};
  mk(&b, b3, 1);
  CHECK(rpol_pseudo_div(&q, &r, &a, &b) == RABIN_SUCCESS, "pseudo 3");
  i64 q3[] = {-13, -6};
  mk(&eq, q3, 1);
  CHECK(rpol_equal(&q, &eq), "pseudo q == -6x - 13");
  i64 r3[] = {21};
  mk(&eq, r3, 0);
  CHECK(rpol_equal(&r, &eq), "pseudo r == 21");
  check_pseudo_identity(&a, &b, &q, &r);

  /* deg a < deg b -> q = 0, r = a */
  i64 small[] = {5};
  mk(&a, small, 0);
  mk(&b, bv, 1);
  CHECK(rpol_pseudo_div(&q, &r, &a, &b) == RABIN_SUCCESS, "pseudo deg a < b");
  CHECK(poly_is_zero(&q), "q == 0 when deg a < deg b");
  CHECK(rpol_equal(&r, &a), "r == a when deg a < deg b");
  check_pseudo_identity(&a, &b, &q, &r);

  /* a = 0 -> q = 0, r = 0 */
  i64 zero[] = {0};
  mk(&a, zero, 0);
  CHECK(rpol_pseudo_div(&q, &r, &a, &b) == RABIN_SUCCESS, "pseudo a == 0");
  CHECK(poly_is_zero(&q) && poly_is_zero(&r), "0 divided leaves 0, 0");

  /* b = 0 -> DIV_BY_ZERO */
  mk(&a, av, 2);
  mk(&b, zero, 0);
  CHECK(rpol_pseudo_div(&q, &r, &a, &b) == RABIN_ERR_DIV_BY_ZERO,
        "division by zero poly");

  /* NULL args */
  CHECK(rpol_pseudo_div(NULL, &r, &a, &b) == RABIN_ERR_NULL_PTR,
        "pseudo NULL q");
  CHECK(rpol_pseudo_div(&q, NULL, &a, &b) == RABIN_ERR_NULL_PTR,
        "pseudo NULL r");
  CHECK(rpol_pseudo_div(&q, &r, NULL, &b) == RABIN_ERR_NULL_PTR,
        "pseudo NULL a");
  CHECK(rpol_pseudo_div(&q, &r, &a, NULL) == RABIN_ERR_NULL_PTR,
        "pseudo NULL b");

  rpol_clear(&a);
  rpol_clear(&b);
  rpol_clear(&q);
  rpol_clear(&r);
  rpol_clear(&eq);
}

static void test_pseudo_div_random(int rounds)
{
  rpol_t a = {0}, b = {0}, q = {0}, r = {0};
  rpol_init(&a);
  rpol_init(&b);
  rpol_init(&q);
  rpol_init(&r);

  for (int t = 0; t < rounds; t++) {
    rnd_poly(&a, 5, 9);
    rnd_poly(&b, 3, 9);

    CHECK(rpol_pseudo_div(&q, &r, &a, &b) == RABIN_SUCCESS, "pseudo random");

    /* remainder must shrink: r == 0 or deg r < deg b */
    CHECK(poly_is_zero(&r) || r.deg < b.deg, "pseudo remainder degree");

    /* quotient degree bound */
    if (a.deg >= b.deg) {
      CHECK(poly_is_zero(&q) || q.deg <= a.deg - b.deg,
            "pseudo quotient degree");
    } else {
      rpol_t a_t = {0};
      trim_copy(&a_t, &a);
      CHECK(poly_is_zero(&q), "q == 0 when deg a < deg b");
      CHECK(rpol_equal(&r, &a_t), "r == a when deg a < deg b (trimmed)");
      rpol_clear(&a_t);
    }

    check_pseudo_identity(&a, &b, &q, &r);
  }

  rpol_clear(&a);
  rpol_clear(&b);
  rpol_clear(&q);
  rpol_clear(&r);
}

/* ------------------------------------------------------------------ */
/* exact division                                                      */
/* ------------------------------------------------------------------ */

static void test_div_exact_hand(void)
{
  rpol_t a = {0}, b = {0}, q = {0}, eq = {0};
  rpol_init(&a);
  rpol_init(&b);
  rpol_init(&q);
  rpol_init(&eq);

  /* (3x^2 + 5x + 2) / (x + 1) = 3x + 2 */
  i64 av[] = {2, 5, 3};
  i64 bv[] = {1, 1};
  mk(&a, av, 2);
  mk(&b, bv, 1);
  CHECK(rpol_div_exact(&q, &a, &b) == RABIN_SUCCESS, "div exact 1");
  i64 qv[] = {2, 3};
  mk(&eq, qv, 1);
  CHECK(rpol_equal(&q, &eq), "(3x^2+5x+2)/(x+1) == 3x+2");

  /* (2x^2 - 5x - 3) / (2x + 1) = x - 3 */
  i64 a2[] = {-3, -5, 2};
  i64 b2[] = {1, 2};
  mk(&a, a2, 2);
  mk(&b, b2, 1);
  CHECK(rpol_div_exact(&q, &a, &b) == RABIN_SUCCESS, "div exact 2");
  i64 q2[] = {-3, 1};
  mk(&eq, q2, 1);
  CHECK(rpol_equal(&q, &eq), "(2x^2-5x-3)/(2x+1) == x-3");

  /* a = 0 -> q = 0 for any nonzero b */
  i64 zero[] = {0};
  mk(&a, zero, 0);
  CHECK(rpol_div_exact(&q, &a, &b) == RABIN_SUCCESS, "div exact a=0");
  CHECK(poly_is_zero(&q), "0 / b == 0");

  /* not divisible: lc not divisible -> INVALID_ARG, q untouched */
  i64 sent[] = {9, 9, 9};
  mk(&q, sent, 2);
  i64 a3[] = {1, 0, 1}; /* x^2 + 1 */
  i64 b3[] = {1, 2};    /* 2x + 1 */
  mk(&a, a3, 2);
  mk(&b, b3, 1);
  CHECK(rpol_div_exact(&q, &a, &b) == RABIN_ERR_INVALID_ARG,
        "non-divisible -> INVALID_ARG");
  mk(&eq, sent, 2);
  CHECK(rpol_equal(&q, &eq), "q untouched after failure");

  /* divisible coefficients but nonzero remainder: (2x^2+x+1)/(2x+1) */
  i64 a4[] = {1, 1, 2};
  mk(&a, a4, 2);
  CHECK(rpol_div_exact(&q, &a, &b) == RABIN_ERR_INVALID_ARG,
        "nonzero remainder -> INVALID_ARG");

  /* deg a < deg b with a != 0 is never exactly divisible */
  i64 a5[] = {7};
  mk(&a, a5, 0);
  CHECK(rpol_div_exact(&q, &a, &b) == RABIN_ERR_INVALID_ARG,
        "deg a < deg b -> INVALID_ARG");

  /* b = 0 */
  i64 zb[] = {0};
  mk(&b, zb, 0);
  CHECK(rpol_div_exact(&q, &a, &b) == RABIN_ERR_DIV_BY_ZERO,
        "div exact by zero poly");

  /* NULL args */
  mk(&b, b3, 1);
  CHECK(rpol_div_exact(NULL, &a, &b) == RABIN_ERR_NULL_PTR, "div exact NULL q");
  CHECK(rpol_div_exact(&q, NULL, &b) == RABIN_ERR_NULL_PTR, "div exact NULL a");
  CHECK(rpol_div_exact(&q, &a, NULL) == RABIN_ERR_NULL_PTR, "div exact NULL b");

  rpol_clear(&a);
  rpol_clear(&b);
  rpol_clear(&q);
  rpol_clear(&eq);
}

static void test_div_exact_roundtrip(int rounds)
{
  rpol_t f = {0}, g = {0}, neg_g = {0}, neg_f = {0}, a = {0}, q = {0}, s = {0},
         bad = {0};
  rpol_init(&f);
  rpol_init(&g);
  rpol_init(&neg_g);
  rpol_init(&neg_f);
  rpol_init(&a);
  rpol_init(&q);
  rpol_init(&s);
  rpol_init(&bad);
  i64 one[] = {1};
  mk(&s, one, 0);

  for (int t = 0; t < rounds; t++) {
    rnd_poly(&f, 4, 9);
    rnd_poly(&g, 3, 9);

    CHECK(rpol_mul(&a, &f, &g) == RABIN_SUCCESS, "roundtrip mul");
    CHECK(rpol_div_exact(&q, &a, &g) == RABIN_SUCCESS, "(f*g)/g succeeds");
    rpol_t f_t = {0};
    trim_copy(&f_t, &f);
    CHECK(rpol_equal(&q, &f_t), "(f*g)/g == trimmed f");
    if (!rpol_equal(&q, &f_t)) {
      printf("  f: ");
      rpol_print(&f_t);
      printf("  g: ");
      rpol_print(&g);
      printf("  q: ");
      rpol_print(&q);
    }
    rpol_clear(&f_t);

    /* negative divisor must work too: (f*g)/(-g) == -f */
    poly_scale_by(&neg_g, &g, -1);
    poly_scale_by(&neg_f, &f, -1);
    CHECK(rpol_div_exact(&q, &a, &neg_g) == RABIN_SUCCESS, "(f*g)/(-g) ok");
    CHECK(rpol_equal(&q, &neg_f), "(f*g)/(-g) == -f");

    /* f*g + 1 is never exactly divisible by g unless g is a unit */
    if (!poly_is_unit(&g)) {
      CHECK(rpol_add(&bad, &a, &s) == RABIN_SUCCESS, "bad poly build");
      CHECK(rpol_div_exact(&q, &bad, &g) == RABIN_ERR_INVALID_ARG,
            "f*g + 1 not divisible by g");
    }
  }

  rpol_clear(&f);
  rpol_clear(&g);
  rpol_clear(&neg_g);
  rpol_clear(&neg_f);
  rpol_clear(&a);
  rpol_clear(&q);
  rpol_clear(&s);
  rpol_clear(&bad);
}

/* ------------------------------------------------------------------ */
/* content and primitive part                                          */
/* ------------------------------------------------------------------ */

static void test_content_pp(void)
{
  rpol_t p = {0}, pp = {0}, eq = {0};
  rz_t c, zero;
  rpol_init(&p);
  rpol_init(&pp);
  rpol_init(&eq);
  rz_init_multi(&c, &zero, NULL);
  rz_set_u64(&zero, 0);

  /* content(6x^3 + 9x^2 - 15x + 3) = 3 */
  i64 v1[] = {3, -15, 9, 6};
  mk(&p, v1, 3);
  CHECK(rpol_content(&c, &p) == RABIN_SUCCESS, "content 1");
  CHECK(rz_is_eq_i64(&c, 3), "content == 3");

  /* content(0) = 0 */
  i64 zero_c[] = {0};
  mk(&p, zero_c, 0);
  CHECK(rpol_content(&c, &p) == RABIN_SUCCESS, "content 0");
  CHECK(rz_cmp(&c, &zero) == 0, "content(0) == 0");

  /* content(6x - 4) = 2 */
  i64 v2[] = {-4, 6};
  mk(&p, v2, 1);
  CHECK(rpol_content(&c, &p) == RABIN_SUCCESS, "content 2");
  CHECK(rz_is_eq_i64(&c, 2), "content(6x-4) == 2");

  /* content(3x + 2) = 1 */
  i64 v3[] = {2, 3};
  mk(&p, v3, 1);
  CHECK(rpol_content(&c, &p) == RABIN_SUCCESS, "content 3");
  CHECK(rz_is_eq_i64(&c, 1), "content(3x+2) == 1");

  /* pp(6x^3 + 9x^2 - 15x + 3) = 2x^3 + 3x^2 - 5x + 1 */
  mk(&p, v1, 3);
  CHECK(rpol_primitive_part(&pp, &p) == RABIN_SUCCESS, "pp 1");
  i64 e1[] = {1, -5, 3, 2};
  mk(&eq, e1, 3);
  CHECK(rpol_equal(&pp, &eq), "pp == 2x^3+3x^2-5x+1");

  /* pp(-4x^2 + 6x - 2) = 2x^2 - 3x + 1 (leading coeff flipped to > 0) */
  i64 v4[] = {-2, 6, -4};
  mk(&p, v4, 2);
  CHECK(rpol_primitive_part(&pp, &p) == RABIN_SUCCESS, "pp 2");
  i64 e2[] = {1, -3, 2};
  mk(&eq, e2, 2);
  CHECK(rpol_equal(&pp, &eq), "pp == 2x^2-3x+1");

  /* pp(0) = 0 */
  mk(&p, zero_c, 0);
  CHECK(rpol_primitive_part(&pp, &p) == RABIN_SUCCESS, "pp 0");
  CHECK(poly_is_zero(&pp), "pp(0) == 0");

  rpol_clear(&eq);
  rpol_clear(&p);
  rpol_clear(&pp);
  rz_clear_multi(&c, &zero, NULL);
}

static void test_pp_invariants(int rounds)
{
  rpol_t p = {0}, src = {0}, pp = {0}, scaled = {0}, negp = {0}, sp = {0};
  rz_t c;
  rpol_init(&p);
  rpol_init(&src);
  rpol_init(&pp);
  rpol_init(&scaled);
  rpol_init(&negp);
  rpol_init(&sp);
  rz_init(&c);

  for (int t = 0; t < rounds; t++) {
    rnd_poly(&src, 4, 7);
    /* boost chance of nontrivial content: double all coefficients */
    poly_scale_by(&p, &src, 2);

    CHECK(rpol_content(&c, &p) == RABIN_SUCCESS, "inv content");
    CHECK(rpol_primitive_part(&pp, &p) == RABIN_SUCCESS, "inv pp");

    /* pp is primitive with positive leading coefficient */
    rz_t c2;
    rz_init(&c2);
    CHECK(rpol_content(&c2, &pp) == RABIN_SUCCESS, "inv content pp");
    CHECK(rz_is_eq_i64(&c2, 1), "pp is primitive");
    rz_clear(&c2);

    CHECK(poly_lc_positive(&pp), "pp leading coeff > 0");

    /* p = +-content(p) * pp(p) */
    CHECK(rpol_set(&sp, &c, 0) == RABIN_SUCCESS, "inv scale set");
    CHECK(rpol_mul(&scaled, &sp, &pp) == RABIN_SUCCESS, "inv scale mul");
    poly_scale_by(&negp, &p, -1);
    bool direct = rpol_equal(&scaled, &p);
    bool flipped = rpol_equal(&scaled, &negp);
    CHECK(direct || flipped, "content * pp == +-p");
    if (!direct && !flipped) {
      printf("  case: ");
      rpol_print(&p);
      printf("  c=");
      rz_print(&c);
      printf("  c*pp=");
      rpol_print(&scaled);
    }
  }

  rpol_clear(&p);
  rpol_clear(&src);
  rpol_clear(&pp);
  rpol_clear(&scaled);
  rpol_clear(&negp);
  rpol_clear(&sp);
  rz_clear(&c);
}

/* ------------------------------------------------------------------ */
/* gcd: independent small __int128 pseudo-remainder sequence oracle    */
/* ------------------------------------------------------------------ */

#define IP_CAP 24

typedef struct {
  __int128 c[IP_CAP];
} ip_t;

static int ip_deg(const ip_t* p)
{
  int d = IP_CAP - 2;
  while (d > 0 && p->c[d] == 0) d--;
  return d;
}

static bool ip_is_zero(const ip_t* p) { return ip_deg(p) == 0 && p->c[0] == 0; }

static __int128 ip_gcd128(__int128 a, __int128 b)
{
  if (a < 0) a = -a;
  if (b < 0) b = -b;
  while (b != 0) {
    __int128 t = a % b;
    a = b;
    b = t;
  }
  return a;
}

static __int128 ip_content(const ip_t* p)
{
  __int128 g = 0;
  for (int i = 0; i <= ip_deg(p); i++) g = ip_gcd128(g, p->c[i]);
  return g;
}

/* divide out the content and normalize lc > 0; zero poly maps to itself */
static ip_t ip_make_prim(ip_t p)
{
  if (ip_is_zero(&p)) return p;

  __int128 ct = ip_content(&p);
  if (ct > 1) {
    for (int i = 0; i < IP_CAP; i++) p.c[i] /= ct;
  }
  if (p.c[ip_deg(&p)] < 0) {
    for (int i = 0; i < IP_CAP; i++) p.c[i] = -p.c[i];
  }
  return p;
}

static ip_t ip_sign_normalize(ip_t p)
{
  if (!ip_is_zero(&p) && p.c[ip_deg(&p)] < 0) {
    for (int i = 0; i < IP_CAP; i++) p.c[i] = -p.c[i];
  }
  return p;
}

/* pseudo-remainder: scale a by lc(b)^(deg a - deg b + 1), then divide */
static ip_t ip_prem(ip_t a, ip_t b)
{
  int db = ip_deg(&b);
  int da = ip_deg(&a);
  if (da < db) return a;

  __int128 lc = b.c[db];
  for (int e = 0; e <= da - db; e++) {
    for (int i = 0; i <= da; i++) a.c[i] *= lc;
  }

  /* ip_deg() reports 0 for the zero polynomial, so a plain da >= db test
   * would spin forever once the remainder is exhausted with db == 0 */
  while (!ip_is_zero(&a) && (da = ip_deg(&a)) >= db) {
    __int128 qk = a.c[da] / lc; /* exact by pseudo-division theory */
    int k = da - db;
    for (int j = 0; j <= db; j++) a.c[k + j] -= qk * b.c[j];
  }
  return a;
}

/* reduced pseudo-remainder sequence gcd over Z[x] (oracle) */
static ip_t ip_gcd(ip_t a, ip_t b)
{
  if (ip_is_zero(&a) && ip_is_zero(&b)) return (ip_t){{0}};
  if (ip_is_zero(&a)) return ip_sign_normalize(b);
  if (ip_is_zero(&b)) return ip_sign_normalize(a);

  __int128 ct = ip_gcd128(ip_content(&a), ip_content(&b));
  a = ip_make_prim(a);
  b = ip_make_prim(b);

  if (ip_deg(&a) < ip_deg(&b)) {
    ip_t t = a;
    a = b;
    b = t;
  }

  while (!ip_is_zero(&b)) {
    ip_t r = ip_prem(a, b);
    a = b; /* a is the last nonzero remainder once the loop exits */
    if (ip_is_zero(&r)) break;
    r = ip_make_prim(r);
    b = r;
  }

  if (ct > 1) {
    for (int i = 0; i < IP_CAP; i++) a.c[i] *= ct;
  }
  return a;
}

static void ip_to_rpol(const ip_t* p, rpol_t* out)
{
  i64 vals[IP_CAP];
  int deg = ip_deg(p);
  for (int i = 0; i <= deg; i++) {
    __int128 v = p->c[i];
    CHECK(v >= (__int128)INT64_MIN && v <= (__int128)INT64_MAX,
          "oracle value fits i64");
    vals[i] = (i64)v;
  }
  CHECK(rpol_set_i64(out, vals, (u64)deg) == RABIN_SUCCESS, "oracle set");
}

static void test_gcd_hand(void)
{
  rpol_t a = {0}, b = {0}, g = {0}, eq = {0}, neg_a = {0}, g2 = {0};
  rpol_init(&a);
  rpol_init(&b);
  rpol_init(&g);
  rpol_init(&eq);
  rpol_init(&neg_a);
  rpol_init(&g2);

  /* gcd(6x^2+7x+2, 6x^2-x-2) = gcd((2x+1)(3x+2), (2x+1)(3x-2)) = 2x+1 */
  i64 av[] = {2, 7, 6};
  i64 bv[] = {-2, -1, 6};
  mk(&a, av, 2);
  mk(&b, bv, 2);
  CHECK(rpol_gcd(&g, &a, &b) == RABIN_SUCCESS, "gcd hand 1");
  i64 e1[] = {1, 2};
  mk(&eq, e1, 1);
  CHECK(rpol_equal(&g, &eq), "gcd == 2x+1");

  /* gcd(12x^3, 8x^2) = 4x^2 */
  i64 a2[] = {0, 0, 0, 12};
  i64 b2[] = {0, 0, 8};
  mk(&a, a2, 3);
  mk(&b, b2, 2);
  CHECK(rpol_gcd(&g, &a, &b) == RABIN_SUCCESS, "gcd hand 2");
  i64 e2[] = {0, 0, 4};
  mk(&eq, e2, 2);
  CHECK(rpol_equal(&g, &eq), "gcd == 4x^2");

  /* gcd(0, 0) = 0 */
  i64 z[] = {0};
  mk(&a, z, 0);
  mk(&b, z, 0);
  CHECK(rpol_gcd(&g, &a, &b) == RABIN_SUCCESS, "gcd 0,0");
  CHECK(poly_is_zero(&g), "gcd(0,0) == 0");

  /* gcd(0, -5x+10) = 5x-10 (lc normalized > 0) */
  i64 nz[] = {10, -5};
  mk(&b, nz, 1);
  CHECK(rpol_gcd(&g, &a, &b) == RABIN_SUCCESS, "gcd 0,b");
  i64 e3[] = {-10, 5};
  mk(&eq, e3, 1);
  CHECK(rpol_equal(&g, &eq), "gcd(0,-5x+10) == 5x-10");

  /* gcd(x^2 - 1, x^3 - 1) = x - 1 */
  i64 a4[] = {-1, 0, 1};
  i64 b4[] = {-1, 0, 0, 1};
  mk(&a, a4, 2);
  mk(&b, b4, 3);
  CHECK(rpol_gcd(&g, &a, &b) == RABIN_SUCCESS, "gcd hand 3");
  i64 e4[] = {-1, 1};
  mk(&eq, e4, 1);
  CHECK(rpol_equal(&g, &eq), "gcd(x^2-1, x^3-1) == x-1");

  /* coprime constants: gcd(4, 6) = 2 as constant poly */
  i64 a5[] = {4};
  i64 b5[] = {6};
  mk(&a, a5, 0);
  mk(&b, b5, 0);
  CHECK(rpol_gcd(&g, &a, &b) == RABIN_SUCCESS, "gcd const, const");
  CHECK(rz_is_eq_i64(&g.coeff[0], 2), "gcd(4,6) == 2");

  /* sign and order invariance for the first case */
  mk(&a, av, 2);
  mk(&b, bv, 2);
  poly_scale_by(&neg_a, &a, -1);
  CHECK(rpol_gcd(&g, &a, &b) == RABIN_SUCCESS, "gcd(a,b)");
  CHECK(rpol_gcd(&g2, &b, &a) == RABIN_SUCCESS, "gcd(b,a)");
  CHECK(rpol_equal(&g, &g2), "gcd symmetric");
  CHECK(rpol_gcd(&g2, &neg_a, &b) == RABIN_SUCCESS, "gcd(-a,b)");
  CHECK(rpol_equal(&g, &g2), "gcd sign invariant");

  /* NULL args */
  CHECK(rpol_gcd(NULL, &a, &b) == RABIN_ERR_NULL_PTR, "gcd NULL g");
  CHECK(rpol_gcd(&g, NULL, &b) == RABIN_ERR_NULL_PTR, "gcd NULL a");
  CHECK(rpol_gcd(&g, &a, NULL) == RABIN_ERR_NULL_PTR, "gcd NULL b");

  rpol_clear(&a);
  rpol_clear(&b);
  rpol_clear(&g);
  rpol_clear(&eq);
  rpol_clear(&neg_a);
  rpol_clear(&g2);
}

static void test_gcd_oracle_small(int rounds)
{
  /* tiny coefficients keep the __int128 oracle comfortably in range */
  rpol_t a = {0}, b = {0}, g = {0}, expect = {0}, q = {0};
  rpol_init(&a);
  rpol_init(&b);
  rpol_init(&g);
  rpol_init(&expect);
  rpol_init(&q);

  for (int t = 0; t < rounds; t++) {
    ip_t ia = {{0}}, ib = {{0}};
    u64 da = (u64)(rand() % 4), db = (u64)(rand() % 4);
    for (u64 i = 0; i <= da; i++) ia.c[i] = rnd_coef(3);
    for (u64 i = 0; i <= db; i++) ib.c[i] = rnd_coef(3);

    ip_to_rpol(&ia, &a);
    ip_to_rpol(&ib, &b);

    CHECK(rpol_gcd(&g, &a, &b) == RABIN_SUCCESS, "gcd oracle run");

    ip_t ig = ip_gcd(ia, ib);
    ip_to_rpol(&ig, &expect);
    CHECK(rpol_equal(&g, &expect), "gcd matches naive PRS oracle");
    if (!rpol_equal(&g, &expect)) {
      printf("  a: ");
      rpol_print(&a);
      printf("  b: ");
      rpol_print(&b);
      printf("  g: ");
      rpol_print(&g);
      printf("  expected: ");
      rpol_print(&expect);
    }

    /* the gcd must divide both inputs exactly */
    if (!poly_is_zero(&g)) {
      if (!poly_is_zero(&a))
        CHECK(rpol_div_exact(&q, &a, &g) == RABIN_SUCCESS, "gcd divides a");
      if (!poly_is_zero(&b))
        CHECK(rpol_div_exact(&q, &b, &g) == RABIN_SUCCESS, "gcd divides b");
    }
  }

  rpol_clear(&a);
  rpol_clear(&b);
  rpol_clear(&g);
  rpol_clear(&expect);
  rpol_clear(&q);
}

static void test_gcd_constructed(int rounds)
{
  /* a = f*u*k1, b = f*v*k2: the gcd must be divisible by f and divide
   * both inputs */
  rpol_t f = {0}, u = {0}, v = {0}, a = {0}, b = {0}, g = {0}, t1 = {0},
         t2 = {0}, q = {0};
  rpol_init(&f);
  rpol_init(&u);
  rpol_init(&v);
  rpol_init(&a);
  rpol_init(&b);
  rpol_init(&g);
  rpol_init(&t1);
  rpol_init(&t2);
  rpol_init(&q);

  for (int t = 0; t < rounds; t++) {
    rnd_poly(&f, 3, 4);
    rnd_poly(&u, 2, 4);
    rnd_poly(&v, 2, 4);

    int k1 = 1 + rand() % 3, k2 = 1 + rand() % 3;

    CHECK(rpol_mul(&t1, &f, &u) == RABIN_SUCCESS, "construct t1");
    poly_scale_by(&a, &t1, k1);
    CHECK(rpol_mul(&t2, &f, &v) == RABIN_SUCCESS, "construct t2");
    poly_scale_by(&b, &t2, k2);

    CHECK(rpol_gcd(&g, &a, &b) == RABIN_SUCCESS, "gcd constructed");
    CHECK(!poly_is_zero(&g), "constructed gcd nonzero");
    CHECK(poly_lc_positive(&g), "constructed gcd lc > 0");

    /* the gcd divides both inputs */
    CHECK(rpol_div_exact(&q, &a, &g) == RABIN_SUCCESS, "constructed g | a");
    CHECK(rpol_div_exact(&q, &b, &g) == RABIN_SUCCESS, "constructed g | b");

    /* f divides the gcd */
    CHECK(rpol_div_exact(&q, &g, &f) == RABIN_SUCCESS, "f | constructed g");
    if (rpol_div_exact(&q, &g, &f) != RABIN_SUCCESS) {
      printf("  f: ");
      rpol_print(&f);
      printf("  a: ");
      rpol_print(&a);
      printf("  b: ");
      rpol_print(&b);
      printf("  g: ");
      rpol_print(&g);
    }
  }

  rpol_clear(&f);
  rpol_clear(&u);
  rpol_clear(&v);
  rpol_clear(&a);
  rpol_clear(&b);
  rpol_clear(&g);
  rpol_clear(&t1);
  rpol_clear(&t2);
  rpol_clear(&q);
}

#ifdef HAVE_FLINT
static i64 rpol_coef_si(const rpol_t* p, u64 i)
{
  if (i >= p->size) return 0;
  return (i64)p->coeff[i].limbs[0] * (p->coeff[i].is_neg ? -1 : 1);
}

static void test_gcd_flint(int rounds)
{
  rpol_t a = {0}, b = {0}, g = {0}, expect = {0};
  rpol_init(&a);
  rpol_init(&b);
  rpol_init(&g);
  rpol_init(&expect);

  fmpz_poly_t fa, fb, fg;
  fmpz_poly_init(fa);
  fmpz_poly_init(fb);
  fmpz_poly_init(fg);

  int checked = 0;
  for (int t = 0; t < rounds; t++) {
    rnd_poly(&a, 5, 9);
    rnd_poly(&b, 5, 9);

    fmpz_poly_zero(fa);
    fmpz_poly_zero(fb);
    for (u64 i = 0; i <= a.deg; i++)
      fmpz_poly_set_coeff_si(fa, (slong)i, rpol_coef_si(&a, i));
    for (u64 i = 0; i <= b.deg; i++)
      fmpz_poly_set_coeff_si(fb, (slong)i, rpol_coef_si(&b, i));

    fmpz_poly_gcd(fg, fa, fb);

    /* FLINT normalizes to positive leading coefficient, same as rpol_gcd */
    slong fdeg = fmpz_poly_degree(fg);
    CHECK(rpol_gcd(&g, &a, &b) == RABIN_SUCCESS, "flint gcd run");
    if (fdeg < 0) {
      CHECK(poly_is_zero(&g), "flint: zero gcd");
      continue;
    }

    i64 ev[64];
    for (slong i = 0; i <= fdeg; i++) {
      fmpz_t fc;
      fmpz_init(fc);
      fmpz_poly_get_coeff_fmpz(fc, fg, i);
      ev[i] = fmpz_get_si(fc);
      fmpz_clear(fc);
    }
    mk(&expect, ev, (u64)fdeg);
    CHECK(rpol_equal(&g, &expect), "gcd matches FLINT fmpz_poly_gcd");
    if (!rpol_equal(&g, &expect)) {
      printf("  a: ");
      rpol_print(&a);
      printf("  b: ");
      rpol_print(&b);
      printf("  rabin gcd: ");
      rpol_print(&g);
      printf("  flint gcd: ");
      rpol_print(&expect);
    }
    checked++;
  }

  printf("test_poly_div: FLINT cross-check covered %d random pairs\n", checked);

  fmpz_poly_clear(fa);
  fmpz_poly_clear(fb);
  fmpz_poly_clear(fg);
  rpol_clear(&a);
  rpol_clear(&b);
  rpol_clear(&g);
  rpol_clear(&expect);
}
#endif

/* ------------------------------------------------------------------ */
/* evaluation and derivation                                           */
/* ------------------------------------------------------------------ */

/* independent O(d^2) evaluation: sum c_i * x^i, no Horner */
static void eval_direct(rz_t* out, const rpol_t* p, const rz_t* x)
{
  rz_t acc, term, pow_x;
  rz_init_multi(&acc, &term, &pow_x, NULL);
  rz_set_u64(&acc, 0);
  rz_set_u64(&pow_x, 1);

  for (u64 i = 0; i <= p->deg; i++) {
    CHECK(rz_mul(&term, &p->coeff[i], &pow_x) == RABIN_SUCCESS, "direct term");
    CHECK(rz_add(&acc, &acc, &term) == RABIN_SUCCESS, "direct add");
    if (i < p->deg) CHECK(rz_mul(&pow_x, &pow_x, x) == RABIN_SUCCESS, "pow");
  }

  CHECK(rz_copy(out, &acc) == RABIN_SUCCESS, "direct copy");
  rz_clear_multi(&acc, &term, &pow_x, NULL);
}

static void test_eval_hand(void)
{
  rpol_t p = {0};
  rz_t x, val, direct;
  rpol_init(&p);
  rz_init_multi(&x, &val, &direct, NULL);

  /* p = 2x^3 - 3x + 5 */
  i64 pv[] = {5, -3, 0, 2};
  mk(&p, pv, 3);

  struct {
    i64 x;
    i64 expect;
  } cases[] = {{0, 5}, {1, 4}, {2, 15}, {-2, -5}, {10, 1975}, {-7, -660}};

  for (size_t t = 0; t < sizeof(cases) / sizeof(cases[0]); t++) {
    CHECK(rz_set_i64(&x, cases[t].x) == RABIN_SUCCESS, "eval x set");
    CHECK(rpol_eval(&val, &p, &x) == RABIN_SUCCESS, "eval hand");
    CHECK(rz_is_eq_i64(&val, cases[t].expect), "eval hand value");
    if (!rz_is_eq_i64(&val, cases[t].expect)) {
      printf("  x=%lld expect=%lld\n", (long long)cases[t].x,
             (long long)cases[t].expect);
      printf("  got: ");
      rz_print(&val);
    }

    /* Horner vs the independent direct-sum oracle */
    eval_direct(&direct, &p, &x);
    CHECK(rz_cmp(&val, &direct) == 0, "Horner == direct sum");
  }

  /* large evaluation point: x = 2^40 + 1 */
  CHECK(rz_lshift(&x, &RZ_ONE, 40) == RABIN_SUCCESS, "big x shift");
  CHECK(rz_add(&x, &x, &RZ_ONE) == RABIN_SUCCESS, "big x inc");
  CHECK(rpol_eval(&val, &p, &x) == RABIN_SUCCESS, "eval big x");
  eval_direct(&direct, &p, &x);
  CHECK(rz_cmp(&val, &direct) == 0, "Horner == direct sum (big x)");

  /* zero polynomial evaluates to 0 everywhere */
  i64 z[] = {0};
  mk(&p, z, 0);
  CHECK(rz_set_i64(&x, 12345) == RABIN_SUCCESS, "zero eval x");
  CHECK(rpol_eval(&val, &p, &x) == RABIN_SUCCESS, "eval zero poly");
  CHECK(rz_is_eq_i64(&val, 0), "zero poly evaluates to 0");

  /* NULL args */
  CHECK(rpol_eval(NULL, &p, &x) == RABIN_ERR_NULL_PTR, "eval NULL r");
  CHECK(rpol_eval(&val, NULL, &x) == RABIN_ERR_NULL_PTR, "eval NULL p");
  CHECK(rpol_eval(&val, &p, NULL) == RABIN_ERR_NULL_PTR, "eval NULL x");

  rpol_clear(&p);
  rz_clear_multi(&x, &val, &direct, NULL);
}

static void test_eval_random(int rounds)
{
  rpol_t p = {0}, q = {0}, prod = {0};
  rz_t x, vp, vq, vpq, expect, direct;
  rpol_init(&p);
  rpol_init(&q);
  rpol_init(&prod);
  rz_init_multi(&x, &vp, &vq, &vpq, &expect, &direct, NULL);

  for (int t = 0; t < rounds; t++) {
    rnd_poly(&p, 5, 9);
    rnd_poly(&q, 4, 9);
    CHECK(rz_set_i64(&x, rnd_coef(11)) == RABIN_SUCCESS, "rand eval x");

    CHECK(rpol_eval(&vp, &p, &x) == RABIN_SUCCESS, "rand eval p");
    eval_direct(&direct, &p, &x);
    CHECK(rz_cmp(&vp, &direct) == 0, "rand Horner == direct sum");

    /* product rule: (p*q)(x) == p(x) * q(x) */
    CHECK(rpol_mul(&prod, &p, &q) == RABIN_SUCCESS, "rand eval prod");
    CHECK(rpol_eval(&vpq, &prod, &x) == RABIN_SUCCESS, "rand eval prod x");
    CHECK(rpol_eval(&vq, &q, &x) == RABIN_SUCCESS, "rand eval q");
    CHECK(rz_mul(&expect, &vp, &vq) == RABIN_SUCCESS, "rand eval vp*vq");
    CHECK(rz_cmp(&vpq, &expect) == 0, "eval product rule");
  }

  rpol_clear(&p);
  rpol_clear(&q);
  rpol_clear(&prod);
  rz_clear_multi(&x, &vp, &vq, &vpq, &expect, &direct, NULL);
}

static void test_derive_hand(void)
{
  rpol_t p = {0}, dp = {0}, eq = {0};
  rpol_init(&p);
  rpol_init(&dp);
  rpol_init(&eq);

  /* d/dx (2x^3 - 3x + 5) = 6x^2 - 3 */
  i64 pv[] = {5, -3, 0, 2};
  mk(&p, pv, 3);
  CHECK(rpol_derive(&dp, &p) == RABIN_SUCCESS, "derive 1");
  i64 e1[] = {-3, 0, 6};
  mk(&eq, e1, 2);
  CHECK(rpol_equal(&dp, &eq), "d/dx(2x^3-3x+5) == 6x^2-3");

  /* d/dx (x^5) = 5x^4 */
  i64 pv2[] = {0, 0, 0, 0, 0, 1};
  mk(&p, pv2, 5);
  CHECK(rpol_derive(&dp, &p) == RABIN_SUCCESS, "derive 2");
  i64 e2[] = {0, 0, 0, 0, 5};
  mk(&eq, e2, 4);
  CHECK(rpol_equal(&dp, &eq), "d/dx(x^5) == 5x^4");

  /* d/dx (3x^4 + 2x^2 + 1) = 12x^3 + 4x */
  i64 pv3[] = {1, 0, 2, 0, 3};
  mk(&p, pv3, 4);
  CHECK(rpol_derive(&dp, &p) == RABIN_SUCCESS, "derive 3");
  i64 e3[] = {0, 4, 0, 12};
  mk(&eq, e3, 3);
  CHECK(rpol_equal(&dp, &eq), "d/dx(3x^4+2x^2+1) == 12x^3+4x");

  /* derivative of constants and of 0 is 0 */
  i64 c[] = {7};
  mk(&p, c, 0);
  CHECK(rpol_derive(&dp, &p) == RABIN_SUCCESS, "derive const");
  CHECK(poly_is_zero(&dp), "d/dx(7) == 0");

  i64 z[] = {0};
  mk(&p, z, 0);
  CHECK(rpol_derive(&dp, &p) == RABIN_SUCCESS, "derive 0");
  CHECK(poly_is_zero(&dp), "d/dx(0) == 0");

  /* NULL args */
  CHECK(rpol_derive(NULL, &p) == RABIN_ERR_NULL_PTR, "derive NULL dp");
  CHECK(rpol_derive(&dp, NULL) == RABIN_ERR_NULL_PTR, "derive NULL p");

  rpol_clear(&p);
  rpol_clear(&dp);
  rpol_clear(&eq);
}

static void test_derive_random(int rounds)
{
  /* product rule: (p*q)' == p'*q + p*q' */
  rpol_t p = {0}, q = {0}, pq = {0}, dpq = {0}, dp = {0}, dq = {0}, t1 = {0},
         t2 = {0}, rhs = {0};
  rpol_init(&p);
  rpol_init(&q);
  rpol_init(&pq);
  rpol_init(&dpq);
  rpol_init(&dp);
  rpol_init(&dq);
  rpol_init(&t1);
  rpol_init(&t2);
  rpol_init(&rhs);

  for (int t = 0; t < rounds; t++) {
    rnd_poly(&p, 4, 7);
    rnd_poly(&q, 4, 7);

    CHECK(rpol_mul(&pq, &p, &q) == RABIN_SUCCESS, "prod rule mul");
    CHECK(rpol_derive(&dpq, &pq) == RABIN_SUCCESS, "prod rule derive");
    CHECK(rpol_derive(&dp, &p) == RABIN_SUCCESS, "prod rule derive p");
    CHECK(rpol_derive(&dq, &q) == RABIN_SUCCESS, "prod rule derive q");

    CHECK(rpol_mul(&t1, &dp, &q) == RABIN_SUCCESS, "prod rule p'q");
    CHECK(rpol_mul(&t2, &p, &dq) == RABIN_SUCCESS, "prod rule pq'");
    CHECK(rpol_add(&rhs, &t1, &t2) == RABIN_SUCCESS, "prod rule add");

    CHECK(rpol_equal(&dpq, &rhs), "(pq)' == p'q + pq'");
    if (!rpol_equal(&dpq, &rhs)) {
      printf("  p: ");
      rpol_print(&p);
      printf("  q: ");
      rpol_print(&q);
      printf("  (pq)': ");
      rpol_print(&dpq);
      printf("  p'q+pq': ");
      rpol_print(&rhs);
    }
  }

  rpol_clear(&p);
  rpol_clear(&q);
  rpol_clear(&pq);
  rpol_clear(&dpq);
  rpol_clear(&dp);
  rpol_clear(&dq);
  rpol_clear(&t1);
  rpol_clear(&t2);
  rpol_clear(&rhs);
}

int main(void)
{
  srand(20261030);
  clock_t start = clock();

  test_pseudo_div_hand();
  test_pseudo_div_random(60);
  test_div_exact_hand();
  test_div_exact_roundtrip(50);
  test_content_pp();
  test_pp_invariants(40);
  test_gcd_hand();
  test_gcd_oracle_small(80);
  test_gcd_constructed(40);
#ifdef HAVE_FLINT
  test_gcd_flint(60);
#endif
  test_eval_hand();
  test_eval_random(40);
  test_derive_hand();
  test_derive_random(40);

  double secs = (double)(clock() - start) / CLOCKS_PER_SEC;
  printf("test_poly_div: %.2fs elapsed\n", secs);

  if (failures == 0) {
    printf("test_poly_div: ALL CHECKS PASSED\n");
  } else {
    printf("test_poly_div: %d CHECK(S) FAILED\n", failures);
  }

  rz_clear_constants();
  return failures == 0 ? 0 : 1;
}
