/*
 * test_serialization.c
 *
 * Round-trip suite for the string grammars of rz_t, rq_t, rvec_t,
 * rmat_t and rpol_t: exact-format tables, malformed-input rejection,
 * failure atomicity, capacity enforcement, and randomized round-trips
 * with independently computed checks.
 *
 * Copyright (C) 2026 Diego Strebel
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <stdlib.h>
#include <string.h>

#include "../include/rmat.h"
#include "../include/rpol.h"
#include "../include/rq.h"
#include "../include/rvec.h"
#include "../include/rzrand.h"
#include "test_util.h"

/* A 40-digit value that exceeds every machine width. */
#define BIGSTR "1234567890123456789012345678901234567890"

/* rpol_equal() dereferences coeff[] even at the empty state; guard it. */
static bool rpol_same(const rpol_t* a, const rpol_t* b)
{
  bool ea = (a->coeff == NULL || a->size == 0);
  bool eb = (b->coeff == NULL || b->size == 0);
  if (ea || eb) return ea && eb;
  return rpol_equal(a, b);
}

/* Random signed decimal string of the given bit size into buf. */
static void rand_dec_string(char* buf, size_t cap, unsigned bits)
{
  rz_t t;
  rz_init(&t);
  rz_gen_random(&t, bits > 0 ? bits : 1);
  if (rand() & 1) rz_neg(&t, &t);

  char* s = rz_to_string(&t);
  snprintf(buf, cap, "%s", s ? s : "0");
  free(s);
  rz_clear(&t);
}

/* Write cap-filling sentinel byte pattern for untouched-buffer checks. */
static void sentinel(char* buf, size_t cap)
{
  memset(buf, '#', cap);
  buf[cap - 1] = '\0';
}

static bool buffer_untouched(const char* buf, size_t cap)
{
  for (size_t i = 0; i + 1 < cap; i++)
    if (buf[i] != '#') return false;
  return buf[cap - 1] == '\0';
}

TEST(rz_round_trip)
{
  static const struct {
    const char* in;
    const char* out;
  } cases[] = {
      {"0", "0"},
      {"-0", "0"},
      {"7", "7"},
      {"-7", "-7"},
      {"007", "7"},
      {"-007", "-7"},
      {"18446744073709551616", "18446744073709551616"}, /* 2^64 */
      {"-18446744073709551616", "-18446744073709551616"},
      {BIGSTR, BIGSTR},
      {"-" BIGSTR, "-" BIGSTR},
  };

  rz_t a, b;
  rz_init_multi(&a, &b, NULL);

  for (size_t i = 0; i < sizeof(cases) / sizeof(cases[0]); i++) {
    TEST_ASSERT(rz_init_val(&a, cases[i].in) == RABIN_SUCCESS, "rz parse");
    char* s = rz_to_string(&a);
    TEST_ASSERT(s != NULL, "rz print");
    if (s != NULL) {
      TEST_ASSERT(strcmp(s, cases[i].out) == 0, "rz exact format");
      if (strcmp(s, cases[i].out) != 0)
        printf("  (\"%s\" -> \"%s\", expected \"%s\")\n", cases[i].in, s,
               cases[i].out);
      TEST_ASSERT(rz_init_val(&b, s) == RABIN_SUCCESS, "rz reparse");
      TEST_ASSERT(rz_cmp(&a, &b) == 0, "rz round trip");
      free(s);
    }
  }

  srand(7);
  for (int i = 0; i < 80; i++) {
    char buf[512];
    rand_dec_string(buf, sizeof buf, 1 + (unsigned)(rand() % 700));
    TEST_ASSERT(rz_init_val(&a, buf) == RABIN_SUCCESS, "rz random parse");
    char* s = rz_to_string(&a);
    TEST_ASSERT(s != NULL && strcmp(s, buf) == 0, "rz random round trip");
    if (s != NULL && strcmp(s, buf) != 0)
      printf("  (random: \"%s\" -> \"%s\")\n", buf, s);
    free(s);
  }

  rz_clear_multi(&a, &b, NULL);
}

TEST(rq_round_trip)
{
  static const struct {
    const char* in;
    const char* out;
  } cases[] = {
      {"3/4", "3/4"},
      {"-3/4", "-3/4"},
      {"6/8", "3/4"},
      {"7/1", "7"},
      {"0/9", "0"},
      {"-17", "-17"},
      {"22/7", "22/7"},
      {BIGSTR "/2",
       "617283945061728394506172839450617283945"}, /* reduces to int */
      {"-" BIGSTR "/3", "-411522630041152263004115226300411522630"},
      {"1234567890123456789012345678901234567891/7",
       "1234567890123456789012345678901234567891/7"}, /* stays a fraction */
  };

  rq_t a, b;
  rq_init_multi(&a, &b, NULL);

  for (size_t i = 0; i < sizeof(cases) / sizeof(cases[0]); i++) {
    TEST_ASSERT(rq_init_val(&a, cases[i].in) == RABIN_SUCCESS, "rq parse");
    size_t need = rq_to_string_len(&a);
    TEST_ASSERT(need == strlen(cases[i].out) + 1, "rq len exact");

    char buf[128];
    TEST_ASSERT(rq_to_string(&a, buf, sizeof buf) == RABIN_SUCCESS, "rq print");
    TEST_ASSERT(strcmp(buf, cases[i].out) == 0, "rq exact format");
    if (strcmp(buf, cases[i].out) != 0)
      printf("  (\"%s\" -> \"%s\", expected \"%s\")\n", cases[i].in, buf,
             cases[i].out);

    TEST_ASSERT(rq_init_val(&b, buf) == RABIN_SUCCESS, "rq reparse");
    TEST_ASSERT(rq_cmp(&a, &b) == 0, "rq round trip");
    /* Reprint must produce the identical canonical string. */
    char buf2[128];
    TEST_ASSERT(rq_to_string(&b, buf2, sizeof buf2) == RABIN_SUCCESS,
                "rq reprint");
    TEST_ASSERT(strcmp(buf, buf2) == 0, "rq canonical stable");
  }

  rq_clear_multi(&a, &b, NULL);
}

TEST(rvec_round_trip)
{
  static const struct {
    const char* in;
    u64 size;
    const char* out; /* NULL: expected to fail */
    rabin_err_t fail;
  } cases[] = {
      {"[]", 0, "[]", 0},
      {"[0]", 1, "[0]", 0},
      {"[-10,+20, 30 ]", 3, "[-10,20,30]", 0},
      {" [ 4 , -5 ] ", 2, "[4,-5]", 0},
      {"[" BIGSTR ",-" BIGSTR "]", 2, NULL, 0}, /* format check skipped */
      {"[1,2]", 3, NULL, RABIN_ERR_INVALID_ARG},
      {"[1,2,3,4]", 3, NULL, RABIN_ERR_INVALID_ARG},
      {"[1,,2]", 3, NULL, RABIN_ERR_INVALID_ARG},
      {"[1,2", 2, NULL, RABIN_ERR_INVALID_ARG},
      {"1,2", 2, NULL, RABIN_ERR_INVALID_ARG},
      {"[abc]", 3, NULL, RABIN_ERR_INVALID_ARG},
      {"[1,2]x", 2, NULL, RABIN_ERR_INVALID_ARG},
      {"[,1]", 2, NULL, RABIN_ERR_INVALID_ARG},
      {"[]", 2, NULL, RABIN_ERR_INVALID_ARG},
  };

  for (size_t i = 0; i < sizeof(cases) / sizeof(cases[0]); i++) {
    rvec_t v = {0};
    rvec_init(&v, cases[i].size);

    rabin_err_t err = rvec_init_val(&v, cases[i].in);
    if (cases[i].fail != RABIN_SUCCESS) {
      TEST_ASSERT(err == cases[i].fail, "rvec rejection");
      if (err != cases[i].fail)
        printf("  (\"%s\" size %llu -> %d)\n", cases[i].in,
               (unsigned long long)cases[i].size, err);
      rvec_clear(&v);
      continue;
    }
    TEST_ASSERT(err == RABIN_SUCCESS, "rvec parse");
    if (err != RABIN_SUCCESS) {
      rvec_clear(&v);
      continue;
    }

    if (cases[i].out != NULL) {
      char buf[128];
      TEST_ASSERT(rvec_to_string(&v, buf, sizeof buf) == RABIN_SUCCESS,
                  "rvec print");
      TEST_ASSERT(strcmp(buf, cases[i].out) == 0, "rvec exact format");
      if (strcmp(buf, cases[i].out) != 0)
        printf("  (\"%s\" -> \"%s\", expected \"%s\")\n", cases[i].in, buf,
               cases[i].out);

      /* Reparse into a fresh vector and compare elementwise. */
      rvec_t w = {0};
      rvec_init(&w, v.size);
      TEST_ASSERT(rvec_init_val(&w, buf) == RABIN_SUCCESS, "rvec reparse");
      for (u64 k = 0; k < v.size; k++)
        TEST_ASSERT(rz_cmp(&v.data[k], &w.data[k]) == 0, "rvec element equal");
      rvec_clear(&w);
    }
    rvec_clear(&v);
  }

  /* Capacity enforcement: too-small cap must leave out untouched. */
  rvec_t v = {0};
  rvec_init(&v, 2);
  TEST_ASSERT(rvec_init_val(&v, "[100,-200]") == RABIN_SUCCESS,
              "rvec cap case");
  size_t need = rvec_to_string_len(&v);
  TEST_ASSERT(need == strlen("[100,-200]") + 1, "rvec len exact");
  char buf[32];
  sentinel(buf, sizeof buf);
  TEST_ASSERT(rvec_to_string(&v, buf, need - 1) == RABIN_ERR_INVALID_ARG,
              "rvec small cap rejected");
  TEST_ASSERT(buffer_untouched(buf, sizeof buf), "rvec small cap leaves out");
  TEST_ASSERT(rvec_to_string(&v, buf, need) == RABIN_SUCCESS, "rvec exact cap");
  TEST_ASSERT(strcmp(buf, "[100,-200]") == 0, "rvec exact cap content");
  rvec_clear(&v);

  /* Failure atomicity: a failed parse leaves the vector unchanged. */
  rvec_t a = {0};
  rvec_init(&a, 3);
  TEST_ASSERT(rvec_init_val(&a, "[1,2,3]") == RABIN_SUCCESS,
              "rvec atomic seed");
  TEST_ASSERT(rvec_init_val(&a, "[9,,9]") != RABIN_SUCCESS, "rvec bad parse");
  TEST_ASSERT(rvec_init_val(&a, "[9,9,9,9]") != RABIN_SUCCESS,
              "rvec bad count parse");
  for (u64 k = 0; k < 3; k++)
    TEST_ASSERT(rz_is_eq_i64(&a.data[k], (i64)(k + 1)),
                "rvec unchanged on fail");
  rvec_clear(&a);

  /* NULL guards. */
  TEST_ASSERT(rvec_to_string_len(NULL) == 0, "rvec NULL len");
  TEST_ASSERT(rvec_to_string(NULL, buf, sizeof buf) == RABIN_ERR_NULL_PTR,
              "rvec NULL print obj");
  TEST_ASSERT(rvec_init_val(NULL, "[1]") == RABIN_ERR_NULL_PTR,
              "rvec NULL parse obj");

  srand(11);
  for (int i = 0; i < 40; i++) {
    u64 size = (u64)(rand() % 9);
    char joined[4096];
    size_t w = 0;
    joined[w++] = '[';
    for (u64 k = 0; k < size; k++) {
      char part[256];
      rand_dec_string(part, sizeof part, 1 + (unsigned)(rand() % 400));
      if (k > 0) joined[w++] = ',';
      w += (size_t)snprintf(joined + w, sizeof joined - w, "%s", part);
    }
    joined[w++] = ']';
    joined[w] = '\0';

    rvec_t v = {0}, wv = {0};
    rvec_init(&v, size);
    rvec_init(&wv, size);
    TEST_ASSERT(rvec_init_val(&v, joined) == RABIN_SUCCESS, "rvec rand parse");
    char out[4096];
    size_t need = rvec_to_string_len(&v);
    TEST_ASSERT(need <= sizeof out, "rvec rand len fits");
    TEST_ASSERT(rvec_to_string(&v, out, sizeof out) == RABIN_SUCCESS,
                "rvec rand print");
    TEST_ASSERT(strcmp(out, joined) == 0, "rvec rand string stable");
    if (strcmp(out, joined) != 0)
      printf("  (rvec rand mismatch round %d)\n", i);
    TEST_ASSERT(rvec_init_val(&wv, out) == RABIN_SUCCESS, "rvec rand reparse");
    for (u64 k = 0; k < size; k++)
      TEST_ASSERT(rz_cmp(&v.data[k], &wv.data[k]) == 0, "rvec rand elements");
    rvec_clear(&v);
    rvec_clear(&wv);
  }
}

TEST(rmat_round_trip)
{
  static const struct {
    const char* in;
    u64 rows;
    u64 cols;
    const char* out;
    rabin_err_t fail;
  } cases[] = {
      {"[[1,2,3],[-4,0,5]]", 2, 3, "[[1,2,3],[-4,0,5]]", 0},
      {" [ [ -1 ] , [ +2 ] ] ", 2, 1, "[[-1],[2]]", 0},
      {"[[]]", 1, 0, "[[]]", 0},
      {"[[],[]]", 2, 0, "[[],[]]", 0},
      {"[[],[]]", 1, 0, NULL, RABIN_ERR_MATRIX_DIM},
      {"[[1,2],[3,4]]", 3, 2, NULL, RABIN_ERR_MATRIX_DIM},
      {"[[1,2,3],[4,5]]", 2, 3, NULL, RABIN_ERR_MATRIX_DIM},
      {"[[1,2,3][4,5,6]]", 2, 3, NULL, RABIN_ERR_INVALID_ARG},
      {"[[1,2,3],[4,5,6]", 2, 3, NULL, RABIN_ERR_INVALID_ARG},
      {"[]", 2, 2, NULL, RABIN_ERR_MATRIX_DIM},
      {"[[1,2,3],[4,5,6]]x", 2, 3, NULL, RABIN_ERR_INVALID_ARG},
      {"[[1,2,3],[4,5,+]]", 2, 3, NULL, RABIN_ERR_INVALID_ARG},
  };

  for (size_t i = 0; i < sizeof(cases) / sizeof(cases[0]); i++) {
    rmat_t A = {0};
    rmat_init(&A, cases[i].rows, cases[i].cols);

    rabin_err_t err = rmat_init_val(&A, cases[i].in);
    if (cases[i].fail != RABIN_SUCCESS) {
      TEST_ASSERT(err == cases[i].fail, "rmat rejection");
      if (err != cases[i].fail)
        printf("  (\"%s\" %llux%llu -> %d)\n", cases[i].in,
               (unsigned long long)cases[i].rows,
               (unsigned long long)cases[i].cols, err);
      rmat_clear(&A);
      continue;
    }
    TEST_ASSERT(err == RABIN_SUCCESS, "rmat parse");
    if (err == RABIN_SUCCESS) {
      char buf[128];
      TEST_ASSERT(rmat_to_string(&A, buf, sizeof buf) == RABIN_SUCCESS,
                  "rmat print");
      TEST_ASSERT(strcmp(buf, cases[i].out) == 0, "rmat exact format");
      if (strcmp(buf, cases[i].out) != 0)
        printf("  (\"%s\" -> \"%s\", expected \"%s\")\n", cases[i].in, buf,
               cases[i].out);

      rmat_t B = {0};
      rmat_init(&B, cases[i].rows, cases[i].cols);
      TEST_ASSERT(rmat_init_val(&B, buf) == RABIN_SUCCESS, "rmat reparse");
      for (u64 k = 0; k < cases[i].rows * cases[i].cols; k++)
        TEST_ASSERT(rz_cmp(&A.data[k], &B.data[k]) == 0, "rmat entry equal");
      rmat_clear(&B);
    }
    rmat_clear(&A);
  }

  /* Capacity enforcement. */
  rmat_t A = {0};
  rmat_init(&A, 1, 2);
  TEST_ASSERT(rmat_init_val(&A, "[[100,-200]]") == RABIN_SUCCESS,
              "rmat cap seed");
  size_t need = rmat_to_string_len(&A);
  TEST_ASSERT(need == strlen("[[100,-200]]") + 1, "rmat len exact");
  char buf[32];
  sentinel(buf, sizeof buf);
  TEST_ASSERT(rmat_to_string(&A, buf, need - 1) == RABIN_ERR_INVALID_ARG,
              "rmat small cap rejected");
  TEST_ASSERT(buffer_untouched(buf, sizeof buf), "rmat small cap leaves out");
  rmat_clear(&A);

  /* Failure atomicity. */
  rmat_t S = {0};
  rmat_init(&S, 2, 2);
  TEST_ASSERT(rmat_init_val(&S, "[[1,2],[3,4]]") == RABIN_SUCCESS,
              "rmat atomic seed");
  TEST_ASSERT(rmat_init_val(&S, "[[1,2],[3,4,5]]") != RABIN_SUCCESS,
              "rmat bad parse");
  TEST_ASSERT(rmat_init_val(&S, "[[1,2],[3,4]x") != RABIN_SUCCESS,
              "rmat trailing parse");
  for (u64 k = 0; k < 4; k++)
    TEST_ASSERT(rz_is_eq_i64(&S.data[k], (i64)(k + 1)),
                "rmat unchanged on fail");
  rmat_clear(&S);

  srand(19);
  for (int i = 0; i < 30; i++) {
    u64 rows = 1 + (u64)(rand() % 5), cols = 1 + (u64)(rand() % 5);
    char joined[8192];
    size_t w = 0;
    joined[w++] = '[';
    for (u64 r = 0; r < rows; r++) {
      if (r > 0) joined[w++] = ',';
      joined[w++] = '[';
      for (u64 c = 0; c < cols; c++) {
        char part[256];
        rand_dec_string(part, sizeof part, 1 + (unsigned)(rand() % 400));
        if (c > 0) joined[w++] = ',';
        w += (size_t)snprintf(joined + w, sizeof joined - w, "%s", part);
      }
      joined[w++] = ']';
    }
    joined[w++] = ']';
    joined[w] = '\0';

    rmat_t A = {0}, B = {0};
    rmat_init(&A, rows, cols);
    rmat_init(&B, rows, cols);
    TEST_ASSERT(rmat_init_val(&A, joined) == RABIN_SUCCESS, "rmat rand parse");
    char out[8192];
    TEST_ASSERT(rmat_to_string(&A, out, sizeof out) == RABIN_SUCCESS,
                "rmat rand print");
    TEST_ASSERT(strcmp(out, joined) == 0, "rmat rand string stable");
    if (strcmp(out, joined) != 0)
      printf("  (rmat rand mismatch round %d)\n", i);
    TEST_ASSERT(rmat_init_val(&B, out) == RABIN_SUCCESS, "rmat rand reparse");
    for (u64 k = 0; k < rows * cols; k++)
      TEST_ASSERT(rz_cmp(&A.data[k], &B.data[k]) == 0, "rmat rand entries");
    rmat_clear(&A);
    rmat_clear(&B);
  }
}

TEST(rpol_round_trip)
{
  static const struct {
    const char* in;
    const char* out; /* NULL: must fail */
    rabin_err_t fail;
  } cases[] = {
      {"[1,-2,3]", "[1,-2,3]", 0},
      {"[ 5 , 0 , 0 ]", "[5]", 0}, /* high-degree zeros trimmed */
      {"[]", "[]", 0},
      {"[0]", "[0]", 0},
      {"[0,0,0]", "[0]", 0},
      {"[-7]", "[-7]", 0},
      {"[" BIGSTR "]", "[" BIGSTR "]", 0},
      {"[1,2]", "[1,2]", 0},
  };

  for (size_t i = 0; i < sizeof(cases) / sizeof(cases[0]); i++) {
    rpol_t p = {0}, q = {0};
    rpol_init(&p);
    rpol_init(&q);

    rabin_err_t err = rpol_init_val(&p, cases[i].in);
    if (cases[i].fail != RABIN_SUCCESS) {
      TEST_ASSERT(err == cases[i].fail, "rpol rejection");
      rpol_clear(&p);
      rpol_clear(&q);
      continue;
    }
    TEST_ASSERT(err == RABIN_SUCCESS, "rpol parse");
    if (err == RABIN_SUCCESS) {
      char buf[128];
      size_t need = rpol_to_string_len(&p);
      TEST_ASSERT(need == strlen(cases[i].out) + 1, "rpol len exact");
      TEST_ASSERT(rpol_to_string(&p, buf, sizeof buf) == RABIN_SUCCESS,
                  "rpol print");
      TEST_ASSERT(strcmp(buf, cases[i].out) == 0, "rpol exact format");
      if (strcmp(buf, cases[i].out) != 0)
        printf("  (\"%s\" -> \"%s\", expected \"%s\")\n", cases[i].in, buf,
               cases[i].out);
      TEST_ASSERT(rpol_init_val(&q, buf) == RABIN_SUCCESS, "rpol reparse");
      TEST_ASSERT(rpol_same(&p, &q), "rpol round trip");
    }
    rpol_clear(&p);
    rpol_clear(&q);
  }

  /* Malformed inputs (the placeholder row above). */
  static const char* bad[] = {"[1,2",   "[1,,2]", "1,2,3", "[abc]",
                              "[1,2]x", "[,1]",   "[1,,]", "[]x"};
  for (size_t i = 0; i < sizeof(bad) / sizeof(bad[0]); i++) {
    rpol_t p = {0};
    rpol_init(&p);
    TEST_ASSERT(rpol_init_val(&p, bad[i]) == RABIN_ERR_INVALID_ARG,
                "rpol malformed rejected");
    rpol_clear(&p);
  }

  /* Trim semantics: "[5,0,0]" stores degree 0. */
  rpol_t p = {0};
  rpol_init(&p);
  TEST_ASSERT(rpol_init_val(&p, "[5,0,0]") == RABIN_SUCCESS, "rpol trim seed");
  TEST_ASSERT(p.deg == 0, "rpol trimmed degree");
  TEST_ASSERT(rz_is_eq_i64(&p.coeff[0], 5), "rpol trimmed coefficient");

  /* Failure atomicity: bad parse leaves p unchanged. */
  TEST_ASSERT(rpol_init_val(&p, "[1,2,") == RABIN_ERR_INVALID_ARG,
              "rpol bad parse");
  TEST_ASSERT(p.deg == 0 && rz_is_eq_i64(&p.coeff[0], 5),
              "rpol unchanged on fail");

  /* Capacity enforcement. */
  size_t need = rpol_to_string_len(&p);
  char buf[16];
  sentinel(buf, sizeof buf);
  TEST_ASSERT(rpol_to_string(&p, buf, need - 1) == RABIN_ERR_INVALID_ARG,
              "rpol small cap rejected");
  TEST_ASSERT(buffer_untouched(buf, sizeof buf), "rpol small cap leaves out");
  rpol_clear(&p);

  /* Independent numeric check: 1 - 2x + 3x^2 at x = 3 is 22. */
  rpol_init(&p);
  TEST_ASSERT(rpol_init_val(&p, "[1,-2,3]") == RABIN_SUCCESS, "rpol eval seed");
  rz_t x, y;
  rz_init_multi(&x, &y, NULL);
  rz_set_u64(&x, 3);
  TEST_ASSERT(rpol_eval(&y, &p, &x) == RABIN_SUCCESS, "rpol eval");
  TEST_ASSERT(rz_is_eq_i64(&y, 22), "rpol eval oracle");
  rz_clear_multi(&x, &y, NULL);
  rpol_clear(&p);

  srand(23);
  for (int i = 0; i < 40; i++) {
    u64 count = (u64)(rand() % 9); /* 0 -> empty polynomial */
    char joined[4096];
    size_t w = 0;
    joined[w++] = '[';
    i64 last_nonzero = -1;
    i64 coeffs[9] = {0};
    for (u64 k = 0; k < count; k++) {
      i64 c = (i64)(rand() % 2001) - 1000;
      if (rand() % 4 == 0) c = 0; /* exercise trimming often */
      coeffs[k] = c;
      if (c != 0) last_nonzero = (i64)k;
      if (k > 0) joined[w++] = ',';
      w +=
          (size_t)snprintf(joined + w, sizeof joined - w, "%lld", (long long)c);
    }
    joined[w++] = ']';
    joined[w] = '\0';

    rpol_t p = {0}, q = {0};
    rpol_init(&p);
    rpol_init(&q);
    TEST_ASSERT(rpol_init_val(&p, joined) == RABIN_SUCCESS, "rpol rand parse");

    if (count == 0) {
      TEST_ASSERT(p.coeff == NULL && p.size == 0, "rpol empty state");
    } else {
      u64 expect_deg = last_nonzero >= 0 ? (u64)last_nonzero : 0;
      TEST_ASSERT(p.deg == expect_deg, "rpol rand trimmed degree");
      for (u64 k = 0; k <= p.deg; k++)
        TEST_ASSERT(rz_is_eq_i64(&p.coeff[k], coeffs[k]),
                    "rpol rand coefficient kept");
    }

    char out[4096];
    TEST_ASSERT(rpol_to_string(&p, out, sizeof out) == RABIN_SUCCESS,
                "rpol rand print");
    TEST_ASSERT(rpol_init_val(&q, out) == RABIN_SUCCESS, "rpol rand reparse");
    TEST_ASSERT(rpol_same(&p, &q), "rpol rand round trip");
    char out2[4096];
    TEST_ASSERT(rpol_to_string(&q, out2, sizeof out2) == RABIN_SUCCESS,
                "rpol rand reprint");
    TEST_ASSERT(strcmp(out, out2) == 0, "rpol rand string stable");
    rpol_clear(&p);
    rpol_clear(&q);
  }
}

int main(void)
{
  rz_init_constants();

  RUN_TEST(rz_round_trip);
  RUN_TEST(rq_round_trip);
  RUN_TEST(rvec_round_trip);
  RUN_TEST(rmat_round_trip);
  RUN_TEST(rpol_round_trip);

  rz_clear_constants();
  return TEST_REPORT();
}
