#include "../../include/rq.h"

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

rabin_err_t rq_add(rq_t* r, const rq_t* p, const rq_t* q)
{
  if (r == NULL || p == NULL || q == NULL) return RABIN_ERR_NULL_PTR;

  /* p = a / c, q = b /d
   * p + q = ad + bc / cd
   * r may alias p or q, so read all of p and q before writing r.
   */

  rz_t t1, t2, td;
  rz_init_multi(&t1, &t2, &td, NULL);

  rabin_err_t err = rz_mul(&td, &p->den, &q->den);
  if (err == RABIN_SUCCESS) err = rz_mul(&t1, &p->num, &q->den);
  if (err == RABIN_SUCCESS) err = rz_mul(&t2, &q->num, &p->den);
  if (err == RABIN_SUCCESS) err = rz_add(&r->num, &t1, &t2);
  if (err == RABIN_SUCCESS) err = rz_copy(&r->den, &td);

  if (err == RABIN_SUCCESS) {
    err = rq_normalize(r);
  }

  rz_clear_multi(&t1, &t2, &td, NULL);
  return err;
}

rabin_err_t rq_sub(rq_t* r, const rq_t* p, const rq_t* q)
{
  if (r == NULL || p == NULL || q == NULL) return RABIN_ERR_NULL_PTR;

  /* p = a / c, q = b /d
   * p - q = ad - bc / cd
   * r may alias p or q, so read all of p and q before writing r.
   */

  rz_t t1, t2, td;
  rz_init_multi(&t1, &t2, &td, NULL);

  rabin_err_t err = rz_mul(&td, &p->den, &q->den);
  if (err == RABIN_SUCCESS) err = rz_mul(&t1, &p->num, &q->den);
  if (err == RABIN_SUCCESS) err = rz_mul(&t2, &q->num, &p->den);
  if (err == RABIN_SUCCESS) err = rz_sub(&r->num, &t1, &t2);
  if (err == RABIN_SUCCESS) err = rz_copy(&r->den, &td);

  if (err == RABIN_SUCCESS) {
    err = rq_normalize(r);
  }

  rz_clear_multi(&t1, &t2, &td, NULL);
  return err;
}

rabin_err_t rq_print(const rq_t* r)
{
  if (r == NULL) return RABIN_ERR_NULL_PTR;

  if (rz_is_zero(&r->num)) {
    printf("0\n");
  } else if (rz_is_eq_i64(&r->den, 1)) {
    rz_println(&r->num);
  } else {
    rz_print(&r->num);
    printf("/");
    rz_println(&r->den);
  }
  return RABIN_SUCCESS;
}

rabin_err_t rq_normalize(rq_t* r)
{
  if (r == NULL) return RABIN_ERR_NULL_PTR;

  // 0 case
  if (rz_is_zero(&r->num)) {
    return rz_set_u64(&r->den, 1);
  }

  // den can't be 0
  if (rz_is_zero(&r->den)) {
    return RABIN_ERR_INVALID_ARG;
  }

  // move negative sign
  if (r->den.is_neg) {
    rabin_err_t err = rz_neg(&r->num, &r->num);
    if (err != RABIN_SUCCESS) return err;
    if ((err = rz_neg(&r->den, &r->den)) != RABIN_SUCCESS) return err;
  }

  rz_t d, num_a;
  if (rz_init_multi(&d, &num_a, NULL) != RABIN_SUCCESS)
    return RABIN_ERR_OUT_OF_MEMORY;

  rabin_err_t err = rz_copy(&num_a, &r->num);
  if (err == RABIN_SUCCESS) num_a.is_neg = false;
  if (err == RABIN_SUCCESS) err = rz_gcd_lehmer(&d, &num_a, &r->den);
  if (err == RABIN_SUCCESS) err = rz_div(&r->den, &r->den, &d);
  if (err == RABIN_SUCCESS) err = rz_div(&r->num, &r->num, &d);

  rz_clear_multi(&d, &num_a, NULL);

  return err;
}

rabin_err_t rq_neg(rq_t* r)
{
  if (r == NULL) return RABIN_ERR_NULL_PTR;

  return rz_neg(&r->num, &r->num);
}

rabin_err_t rq_mul(rq_t* r, const rq_t* p, const rq_t* q)
{
  if (r == NULL || p == NULL || q == NULL) return RABIN_ERR_NULL_PTR;

  rz_t g_1, g_2, a_n, a_d, b_n, b_d;
  rz_init_multi(&g_1, &g_2, &a_n, &a_d, &b_n, &b_d, NULL);

  // p = a / c, q = b / d
  // g_1 = gcd(a, d)
  // g_2 = gcd(c, b)
  rabin_err_t err = rz_gcd_lehmer(&g_1, &p->num, &q->den);
  if (err == RABIN_SUCCESS) err = rz_gcd_lehmer(&g_2, &p->den, &q->num);
  if (err == RABIN_SUCCESS) err = rz_div(&a_n, &p->num, &g_1);
  if (err == RABIN_SUCCESS) err = rz_div(&a_d, &q->den, &g_1);
  if (err == RABIN_SUCCESS) err = rz_div(&b_n, &q->num, &g_2);
  if (err == RABIN_SUCCESS) err = rz_div(&b_d, &p->den, &g_2);

  // ((a/g1)(c/g2) / (d/g1)(b/g2))
  if (err == RABIN_SUCCESS) err = rz_mul(&r->num, &a_n, &b_n);
  if (err == RABIN_SUCCESS) err = rz_mul(&r->den, &a_d, &b_d);

  if (err == RABIN_SUCCESS) {
    err = rq_normalize(r);
  }

  rz_clear_multi(&g_1, &g_2, &a_n, &a_d, &b_n, &b_d, NULL);

  return err;
}

rabin_err_t rq_div(rq_t* r, const rq_t* p, const rq_t* q)
{
  if (r == NULL || p == NULL || q == NULL) return RABIN_ERR_NULL_PTR;
  if (rz_is_zero(&q->num)) return RABIN_ERR_DIV_BY_ZERO;

  /* p = a/c, q = b/d
   * p / q = a/c * d/b
   * g_1 = gcd(a, b)
   * g_2 = gcd(c, d)
   * p / q = (a/g_1) / (c/g_2) * (d/g_2) / (b/g_1)
   */

  rz_t g_1, g_2, a_n, b_n, c_n, d_n;
  rz_init_multi(&g_1, &g_2, &a_n, &b_n, &c_n, &d_n, NULL);

  rabin_err_t err = rz_gcd_lehmer(&g_1, &p->num, &q->num);
  if (err == RABIN_SUCCESS) err = rz_gcd_lehmer(&g_2, &p->den, &q->den);
  if (err == RABIN_SUCCESS) err = rz_div(&a_n, &p->num, &g_1);
  if (err == RABIN_SUCCESS) err = rz_div(&b_n, &q->num, &g_1);
  if (err == RABIN_SUCCESS) err = rz_div(&c_n, &p->den, &g_2);
  if (err == RABIN_SUCCESS) err = rz_div(&d_n, &q->den, &g_2);

  if (err == RABIN_SUCCESS) err = rz_mul(&r->num, &a_n, &d_n);
  if (err == RABIN_SUCCESS) err = rz_mul(&r->den, &c_n, &b_n);

  if (err == RABIN_SUCCESS) {
    err = rq_normalize(r);
  }

  rz_clear_multi(&g_1, &g_2, &a_n, &b_n, &c_n, &d_n, NULL);
  return err;
}

rabin_err_t rq_init(rq_t* r)
{
  if (r == NULL) return RABIN_ERR_NULL_PTR;

  rabin_err_t err = rz_init(&r->den);
  if (err == RABIN_SUCCESS) err = rz_init(&r->num);
  return err;
}

rabin_err_t rq_init_multi(rq_t* r, ...)
{
  if (r == NULL) return RABIN_ERR_NULL_PTR;

  rabin_err_t err = rq_init(r);

  va_list arg;
  va_start(arg, r);

  rq_t* next;
  while (err == RABIN_SUCCESS && (next = va_arg(arg, rq_t*)) != NULL) {
    err = rq_init(next);
  }

  va_end(arg);
  return err;
}

rabin_err_t rq_clear(rq_t* r)
{
  if (r == NULL) return RABIN_ERR_NULL_PTR;

  rabin_err_t err = rz_clear(&r->den);
  if (err == RABIN_SUCCESS) err = rz_clear(&r->num);
  return err;
}

rabin_err_t rq_clear_multi(rq_t* r, ...)
{
  if (r == NULL) return RABIN_ERR_NULL_PTR;

  rabin_err_t err = rq_clear(r);

  va_list arg;
  va_start(arg, r);

  rq_t* next;
  while (err == RABIN_SUCCESS && (next = va_arg(arg, rq_t*)) != NULL) {
    err = rq_clear(next);
  }

  va_end(arg);
  return err;
}

rabin_err_t rq_set(rq_t* r, const rz_t* a)
{
  if (r == NULL || a == NULL) return RABIN_ERR_NULL_PTR;

  rabin_err_t err = rz_copy(&r->num, a);
  if (err == RABIN_SUCCESS) err = rz_set_u64(&r->den, 1);
  return err;
}

bool rq_is_zero(const rq_t* r)
{
  if (r == NULL) return false;

  return rz_is_zero(&r->num);
}

int rq_cmp(const rq_t* a, const rq_t* b)
{
  if (a == NULL || b == NULL) return 0;

  rz_t t1, t2;
  if (rz_init_multi(&t1, &t2, NULL) != RABIN_SUCCESS) return 0;

  rz_mul(&t1, &a->num, &b->den);
  rz_mul(&t2, &b->num, &a->den);
  int res = rz_cmp(&t1, &t2);

  rz_clear_multi(&t1, &t2, NULL);
  return res;
}

rabin_err_t rq_copy(rq_t* r, const rq_t* a)
{
  if (r == NULL || a == NULL) return RABIN_ERR_NULL_PTR;
  if (r == a) return RABIN_SUCCESS;

  rabin_err_t err = rz_copy(&r->num, &a->num);
  if (err == RABIN_SUCCESS) err = rz_copy(&r->den, &a->den);
  if (err != RABIN_SUCCESS) return err;

  return rq_normalize(r);
}

rabin_err_t rq_sqr(rq_t* r, const rq_t* a)
{
  if (r == NULL || a == NULL) return RABIN_ERR_NULL_PTR;

  return rq_mul(r, a, a);
}

rabin_err_t rq_set_u64(rq_t* r, u64 a)
{
  if (r == NULL) return RABIN_ERR_NULL_PTR;

  rabin_err_t err = rz_set_u64(&r->num, a);
  if (err == RABIN_SUCCESS) err = rz_set_u64(&r->den, 1);
  return err;
}

static int rq_char_is_ws(char c)
{
  return c == ' ' || c == '\t' || c == '\n' || c == '\r';
}

/* rz_init_val on the digit run [i, j): copy out, NUL-terminate, parse. */
static rabin_err_t rq_parse_uint(rz_t* out, const char* s, size_t i, size_t j)
{
  char* tmp = strndup(s + i, j - i);
  if (tmp == NULL) return RABIN_ERR_OUT_OF_MEMORY;

  rabin_err_t err = rz_init_val(out, tmp);
  free(tmp);
  return err;
}

rabin_err_t rq_init_val(rq_t* r, const char* str)
{
  if (r == NULL || str == NULL) return RABIN_ERR_NULL_PTR;

  size_t n = strlen(str);
  size_t i = 0, j = n;
  while (i < j && rq_char_is_ws(str[i])) i++;
  while (j > i && rq_char_is_ws(str[j - 1])) j--;

  int neg = 0;
  if (i < j && (str[i] == '+' || str[i] == '-')) {
    neg = (str[i] == '-');
    i++;
  }

  const size_t ns = i;
  while (i < j && str[i] >= '0' && str[i] <= '9') i++;
  const size_t ne = i;
  if (ne == ns) {
    rz_set_u64(&r->num, 0);
    rz_set_u64(&r->den, 0);
    return RABIN_ERR_INVALID_ARG;
  }

  int have_den = 0;
  size_t den_start = 0, den_end = 0;
  if (i < j && str[i] == '/') {
    i++;
    den_start = i;
    while (i < j && str[i] >= '0' && str[i] <= '9') i++;
    den_end = i;
    if (den_end == den_start) {
      rz_set_u64(&r->num, 0);
      rz_set_u64(&r->den, 0);
      return RABIN_ERR_INVALID_ARG;
    }
    have_den = 1;
  }
  if (i != j) {
    rz_set_u64(&r->num, 0);
    rz_set_u64(&r->den, 0);
    return RABIN_ERR_INVALID_ARG;
  }

  rabin_err_t err = rq_parse_uint(&r->num, str, ns, ne);
  if (err == RABIN_SUCCESS) {
    if (have_den)
      err = rq_parse_uint(&r->den, str, den_start, den_end);
    else
      err = rz_set_u64(&r->den, 1);
  }
  if (err != RABIN_SUCCESS) {
    rz_set_u64(&r->num, 0);
    rz_set_u64(&r->den, 0);
    return err;
  }

  if (rz_is_zero(&r->den)) {
    rz_set_u64(&r->num, 0);
    rz_set_u64(&r->den, 0);
    return RABIN_ERR_DIV_BY_ZERO;
  }

  if (neg) rz_neg(&r->num, &r->num);
  return rq_normalize(r);
}

size_t rq_to_string_len(const rq_t* r)
{
  if (r == NULL) return 0;

  char* ns = rz_to_string(&r->num);
  if (ns == NULL) return 0;
  size_t len = strlen(ns) + 1; /* digits + NUL */

  if (!rz_is_one(&r->den)) {
    char* ds = rz_to_string(&r->den);
    if (ds == NULL) {
      free(ns);
      return 0;
    }
    len += strlen(ds) + 1; /* '/' + digits */
    free(ds);
  }

  free(ns);
  return len;
}

rabin_err_t rq_to_string(const rq_t* r, char* out, size_t cap)
{
  if (r == NULL || out == NULL) return RABIN_ERR_NULL_PTR;
  if (cap == 0) return RABIN_ERR_INVALID_ARG;

  char* ns = rz_to_string(&r->num);
  if (ns == NULL) return RABIN_ERR_OUT_OF_MEMORY;

  rabin_err_t err = RABIN_SUCCESS;
  if (rz_is_one(&r->den)) {
    if (strlen(ns) + 1 > cap)
      err = RABIN_ERR_INVALID_ARG;
    else
      snprintf(out, cap, "%s", ns);
  } else {
    char* ds = rz_to_string(&r->den);
    if (ds == NULL) {
      free(ns);
      return RABIN_ERR_OUT_OF_MEMORY;
    }
    if (strlen(ns) + strlen(ds) + 2 > cap)
      err = RABIN_ERR_INVALID_ARG;
    else
      snprintf(out, cap, "%s/%s", ns, ds);
    free(ds);
  }

  free(ns);
  return err;
}

rabin_err_t rq_set_i64(rq_t* r, i64 a)
{
  if (r == NULL) return RABIN_ERR_NULL_PTR;

  rabin_err_t err = rz_set_i64(&r->num, a);
  if (err == RABIN_SUCCESS) err = rz_set_u64(&r->den, 1);
  return err;
}

rabin_err_t rq_inv(rq_t* r, const rq_t* a)
{
  if (r == NULL || a == NULL) return RABIN_ERR_NULL_PTR;
  if (rz_is_zero(&a->num)) return RABIN_ERR_DIV_BY_ZERO;
  if (rz_is_zero(&a->den)) return RABIN_ERR_INVALID_ARG;

  /* r may alias a: stage both fields first. */
  rz_t n, d;
  rz_init_multi(&n, &d, NULL);

  rabin_err_t err = rz_copy(&n, &a->den);
  if (err == RABIN_SUCCESS) err = rz_copy(&d, &a->num);
  if (err == RABIN_SUCCESS) err = rz_copy(&r->num, &n);
  if (err == RABIN_SUCCESS) err = rz_copy(&r->den, &d);
  if (err == RABIN_SUCCESS) err = rq_normalize(r);

  rz_clear_multi(&n, &d, NULL);
  return err;
}

rabin_err_t rq_get_num(rz_t* out, const rq_t* r)
{
  if (out == NULL || r == NULL) return RABIN_ERR_NULL_PTR;
  return rz_copy(out, &r->num);
}

rabin_err_t rq_get_den(rz_t* out, const rq_t* r)
{
  if (out == NULL || r == NULL) return RABIN_ERR_NULL_PTR;
  return rz_copy(out, &r->den);
}

int rq_cmp_i64(const rq_t* r, i64 a)
{
  if (r == NULL) return 0;
  if (rz_is_zero(&r->den)) return 0;

  rz_t v;
  if (rz_init(&v) != RABIN_SUCCESS) return 0;

  int res = 0;
  if (rz_set_i64(&v, a) == RABIN_SUCCESS &&
      rz_mul(&v, &v, &r->den) == RABIN_SUCCESS) {
    res = rz_cmp(&r->num, &v);
  }

  rz_clear(&v);
  return res;
}
