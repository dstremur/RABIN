#include "../../include/bigrat.h"

#include <stdarg.h>

void br_add(bigrat* r, const bigrat* a, const bigrat* b) {}

void br_print(const bigrat* r)
{
  bn_print(&r->num);
  printf("/");
  bn_println(&r->den);
}

bool br_normalize(bigrat* r)
{
  // 0 case
  if (bn_is_zero(&r->num)) {
    bn_set_u64(&r->den, 1);
    return true;
  }

  // den can't be 0
  if (bn_is_zero(&r->den)) {
    return false;
  }

  // move negative sign
  if (r->den.is_neg) {
    bn_neg(&r->num, &r->num);
    bn_neg(&r->den, &r->den);
  }

  bignum d, num_a;
  bn_init_multi(&d, &num_a);
  bn_copy(&num_a, &r->num);
  num_a.is_neg = false;

  bn_gcd_lehmer(&d, &num_a, &r->den);

  bn_div(&r->den, &r->den, &d);
  bn_div(&r->num, &r->num, &d);

  bn_free_multi(&d, &num_a);

  return true;
}

void br_neg(bigrat* r) { bn_neg(&r->num, &r->num); }

void br_mul(bigrat* r, const bigrat* p, const bigrat* q)
{
  bignum g_1, g_2, a_n, a_d, b_n, b_d;
  bn_init_multi(&g_1, &g_2, &a_n, &a_d, &b_n, &b_d, NULL);

  // p = a / c, q = b / d
  // g_1 = gcd(a, d)
  // g_2 = gcd(c, b)
  bn_gcd_lehmer(&g_1, &p->num, &q->den);

  bn_gcd_lehmer(&g_2, &p->den, &q->num);

  bn_div(&a_n, &p->num, &g_1);
  bn_div(&a_d, &q->den, &g_1);

  bn_div(&b_n, &q->num, &g_2);
  bn_div(&b_d, &p->den, &g_2);

  // ((a/g1)(c/g2) / (d/g1)(b/g2))
  bn_mul(&r->num, &a_n, &b_n);
  bn_mul(&r->den, &a_d, &b_d);

  br_normalize(r);

  bn_free_multi(&g_1, &g_2, &a_n, &a_d, &b_n, &b_d, NULL);

  return;
}

void br_div(bigrat* r, bigrat* p, bigrat* q)
{
  /* p = a/c, q = b/d
   * p / q = a/c * d/b
   * g_1 = gcd(a, b)
   * g_2 = gcd(c, d)
   * p / q = (a/g_1) / (c/g_2) * (d/g_2) / (b/g_1)
   */

  bignum g_1, g_2, a_n, b_n, c_n, d_n;
  bn_init_multi(&g_1, &g_2, &a_n, &b_n, &c_n, &d_n, NULL);

  bn_gcd_lehmer(&g_1, &p->num, &q->num);
  bn_gcd_lehmer(&g_2, &p->den, &q->den);

  bn_div(&a_n, &p->num, &g_1);
  bn_div(&b_n, &q->num, &g_1);
  bn_div(&c_n, &p->den, &g_2);
  bn_div(&d_n, &q->den, &g_2);

  bn_mul(&r->num, &a_n, &d_n);
  bn_mul(&r->den, &c_n, &b_n);

  bn_free_multi(&g_1, &g_2, &a_n, &b_n, &c_n, &d_n, NULL);
}

void br_init(bigrat* r)
{
  bn_init(&r->den);
  bn_init(&r->num);
}

void br_init_multi(bigrat* r, ...)
{
  if (r == NULL) return;

  br_init(r);

  va_list arg;
  va_start(arg, r);

  bigrat* next;
  while ((next = va_arg(arg, bigrat*)) != NULL) {
    br_init(next);
  }

  va_end(arg);
}

void br_free(bigrat* r)
{
  bn_free(&r->den);
  bn_free(&r->num);
}

void br_free_multi(bigrat* r, ...)
{
  if (r == NULL) return;

  br_free(r);

  va_list arg;
  va_start(arg, r);

  bigrat* next;
  while ((next = va_arg(arg, bigrat*)) != NULL) {
    br_free(next);
  }

  va_end(arg);
}
