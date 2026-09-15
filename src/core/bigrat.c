#include "../../include/bigrat.h"

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

void br_mul(bigrat* r, const bigrat* a, const bigrat* b)
{
  bignum g_1, g_2, a_n, a_d, b_n, b_d;
  bn_init_multi(&g_1, &g_2, &a_n, &a_d, &b_n, &b_d, NULL);

  bn_gcd_lehmer(&g_1, &a->num, &b->den);

  bn_gcd_lehmer(&g_2, &a->den, &b->num);

  bn_div(&a_n, &a->num, &g_1);
  bn_div(&a_d, &b->den, &g_1);

  bn_div(&b_n, &b->num, &g_2);
  bn_div(&b_d, &a->den, &g_2);

  bn_mul(&r->num, &a_n, &b_n);
  bn_mul(&r->den, &a_d, &b_d);

  br_normalize(r);

  bn_free_multi(&g_1, &g_2, &a_n, &a_d, &b_n, &b_d, NULL);

  return;
}

void br_init(bigrat* r)
{
  bn_init(&r->den);
  bn_init(&r->num);
}

void br_free(bigrat* r)
{
  bn_free(&r->den);
  bn_free(&r->num);
}
