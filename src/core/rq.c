#include "../../include/rq.h"

#include <stdarg.h>

rabin_err_t rq_add(rq_t* r, const rq_t* p, const rq_t* q)
{
  if (r == NULL || p == NULL || q == NULL) return RABIN_ERR_NULL_PTR;

  /* p = a / c, q = b /d
   * p + q = ad + bc / cd
   */

  rz_t t1, t2;
  rz_init_multi(&t1, &t2, NULL);

  rabin_err_t err = rz_mul(&r->den, &p->den, &q->den);
  if (err == RABIN_SUCCESS) err = rz_mul(&t1, &p->num, &q->den);
  if (err == RABIN_SUCCESS) err = rz_mul(&t2, &q->num, &p->den);
  if (err == RABIN_SUCCESS) err = rz_add(&r->num, &t1, &t2);

  if (err == RABIN_SUCCESS) {
    err = rq_normalize(r);
  }

  rz_clear_multi(&t1, &t2, NULL);
  return err;
}

rabin_err_t rq_sub(rq_t* r, const rq_t* p, const rq_t* q)
{
  if (r == NULL || p == NULL || q == NULL) return RABIN_ERR_NULL_PTR;

  /* p = a / c, q = b /d
   * p - q = ad - bc / cd
   */

  rz_t t1, t2;
  rz_init_multi(&t1, &t2, NULL);

  rabin_err_t err = rz_mul(&r->den, &p->den, &q->den);
  if (err == RABIN_SUCCESS) err = rz_mul(&t1, &p->num, &q->den);
  if (err == RABIN_SUCCESS) err = rz_mul(&t2, &q->num, &p->den);
  if (err == RABIN_SUCCESS) err = rz_sub(&r->num, &t1, &t2);

  if (err == RABIN_SUCCESS) {
    err = rq_normalize(r);
  }

  rz_clear_multi(&t1, &t2, NULL);
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
