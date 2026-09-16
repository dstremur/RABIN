
#include "../include/bigrat.h"

int main()
{
  bigrat r, q, p;
  br_init(&r);
  br_init(&q);
  br_init(&p);

  bn_set_i64(&r.den, 8905493208);
  bn_set_i64(&r.num, 3);

  bn_set_i64(&q.den, 23534563456);
  bn_set_i64(&q.num, 3);

  br_print(&r);

  br_mul(&p, &r, &q);
  br_mul(&r, &p, &q);

  br_print(&p);
  br_print(&r);

  bn_set_i64(&r.den, 3);
  bn_set_i64(&r.num, 1);

  bn_set_i64(&q.den, 4);
  bn_set_i64(&q.num, 1);

  br_div(&p, &r, &q);

  br_print(&p);

  br_free(&r);
  br_free(&p);
  br_free(&q);
}
