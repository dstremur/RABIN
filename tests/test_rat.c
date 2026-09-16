
#include "../include/bigrat.h"

int main()
{
  bigrat r, q, p, k, l;
  br_init(&r);
  br_init(&q);
  br_init(&p);
  br_init_multi(&k, &l, NULL);

  bn_set_i64(&r.den, 8905493208);
  bn_set_i64(&r.num, 3);

  bn_set_i64(&q.den, 23534563456);
  bn_set_i64(&q.num, 3);

  br_print(&r);

  br_mul(&p, &r, &q);
  br_mul(&r, &p, &q);

  br_print(&p);
  br_print(&r);

  // r = 1 / 3
  bn_set_i64(&r.den, 4);
  bn_set_i64(&r.num, -1);

  // q = 1 / 4
  bn_set_i64(&q.den, 4);
  bn_set_i64(&q.num, 1);

  br_div(&p, &r, &q);

  br_print(&p);

  br_add(&p, &r, &q);

  br_print(&p);

  /*
  for (u64 x = 0; x < 100; x++) {
    bn_gen_random(&r.den, 100);
    bn_gen_random(&r.num, 100);
    bn_gen_random(&q.den, 100);
    bn_gen_random(&q.num, 100);

    br_div(&k, &r, &q);

        br_print(&k);
    br_div(&l, &q, &r);

    br_mul(&r, &k, &l);

  //  br_print(&r);
  }
*/
  br_free(&r);
  br_free(&p);
  br_free(&q);
  br_free_multi(&k, &l, NULL);
}
