#include "../../include/bignum.h"

void bn_isqrt(bignum* r, bignum* a)
{
  if (a->is_neg) return;

  if (bn_is_zero(a)) {
    bn_set_u64(r, 0);
    return;
  }

  bignum xn, xnext, tmp;
  bn_init_multi(&xn, &xnext, &tmp, NULL);

  bn_copy(&xn, a);

  while (true) {
    bn_div(&tmp, a, &xn);

    bn_add(&xnext, &xn, &tmp);

    bn_rshift1(&xnext);
    if (bn_cmp(&xnext, &xn) >= 0) {
      bn_copy(r, &xn);
      break;
    }

    bn_copy(&xn, &xnext);
  }

  bn_free_multi(&xn, &xnext, &tmp, NULL);
}
