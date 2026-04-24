#include "../../include/bignum.h"

// calculates the integer square root of a using Heron's method
void bn_isqrt_heron(bignum* r, bignum* a)
{
  if (a->is_neg) return;

  if (bn_is_zero(a)) {
    bn_set_u64(r, 0);
    return;
  }

  bignum xn, xnext, tmp;
  bn_init_multi(&xn, &xnext, &tmp, NULL);

  // set x0 = 2^{log_2(a) / 2 + 1}
  u64 k = bn_bit_length(a);
  bn_set_u64(&xn, 2);

  bn_set_u64(&tmp, k + 2);
  bn_rshift1(&tmp);
  bn_pow(&xn, &xn, &tmp);

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

void bn_isqrt(bignum* r, bignum* a) { bn_isqrt_heron(r, a); }
