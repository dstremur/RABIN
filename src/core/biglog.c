#include "../../include/bignum.h"

// calculates the integer natural logarithm of the largest power of 2 below a
void bn_ln(bignum* r, bignum* a)
{
  if (bn_is_zero(a)) {
    return;
  }

  // ln(2) = 0.6931471805599453094172321214

  bignum log2, tmp, tmp2;
  bn_init_multi(&log2, &tmp, &tmp2, NULL);

  bn_log_2(&log2, a);

  bn_set_u64(&tmp, 69314718);
  bn_set_u64(&tmp2, 100000000);

  // r = (k * 69314718) / 100000000
  bn_mul(r, &log2, &tmp);
  bn_div(r, r, &tmp2);

  bn_free_multi(&log2, &tmp, &tmp2, NULL);
}

// calculates the integer logarithm base 2 of a
void bn_log_2(bignum* r, bignum* a)
{
  if (bn_is_zero(a)) {
    return;
  }

  u64 k = bn_bit_length(a) - 1;

  bn_set_u64(r, k);
}
