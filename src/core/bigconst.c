#include "../../include/bignum.h"

bignum BN_ZERO;
bignum BN_ONE;
bignum BN_TWO;

void bn_init_constants()
{
  bn_init(&BN_ZERO);
  bn_set_u64(&BN_ZERO, 0);

  bn_init(&BN_ONE);
  bn_set_u64(&BN_ONE, 1);

  bn_init(&BN_TWO);
  bn_set_u64(&BN_TWO, 2);
}

void bn_free_constants() { bn_free_multi(&BN_ZERO, &BN_ONE, &BN_TWO, NULL); }
