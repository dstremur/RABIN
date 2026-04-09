#include "../include/bignum.h"
#include "stdio.h"
bool bn_pollard_rho(bignum* f, const bignum* n)
{
  if (bn_is_even(n)) {
    bn_set_u64(f, 2);
    return true;
  }

  bignum a, b, d, tmp;
  bn_init_multi(&a, &b, &d, &tmp, NULL);
  bn_set_u64(&a, 2);
  bn_set_u64(&b, 2);

  for (u64 i = 0; i < 100000; i++) {
    bn_mul(&a, &a, &a);
    bn_add_u64(&a, &a, 1);
    bn_mod(&a, &a, n);

    bn_mul(&b, &b, &b);
    bn_add_u64(&b, &b, 1);
    bn_mod(&b, &b, n);

    bn_mul(&b, &b, &b);
    bn_add_u64(&b, &b, 1);
    bn_mod(&b, &b, n);

    if (bn_cmp(&a, &b) >= 0) {
      bn_sub(&tmp, &a, &b);
    } else {
      bn_sub(&tmp, &b, &a);
    }

    bn_gcd(&d, &tmp, n);

    bn_set_u64(&tmp, 1);
    if (bn_cmp(&tmp, &d) == -1 && bn_cmp(&d, n) == -1) {
      bn_copy(f, &d);
      bn_free_multi(&a, &b, &d, &tmp, NULL);
      return true;
    }

    if (bn_cmp(&d, n) == 0) {
      printf("Fail\n");
      bn_free_multi(&a, &b, &d, &tmp, NULL);
      return false;
    }
  }

  bn_free_multi(&a, &b, &d, &tmp, NULL);

  return false;
}

void bn_pollard_p_minus_one(bignum* f, bignum* n, bignum* B) {}
