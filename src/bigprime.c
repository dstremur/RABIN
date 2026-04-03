#include "../include/bignum.h"

bool bpsw(bignum* n) {
  // TODO trial division

  // 1. run miller rabin base 2
  bignum two;
  bn_init(&two);
  bn_set_u64(&two, 2);
  if (!bn_rabin(n, &two)) {
    bn_free(&two);
    return false;
  }

  // check if perfect square
  bignum D, magnitude;
  bn_init(&D);
  bn_init(&magnitude);
  bn_set_i64(&magnitude, 5);

  bool negative = false;

  while (1) {
    bn_copy(&D, &magnitude);
    D.is_neg = negative;

    i64 jacobi = bn_jacobi(&D, n);
    if (jacobi == -1) {
      break;
    }

    if (jacobi == 0) {
      bn_free(&D);
      bn_free(&magnitude);
      return false;
    }

    bn_add_u64(&magnitude, &magnitude, 2);
    negative = !negative;
  }

  bn_free(&D);
  bn_free(&magnitude);
  return true;
}
