#include <ctype.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>

#include "../include/bignum.h"

i64 bn_jacobi(bignum* a, bignum* m) {
  if (bn_is_zero(m) || bn_is_even(m)) {
    return 0;
  }

  bignum temp_a, temp_m;
  bn_init_multi(&temp_a, &temp_m);
  bn_mod(&temp_a, a, m);
  bn_copy(&temp_m, m);

  int t = 1;

  while (!bn_is_zero(&temp_a)) {
    u64 k = bn_cnt_trailing_zeros(&temp_a);
    bn_rshift(&temp_a, &temp_a, k);

    if (k & 1) {
      u64 m_mod8 = temp_m.limbs[0] & 7;
      if (m_mod8 == 3 || m_mod8 == 5) {
        t = -t;
      }
    }

    // 3. Quadratic Reciprocity swap: (a/m) -> (m/a) * (-1 if both % 4 == 3)
    u64 a_mod4 = temp_a.limbs[0] & 3;
    u64 m_mod4 = temp_m.limbs[0] & 3;
    if (a_mod4 == 3 && m_mod4 == 3) {
      t = -t;
    }

    // Swap and Mod
    // Reuse temp_a as a placeholder to avoid more inits
    bignum swap_tmp = temp_a;
    temp_a = temp_m;
    temp_m = swap_tmp;

    bn_mod(&temp_a, &temp_a, &temp_m);
  }

  int result = t;
  bignum one;
  bn_init(&one);
  bn_set_u64(&one, 1);
  if (bn_cmp(&temp_m, &one) == 0) {
    result = 0;
  }

  bn_free(&temp_a);
  bn_free(&temp_m);
  return result;
}
