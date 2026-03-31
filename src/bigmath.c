#include <ctype.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>

#include "../include/bignum.h"

// based on wikipedia implementation
// return the Jacobi symbol
i64 bn_jacobi(bignum* a, bignum* m) {
  if (bn_is_zero(m) || bn_is_even(m)) {
    printf("m must be positive and odd\n");
    return 0;
  }

  bignum A, M, R, one;
  bn_init_multi(&A, &M, &R, &one);

  bn_set_u64(&one, 1);
  bn_copy(&M, m);

  bn_mod(&A, a, &M);
  if (a->is_neg && !bn_is_zero(&A)) {
    bn_sub(&A, &M, &A);
  }

  int t = 0;

  while (!bn_is_zero(&A)) {
    u64 z = bn_cnt_trailing_zeros(&A);

    if (z > 0) {
      bn_rshift(&A, &A, z);

      if (z % 2 == 1) {
        u64 m_val = M.size > 0 ? M.limbs[0] : 0;
        t ^= (m_val & 7);
      }
    }

    u64 a_val = A.size > 0 ? A.limbs[0] : 0;
    u64 m_val = M.size > 0 ? M.limbs[0] : 0;
    t ^= (a_val & m_val & 2);

    bn_mod(&R, &M, &A);
    bn_copy(&M, &A);
    bn_copy(&A, &R);
  }

  i64 res;

  if (bn_cmp(&M, &one) != 0) {
    res = 0;
  } else if ((t + 2) & 4) {
    res = -1;
  } else {
    res = 1;
  }

  bn_free(&A);
  bn_free(&M);
  bn_free(&R);
  bn_free(&one);

  return res;
}

// return the euler euler criterion
// p must be an odd prime and a comprime to p
i64 euler_criterion(bignum* a, bignum* p) {
  bignum exp, r, one, p_minus_one;

  bn_init_multi(&exp, &r, &one, &p_minus_one);

  bn_set_u64(&one, 1);

  // Step 1: Handle the a = 0 mod p case
  bn_mod(&r, a, p);
  if (bn_is_zero(&r)) {
    bn_free(&exp);
    bn_free(&r);
    bn_free(&one);
    bn_free(&p_minus_one);
    return 0;
  }

  // Step 2: Calculate p - 1
  bn_sub(&p_minus_one, p, &one);

  // Step 3: Calculate exp = (p - 1) / 2
  bn_copy(&exp, &p_minus_one);
  bn_rshift1(&exp);  // In-place divide by 2

  // Step 4: Calculate r = a^((p-1)/2) mod p
  bn_mod_exp(&r, a, &exp, p);

  i64 result;
  if (bn_cmp(&r, &one) == 0) {
    result = 1;
  } else if (bn_cmp(&r, &p_minus_one) == 0) {
    // In modular arithmetic, p - 1 is equivalent to -1
    result = -1;
  } else {
    // If we get here, p is definitely not a prime.
    result = -2;
  }

  bn_free(&exp);
  bn_free(&r);
  bn_free(&one);
  bn_free(&p_minus_one);

  return result;
}
/* Inputs:
p, a prime
n, an element of Z / p Z such that solutions to the congruence r2 = n exist;
when this is so we say that n is a quadratic residue mod p.
*/
void tonelli_shanks(bignum* r, bignum* p, bignum* n) {
  bignum p_minus_one, one, Q;
  bn_init_multi(&p_minus_one, &one, &Q);
  bn_copy(&p_minus_one, p);
  bn_sub(&p_minus_one, &p_minus_one, &one);
  bn_copy(&Q, &p_minus_one);
  u64 S = 0;
  while (bn_is_even(&Q)) {
    bn_rshift1(&Q);
    S++;
  }
}
