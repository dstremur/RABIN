#include "../../include/bignum.h"
#include "../../include/primes.h"
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

/*
 1. B = 1000, 3 rounds

 2. B = 10000, 6 rounds

 3. B = 100000, 10 rounds
*/
bool bn_pollard_p_minus_one_stage_1(bignum* f, bignum* n, u64 B, u64 iterations)
{
  bignum M, a, tmp, n_min1, q, ln_q, ln_n, l;
  bn_init_multi(&M, &a, &tmp, &n_min1, &q, &ln_q, &ln_n, &l, NULL);

  bool res = false;
  bn_sub(&n_min1, n, &BN_ONE);

  for (u64 attempt = 0; attempt < iterations; attempt++) {
    // set a to random number in 2 <= a <= a - 1
    bn_gen_random_range(&a, &tmp, &n_min1);

    bn_gcd(&tmp, &a, n);

    // if gcd(a,n) >= 2 found a factor
    if (bn_cmp(&tmp, &BN_ONE) > 0) {
      bn_copy(f, &tmp);
      res = true;
      goto cleanup;
    }

    int i = 0;

    // go through all primes upto B
    while (primes[i] <= B) {
      u64 q_val = primes[i];
      u64 q_pow = q_val;

      while (q_pow <= B / q_val) {
        q_pow *= q_val;
      }

      bn_set_u64(&tmp, q_pow);
      bn_mod_exp(&a, &a, &tmp, n);

      i++;
    }

    // nmin1 = a - 1
    bn_sub(&n_min1, &a, &M);

    bn_gcd(&tmp, &n_min1, n);

    if (bn_cmp(&tmp, &BN_ONE) > 0 && bn_cmp(&tmp, n) < 0) {
      bn_copy(f, &tmp);
      res = true;
      goto cleanup;
    }
  }

cleanup:
  bn_free_multi(&M, &a, &tmp, &n_min1, &q, &ln_q, &ln_n, &l, NULL);
  return res;
}

bool bn_pollard_p_minus_one(bignum* f, bignum* n)
{
  if (bn_pollard_p_minus_one_stage_1(f, n, 1000, 1000)) {
    return true;
  }

  if (bn_pollard_p_minus_one_stage_1(f, n, 10000, 100)) {
    return true;
  }

  return bn_pollard_p_minus_one_stage_1(f, n, 100000, 100);
}
