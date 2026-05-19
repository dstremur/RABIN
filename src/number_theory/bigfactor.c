#include "../../include/bignum.h"
#include "../../include/primes.h"
#include "../../include/bigvector.h" 
#include "stdio.h"

// use randomized pollard rho
bool bn_pollard_rho(bignum* f, const bignum* n)
{
  if (bn_is_even(n)) {
    bn_set_u64(f, 2);
    return true;
  }

  bignum a, b, c, d, tmp;
  bn_init_multi(&a, &b, &c, &d, &tmp, NULL);
  bn_gen_random_range(&a, &BN_TWO, n);
  bn_copy(&b, &a);

  bn_gen_random_range(&c, &BN_ONE, n);


  for (u64 i = 0; i < 100000; i++) {
    bn_mul(&a, &a, &a);
    bn_add(&a, &a, &c);
    bn_mod(&a, &a, n);

    bn_mul(&b, &b, &b);
    bn_add(&b, &b, &c);
    bn_mod(&b, &b, n);

    bn_mul(&b, &b, &b);
    bn_add(&b, &b, &c);
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
      bn_free_multi(&a, &b, &c, &d, &tmp, NULL);
      return true;
    }

    if (bn_cmp(&d, n) == 0) {
      printf("Fail\n");
      bn_free_multi(&a, &b, &c, &d, &tmp, NULL);
      return false;
    }
  }

  bn_free_multi(&a, &b, &c, &d, &tmp, NULL);

  return false;
}

/*
 1. B = 1000, 3 rounds

 2. B = 10000, 6 rounds

 3. B = 100000, 10 rounds
*/
void bn_factorize(bigvector* v, bignum* n)
{	
	if (bn_cmp(n, &BN_ONE) == 0 || bn_is_zero(n)) return;

	if (bn_cmp(n, &BN_TWO) == 0) {
		printf("factor 2\n");
		bigvector_append(v, &BN_TWO);
		return;
	}

	if (!bn_is_even(n) && bn_bpsw(n)){
		bn_println(n);
		bigvector_append(v, n);
		return; 
	}

	bignum f, n1;
	bn_init_multi(&f, &n1, NULL);

	bn_pollard_rho(&f, n);

	if (bn_is_zero(&f) ||
    bn_cmp(&f, &BN_ONE) == 0 ||
    bn_cmp(&f, n) == 0)
{
    fprintf(stderr, "pollard rho failed\n");
		
	bn_pollard_p_minus_one(&f, n);
}

	printf("f: ");
	bn_println(&f);
	bn_div(&n1, n, &f);

	bn_factorize(v, &f);
	bn_factorize(v, &n1);
	
	bn_free_multi(&f, &n1, NULL); 

	

}


// completely factorize n using recursive applications of pollard_rho, and stores in array factors

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
    bn_sub(&n_min1, &a, &BN_ONE);

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
