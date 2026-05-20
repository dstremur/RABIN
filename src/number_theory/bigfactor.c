#include "../../include/bignum.h"
#include "../../include/bigvector.h"
#include "../../include/primes.h"
#include "stdio.h"

bool trialdiv(bignum* n, u64 g)
{
  u64 i = 0;

  while (i < 50000 && primes[i] < g) {
    if (bn_mod_u64(n, primes[i]) == 0) {
      return false;
    }
    i++;
  }

  return true;
}

bool trialdiv_factor(bignum* f, bignum* n, u64 g)
{
  u64 i = 0;

  while (i < 50000 && primes[i] < g) {
    if (bn_mod_u64(n, primes[i]) == 0) {
      bn_set_u64(f, primes[i]);
      return false;
    }
    i++;
  }

  return true;
}

// use randomized pollard rho
bool bn_pollard_rho_inner(bignum* f, const bignum* n)
{
  if (bn_is_even(n)) {
    bn_set_u64(f, 2);
    return true;
  }

  bignum a, b, c, d, tmp;
  bn_init_multi(&a, &b, &c, &d, &tmp, NULL);

  bn_set_u64(&a, 2);
  bn_set_u64(&b, 2);
  bn_gen_random_range(&c, &BN_TWO, n);

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
      bn_free_multi(&a, &b, &c, &d, &tmp, NULL);
      return false;
    }
  }
  bn_free_multi(&a, &b, &c, &d, &tmp, NULL);

  return false;
}

bool bn_pollard_rho(bignum* f, const bignum* n)
{
  if (bn_is_even(n)) {
    bn_set_u64(f, 2);
    return true;
  }

  // try several times
  for (u64 i = 0; i < 64; i++) {
    if (bn_pollard_rho_inner(f, n)) {
      return true;
    }
  }

  return false;
}

static void bigvector_append_distinct(bigvector* v, const bignum* f)
{
  // 1. Compare values deeply using bn_cmp (never compare raw struct memory
  // blocks)
  if (v->size > 0) {
    if (bn_cmp(&v->data[v->size - 1], f) == 0) {
      return;  // Element matches the last factor found, exit out to stay
               // distinct
    }
  }

  // 2. Safely hand off to your append function, which executes its own safe
  // bn_copy
  bigvector_append(v, (bignum*)f);
}

// completely factorize n using recursive applications of pollard_rho, and
// stores in array factors
void bn_factorize(bigvector* v, bignum* n)
{
  if (bn_cmp(n, &BN_ONE) == 0 || bn_is_zero(n)) {
    return;
  }

  if (bn_is_eq_i64(n, 2) || bn_is_eq_i64(n, 3)) {
    bigvector_append_distinct(v, n);
    return;
  }

  if (bn_bpsw(n)) {
    bigvector_append_distinct(v, n);
    return;
  }

  bignum tmp, f, n1, rem;
  bn_init_multi(&tmp, &f, &n1, &rem, NULL);
  bn_copy(&tmp, n);

  // even case
  if (bn_is_even(&tmp)) {
    bn_set_u64(&f, 2);
    bigvector_append_distinct(v, &f);

    while (bn_is_even(&tmp)) {
      bn_div(&tmp, &tmp, &f);
    }

    bn_factorize(v, &tmp);
    goto cleanup;
  }

  // first trialdiv
  if (!trialdiv_factor(&f, &tmp, 10000)) {
    bigvector_append_distinct(v, &f);

    bn_div(&tmp, &tmp, &f);
    bn_mod(&rem, &tmp, &f);
    while (bn_is_zero(&rem)) {
      bn_div(&tmp, &tmp, &f);
      bn_mod(&rem, &tmp, &f);
    }

    bn_factorize(v, &tmp);
    goto cleanup;
  }

  bool res = false;

  // pollard rho
  for (u64 i = 0; i < 16 && !res; i++) {
    if (bn_pollard_rho(&f, &tmp)) {
      res = true;
    }
  }

  // p - 1
  if (!res) {
    printf("trying p-1 \n");
    if (bn_pollard_p_minus_one(&f, &tmp)) {
      res = true;
    }
  }

  if (!res) {
    printf("failed to split composite branch completely!\n");
    goto cleanup;
  }

  // bn_println(&f);

  if (trialdiv(&f, 10000) && (&f)) {
    bigvector_append_distinct(v, &f);

    bn_div(&tmp, &tmp, &f);
    bn_mod(&rem, &tmp, &f);
    while (bn_is_zero(&rem)) {
      bn_div(&tmp, &tmp, &f);
      bn_mod(&rem, &tmp, &f);
    }

    bn_factorize(v, &tmp);
  } else {
    bn_div(&n1, &tmp, &f);
    bn_factorize(v, &f);
    bn_factorize(v, &n1);
  }

cleanup:
  bn_free_multi(&tmp, &f, &n1, &rem, NULL);
}

bool bn_pollard_p_minus_one_stage_1(bignum* f, bignum* n, u64 B, u64 iterations)
{
  bignum M, a, tmp, n_min1, a_min_1, q, ln_q, ln_n, l;
  bn_init_multi(&M, &a, &tmp, &n_min1, &a_min_1, &q, &ln_q, &ln_n, &l, NULL);

  bool res = false;
  bn_sub(&n_min1, n, &BN_ONE);

  for (u64 attempt = 0; attempt < iterations; attempt++) {
    // set a to random number in 2 <= a <= a - 1
    bn_gen_random_range(&a, &BN_TWO, &n_min1);

    bn_gcd(&tmp, &a, n);

    // if gcd(a,n) >= 2 found a factor
    if (bn_cmp(&tmp, &BN_ONE) > 0 && bn_cmp(&tmp, n) < 0) {
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
    bn_sub(&a_min_1, &a, &BN_ONE);

    bn_gcd(&tmp, &a_min_1, n);

    if (bn_cmp(&tmp, &BN_ONE) > 0 && bn_cmp(&tmp, n) < 0) {
      bn_copy(f, &tmp);
      res = true;
      goto cleanup;
    }
  }

cleanup:
  bn_free_multi(&M, &a, &tmp, &n_min1, &a_min_1, &q, &ln_q, &ln_n, &l, NULL);
  return res;
}

/*
 1. B = 1000, 3 rounds

 2. B = 10000, 6 rounds

 3. B = 100000, 10 rounds
*/

bool bn_pollard_p_minus_one(bignum* f, bignum* n)
{
  if (bn_pollard_p_minus_one_stage_1(f, n, 1000, 1000)) {
    return true;
  }

  if (bn_pollard_p_minus_one_stage_1(f, n, 10000, 100)) {
    return true;
  }

  return bn_pollard_p_minus_one_stage_1(f, n, 50000, 100);
}
