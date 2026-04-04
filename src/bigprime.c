#include "../include/bignum.h"

bool bn_is_perfect_square(bignum* n)
{
  // 0 is a perfect square
  if (bn_is_zero(n)) return true;

  // Negative numbers cannot be perfect squares
  if (n->is_neg) return false;

  // --- STEP 2: Newton's Method Calculation ---
  bignum x, y, tmp, rem;
  bn_init_multi(&x, &y, &tmp, &rem, NULL);

  // Initial guess: x = 2^(bits/2)
  int bits = bn_bit_length(n);
  bn_set_u64(&x, 1);
  bn_lshift(&x, &x, (bits + 1) / 2);

  // Newton-Raphson iteration
  while (true) {
    // y = (x + n/x) / 2
    bn_div(&tmp, n, &x);   // tmp = n / x
    bn_add(&y, &x, &tmp);  // y = x + n/x
    bn_rshift1(&y);        // y = y / 2

    if (bn_cmp(&y, &x) >= 0) {
      break;
    }
    bn_copy(&x, &y);
  }

  // Check if x * x == n
  bn_mul(&tmp, &x, &x);
  bool is_square = (bn_cmp(&tmp, n) == 0);

  bn_free(&x);
  bn_free(&y);
  bn_free(&tmp);
  bn_free(&rem);

  return is_square;
}

bool bn_stronglucas(bignum* n, bignum* P, bignum* Q)
{
  bignum p, q, d, u, v, qn, tmp, n_plus_1;
  bn_init_multi(&p, &q, &d, &u, &v, &qn, &tmp, &n_plus_1, NULL);

  bool prime = false;

  int s = 0;

  bn_copy(&n_plus_1, n);
  bn_add_u64(&n_plus_1, &n_plus_1, 1);
  bn_copy(&d, &n_plus_1);

  while (bn_is_even(&d)) {
    bn_rshift1(&d);
    s++;
  }

  bn_lucas_solve_mod(&u, &v, P, Q, &qn, &d, n);
  // 1. condition
  if (bn_is_zero(&u)) {
    prime = true;
    goto cleanup;
  }
  // 2. condition
  if (bn_is_zero(&v)) {
    prime = true;
    goto cleanup;
  }
  // check 2. condition for all d * 2^r
  for (int r = 1; r < s; r++) {
    bn_mul(&tmp, &v, &v);
    bn_mod(&tmp, &tmp, n);

    bn_add(&tmp, &tmp, n);
    bn_sub(&tmp, &tmp, &qn);
    bn_add(&tmp, &tmp, n);
    bn_sub(&tmp, &tmp, &qn);
    bn_mod(&v, &tmp, n);

    bn_mul(&qn, &qn, &qn);
    bn_mod(&qn, &qn, n);

    if (bn_is_eq_i64(&v, 0)) {
      prime = true;
      goto cleanup;
    }
  }

cleanup:
  bn_free(&p);
  bn_free(&q);
  bn_free(&d);
  bn_free(&u);
  bn_free(&v);
  bn_free(&qn);
  bn_free(&tmp);
  bn_free(&n_plus_1);
  return prime;
}

bool bn_bpsw(bignum* n)
{
  // TODO trial division
  bool prime = false;
  if (bn_is_perfect_square(n)) {
    return false;
  }
  // 1. run miller rabin base 2
  bignum two;
  bn_init(&two);
  bn_set_u64(&two, 2);
  if (!bn_rabin_mont(n, &two)) {
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

  bignum P, Q;
  bn_init_multi(&P, &Q, NULL);
  bn_set_u64(&P, 1);
  bn_set_u64(&Q, 1);

  if (negative) {
    bn_copy(&Q, &magnitude);
    bn_add_u64(&Q, &Q, 1);
    bn_divmod_u64(&Q, &Q, 4);
  } else {
    bn_copy(&Q, &magnitude);
    bn_sub(&Q, &Q, &P);
    bn_divmod_u64(&Q, &Q, 4);
    bn_sub(&Q, n, &Q);
  }

  if (bn_is_eq_i64(&D, 5)) {
    bn_set_u64(&P, 5);
    bn_set_u64(&Q, 5);
  }

  bool res = bn_stronglucas(n, &P, &Q);

  bn_free(&P);
  bn_free(&Q);
  bn_free(&D);
  bn_free(&magnitude);
  return res;
}
