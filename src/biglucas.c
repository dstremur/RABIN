#include <assert.h>
#include <ctype.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>

#include "../include/bignum.h"
#include "alloca.h"

/*
 Lucas sequence relations

U_2n = U_n V_n
V_2n = V_n^2 - 2Q^n
U_2n+1 = U_n+1V_n - Q^n
V_2n+1 = V_n+1V_n - P(Q)^n

*/

void bn_lucas_solve(bignum* u, bignum* v, bignum* p, bignum* q, bignum* qn,
                    bignum* n) {
  if (bn_is_eq_i64(n, 0)) {
    bn_set_u64(u, 0);
    bn_set_u64(v, 2);
    bn_set_u64(qn, 1);
    return;
  }

  if (bn_is_eq_i64(n, 1)) {
    bn_set_u64(u, 1);
    bn_copy(v, p);
    bn_copy(qn, q);
    return;
  }

  bignum n_half, u_n, v_n, q_n, a, b, tmp;
  bn_init_multi(&n_half, &u_n, &v_n, &tmp, &q_n, &a, &b, NULL);

  bn_divmod_u64(&n_half, n, 2);

  bn_lucas_solve(&u_n, &v_n, p, q, &q_n, &n_half);

  // Precompute used values
  bn_mul(&a, &u_n, &v_n);
  bn_mul(&b, &v_n, &v_n);

  if (bn_is_even(n)) {
    // U_2n = a
    bn_copy(u, &a);
    // V_2n = b - 2 * Q^k
    bn_copy(v, &b);
    bn_sub(v, v, &q_n);
    bn_sub(v, v, &q_n);
    // Q^n = (Q^n/2)^2
    bn_mul(qn, &q_n, &q_n);
  } else {
    bn_mul(&tmp, p, &a);
    bn_add(&tmp, &tmp, &b);
    bn_rshift1(&tmp);
    bn_sub(u, &tmp, &q_n);

    bn_mul(v, p, u);
    bn_mul(&tmp, q, &a);
    bn_sub(v, v, &tmp);
    bn_sub(v, v, &tmp);

    bn_mul(qn, &q_n, &q_n);
    bn_mul(qn, qn, q);
  }

  bn_free(&n_half);
  bn_free(&u_n);
  bn_free(&v_n);
  bn_free(&q_n);
  bn_free(&a);
  bn_free(&b);
  bn_free(&tmp);
}

void bn_lucas(bignum* u, bignum* v, bignum* p, bignum* q, bignum* n) {
  bignum qn;

  bn_init(&qn);
  bn_lucas_solve(u, v, p, q, &qn, n);
  bn_free(&qn);
}
