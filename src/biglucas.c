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
                    bignum* n)
{
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

void bn_lucas_solve_mod_old(bignum* u, bignum* v, bignum* p, bignum* q,
                            bignum* qn, bignum* n, bignum* m)
{
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

    bn_mod(v, v, m);
    bn_mod(qn, qn, m);
    return;
  }

  bignum n_half, u_n, v_n, q_n, a, b, tmp;
  bn_init_multi(&n_half, &u_n, &v_n, &tmp, &q_n, &a, &b, NULL);

  bn_divmod_u64(&n_half, n, 2);

  bn_lucas_solve_mod(&u_n, &v_n, p, q, &q_n, &n_half, m);

  // Precompute used values
  bn_mul(&a, &u_n, &v_n);
  bn_mul(&b, &v_n, &v_n);

  bn_mod(&a, &a, m);
  bn_mod(&b, &b, m);

  if (bn_is_even(n)) {
    // U_2n = a
    bn_copy(u, &a);
    // V_2n = b - 2 * Q^k
    bn_copy(v, &b);
    bn_add(v, v, m);
    bn_sub(v, v, &q_n);
    bn_add(v, v, m);
    bn_sub(v, v, &q_n);
    bn_mod(v, v, m);
    // Q^n = (Q^n/2)^2
    bn_mul(qn, &q_n, &q_n);
    bn_mod(qn, qn, m);
  } else {
    bn_mul(&tmp, p, &a);
    bn_add(&tmp, &tmp, &b);
    // fix for modular arethmetic
    if (!bn_is_even(&tmp)) {
      bn_add(&tmp, &tmp, m);
    }
    bn_rshift1(&tmp);

    bn_add(&tmp, &tmp, m);
    bn_sub(u, &tmp, &q_n);
    bn_mod(u, u, m);

    bn_mul(v, p, u);
    bn_mod(v, v, m);

    bn_mul(&tmp, q, &a);
    bn_mod(&tmp, &tmp, m);

    bn_add(v, v, m);
    bn_sub(v, v, &tmp);
    bn_mod(v, v, m);
    bn_add(v, v, m);
    bn_sub(v, v, &tmp);
    bn_mod(v, v, m);

    bn_mul(qn, &q_n, &q_n);
    bn_mod(qn, qn, m);
    bn_mul(qn, qn, q);
    bn_mod(qn, qn, m);
  }

  bn_free(&n_half);
  bn_free(&u_n);
  bn_free(&v_n);
  bn_free(&q_n);
  bn_free(&a);
  bn_free(&b);
  bn_free(&tmp);
}

void bn_lucas_solve_mod(bignum* u, bignum* v, bignum* p, bignum* q, bignum* qn,
                        bignum* n, bignum* m)
{
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

    bn_mod(v, v, m);
    bn_mod(qn, qn, m);
    return;
  }

  bn_mont_ctx ctx;
  bn_mont_ctx_init(&ctx, m);
  bignum u_bar, v_bar, qn_bar, p_bar, q_bar, one_bar;
  bn_init_multi(&u_bar, &v_bar, &qn_bar, &q_bar, &p_bar, &one_bar, NULL);

  // Transform constants into Montgomery space
  bn_mont_in(&p_bar, p, &ctx);
  bn_mont_in(&q_bar, q, &ctx);
  bn_copy(&one_bar, &ctx.one_mont);

  // Initial state for n=1: U=1, V=P, Q^n=Q
  bn_copy(&u_bar, &one_bar);
  bn_copy(&v_bar, &p_bar);
  bn_copy(&qn_bar, &q_bar);

  bignum a, b, tmp;
  bn_init_multi(&a, &b, &tmp, NULL);

  u64 len = bn_bit_length(n);

  for (i64 i = len - 2; i >= 0; i--) {
    // Precompute used values
    bn_mont_mul(&a, &u_bar, &v_bar, &ctx);

    bn_mont_mul(&b, &v_bar, &v_bar, &ctx);

    if (bn_get_bit(n, i) == 0) {
      // U_2n = a
      bn_copy(&u_bar, &a);
      // V_2n = b - 2 * Q^k
      bn_sub(&v_bar, &b, &qn_bar);
      if (v_bar.is_neg) bn_add(&v_bar, &v_bar, &ctx.n);
      bn_sub(&v_bar, &v_bar, &qn_bar);
      if (v_bar.is_neg) bn_add(&v_bar, &v_bar, &ctx.n);
      // Q^n = (Q^n/2)^2
      bn_mont_mul(&qn_bar, &qn_bar, &qn_bar, &ctx);
    } else {
      bn_mont_mul(&tmp, &p_bar, &a, &ctx);
      bn_add(&tmp, &tmp, &b);
      // fix for modular arethmetic
      if (!bn_is_even(&tmp)) {
        bn_add(&tmp, &tmp, &ctx.n);
      }
      bn_rshift1(&tmp);

      // bn_sub(&u_bar, &tmp, &qn_bar);
      bn_copy(&u_bar, &tmp);
      if (u_bar.is_neg) bn_add(&u_bar, &u_bar, &ctx.n);

      // v = P*u - Q*a
      bn_mont_mul(&v_bar, &p_bar, &u_bar, &ctx);
      bn_mont_mul(&tmp, &q_bar, &a, &ctx);

      // Subtract Q*a twice (matches your original logic)
      for (int j = 0; j < 2; j++) {
        bn_sub(&v_bar, &v_bar, &tmp);
        if (v_bar.is_neg) bn_add(&v_bar, &v_bar, &ctx.n);
      }

      // Qn = Qn^2 * Q
      bn_mont_mul(&qn_bar, &qn_bar, &qn_bar, &ctx);
      bn_mont_mul(&qn_bar, &qn_bar, &q_bar, &ctx);
    }
  }

  bn_mont_out(u, &u_bar, &ctx);
  bn_mont_out(v, &v_bar, &ctx);
  bn_mont_out(qn, &qn_bar, &ctx);

  bn_free(&a);
  bn_free(&b);
  bn_free(&tmp);
  bn_free(&u_bar);
  bn_free(&v_bar);
  bn_free(&qn_bar);
  bn_free(&p_bar);
  bn_free(&q_bar);
  bn_free(&one_bar);
  bn_mont_ctx_free(&ctx);
}

void bn_lucas(bignum* u, bignum* v, bignum* p, bignum* q, bignum* n)
{
  bignum qn;

  bn_init(&qn);
  bn_lucas_solve(u, v, p, q, &qn, n);
  bn_free(&qn);
}

void bn_lucas_mod(bignum* u, bignum* v, bignum* p, bignum* q, bignum* n,
                  bignum* m)
{
  bignum qn;

  bn_init(&qn);
  bn_lucas_solve_mod(u, v, p, q, &qn, n, m);
  bn_free(&qn);
}
