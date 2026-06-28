#include "../../include/bigpoly.h"

#include <inttypes.h>
#include <stdio.h>

#include "../../include/bigntt.h"
#include "../../include/u64.h"

void bigpoly_init(bigpoly* p)
{
  p->coeff = NULL;
  p->deg = 0;
  p->size = 0;
}

void bn_and(bignum* r, const bignum* a, const bignum* mask)
{
  bn_copy(r, a);

  u64 min = MIN(r->size, mask->size);

  for (u64 i = 0; i < min; i++) {
    r->limbs[i] &= mask->limbs[i];
  }

  // zero the rest
  for (u64 i = min; i < r->size; i++) {
    r->limbs[i] = 0;
  }

  bn_trim(r);
}

// sets the value of a polynomial, takes an inputs coefficients array
void bigpoly_set(bigpoly* p, bignum* coeff, u64 deg)
{
  bigpoly_alloc(p, deg + 1);

  for (u64 i = 0; i <= deg; i++) {
    bn_copy(&p->coeff[i], &coeff[i]);
  }

  p->deg = deg;
}

void bigpoly_set_i64(bigpoly* p, i64* coeff, u64 deg)
{
  bigpoly_alloc(p, deg);

  for (u64 i = 0; i <= deg; i++) {
    bn_set_i64(&p->coeff[i], coeff[i]);
  }

  p->deg = deg;
  bigpoly_trim(p);
}

void bigpoly_copy(bigpoly* p, bigpoly* q)
{
  for (u64 i = 0; i <= q->deg; i++) {
    bn_copy(&p->coeff[i], &q->coeff[i]);
  }
}

bool bigpoly_equal(bigpoly* a, bigpoly* b)
{
  if (a->deg != b->deg) return false;
  for (u64 i = 0; i <= a->deg; i++) {
    if (bn_cmp(&a->coeff[i], &b->coeff[i]) != 0) {
      return false;
    }
  }
  return true;
}

// Frees the polynomial struct
void bigpoly_free(bigpoly* p)
{
  if (!p->coeff) return;

  for (u64 i = 0; i < p->size; i++) {
    bn_free(&p->coeff[i]);
  }

  free(p->coeff);
  p->coeff = NULL;
  p->deg = 0;
  p->size = 0;
}

bool bigpoly_alloc(bigpoly* p, u64 deg)
{
  if (p->size >= deg + 1) return true;

  u64 new_cap = (p->size == 0) ? deg + 1 : (p->size * 2);
  if (new_cap < deg + 1) {
    new_cap = deg + 1;
  }

  bignum* new = realloc(p->coeff, new_cap * sizeof(bignum));
  if (!new) return false;

  p->coeff = new;

  for (u64 i = p->size; i < new_cap; i++) {
    bn_init(&p->coeff[i]);
  }

  p->size = new_cap;
  return true;
}
// removes excess zero terms
void bigpoly_trim(bigpoly* p)
{
  while (p->deg > 0 && bn_is_zero(&p->coeff[p->deg])) {
    p->deg--;
  }
}

void bigpoly_print(const bigpoly* p)
{
  for (i64 i = p->deg; i >= 0; i--) {
    bn_print(&p->coeff[i]);
    printf("x^%" PRId64 "", i);
    if (i != 0) {
      printf("+");
    }
  }
  printf("\n");
}

// adds two polynomials
void bigpoly_add(bigpoly* r, const bigpoly* p, const bigpoly* q)
{
  u64 min = MIN(p->deg, q->deg);
  u64 max = MAX(p->deg, q->deg);

  if (!bigpoly_alloc(r, max)) {
    // TODO: error handling
    return;
  }

  // add upto the minimum degree
  for (u64 i = 0; i <= min; i++) {
    bn_add(&r->coeff[i], &p->coeff[i], &q->coeff[i]);
  }

  const bigpoly* longer = (p->deg > q->deg) ? p : q;
  for (u64 i = min + 1; i <= max; i++) {
    bn_copy(&r->coeff[i], &longer->coeff[i]);
  }

  r->deg = max;
  bigpoly_trim(r);
}

// p - q = r
void bigpoly_sub(bigpoly* r, const bigpoly* p, const bigpoly* q)
{
  u64 min = MIN(p->deg, q->deg);
  u64 max = MAX(p->deg, q->deg);

  if (!bigpoly_alloc(r, max)) {
    // TODO: error handling
    return;
  }

  // add upto the minimum degree
  for (u64 i = 0; i <= min; i++) {
    bn_sub(&r->coeff[i], &p->coeff[i], &q->coeff[i]);
  }

  const bigpoly* longer = (p->deg > q->deg) ? p : q;
  for (u64 i = min + 1; i <= max; i++) {
    bn_copy(&r->coeff[i], &longer->coeff[i]);
    r->coeff[i].is_neg = true;
  }

  r->deg = max;
  bigpoly_trim(r);
}

void bigpoly_mul_digit(bigpoly* r, const bigpoly* p, const bigpoly* q,
                       const bignum* m)
{
  u64 max = MAX(p->deg, q->deg);

  bigpoly_alloc(r, max);

  for (u64 i = 0; i <= max; i++) {
    bn_mul(&r->coeff[i], &p->coeff[i], &q->coeff[i]);
    bn_mod(&r->coeff[i], &r->coeff[i], m);
  }

  r->deg = max;
}

// multiply two polynomials using ntt only works for coeffcient below u64
void bigpoly_mul_ntt(bigpoly* r, const bigpoly* p, const bigpoly* q)
{
  ntt_ctx_u64 ctx;

  u64 required_len = p->deg + q->deg + 1;
  u64 ntt_size = 1;
  u64 k = 0;

  // Pad to next power of 2
  while (ntt_size < required_len) {
    ntt_size <<= 1;
    k++;
  }

  if (!ntt_ctx_u64_init_golden(&ctx, k)) {
    return;
  }

  // 1. Allocate flat u64 arrays for the inputs and outputs
  u64* p_arr = calloc(ntt_size, sizeof(u64));
  u64* q_arr = calloc(ntt_size, sizeof(u64));
  u64* p_hat = calloc(ntt_size, sizeof(u64));
  u64* q_hat = calloc(ntt_size, sizeof(u64));
  u64* r_hat = calloc(ntt_size, sizeof(u64));
  u64* r_arr = calloc(ntt_size, sizeof(u64));

  // 2. Extract u64 values from the bignum polynomials
  for (u64 i = 0; i <= p->deg; i++) {
    p_arr[i] = p->coeff[i].limbs[0];
  }
  for (u64 i = 0; i <= q->deg; i++) {
    q_arr[i] = q->coeff[i].limbs[0];
  }

  // 3. Perform Forward NTTs
  ntt_u64_cyclic_forward(p_hat, p_arr, &ctx);
  ntt_u64_cyclic_forward(q_hat, q_arr, &ctx);

  // 4. Pointwise Multiplication
  for (u64 i = 0; i < ntt_size; i++) {
    // Both p_hat and q_hat are currently in Montgomery form.
    // mont_mul computes (P * Q * R^-1), yielding the product still in
    // Montgomery form.
    u64 r_mont = mont_mul(p_hat[i], q_hat[i], &ctx.mctx);

    // ntt_u64_cyclic_inverse expects inputs in standard form (it calls mont_in
    // internally). Therefore, we must convert r_mont OUT of Montgomery space
    // before passing it.
    r_hat[i] = mont_out(r_mont, &ctx.mctx);
  }

  // 5. Perform Inverse NTT
  ntt_u64_cyclic_inverse(r_arr, r_hat, &ctx);

  // 6. Pack the result back into the bigpoly structure
  // Note: Your original code freed 'r' before init. Ensure bigpoly_free is safe
  // on uninitialized 'r'.
  bigpoly_free(r);
  bigpoly_init(r);
  bigpoly_alloc(r, required_len);

  for (u64 i = 0; i < required_len; i++) {
    bn_set_u64(&r->coeff[i], r_arr[i]);
  }

  r->deg = required_len - 1;
  bigpoly_trim(r);

  // 7. Cleanup
  free(p_arr);
  free(q_arr);
  free(p_hat);
  free(q_hat);
  free(r_hat);
  free(r_arr);
  ntt_ctx_u64_free(&ctx);
}

// multiplies two polynomials
void bigpoly_mul_school(bigpoly* r, const bigpoly* p, const bigpoly* q)
{
  u64 r_deg = p->deg + q->deg;

  bigpoly temp;
  bigpoly_init(&temp);
  bigpoly_alloc(&temp, r_deg);

  temp.deg = r_deg;

  bignum t;
  bn_init(&t);

  for (u64 i = 0; i <= p->deg; i++) {
    for (u64 j = 0; j <= q->deg; j++) {
      bn_mul(&t, &p->coeff[i], &q->coeff[j]);

      bn_add(&temp.coeff[i + j], &temp.coeff[i + j], &t);
    }
  }

  bn_free(&t);

  bigpoly_trim(&temp);

  if (bigpoly_alloc(r, temp.deg)) {
    for (u64 i = 0; i <= temp.deg; i++) {
      bn_copy(&r->coeff[i], &temp.coeff[i]);
    }
    r->deg = temp.deg;
  }

  bigpoly_free(&temp);
}

void bigpoly_mul(bigpoly* r, const bigpoly* p, const bigpoly* q)
{
  bigpoly_mul_school(r, p, q);
}

/**
 * @brief Slices a bignum into chunks and stores them as polynomial
 * coefficients.
 * @param p The output polynomial.
 * @param n The source bignum.
 * @param bit_width How many bits of 'n' to store in each coefficient (e.g.,
 * 16).
 */
void bn_decompose(bigpoly* r, const bignum* n, u64 width)
{
  bignum tmp, mask, digit;
  bn_init_multi(&tmp, &mask, &digit);
  bn_copy(&tmp, n);

  // mask = 0xFFFF..
  bn_lshift(&mask, &BN_ONE, width);
  bn_sub(&mask, &mask, &BN_ONE);

  u64 i = 0;
  while (!bn_is_zero(&tmp)) {
    bn_and(&digit, &tmp, &mask);

    bigpoly_alloc(r, i + 1);
    bn_copy(&r->coeff[i], &digit);

    bn_rshift(&tmp, &tmp, width);

    i++;
  }

  r->deg = (i > 0) ? i - 1 : 0;

  bn_free_multi(&tmp, &mask, &digit);
}

void poly_carry_propagation(bigpoly* r, u64 bit_width)
{
  bignum carry, base, mask, total;
  bn_init_multi(&carry, &base, &mask, &total, NULL);

  bn_lshift(&base, &BN_ONE, bit_width);
  bn_sub(&mask, &base, &BN_ONE);
  bn_set_u64(&carry, 0);

  u64 i = 0;
  // Iterate through all coefficients plus any remaining carries
  while (i <= r->deg || !bn_is_zero(&carry)) {
    if (i > r->deg) {
      bigpoly_alloc(r, i + 1);  // Expand poly if carry exceeds current deg
      r->deg = i;
    }

    // total = coeff[i] + carry
    bn_add(&total, &r->coeff[i], &carry);

    // carry = total >> bit_width
    bn_rshift(&carry, &total, bit_width);

    // coeff[i] = total & mask
    bn_and(&r->coeff[i], &total, &mask);

    i++;
  }

  bn_free_multi(&carry, &base, &mask, &total, NULL);
}

void bn_recompose(bignum* n, const bigpoly* p, u64 bit_width)
{
  bn_set_u64(n, 0);
  bignum term;
  bn_init(&term);

  for (u64 i = 0; i <= p->deg; i++) {
    // term = coeff[i] << (i * bit_width)
    bn_lshift(&term, &p->coeff[i], i * bit_width);
    // n += term
    bn_add(n, n, &term);
  }

  bn_free(&term);
}

void bigpoly_test()
{
  printf("--- BigPoly Test --- \n");

  bigpoly p, q, r;
  bigpoly_init(&p);
  bigpoly_init(&q);
  bigpoly_init(&r);

  // Let p(x) = 2x^2 + 3x + 1
  i64 p_vals[] = {1, 3, 2};
  bigpoly_set_i64(&p, p_vals, 2);

  // Let q(x) = 4x + 5
  i64 q_vals[] = {5, 4};
  bigpoly_set_i64(&q, q_vals, 1);

  printf("\nPolynomial P(x):\n");
  bigpoly_print(&p);

  printf("\nPolynomial Q(x):\n");
  bigpoly_print(&q);

  // Test Addition: r = p + q = 2x^2 + 7x + 6
  printf("\n--- Test: Addition (P + Q) ---\n");
  bigpoly_add(&r, &p, &q);
  bigpoly_print(&r);

  // Test Multiplication: r = p * q = 8x^3 + 22x^2 + 19x + 5
  printf("\n--- Test: Multiplication (P * Q) ---\n");
  bigpoly_mul_school(&r, &p, &q);
  bigpoly_print(&r);

  bigpoly_free(&r);
  bigpoly_init(&r);

  bigpoly_mul_ntt(&r, &p, &q);
  bigpoly_print(&r);

  // Cleanup
  bigpoly_free(&p);
  bigpoly_free(&q);
  bigpoly_free(&r);

  printf("\nTests complete.\n");
}
