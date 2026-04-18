#include "../include/bigpoly.h"

#include <inttypes.h>
#include <stdio.h>

void bigpoly_init(bigpoly* p)
{
  p->coeff = NULL;
  p->deg = 0;
  p->size = 0;
}

// sets the value of a polynomial, takes an inputs coefficients array
void bigpoly_set(bigpoly* p, bignum* coeff, u64 deg)
{
  bigpoly_alloc(p, deg + 1);

  for (u64 i = 0; i < deg; i++) {
    bn_copy(&p->coeff[i], &coeff[i]);
  }

  p->deg = deg;
  p->size = deg;
}

void bigpoly_set_i64(bigpoly* p, i64* coeff, u64 deg)
{
  bigpoly_alloc(p, deg + 1);

  for (u64 i = 0; i <= deg; i++) {
    bn_set_i64(&p->coeff[i], coeff[i]);
  }

  p->deg = deg;
  bigpoly_trim(p);
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
  p->deg = 0;
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

void bigpoly_print(bigpoly* p)
{
  for (i64 i = p->deg; i >= 0; i--) {
    bn_print(&p->coeff[i]);
    printf("x^%" PRId64 "\n", i);
  }
}

// adds two polynomials
void bigpoly_add(bigpoly* r, bigpoly* p, bigpoly* q)
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

  bigpoly* longer = (p->deg > q->deg) ? p : q;
  for (u64 i = min + 1; i <= max; i++) {
    bn_copy(&r->coeff[i], &longer->coeff[i]);
  }

  r->deg = max;
  bigpoly_trim(r);
}

// p - q = r
void bigpoly_sub(bigpoly* r, bigpoly* p, bigpoly* q)
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

  bigpoly* longer = (p->deg > q->deg) ? p : q;
  for (u64 i = min + 1; i <= max; i++) {
    bn_copy(&r->coeff[i], &longer->coeff[i]);
    r->coeff[i].is_neg = true;
  }

  r->deg = max;
  bigpoly_trim(r);
}
// multiplies two polynomials
void bigpoly_mul_school(bigpoly* r, bigpoly* p, bigpoly* q)
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

  // Cleanup
  bigpoly_free(&p);
  bigpoly_free(&q);
  bigpoly_free(&r);

  printf("\nTests complete.\n");
}
