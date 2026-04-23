
#include "../../include/bigvector.h"

#include <stdio.h>

#include "../../include/bignum.h"
void bigvector_init(bigvector* a, u64 d)
{
  a->size = d;

  a->data = malloc(d * sizeof(bignum));

  for (u64 i = 0; i < d; i++) {
    bn_init(&a->data[i]);
  }
}

void bigvector_free(bigvector* a)
{
  for (u64 i = 0; i < a->size; i++) {
    bn_free(&a->data[i]);
  }

  free(a->data);
  a->size = 0;
}

void bigvector_set(bigvector* a, bignum* v, u64 i)
{
  if (a->size < i) return;

  bn_copy(&a->data[i], v);
}
void bigvector_add(bigvector* r, const bigvector* a, const bigvector* b)
{
  if (a->size != b->size) return;

  for (u64 i = 0; i < a->size; i++) {
    bn_add(&r->data[i], &a->data[i], &b->data[i]);
  }
}

void bigvector_sub(bigvector* r, const bigvector* a, const bigvector* b)
{
  if (a->size != b->size) return;

  for (u64 i = 0; i < a->size; i++) {
    bn_sub(&r->data[i], &a->data[i], &b->data[i]);
  }
}

void bigvector_dot(bignum* r, const bigvector* a, const bigvector* b)
{
  if (a->size != b->size) return;
  bignum temp;
  bn_init(&temp);
  bn_init(r);

  for (u64 i = 0; i < a->size; i++) {
    bn_mul(&temp, &a->data[i], &b->data[i]);
    bn_add(r, r, &temp);
  }

  bn_free(&temp);
}

void bigvector_print(bigvector* a)
{
  printf("[ ");

  for (u64 i = 0; i < a->size; i++) {
    bn_print(&a->data[i]);
    printf(", ");
  }

  printf("]");
}

void bigvector_println(bigvector* a)
{
  bigvector_print(a);
  printf("\n");
}
