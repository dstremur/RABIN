
#include "../../include/bigvector.h"

#include "../../include/bignum.h"

void bigvector_init(bigvector* a)
{
  a->data = NULL;
  a->size = 0;
}

void bigvector_free(bigvector* a)
{
  for (u64 i = 0; i < a->size; i++) {
    bn_free(&a->data[i]);
  }

  free(a->data);
  a->size = 0;
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
