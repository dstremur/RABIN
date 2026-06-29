
#include "bignum.h"

void bigvector_init(bigvector* a, u64 d);
void bigvector_init_dynamic(bigvector* a);

// only for dynamic vectors
void bigvector_append(bigvector* v, bignum* a);
void bigvector_remove(bigvector* v, u64 index);

void bigvector_set(bigvector* a, bignum* v, u64 i);
void bigvector_free(bigvector* a);

void bigvector_add(bigvector* r, const bigvector* a, const bigvector* b);

void bigvector_sub(bigvector* r, const bigvector* a, const bigvector* b);
void bigvector_dot(bignum* r, const bigvector* a, const bigvector* b);
void bigvector_cross(bigvector* r, const bigvector* a, const bigvector* b);
void bigvector_norm(bignum* r, const bigvector* a);

void bigvector_print(bigvector* a);
void bigvector_println(bigvector* a);
