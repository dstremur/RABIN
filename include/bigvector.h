
#include "bignum.h"

typedef struct {
	bignum* data;
	u64 size;
} bigvector;

void bigvector_init(bigvector* a);

void bigvector_free(bigvector* a);

void bigvector_add(bigvector* r, const bigvector* a, const bigvector* b); 

void bigvector_sub(bigvector* r, const bigvector* a, const bigvector* b); 
void bigvector_dot(bignum* r, const bigvector* a, const bigvector* b); 
void bigvector_cross(bigvector* r, const bigvector* a, const bigvector* b); 
void bigvector_norm(bigvector* r, const bigvector* a);


