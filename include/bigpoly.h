#include "bignum.h" 

typedef struct {
	bignum* coeff;
	u64 size;
	u64 deg;
} bigpoly; 


void bigpoly_init(bigpoly* p); 
void bigpoly_set(bigpoly* p, bignum* coeff, u64 deg);
void bigpoly_print(const bigpoly* p); 
bool bigpoly_alloc(bigpoly* p, u64 deg);
void bigpoly_trim(bigpoly* p);
void bigpoly_mul(bigpoly* r, const bigpoly* p, const bigpoly* q);
void bigpoly_add(bigpoly* r, const bigpoly* p, const bigpoly* q);
void bigpoly_sub(bigpoly* r, const bigpoly* p, const bigpoly* q);
void bigpoly_test(); 

