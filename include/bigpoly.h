#include "bignum.h" 

typedef struct {
	bignum* coeff;
	u64 size;
	u64 deg;
} bigpoly; 


void bigpoly_init(bigpoly* p); 

void bigpoly_set(bigpoly* p, bignum* coeff, u64 deg);

void bigpoly_print(bigpoly* p); 
bool bigpoly_alloc(bigpoly* p, u64 deg);
void bigpoly_trim(bigpoly* p);
void bigpoly_test(); 

