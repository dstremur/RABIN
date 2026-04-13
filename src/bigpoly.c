#include "../include/bigpoly.h"

void bigpoly_init(bigpoly* p)
{
    bn_init(p->coeff);
	p->deg = 0;
	p->size = 0; 
}

void bigpoly_set(bigpoly* p, bignum* coeff, u64 deg)
{
	bn_copy(p->coeff, coeff); 
	p->deg = deg;
	p->size = deg; 


}


