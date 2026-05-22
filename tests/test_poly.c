#include "../include/bignum.h" 
#include "../include/bigpoly.h" 


int main() 
{
	u64 n = 10;
	bigpoly a, b, r;
	bigpoly_init(&a);
	bigpoly_init(&b);
	bigpoly_init(&r);
	bigpoly_alloc(&a, n);
	bigpoly_alloc(&b, n);
	a.deg = b.deg = n - 1;
	bigpoly_alloc(&r, n);

	bignum* C = malloc(sizeof(bignum) * n); 

	for (u64 i = 0; i < n; i++) {
		bn_init(&C[i]);
		bn_set_u64(&C[i], i); 
	}

	bigpoly_set(&a, C, n - 1); 
	bigpoly_set(&b, C, n - 1); 

	bigpoly_print(&a);
	bigpoly_print(&b);

	bigpoly_mul_ntt(&r, &a, &b);
	bigpoly_print(&r);

	for (u64 i = 0; i < n; i++) {
		bn_free(&C[i]);
	}
	free(C);

}
