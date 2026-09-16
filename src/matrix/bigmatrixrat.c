
#include "../../include/bigmatrixrat.h"

/*
void bigmatrixrat_normalize(bigmatrixrat *m) {
    bignum g;
    bn_init(&g);

    // Compute GCD of m->den and all entries in m->num
    bigmatrix_gcd_all(&g, &m->num);
    bn_gcd(&g, &g, &m->den);

    // If GCD > 1, divide both numerator matrix entries and denominator by g
    if (!bignum_is_one(&g)) {
        bigmatrix_exact_div_scalar(&m->num, &m->num, &g);
        bignum_exact_div(&m->den, &m->den, &g);
    }

    bignum_free(&g);
}

*/
