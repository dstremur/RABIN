#ifndef BIGNUM_H
#define BIGNUM_H

/*
 *  Modules:
 *    bigcore.h    bignum type + core arithmetic      (src/core)
 *    bighelper.h  raw-limb helpers                   (src/utils)
 *    u64.h        64-bit single-limb number theory   (src/u64)
 *    bigmatrix.h  dense bignum matrices              (src/matrix)
 *    bigrns.h     residue number systems             (src/rns)
 *    bigprime.h   prime tests, generation, proving   (src/number_theory)
 *    bigrabin.h   Miller-Rabin tests                 (src/number_theory)
 *    bigrand.h    random bignum generation           (src/number_theory)
 *    bigpseudo.h  strong pseudoprime generation      (src/number_theory)
 *    bigfactor.h  factorization (rho, p-1, trial)    (src/number_theory)
 *    bigmath.h    Jacobi / Tonelli-Shanks / gcd      (src/number_theory)
 *    biglucas.h   Lucas sequences                    (src/number_theory)
 *    bigpoly.h    integer polynomials                (src/poly)
 *    bigntt.h     bignum NTT                         (src/number_theory)
 *    bigfield.h   Z_m and Z_m[x]/(q) arithmetic      (src/number_theory)
 *    bigvector.h  bignum vectors                     (src/vector)
 *    bigcert.h    Pocklington certificates           (src/number_theory)
 */

#include "bigcert.h"
#include "bigcore.h"
#include "bigfactor.h"
#include "bigfield.h"
#include "bighelper.h"
#include "biglucas.h"
#include "bigmath.h"
#include "bigmatrix.h"
#include "bigntt.h"
#include "bigpoly.h"
#include "bigprime.h"
#include "bigpseudo.h"
#include "bigrabin.h"
#include "bigrand.h"
#include "bigrns.h"
#include "bigvector.h"
#include "u64.h"

/* Ideas
 Catalan pseudoprime

 Hensel lifting

 berlenkamp algo

 LLL

 Discrete log problem

 suntherlands algorithm

 smith normal form

  */

#endif
