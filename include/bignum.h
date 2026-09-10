#ifndef BIGNUM_H
#define BIGNUM_H

/**
 *Category	Primary Intrinsics	Assembly / Mapping	Core Bignum Use
Case Bitwise Logic

_mm256_and_si256

_mm256_or_si256

_mm256_xor_si256

_mm256_andnot_si256


vpand, vpor

vpxor, vpandn
        Zero-carry operations: bitwise filtering, masks, and logical
adjustments. Packed Arithmetic

_mm256_add_epi64

_mm256_sub_epi64


vpaddq

vpsubq
        Block-wise addition and subtraction prior to carry resolution.
Multiplication

_mm256_mul_epu32

_mm512_madd52lo_epu64


vpmuludq

vpmadd52luq
        Parallel limb multiplication and high-radix cryptographic arithmetic.
Shifts & Routing

_mm256_slli_epi64

_mm256_permute4x64_epi64


vpsllq

vpermq
        Shifting limbs and routing overflow bits across vector lanes.
Comparisons

_mm256_cmpeq_epi64

_mm256_cmpgt_epi64


vpcmpeqq

vpcmpgtq
        Magnitude ordering, equality tests, and conditional checks.
Memory / NTT

_mm256_loadu_si256

_mm256_shuffle_epi32


vmovdqu

vpshufd
        High-bandwidth data loading and fast polynomial multiplication (NTT).
 */

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
