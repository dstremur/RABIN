#ifndef RABIN_H
#define RABIN_H

/**
 * @brief Mapping of the SIMD instructions used by the library kernels.
 *
 * For each class of operation the table lists the primary SSE/AVX
 * intrinsics, the corresponding assembly mnemonics, and the role the
 * operation plays in the bignum kernels.
 *
 * Bitwise logic
 *   Intrinsics: _mm256_and_si256, _mm256_or_si256, _mm256_xor_si256,
 *               _mm256_andnot_si256
 *   Assembly:   vpand, vpor, vpxor, vpandn
 *   Use:        Zero-carry operations: bitwise filtering, masks, and
 *               logical adjustments.
 *
 * Packed arithmetic
 *   Intrinsics: _mm256_add_epi64, _mm256_sub_epi64
 *   Assembly:   vpaddq, vpsubq
 *   Use:        Block-wise addition and subtraction prior to carry
 *               resolution.
 *
 * Multiplication
 *   Intrinsics: _mm256_mul_epu32, _mm512_madd52lo_epu64
 *   Assembly:   vpmuludq, vpmadd52luq
 *   Use:        Parallel limb multiplication and high-radix
 *               cryptographic arithmetic.
 *
 * Shifts & routing
 *   Intrinsics: _mm256_slli_epi64, _mm256_permute4x64_epi64
 *   Assembly:   vpsllq, vpermq
 *   Use:        Shifting limbs and routing overflow bits across vector
 *               lanes.
 *
 * Comparisons
 *   Intrinsics: _mm256_cmpeq_epi64, _mm256_cmpgt_epi64
 *   Assembly:   vpcmpeqq, vpcmpgtq
 *   Use:        Magnitude ordering, equality tests, and conditional
 *               checks.
 *
 * Memory / NTT
 *   Intrinsics: _mm256_loadu_si256, _mm256_shuffle_epi32
 *   Assembly:   vmovdqu, vpshufd
 *   Use:        High-bandwidth data loading and fast polynomial
 *               multiplication (NTT).
 */

/*
 *  Modules:
 *    rz.h    rz_t type + core arithmetic      (src/core)
 *    rzlimb.h  raw-limb helpers                   (src/utils)
 *    u64.h        64-bit single-limb number theory   (src/u64)
 *    rmat.h  dense rz_t matrices              (src/matrix)
 *    rns.h     residue number systems             (src/rns)
 *    rzprime.h   prime tests, generation, proving   (src/number_theory)
 *    rzrabin.h   Miller-Rabin tests                 (src/number_theory)
 *    rzrand.h    random rz_t generation           (src/number_theory)
 *    rzfactor.h  factorization (rho, p-1, trial)    (src/number_theory)
 *    rzmath.h    Jacobi / Tonelli-Shanks / gcd      (src/number_theory)
 *    rzlucas.h   Lucas sequences                    (src/number_theory)
 *    rpol.h    integer polynomials                (src/poly)
 *    rntt.h     rz_t NTT                         (src/number_theory)
 *    rzfield.h   Z_m and Z_m[x]/(q) arithmetic      (src/number_theory)
 *    rvec.h  rz_t vectors                     (src/vector)
 *    rzcert.h    Pocklington certificates           (src/number_theory)
 */

#include "rabin_errors.h"
#include "rmat.h"
#include "rns.h"
#include "rntt.h"
#include "rpol.h"
#include "rvec.h"
#include "rz.h"
#include "rzcert.h"
#include "rzfactor.h"
#include "rzfield.h"
#include "rzlimb.h"
#include "rzlogic.h"
#include "rzlucas.h"
#include "rzmath.h"
#include "rzprime.h"
#include "rzrabin.h"
#include "rzrand.h"
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
