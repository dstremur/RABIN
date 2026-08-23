/*
 * bigpoly.c
 *
 * Polynomial arithmetic over bignums.
 *
 * This file implements polynomials with bignum coefficients: basic
 * management (init/alloc/copy/free/trim), addition, subtraction,
 * schoolbook and NTT-based multiplication (bignum and u64 variants),
 * and the decompose/carry-propagate/recompose pipeline used by the
 * NTT-based bignum multiplication path.
 *
 * A polynomial is stored as a dynamic array of bignum coefficients in
 * ascending order of degree, with a degree and a capacity.
 *
 * Copyright (C) 2026 Diego Strebel
 * SPDX-License-Identifier: GPL-3.0-or-later
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "../../include/bigpoly.h"

#include <inttypes.h>
#include <stdio.h>

#include "../../include/bigntt.h"
#include "../../include/u64.h"

/**
 * @brief Initialize a polynomial to the empty state (no coefficients).
 *
 * Complexity:
 *   Time: O(1)
 *   Auxiliary memory: O(1)
 *   Output memory: O(1)
 *
 * @param[out] p Polynomial to initialize.
 */
void bigpoly_init(bigpoly* p)
{
  p->coeff = NULL;
  p->deg = 0;
  p->size = 0;
}

/**
 * @brief Bitwise AND of two bignums: r = a & mask.
 *
 * Let n = a->size, measured in 64-bit limbs.
 *
 * Copies a into r and ANDs each limb with the corresponding limb of
 * mask (limbs beyond mask->size are zeroed), then trims. Used to
 * extract fixed-width chunks of a bignum.
 *
 * Complexity:
 *   Time: O(n)
 *   Auxiliary memory: O(1)
 *   Output memory: O(n) limbs
 *
 * @param[out]   r    Result storing a & mask.
 * @param[in]    a    First operand.
 * @param[in] mask Mask.
 */
void bn_and(bignum* r, const bignum* a, const bignum* mask)
{
  bn_copy(r, a);

  u64 min = MIN(r->size, mask->size);

  for (u64 i = 0; i < min; i++) {
    r->limbs[i] &= mask->limbs[i];
  }

  // zero the rest
  for (u64 i = min; i < r->size; i++) {
    r->limbs[i] = 0;
  }

  bn_trim(r);
}

/**
 * @brief Set a polynomial from an array of bignum coefficients.
 *
 * Let d = deg.
 *
 * Allocates d + 1 coefficient slots and deep-copies coeff[0..d] into
 * the polynomial, setting the degree to d.
 *
 * Complexity:
 *   Time: O(d * n) where n is the size of the coefficients in limbs
 *   Auxiliary memory: O(1)
 *   Output memory: O(d) bignums
 *
 * @param[out]    p     Polynomial to set.
 * @param[in] coeff Array of coefficients.
 * @param[in]   deg   Degree of the polynomial.
 */
void bigpoly_set(bigpoly* p, bignum* coeff, u64 deg)
{
  bigpoly_alloc(p, deg + 1);

  for (u64 i = 0; i <= deg; i++) {
    bn_copy(&p->coeff[i], &coeff[i]);
  }

  p->deg = deg;
}

/**
 * @brief Set a polynomial from an array of 64-bit signed coefficients.
 *
 * Let d = deg.
 *
 * Allocates coefficient slots and stores coeff[0..d] as single-limb
 * bignums, then trims leading zero coefficients.
 *
 * Complexity:
 *   Time: O(d)
 *   Auxiliary memory: O(1)
 *   Output memory: O(d) bignums
 *
 * @param[out]    p     Polynomial to set.
 * @param[in] coeff Array of 64-bit signed coefficients.
 * @param[in]   deg   Degree of the polynomial.
 */
void bigpoly_set_i64(bigpoly* p, i64* coeff, u64 deg)
{
  bigpoly_alloc(p, deg);

  for (u64 i = 0; i <= deg; i++) {
    bn_set_i64(&p->coeff[i], coeff[i]);
  }

  p->deg = deg;
  bigpoly_trim(p);
}

/**
 * @brief Copy the coefficients of q into p (up to q->deg).
 *
 * Let d = q->deg.
 *
 * p must have capacity for at least d + 1 coefficients. The degree of
 * p is left unchanged.
 *
 * Complexity:
 *   Time: O(d * n) where n is the size of the coefficients in limbs
 *   Auxiliary memory: O(1)
 *   Output memory: O(d) bignums
 *
 * @param[out] p Destination polynomial.
 * @param[in]  q Source polynomial.
 */
void bigpoly_copy(bigpoly* p, bigpoly* q)
{
  for (u64 i = 0; i <= q->deg; i++) {
    bn_copy(&p->coeff[i], &q->coeff[i]);
  }
}

/**
 * @brief Test whether two polynomials are equal.
 *
 * Let d = max(a->deg, b->deg).
 *
 * Returns true iff the degrees match and all coefficients up to the
 * degree are equal (deep bignum comparison).
 *
 * Complexity:
 *   Time: O(d * n) where n is the size of the coefficients in limbs
 *   Auxiliary memory: O(1)
 *   Output memory: O(1)
 *
 * @param[in] a First polynomial.
 * @param[in] b Second polynomial.
 *
 * @return true  If the polynomials are equal.
 * @return false If the polynomials differ.
 */
bool bigpoly_equal(bigpoly* a, bigpoly* b)
{
  if (a->deg != b->deg) return false;
  for (u64 i = 0; i <= a->deg; i++) {
    if (bn_cmp(&a->coeff[i], &b->coeff[i]) != 0) {
      return false;
    }
  }
  return true;
}

/**
 * @brief Free all coefficient storage of a polynomial.
 *
 * Complexity:
 *   Time: O(size)
 *   Auxiliary memory: O(1)
 *   Output memory: O(1)
 *
 * @param[in,out] p Polynomial to free.
 */
void bigpoly_free(bigpoly* p)
{
  if (!p->coeff) return;

  for (u64 i = 0; i < p->size; i++) {
    bn_free(&p->coeff[i]);
  }

  free(p->coeff);
  p->coeff = NULL;
  p->deg = 0;
  p->size = 0;
}

/**
 * @brief Grow the coefficient array of a polynomial to at least deg + 1
 * slots.
 *
 * If the capacity is already sufficient, nothing happens. Otherwise
 * the array is reallocated with exponential growth and the new slots
 * are initialized to zero bignums.
 *
 * Returns true on success, false if the realloc fails.
 *
 * Complexity:
 *   Time: O(new_cap) for the initialization of the new slots
 *   Auxiliary memory: O(new_cap) during the realloc
 *   Output memory: O(new_cap) bignums
 *
 * @param[in,out] p   Polynomial to grow.
 * @param[in]  deg Minimum number of coefficient slots minus one.
 *
 * @return true  On success.
 * @return false If the realloc fails.
 */
bool bigpoly_alloc(bigpoly* p, u64 deg)
{
  if (p->size >= deg + 1) return true;

  u64 new_cap = (p->size == 0) ? deg + 1 : (p->size * 2);
  if (new_cap < deg + 1) {
    new_cap = deg + 1;
  }

  bignum* new = realloc(p->coeff, new_cap * sizeof(bignum));
  if (!new) return false;

  p->coeff = new;

  for (u64 i = p->size; i < new_cap; i++) {
    bn_init(&p->coeff[i]);
  }

  p->size = new_cap;
  return true;
}

/**
 * @brief Remove leading zero coefficients from a polynomial.
 *
 * Lowers p->deg until the top coefficient is nonzero (or the degree
 * is 0). The capacity is left unchanged.
 *
 * Complexity:
 *   Time: O(number of trimmed coefficients)
 *   Auxiliary memory: O(1)
 *   Output memory: O(1)
 *
 * @param[in,out] p Polynomial to trim.
 */
void bigpoly_trim(bigpoly* p)
{
  while (p->deg > 0 && bn_is_zero(&p->coeff[p->deg])) {
    p->deg--;
  }
}

/**
 * @brief Print a polynomial to stdout in descending powers of x.
 *
 * Complexity:
 *   Time: O(d * n) where d is the degree and n the coefficient size
 *   Auxiliary memory: O(1)
 *   Output memory: O(d * n) characters written
 *
 * @param[in] p Polynomial to print.
 */
void bigpoly_print(const bigpoly* p)
{
  for (i64 i = p->deg; i >= 0; i--) {
    bn_print(&p->coeff[i]);
    printf("x^%" PRId64 "", i);
    if (i != 0) {
      printf("+");
    }
  }
  printf("\n");
}

/**
 * @brief Add two polynomials: r = p + q.
 *
 * Let d = max(p->deg, q->deg).
 *
 * Coefficients up to the smaller degree are added with bn_add(); the
 * remaining coefficients of the longer polynomial are copied. The
 * result is trimmed.
 *
 * Complexity:
 *   Time: O(d * n) where n is the size of the coefficients in limbs
 *   Auxiliary memory: O(1)
 *   Output memory: O(d) bignums
 *
 * @param[out] r Result polynomial.
 * @param[in]  p First polynomial.
 * @param[in]  q Second polynomial.
 */
void bigpoly_add(bigpoly* r, const bigpoly* p, const bigpoly* q)
{
  u64 min = MIN(p->deg, q->deg);
  u64 max = MAX(p->deg, q->deg);

  if (!bigpoly_alloc(r, max)) {
    // TODO: error handling
    return;
  }

  // add upto the minimum degree
  for (u64 i = 0; i <= min; i++) {
    bn_add(&r->coeff[i], &p->coeff[i], &q->coeff[i]);
  }

  const bigpoly* longer = (p->deg > q->deg) ? p : q;
  for (u64 i = min + 1; i <= max; i++) {
    bn_copy(&r->coeff[i], &longer->coeff[i]);
  }

  r->deg = max;
  bigpoly_trim(r);
}

/**
 * @brief Subtract two polynomials: r = p - q.
 *
 * Let d = max(p->deg, q->deg).
 *
 * Coefficients up to the smaller degree are subtracted with bn_sub();
 * the remaining coefficients of the longer polynomial are copied with
 * their sign flipped. The result is trimmed.
 *
 * Complexity:
 *   Time: O(d * n) where n is the size of the coefficients in limbs
 *   Auxiliary memory: O(1)
 *   Output memory: O(d) bignums
 *
 * @param[out] r Result polynomial.
 * @param[in]  p First polynomial.
 * @param[in]  q Second polynomial.
 */
void bigpoly_sub(bigpoly* r, const bigpoly* p, const bigpoly* q)
{
  u64 min = MIN(p->deg, q->deg);
  u64 max = MAX(p->deg, q->deg);

  if (!bigpoly_alloc(r, max)) {
    // TODO: error handling
    return;
  }

  // add upto the minimum degree
  for (u64 i = 0; i <= min; i++) {
    bn_sub(&r->coeff[i], &p->coeff[i], &q->coeff[i]);
  }

  const bigpoly* longer = (p->deg > q->deg) ? p : q;
  for (u64 i = min + 1; i <= max; i++) {
    bn_copy(&r->coeff[i], &longer->coeff[i]);
    r->coeff[i].is_neg = true;
  }

  r->deg = max;
  bigpoly_trim(r);
}

/**
 * @brief Multiply two polynomials with the bignum NTT: r = p * q.
 *
 * Let d = p->deg + q->deg + 1 (length of the true product) and
 * n = 2^k the next power of two >= d.
 *
 * Pads both polynomials to length n, applies the forward bignum NTT
 * (Goldilocks field) to each, multiplies pointwise in the Montgomery
 * domain, and applies the inverse NTT. The first d coefficients of the
 * cyclic convolution equal the true product coefficients (no
 * wrap-around, since the Goldilocks prime is large enough for the
 * coefficient sizes used).
 *
 * A fresh NTT context is initialized per call.
 *
 * Complexity:
 *   Time: O(n log n * k_l^2) where k_l is the size of the modulus in
 *         limbs, plus O(n * k_l^2) for the context tables
 *   Auxiliary memory: O(n) bignums for the transformed copies
 *   Output memory: O(d) bignums
 *
 * @param[out] r Result polynomial.
 * @param[in]  p First polynomial.
 * @param[in]  q Second polynomial.
 */
void bigpoly_mul_ntt(bigpoly* r, const bigpoly* p, const bigpoly* q)
{
  ntt_ctx ctx;

  u64 required_len = p->deg + q->deg + 1;
  u64 ntt_size = 1;
  u64 k = 0;
  // pad to next power of 2
  while (ntt_size < required_len) {
    ntt_size <<= 1;
    k++;
  }

  bigpoly_alloc(r, required_len);

  if (!bigntt_ctx_init_golden(&ctx, k)) {
    return;
  }

  bigpoly p_hat, q_hat, r_hat;
  bigpoly_init(&p_hat);
  bigpoly_init(&q_hat);
  bigpoly_init(&r_hat);

  bigntt_cyclic_forward(&p_hat, p, &ctx);
  bigntt_cyclic_forward(&q_hat, q, &ctx);
  bigpoly_alloc(&r_hat, ntt_size);

  // bigpoly_mul_digit(&r_hat, &p_hat, &q_hat, &ctx.q);

  for (u64 i = 0; i < ntt_size; i++) {
    // MontMul(pR, qR) = (pR * qR * R^-1) mod q = (p*q)R mod q
    bn_mont_mul(&r_hat.coeff[i], &p_hat.coeff[i], &q_hat.coeff[i], &ctx.mctx);
  }
  r_hat.deg = ntt_size - 1;

  bigpoly r_ntt;
  bigpoly_init(&r_ntt);
  bigntt_cyclic_inverse_mont_in(&r_ntt, &r_hat, &ctx);

  bigpoly_free(r);
  bigpoly_init(r);
  bigpoly_alloc(r, required_len);

  for (u64 i = 0; i < required_len; i++) {
    bn_copy(&r->coeff[i], &r_ntt.coeff[i]);
  }

  r->deg = required_len - 1;

  bigpoly_trim(r);

  bigpoly_free(&p_hat);
  bigpoly_free(&q_hat);
  bigpoly_free(&r_hat);
  bigpoly_free(&r_ntt);
  bigntt_ctx_free(&ctx);
}

/**
 * @brief Multiply two polynomials with the u64 NTT: r = p * q.
 *
 * Let d = p->deg + q->deg + 1 (length of the true product) and
 * n = 2^k the next power of two >= d.
 *
 * Only works for coefficients that fit in a single u64 limb: each
 * coefficient is reduced to its low limb, the convolution is computed
 * with the fast u64 Goldilocks NTT (flat arrays, Montgomery
 * butterflies), and the result limbs are lifted back to bignums.
 *
 * The Goldilocks NTT context for the transform size is taken from a
 * shared per-size cache (built once per k, see
 * ntt_ctx_u64_golden_cached), and the flat u64 arrays come from the
 * thread-local scratch arena.
 *
 * Complexity:
 *   Time: O(n log n) for the NTTs
 *   Auxiliary memory: O(n) u64s for the flat arrays (scratch arena)
 *   Output memory: O(d) bignums
 *
 * @param[out] r Result polynomial.
 * @param[in]  p First polynomial.
 * @param[in]  q Second polynomial.
 */
void bigpoly_mul_ntt_u64(bigpoly* r, const bigpoly* p, const bigpoly* q)
{
  u64 required_len = p->deg + q->deg + 1;
  u64 ntt_size = 1;
  u64 k = 0;

  // Pad to next power of 2
  while (ntt_size < required_len) {
    ntt_size <<= 1;
    k++;
  }

  ntt_ctx_u64* ctx = ntt_ctx_u64_golden_cached(k);
  if (!ctx) {
    return;
  }

  // 1. One arena block for all six flat u64 arrays (zeroed, like calloc)
  u64* buf = bn_scratch_get(6 * ntt_size);
  memset(buf, 0, 6 * ntt_size * sizeof(u64));
  u64* p_arr = buf;
  u64* q_arr = buf + ntt_size;
  u64* p_hat = buf + 2 * ntt_size;
  u64* q_hat = buf + 3 * ntt_size;
  u64* r_hat = buf + 4 * ntt_size;
  u64* r_arr = buf + 5 * ntt_size;

  // 2. Extract u64 values from the bignum polynomials
  for (u64 i = 0; i <= p->deg; i++) {
    p_arr[i] = p->coeff[i].limbs[0];
  }
  for (u64 i = 0; i <= q->deg; i++) {
    q_arr[i] = q->coeff[i].limbs[0];
  }

  // 3. Perform Forward NTTs
  ntt_u64_cyclic_forward(p_hat, p_arr, ctx);
  ntt_u64_cyclic_forward(q_hat, q_arr, ctx);

  // 4. Pointwise Multiplication
  for (u64 i = 0; i < ntt_size; i++) {
    r_hat[i] = mont_mul(p_hat[i], q_hat[i], &ctx->mctx);
  }

  // 5. Perform Inverse NTT
  ntt_u64_cyclic_inverse_montgomery_in(r_arr, r_hat, ctx);

  bigpoly_free(r);
  bigpoly_init(r);
  bigpoly_alloc(r, required_len);

  for (u64 i = 0; i < required_len; i++) {
    bn_set_u64(&r->coeff[i], r_arr[i]);
  }

  r->deg = required_len - 1;
  bigpoly_trim(r);

  // 7. Cleanup
  bn_scratch_release();
}

/**
 * @brief Multiply two polynomials with the schoolbook algorithm: r = p * q.
 *
 * Let d = p->deg + q->deg.
 *
 * Each pair of coefficients is multiplied and accumulated into the
 * corresponding slot of a temporary polynomial, which is then copied
 * into r and trimmed.
 *
 * Complexity:
 *   Time: O((p->deg + 1) * (q->deg + 1) * n^2) where n is the size of
 *         the coefficients in limbs
 *   Auxiliary memory: O(d) bignums for the temporary
 *   Output memory: O(d) bignums
 *
 * @param[out] r Result polynomial.
 * @param[in]  p First polynomial.
 * @param[in]  q Second polynomial.
 */
void bigpoly_mul_school(bigpoly* r, const bigpoly* p, const bigpoly* q)
{
  u64 r_deg = p->deg + q->deg;

  bigpoly temp;
  bigpoly_init(&temp);
  bigpoly_alloc(&temp, r_deg);

  temp.deg = r_deg;

  bignum t;
  bn_init(&t);

  for (u64 i = 0; i <= p->deg; i++) {
    for (u64 j = 0; j <= q->deg; j++) {
      bn_mul(&t, &p->coeff[i], &q->coeff[j]);

      bn_add(&temp.coeff[i + j], &temp.coeff[i + j], &t);
    }
  }

  bn_free(&t);

  bigpoly_trim(&temp);

  if (bigpoly_alloc(r, temp.deg)) {
    for (u64 i = 0; i <= temp.deg; i++) {
      bn_copy(&r->coeff[i], &temp.coeff[i]);
    }
    r->deg = temp.deg;
  }

  bigpoly_free(&temp);
}

/**
 * @brief Multiply two polynomials: r = p * q.
 *
 * Currently a thin wrapper around bigpoly_mul_school().
 *
 * Complexity:
 *   Time: see bigpoly_mul_school()
 *   Auxiliary memory: O(d) bignums
 *   Output memory: O(d) bignums
 *
 * @param[out] r Result polynomial.
 * @param[in]  p First polynomial.
 * @param[in]  q Second polynomial.
 */
void bigpoly_mul(bigpoly* r, const bigpoly* p, const bigpoly* q)
{
  bigpoly_mul_school(r, p, q);
}

/**
 * @brief Decompose a bignum into a polynomial of fixed-width chunks.
 *
 * Let n = n->size, measured in 64-bit limbs, and W = width.
 *
 * Slices n into chunks of W bits (base 2^W) and stores them as the
 * coefficients of r in ascending order:
 *
 *   n = sum_i r->coeff[i] * (2^W)^i
 *
 * The degree is ceil(n / W) - 1 (0 for n == 0).
 *
 * Complexity:
 *   Time: O(n * ceil(n/W)) - one O(n) shift per chunk
 *   Auxiliary memory: O(n) limbs for temporaries
 *   Output memory: O(n/W) bignums
 *
 * @param[out]    r    Result polynomial storing the chunks.
 * @param[in]     n    Bignum to decompose.
 * @param[in] width Chunk width in bits.
 */
void bn_decompose(bigpoly* r, const bignum* n, u64 width)
{
  bignum tmp, mask, digit;
  bn_init_multi(&tmp, &mask, &digit);
  bn_copy(&tmp, n);

  // mask = 0xFFFF..
  bn_lshift(&mask, &BN_ONE, width);
  bn_sub(&mask, &mask, &BN_ONE);

  u64 i = 0;
  while (!bn_is_zero(&tmp)) {
    bn_and(&digit, &tmp, &mask);

    bigpoly_alloc(r, i + 1);
    bn_copy(&r->coeff[i], &digit);

    bn_rshift(&tmp, &tmp, width);

    i++;
  }

  r->deg = (i > 0) ? i - 1 : 0;

  bn_free_multi(&tmp, &mask, &digit);
}

/**
 * @brief Propagate carries between the fixed-width coefficient slots of a
 * polynomial, in place.
 *
 * Let d = r->deg and W = bit_width.
 *
 * After an NTT-based convolution the coefficient slots may exceed
 * 2^W - 1. This normalizes them: for each slot, total = coeff + carry,
 * coeff = total mod 2^W, carry = total / 2^W, rippling upward and
 * growing the polynomial if the carry outlives the current degree.
 *
 * Complexity:
 *   Time: O(d * n) where n is the size of the coefficients in limbs
 *   Auxiliary memory: O(n) limbs for temporaries
 *   Output memory: O(d) bignums (possibly grown by the final carry)
 *
 * @param[in,out]      r       Polynomial to normalize (modified in place).
 * @param[in] bit_width Width of each coefficient slot in bits.
 */
void poly_carry_propagation(bigpoly* r, u64 bit_width)
{
  bignum carry, base, mask, total;
  bn_init_multi(&carry, &base, &mask, &total, NULL);

  bn_lshift(&base, &BN_ONE, bit_width);
  bn_sub(&mask, &base, &BN_ONE);
  bn_set_u64(&carry, 0);

  u64 i = 0;
  // Iterate through all coefficients plus any remaining carries
  while (i <= r->deg || !bn_is_zero(&carry)) {
    if (i > r->deg) {
      bigpoly_alloc(r, i + 1);  // Expand poly if carry exceeds current deg
      r->deg = i;
    }

    // total = coeff[i] + carry
    bn_add(&total, &r->coeff[i], &carry);

    // carry = total >> bit_width
    bn_rshift(&carry, &total, bit_width);

    // coeff[i] = total & mask
    bn_and(&r->coeff[i], &total, &mask);

    i++;
  }

  bn_free_multi(&carry, &base, &mask, &total, NULL);
}

/**
 * @brief Recompose a polynomial of fixed-width chunks into a bignum.
 *
 * Let d = p->deg and W = bit_width.
 *
 * The inverse of bn_decompose():
 *
 *   n = sum_i p->coeff[i] << (i * W)
 *
 * Complexity:
 *   Time: O(d * n) where n is the size of the result in limbs
 *   Auxiliary memory: O(n) limbs for the running term
 *   Output memory: O(n) limbs
 *
 * @param[out]        n         Result bignum.
 * @param[in]         p         Polynomial of chunks.
 * @param[in] bit_width Width of each chunk in bits.
 */
void bn_recompose(bignum* n, const bigpoly* p, u64 bit_width)
{
  bn_set_u64(n, 0);
  bignum term;
  bn_init(&term);

  for (u64 i = 0; i <= p->deg; i++) {
    // term = coeff[i] << (i * bit_width)
    bn_lshift(&term, &p->coeff[i], i * bit_width);
    // n += term
    bn_add(n, n, &term);
  }

  bn_free(&term);
}

/**
 * @brief Self-test of the polynomial operations.
 *
 * Builds p(x) = 2x^2 + 3x + 1 and q(x) = 4x + 5, prints them, and
 * prints p + q and p * q (schoolbook and u64-NTT). Intended for
 * manual verification, not part of the library API.
 *
 * Complexity:
 *   Time: O(1) (fixed-size inputs)
 *   Auxiliary memory: O(1)
 *   Output memory: O(1)
 */
void bigpoly_test()
{
  printf("--- BigPoly Test --- \n");

  bigpoly p, q, r;
  bigpoly_init(&p);
  bigpoly_init(&q);
  bigpoly_init(&r);

  // Let p(x) = 2x^2 + 3x + 1
  i64 p_vals[] = {1, 3, 2};
  bigpoly_set_i64(&p, p_vals, 2);

  // Let q(x) = 4x + 5
  i64 q_vals[] = {5, 4};
  bigpoly_set_i64(&q, q_vals, 1);

  printf("\nPolynomial P(x):\n");
  bigpoly_print(&p);

  printf("\nPolynomial Q(x):\n");
  bigpoly_print(&q);

  // Test Addition: r = p + q = 2x^2 + 7x + 6
  printf("\n--- Test: Addition (P + Q) ---\n");
  bigpoly_add(&r, &p, &q);
  bigpoly_print(&r);

  // Test Multiplication: r = p * q = 8x^3 + 22x^2 + 19x + 5
  printf("\n--- Test: Multiplication (P * Q) ---\n");
  bigpoly_mul_school(&r, &p, &q);
  bigpoly_print(&r);

  bigpoly_free(&r);
  bigpoly_init(&r);

  bigpoly_mul_ntt_u64(&r, &p, &q);
  bigpoly_print(&r);

  // Cleanup
  bigpoly_free(&p);
  bigpoly_free(&q);
  bigpoly_free(&r);

  printf("\nTests complete.\n");
}
