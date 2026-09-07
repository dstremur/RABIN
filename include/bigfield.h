#ifndef BIGFIELD_H
#define BIGFIELD_H

#include "bigcore.h"
#include "bigpoly.h"

/*
 * field_ctx: arithmetic in the ring Z_m = Z / m Z.
 *
 * If m is odd and > 1 the elements are stored in the Montgomery domain
 * (a*R mod m) and the operations use the fast REDC path. Otherwise (m
 * even) the elements are stored in plain form [0, m) and the operations
 * use plain multiply + reduce. The representation is an internal detail:
 * every field_* operation is uniform across both paths. Use field_in() /
 * field_out() to convert between the external plain form and the internal
 * representation.
 *
 * Inverses exist exactly for the units of Z_m (gcd(a, m) == 1); for a
 * prime modulus that is every nonzero element.
 */
typedef struct {
  bignum* q;        /* modulus m (owned heap copy)            */
  bool mont;        /* true if the Montgomery path is active  */
  bn_mont_ctx mctx; /* valid only if mont                     */
} field_ctx;

/*
 * poly_ring: the ring of polynomials over Z_m modulo q, i.e. Z_m[x]/(q).
 *
 * Coefficients are field elements of fctx (internal representation). The
 * modulus polynomial q is expected to have a unit leading coefficient
 * (in particular it is usually taken monic) so that reduction is always
 * possible.
 */
typedef struct {
  bigpoly* q;      /* modulus polynomial (owned heap copy)    */
  field_ctx* fctx; /* coefficient ring Z_m (borrowed, not owned) */
} poly_ring;

/**
 * @brief Initialize a Z_m context.
 *
 * Copies the modulus m and selects the arithmetic path: the Montgomery
 * path is used iff m is odd and > 1 (REDC requires an odd modulus), the
 * plain multiply + reduce path otherwise.
 *
 * Complexity:
 *   Time: O(n^2) for the Montgomery setup (one division for R mod m and
 *         one for R^2 mod m), O(n) otherwise, where n = m->size
 *   Auxiliary memory: O(n) limbs
 *   Output memory: O(n) limbs per context field
 *
 * @param[out] ctx Context to initialize.
 * @param[in]  m   Modulus (must be > 1).
 */
void field_ctx_init(field_ctx* ctx, const bignum* m);

/**
 * @brief Free a Z_m context.
 *
 * Complexity:
 *   Time: O(n)
 *   Auxiliary memory: O(1)
 *   Output memory: O(1)
 *
 * @param[in,out] ctx Context to free.
 */
void field_ctx_free(field_ctx* ctx);

/**
 * @brief Convert a plain value into the internal representation.
 *
 * Montgomery path: r = a*R mod m. Plain path: r = a mod m in [0, m).
 *
 * Complexity:
 *   Time: O(n^2) (Montgomery) or O(n^2) (one division, plain)
 *   Auxiliary memory: O(n) limbs
 *   Output memory: O(n) limbs
 *
 * @param[out] r Result in the internal representation.
 * @param[in]  a Plain value.
 * @param[in]  ctx Initialized context.
 */
void field_in(bignum* r, const bignum* a, field_ctx* ctx);

/**
 * @brief Convert an internal value back to the plain form [0, m).
 *
 * Montgomery path: r = a*R^(-1) mod m. Plain path: r = a.
 *
 * Complexity:
 *   Time: O(n^2) (Montgomery) or O(n) (plain copy)
 *   Auxiliary memory: O(n) limbs
 *   Output memory: O(n) limbs
 *
 * @param[out] r Result in the plain form.
 * @param[in]  a Value in the internal representation.
 * @param[in]  ctx Initialized context.
 */
void field_out(bignum* r, const bignum* a, field_ctx* ctx);

/**
 * @brief Set r to the 64-bit value v in the internal representation.
 *
 * Complexity:
 *   Time: O(n^2) (Montgomery) or O(1) (plain)
 *   Auxiliary memory: O(n) limbs
 *   Output memory: O(n) limbs
 *
 * @param[out] r Result.
 * @param[in]  v 64-bit value.
 * @param[in]  ctx Initialized context.
 */
void field_set_u64(bignum* r, u64 v, field_ctx* ctx);

/**
 * @brief r = (a + b) mod m.
 *
 * Both operands are in [0, m), so a + b < 2m and a single conditional
 * subtraction suffices. Addition is representation independent.
 *
 * Complexity:
 *   Time: O(n)
 *   Auxiliary memory: O(1)
 *   Output memory: O(n) limbs
 *
 * @param[out] r Result.
 * @param[in]  a First operand.
 * @param[in]  b Second operand.
 * @param[in]  ctx Initialized context.
 */
void field_add(bignum* r, const bignum* a, const bignum* b, field_ctx* ctx);

/**
 * @brief r = (a - b) mod m.
 *
 * Both operands are in [0, m), so a - b is in (-m, m); if the result is
 * negative, m is added to bring it into [0, m). Subtraction is
 * representation independent.
 *
 * Complexity:
 *   Time: O(n)
 *   Auxiliary memory: O(1)
 *   Output memory: O(n) limbs
 *
 * @param[out] r Result.
 * @param[in]  a First operand.
 * @param[in]  b Second operand.
 * @param[in]  ctx Initialized context.
 */
void field_sub(bignum* r, const bignum* a, const bignum* b, field_ctx* ctx);

/**
 * @brief r = (-a) mod m.
 *
 * Computes m - a and reduces the result (which equals m when a == 0)
 * back into [0, m).
 *
 * Complexity:
 *   Time: O(n)
 *   Auxiliary memory: O(1)
 *   Output memory: O(n) limbs
 *
 * @param[out] r Result.
 * @param[in]  a Operand.
 * @param[in]  ctx Initialized context.
 */
void field_neg(bignum* r, const bignum* a, field_ctx* ctx);

/**
 * @brief r = (a * b) mod m.
 *
 * Montgomery path: one REDC of the product. Plain path: full multiply
 * followed by one reduction.
 *
 * Complexity:
 *   Time: O(n^2)
 *   Auxiliary memory: O(n) limbs
 *   Output memory: O(n) limbs
 *
 * @param[out] r Result.
 * @param[in]  a First operand.
 * @param[in]  b Second operand.
 * @param[in]  ctx Initialized context.
 */
void field_mul(bignum* r, const bignum* a, const bignum* b, field_ctx* ctx);

/**
 * @brief r = a^(-1) mod m.
 *
 * The inverse exists exactly when a is a unit of Z_m (gcd(a, m) == 1);
 * for a prime modulus that is every nonzero element. The Montgomery path
 * (odd m) uses the fast binary extended GCD (bn_mod_inverse) after
 * converting the operand to the plain domain. The plain path (even m)
 * uses the classical extended Euclidean algorithm (mod_inverse_euclid),
 * because the binary GCD is only correct for odd moduli.
 *
 * Complexity:
 *   Time: O(n^2) (Montgomery) or O(n^3) (plain)
 *   Auxiliary memory: O(n) limbs
 *   Output memory: O(n) limbs
 *
 * @param[out] r Receives the inverse if it exists.
 * @param[in]  a Value to invert.
 * @param[in]  ctx Initialized context.
 *
 * @return true  If a is a unit; r is set.
 * @return false If a is zero or not coprime to m; r is left unchanged.
 */
bool field_inv(bignum* r, const bignum* a, field_ctx* ctx);

/**
 * @brief r = (a * b^(-1)) mod m.
 *
 * Complexity:
 *   Time: O(n^2)
 *   Auxiliary memory: O(n) limbs
 *   Output memory: O(n) limbs
 *
 * @param[out] r Receives the quotient if b is a unit.
 * @param[in]  a Numerator.
 * @param[in]  b Denominator.
 * @param[in]  ctx Initialized context.
 *
 * @return true  If b is a unit; r is set.
 * @return false If b is zero or not coprime to m; r is left unchanged.
 */
bool field_div(bignum* r, const bignum* a, const bignum* b, field_ctx* ctx);

/**
 * @brief r = a^e mod m.
 *
 * Montgomery path: binary exponentiation with REDC (the base is already
 * in the Montgomery domain). Plain path: bn_mod_exp, which handles even
 * moduli with its slow multiply + divide loop.
 *
 * Complexity:
 *   Time: O(e * n^2) where e = e->size in limbs
 *   Auxiliary memory: O(n) limbs
 *   Output memory: O(n) limbs
 *
 * @param[out] r Result.
 * @param[in]  a Base.
 * @param[in]  e Exponent (plain, nonnegative).
 * @param[in]  ctx Initialized context.
 */
void field_pow(bignum* r, const bignum* a, const bignum* e, field_ctx* ctx);

/**
 * @brief Test whether a field element is zero.
 *
 * Complexity:
 *   Time: O(1)
 *   Auxiliary memory: O(1)
 *   Output memory: O(1)
 *
 * @param[in] a Field element.
 *
 * @return true If a == 0.
 */
bool field_is_zero(const bignum* a);

/**
 * @brief Test whether two field elements are equal.
 *
 * Complexity:
 *   Time: O(n)
 *   Auxiliary memory: O(1)
 *   Output memory: O(1)
 *
 * @param[in] a First element.
 * @param[in] b Second element.
 *
 * @return true If a == b.
 */
bool field_equal(const bignum* a, const bignum* b);

/**
 * @brief Polynomial long division in Z_m[x]: q = a / b, r = a mod b.
 *
 * Let d_a = a->deg and d_b = b->deg.
 *
 * Standard synthetic division: while deg(r) >= deg(b), the quotient term
 * is factor = r[deg r] * (b[deg b])^(-1) and r -= factor * x^shift * b.
 * The inverse of the leading coefficient of b is computed once.
 *
 * Requires the leading coefficient of b to be a unit of Z_m (always true
 * for a monic b). If b is zero or its leading coefficient is not
 * invertible, returns false with r = a and q = 0.
 *
 * Complexity:
 *   Time: O((d_a - d_b + 1) * d_b * n^2) field multiplications plus one
 *         O(n^2) inversion
 *   Auxiliary memory: O(d_a) bignums
 *   Output memory: O(d_a) bignums
 *
 * @param[out] q Quotient (may be NULL).
 * @param[out] r Remainder (may be NULL).
 * @param[in]  a Dividend.
 * @param[in]  b Divisor (nonzero, unit leading coefficient).
 * @param[in]  fctx Coefficient ring.
 *
 * @return true  On success.
 * @return false If b is zero or its leading coefficient is not a unit.
 */
bool bigpoly_divmod(bigpoly* q, bigpoly* r, const bigpoly* a, const bigpoly* b,
                    field_ctx* fctx);

/**
 * @brief Extended polynomial GCD in Z_m[x].
 *
 * Let d = max(a->deg, b->deg).
 *
 * Iterated bigpoly_divmod() tracking the Bezout coefficients:
 *
 *   (r0, s0, t0) = (a, 1, 0)
 *   (r1, s1, t1) = (b, 0, 1)
 *   while r1 != 0:
 *     (qq, r2) = divmod(r0, r1)
 *     (r0, r1) = (r1, r2)
 *     (s0, s1) = (s1, s0 - qq*s1)
 *     (t0, t1) = (t1, t0 - qq*t1)
 *
 * On success g = r0 is the GCD and x = s0, y = t0 satisfy
 * x*a + y*b = g.
 *
 * Complexity:
 *   Time: O(d^2 * n^2) in the worst case (Euclidean algorithm)
 *   Auxiliary memory: O(d) bignums
 *   Output memory: O(d) bignums
 *
 * @param[out] g GCD (may be NULL).
 * @param[out] x Bezout coefficient of a (may be NULL).
 * @param[out] y Bezout coefficient of b (may be NULL).
 * @param[in]  a First polynomial.
 * @param[in]  b Second polynomial.
 * @param[in]  fctx Coefficient ring.
 *
 * @return true  On success.
 * @return false If a division step fails (non-unit leading coefficient).
 */
bool bigpoly_xgcd(bigpoly* g, bigpoly* x, bigpoly* y, const bigpoly* a,
                  const bigpoly* b, field_ctx* fctx);

/**
 * @brief Initialize a Z_m[x]/(q) ring context.
 *
 * Deep-copies the modulus polynomial and stores a borrowed pointer to the
 * coefficient ring.
 *
 * Complexity:
 *   Time: O(d * n) where d = q->deg
 *   Auxiliary memory: O(1)
 *   Output memory: O(d) bignums
 *
 * @param[out] ctx Ring context to initialize.
 * @param[in]  q   Modulus polynomial (unit leading coefficient).
 * @param[in]  fctx Coefficient ring (borrowed).
 */
void poly_ring_init(poly_ring* ctx, const bigpoly* q, field_ctx* fctx);

/**
 * @brief Free a Z_m[x]/(q) ring context.
 *
 * The coefficient ring fctx is borrowed and NOT freed.
 *
 * Complexity:
 *   Time: O(d)
 *   Auxiliary memory: O(1)
 *   Output memory: O(1)
 *
 * @param[in,out] ctx Ring context to free.
 */
void poly_ring_free(poly_ring* ctx);

/**
 * @brief r = (a + b) mod q in Z_m[x]/(q).
 *
 * Complexity:
 *   Time: O(d * n) where d = max(a->deg, b->deg, q->deg)
 *   Auxiliary memory: O(d) bignums
 *   Output memory: O(q->deg) bignums
 *
 * @param[out] r Result.
 * @param[in]  a First element.
 * @param[in]  b Second element.
 * @param[in]  ctx Ring context.
 */
void poly_ring_add(bigpoly* r, const bigpoly* a, const bigpoly* b,
                   poly_ring* ctx);

/**
 * @brief r = (a - b) mod q in Z_m[x]/(q).
 *
 * Complexity:
 *   Time: O(d * n) where d = max(a->deg, b->deg, q->deg)
 *   Auxiliary memory: O(d) bignums
 *   Output memory: O(q->deg) bignums
 *
 * @param[out] r Result.
 * @param[in]  a First element.
 * @param[in]  b Second element.
 * @param[in]  ctx Ring context.
 */
void poly_ring_sub(bigpoly* r, const bigpoly* a, const bigpoly* b,
                   poly_ring* ctx);

/**
 * @brief r = (-a) mod q in Z_m[x]/(q).
 *
 * Complexity:
 *   Time: O(d * n) where d = max(a->deg, q->deg)
 *   Auxiliary memory: O(d) bignums
 *   Output memory: O(q->deg) bignums
 *
 * @param[out] r Result.
 * @param[in]  a Element.
 * @param[in]  ctx Ring context.
 */
void poly_ring_neg(bigpoly* r, const bigpoly* a, poly_ring* ctx);

/**
 * @brief r = (a * b) mod q in Z_m[x]/(q).
 *
 * Full schoolbook product over Z_m followed by reduction mod q.
 *
 * Complexity:
 *   Time: O((a->deg + 1) * (b->deg + 1) * n^2) for the product plus
 *         O((a->deg + b->deg) * q->deg * n^2) for the reduction
 *   Auxiliary memory: O(a->deg + b->deg) bignums
 *   Output memory: O(q->deg) bignums
 *
 * @param[out] r Result.
 * @param[in]  a First element.
 * @param[in]  b Second element.
 * @param[in]  ctx Ring context.
 */
void poly_ring_mul(bigpoly* r, const bigpoly* a, const bigpoly* b,
                   poly_ring* ctx);

/**
 * @brief r = a^(-1) mod q in Z_m[x]/(q).
 *
 * Uses the extended GCD: x*a + y*q = g. If g is a nonzero constant then
 * a is a unit and its inverse is x * g^(-1) mod q.
 *
 * Complexity:
 *   Time: O(d^2 * n^2) for the xgcd plus O(d * n) for the normalisation
 *   Auxiliary memory: O(d) bignums
 *   Output memory: O(q->deg) bignums
 *
 * @param[out] r Receives the inverse if a is a unit.
 * @param[in]  a Element to invert.
 * @param[in]  ctx Ring context.
 *
 * @return true  If a is a unit; r is set.
 * @return false If a and q are not coprime; r is left unchanged.
 */
bool poly_ring_inv(bigpoly* r, const bigpoly* a, poly_ring* ctx);

/**
 * @brief r = (a * b^(-1)) mod q in Z_m[x]/(q).
 *
 * Complexity:
 *   Time: see poly_ring_inv() plus one poly_ring_mul()
 *   Auxiliary memory: O(d) bignums
 *   Output memory: O(q->deg) bignums
 *
 * @param[out] r Receives the quotient if b is a unit.
 * @param[in]  a Numerator.
 * @param[in]  b Denominator.
 * @param[in]  ctx Ring context.
 *
 * @return true  If b is a unit; r is set.
 * @return false If b is not a unit; r is left unchanged.
 */
bool poly_ring_div(bigpoly* r, const bigpoly* a, const bigpoly* b,
                   poly_ring* ctx);

/**
 * @brief r = a^e mod q in Z_m[x]/(q).
 *
 * Left-to-right binary exponentiation with ring multiplication.
 *
 * Complexity:
 *   Time: O(e_bits * (q->deg^2 * n^2)) where e_bits = bit length of e
 *   Auxiliary memory: O(q->deg) bignums
 *   Output memory: O(q->deg) bignums
 *
 * @param[out] r Result.
 * @param[in]  a Base.
 * @param[in]  e Exponent (plain, nonnegative).
 * @param[in]  ctx Ring context.
 */
void poly_ring_pow(bigpoly* r, const bigpoly* a, const bignum* e,
                   poly_ring* ctx);

/**
 * @brief Test whether a ring element is zero.
 *
 * Complexity:
 *   Time: O(d) where d = a->deg
 *   Auxiliary memory: O(1)
 *   Output memory: O(1)
 *
 * @param[in] a Ring element.
 *
 * @return true If a is the zero element.
 */
bool poly_ring_is_zero(const bigpoly* a);

/**
 * @brief Test whether two ring elements are equal.
 *
 * Compares coefficient by coefficient up to the larger degree, treating
 * missing coefficients as zero, so elements need not be trimmed.
 *
 * Complexity:
 *   Time: O(d * n) where d = max(a->deg, b->deg)
 *   Auxiliary memory: O(1)
 *   Output memory: O(1)
 *
 * @param[in] a First element.
 * @param[in] b Second element.
 *
 * @return true If a and b are equal ring elements.
 */
bool poly_ring_equal(const bigpoly* a, const bigpoly* b);

#endif
