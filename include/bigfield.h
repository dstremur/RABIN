#ifndef BIGFIELD_H
#define BIGFIELD_H

#include "bigcore.h"
#include "bigpoly.h"

typedef struct {
  bignum p;        // The prime modulus
  bignum r_sq;     // R^2 mod p (used for entering Montgomery form)
  uint64_t p_inv;  // -p^-1 mod 2^64 (Montgomery reduction constant)
  size_t limbs;    // Number of machine words p takes up
} fp_ctx;

typedef struct {
  bignum a;      // Curve parameter 'a'
  bignum b;      // Curve parameter 'b'
  fp_ctx field;  // The underlying finite field
} ec_curve_t;

void ec_point_add(ec_point res, ec_point p, ec_point q, const ec_curve* curve);

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
 * @brief Initialize a \f$Z_m\f$ context.
 *
 * Copies the modulus \f$m\f$ and selects the arithmetic path: the Montgomery
 * path is used iff \f$m\f$ is odd and \f$> 1\f$ (REDC requires an odd modulus),
 * the plain multiply + reduce path otherwise.
 *
 * Complexity:
 *   - Time: \f$O(n^2)\f$ for the Montgomery setup (one division for \f$R \bmod
 * m\f$ and one for \f$R^2 \bmod m\f$), \f$O(n)\f$ otherwise, where \f$n =\f$
 * m->size
 *   - Auxiliary memory: \f$O(n)\f$ limbs
 *   - Output memory: \f$O(n)\f$ limbs per context field
 *
 * @param[out] ctx Context to initialize.
 * @param[in]  m   Modulus (must be > 1).
 */
void field_ctx_init(field_ctx* ctx, const bignum* m);

/**
 * @brief Free a \f$Z_m\f$ context.
 *
 * Complexity:
 *   - Time: \f$O(n)\f$
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(1)\f$
 *
 * @param[in,out] ctx Context to free.
 */
void field_ctx_free(field_ctx* ctx);

/**
 * @brief Convert a plain value into the internal representation.
 *
 * Montgomery path: \f$r = a\cdotR \bmod m\f$. Plain path: \f$r = a \bmod m\f$
 * in \f$[0, m)\f$.
 *
 * Complexity:
 *   - Time: \f$O(n^2)\f$ (Montgomery) or \f$O(n^2)\f$ (one division, plain)
 *   - Auxiliary memory: \f$O(n)\f$ limbs
 *   - Output memory: \f$O(n)\f$ limbs
 *
 * @param[out] r Result in the internal representation.
 * @param[in]  a Plain value.
 * @param[in]  ctx Initialized context.
 */
void field_in(bignum* r, const bignum* a, field_ctx* ctx);

/**
 * @brief Convert an internal value back to the plain form \f$[0, m)\f$.
 *
 * Montgomery path: \f$r = a\cdotR^{-1} \bmod m\f$. Plain path: \f$r = a\f$.
 *
 * Complexity:
 *   - Time: \f$O(n^2)\f$ (Montgomery) or \f$O(n)\f$ (plain copy)
 *   - Auxiliary memory: \f$O(n)\f$ limbs
 *   - Output memory: \f$O(n)\f$ limbs
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
 *   - Time: \f$O(n^2)\f$ (Montgomery) or \f$O(1)\f$ (plain)
 *   - Auxiliary memory: \f$O(n)\f$ limbs
 *   - Output memory: \f$O(n)\f$ limbs
 *
 * @param[out] r Result.
 * @param[in]  v 64-bit value.
 * @param[in]  ctx Initialized context.
 */
void field_set_u64(bignum* r, u64 v, field_ctx* ctx);

/**
 * @brief \f$r = (a + b) \bmod m\f$.
 *
 * Both operands are in \f$[0, m)\f$, so \f$a + b < 2m\f$ and a single
 * conditional subtraction suffices. Addition is representation independent.
 *
 * Complexity:
 *   - Time: \f$O(n)\f$
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(n)\f$ limbs
 *
 * @param[out] r Result.
 * @param[in]  a First operand.
 * @param[in]  b Second operand.
 * @param[in]  ctx Initialized context.
 */
void field_add(bignum* r, const bignum* a, const bignum* b, field_ctx* ctx);

/**
 * @brief \f$r = (a - b) \bmod m\f$.
 *
 * Both operands are in \f$[0, m)\f$, so \f$a - b\f$ is in \f$(-m, m)\f$; if the
 * result is negative, \f$m\f$ is added to bring it into \f$[0, m)\f$.
 * Subtraction is representation independent.
 *
 * Complexity:
 *   - Time: \f$O(n)\f$
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(n)\f$ limbs
 *
 * @param[out] r Result.
 * @param[in]  a First operand.
 * @param[in]  b Second operand.
 * @param[in]  ctx Initialized context.
 */
void field_sub(bignum* r, const bignum* a, const bignum* b, field_ctx* ctx);

/**
 * @brief \f$r = (-a) \bmod m\f$.
 *
 * Computes \f$m - a\f$ and reduces the result (which equals \f$m\f$ when \f$a =
 * 0\f$) back into \f$[0, m)\f$.
 *
 * Complexity:
 *   - Time: \f$O(n)\f$
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(n)\f$ limbs
 *
 * @param[out] r Result.
 * @param[in]  a Operand.
 * @param[in]  ctx Initialized context.
 */
void field_neg(bignum* r, const bignum* a, field_ctx* ctx);

/**
 * @brief \f$r = (a \cdot b) \bmod m\f$.
 *
 * Montgomery path: one REDC of the product. Plain path: full multiply
 * followed by one reduction.
 *
 * Complexity:
 *   - Time: \f$O(n^2)\f$
 *   - Auxiliary memory: \f$O(n)\f$ limbs
 *   - Output memory: \f$O(n)\f$ limbs
 *
 * @param[out] r Result.
 * @param[in]  a First operand.
 * @param[in]  b Second operand.
 * @param[in]  ctx Initialized context.
 */
void field_mul(bignum* r, const bignum* a, const bignum* b, field_ctx* ctx);

/**
 * @brief \f$r = a^{-1} \bmod m\f$.
 *
 * The inverse exists exactly when \f$a\f$ is a unit of \f$Z_m\f$ (\f$\gcd(a, m)
 * = 1\f$); for a prime modulus that is every nonzero element. The Montgomery
 * path (odd \f$m\f$) uses the fast binary extended GCD (bn_mod_inverse) after
 * converting the operand to the plain domain. The plain path (even \f$m\f$)
 * uses the classical extended Euclidean algorithm (mod_inverse_euclid),
 * because the binary GCD is only correct for odd moduli.
 *
 * Complexity:
 *   - Time: \f$O(n^2)\f$ (Montgomery) or \f$O(n^3)\f$ (plain)
 *   - Auxiliary memory: \f$O(n)\f$ limbs
 *   - Output memory: \f$O(n)\f$ limbs
 *
 * @param[out] r Receives the inverse if it exists.
 * @param[in]  a Value to invert.
 * @param[in]  ctx Initialized context.
 *
 * @return true  If \f$a\f$ is a unit; \f$r\f$ is set.
 * @return false If \f$a\f$ is zero or not coprime to \f$m\f$; \f$r\f$ is left
 * unchanged.
 */
bool field_inv(bignum* r, const bignum* a, field_ctx* ctx);

/**
 * @brief \f$r = (a \cdot b^{-1}) \bmod m\f$.
 *
 * Complexity:
 *   - Time: \f$O(n^2)\f$
 *   - Auxiliary memory: \f$O(n)\f$ limbs
 *   - Output memory: \f$O(n)\f$ limbs
 *
 * @param[out] r Receives the quotient if \f$b\f$ is a unit.
 * @param[in]  a Numerator.
 * @param[in]  b Denominator.
 * @param[in]  ctx Initialized context.
 *
 * @return true  If \f$b\f$ is a unit; \f$r\f$ is set.
 * @return false If \f$b\f$ is zero or not coprime to \f$m\f$; \f$r\f$ is left
 * unchanged.
 */
bool field_div(bignum* r, const bignum* a, const bignum* b, field_ctx* ctx);

/**
 * @brief \f$r = a^e \bmod m\f$.
 *
 * Montgomery path: binary exponentiation with REDC (the base is already
 * in the Montgomery domain). Plain path: bn_mod_exp, which handles even
 * moduli with its slow multiply + divide loop.
 *
 * Complexity:
 *   - Time: \f$O(e \cdot n^2)\f$ where \f$e =\f$ e->size in limbs
 *   - Auxiliary memory: \f$O(n)\f$ limbs
 *   - Output memory: \f$O(n)\f$ limbs
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
 *   - Time: \f$O(1)\f$
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(1)\f$
 *
 * @param[in] a Field element.
 *
 * @return true If \f$a = 0\f$.
 */
bool field_is_zero(const bignum* a);

/**
 * @brief Test whether two field elements are equal.
 *
 * Complexity:
 *   - Time: \f$O(n)\f$
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(1)\f$
 *
 * @param[in] a First element.
 * @param[in] b Second element.
 *
 * @return true If \f$a = b\f$.
 */
bool field_equal(const bignum* a, const bignum* b);

/**
 * @brief Polynomial long division in \f$Z_m[x]\f$: \f$q = a / b\f$, \f$r = a
 * \bmod b\f$.
 *
 * Let \f$d_a =\f$ a->deg and \f$d_b =\f$ b->deg.
 *
 * Standard synthetic division: while \f$\deg(r) \ge \deg(b)\f$, the quotient
 * term is \f$factor = r[\deg r] \cdot (b[\deg b])^{-1}\f$ and r -= factor *
 * x^shift * b. The inverse of the leading coefficient of \f$b\f$ is computed
 * once.
 *
 * Requires the leading coefficient of \f$b\f$ to be a unit of \f$Z_m\f$ (always
 * true for a monic \f$b\f$). If \f$b\f$ is zero or its leading coefficient is
 * not invertible, returns false with \f$r = a\f$ and \f$q = 0\f$.
 *
 * Complexity:
 *   - Time: \f$O((d_a - d_b + 1) \cdot d_b \cdot n^2)\f$ field multiplications
 * plus one
 *           \f$O(n^2)\f$ inversion
 *   - Auxiliary memory: \f$O(d_a)\f$ bignums
 *   - Output memory: \f$O(d_a)\f$ bignums
 *
 * @param[out] q Quotient (may be NULL).
 * @param[out] r Remainder (may be NULL).
 * @param[in]  a Dividend.
 * @param[in]  b Divisor (nonzero, unit leading coefficient).
 * @param[in]  fctx Coefficient ring.
 *
 * @return true  On success.
 * @return false If \f$b\f$ is zero or its leading coefficient is not a unit.
 */
bool bigpoly_divmod(bigpoly* q, bigpoly* r, const bigpoly* a, const bigpoly* b,
                    field_ctx* fctx);

/**
 * @brief Extended polynomial GCD in \f$Z_m[x]\f$.
 *
 * Let \f$d =\f$ max(a->deg, b->deg).
 *
 * Iterated bigpoly_divmod() tracking the Bezout coefficients:
 *
 *   \f$(r_0, s_0, t_0) = (a, 1, 0)\f$
 *   \f$(r_1, s_1, t_1) = (b, 0, 1)\f$
 *   while \f$r_1 \neq 0\f$:
 *     \f$(qq, r_2) = divmod(r_0, r_1)\f$
 *     \f$(r_0, r_1) = (r_1, r_2)\f$
 *     \f$(s_0, s_1) = (s_1, s_0 - qq\cdots_1)\f$
 *     \f$(t_0, t_1) = (t_1, t_0 - qq\cdott_1)\f$
 *
 * On success \f$g = r_0\f$ is the GCD and \f$x = s_0\f$, \f$y = t_0\f$ satisfy
 * \f$x\cdota + y\cdotb = g\f$.
 *
 * Complexity:
 *   - Time: \f$O(d^2 \cdot n^2)\f$ in the worst case (Euclidean algorithm)
 *   - Auxiliary memory: \f$O(d)\f$ bignums
 *   - Output memory: \f$O(d)\f$ bignums
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
 * @brief Initialize a \f$Z_m[x]/(q)\f$ ring context.
 *
 * Deep-copies the modulus polynomial and stores a borrowed pointer to the
 * coefficient ring.
 *
 * Complexity:
 *   - Time: \f$O(d \cdot n)\f$ where \f$d =\f$ q->deg
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(d)\f$ bignums
 *
 * @param[out] ctx Ring context to initialize.
 * @param[in]  q   Modulus polynomial (unit leading coefficient).
 * @param[in]  fctx Coefficient ring (borrowed).
 */
void poly_ring_init(poly_ring* ctx, const bigpoly* q, field_ctx* fctx);

/**
 * @brief Free a \f$Z_m[x]/(q)\f$ ring context.
 *
 * The coefficient ring fctx is borrowed and NOT freed.
 *
 * Complexity:
 *   - Time: \f$O(d)\f$
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(1)\f$
 *
 * @param[in,out] ctx Ring context to free.
 */
void poly_ring_free(poly_ring* ctx);

/**
 * @brief \f$r = (a + b) \bmod q\f$ in \f$Z_m[x]/(q)\f$.
 *
 * Let \f$d_a =\f$ a->deg, \f$d_b =\f$ b->deg, \f$d_q =\f$ q->deg.
 *
 * Complexity:
 *   - Time: \f$O(d \cdot n)\f$ where \f$d =\f$ max(\f$d_a, d_b, d_q\f$)
 *   - Auxiliary memory: \f$O(d)\f$ bignums
 *   - Output memory: \f$O(d_q)\f$ bignums
 *
 * @param[out] r Result.
 * @param[in]  a First element.
 * @param[in]  b Second element.
 * @param[in]  ctx Ring context.
 */
void poly_ring_add(bigpoly* r, const bigpoly* a, const bigpoly* b,
                   poly_ring* ctx);

/**
 * @brief \f$r = (a - b) \bmod q\f$ in \f$Z_m[x]/(q)\f$.
 *
 * Let \f$d_a =\f$ a->deg, \f$d_b =\f$ b->deg, \f$d_q =\f$ q->deg.
 *
 * Complexity:
 *   - Time: \f$O(d \cdot n)\f$ where \f$d =\f$ max(\f$d_a, d_b, d_q\f$)
 *   - Auxiliary memory: \f$O(d)\f$ bignums
 *   - Output memory: \f$O(d_q)\f$ bignums
 *
 * @param[out] r Result.
 * @param[in]  a First element.
 * @param[in]  b Second element.
 * @param[in]  ctx Ring context.
 */
void poly_ring_sub(bigpoly* r, const bigpoly* a, const bigpoly* b,
                   poly_ring* ctx);

/**
 * @brief \f$r = (-a) \bmod q\f$ in \f$Z_m[x]/(q)\f$.
 *
 * Let \f$d_a =\f$ a->deg and \f$d_q =\f$ q->deg.
 *
 * Complexity:
 *   - Time: \f$O(d \cdot n)\f$ where \f$d =\f$ max(\f$d_a, d_q\f$)
 *   - Auxiliary memory: \f$O(d)\f$ bignums
 *   - Output memory: \f$O(d_q)\f$ bignums
 *
 * @param[out] r Result.
 * @param[in]  a Element.
 * @param[in]  ctx Ring context.
 */
void poly_ring_neg(bigpoly* r, const bigpoly* a, poly_ring* ctx);

/**
 * @brief \f$r = (a \cdot b) \bmod q\f$ in \f$Z_m[x]/(q)\f$.
 *
 * Full schoolbook product over \f$Z_m\f$ followed by reduction \f$\bmod q\f$.
 *
 * Let \f$d_a =\f$ a->deg, \f$d_b =\f$ b->deg, \f$d_q =\f$ q->deg.
 *
 * Complexity:
 *   - Time: \f$O((d_a + 1) \cdot (d_b + 1) \cdot n^2)\f$ for the product plus
 *           \f$O((d_a + d_b) \cdot d_q \cdot n^2)\f$ for the reduction
 *   - Auxiliary memory: \f$O(d_a + d_b)\f$ bignums
 *   - Output memory: \f$O(d_q)\f$ bignums
 *
 * @param[out] r Result.
 * @param[in]  a First element.
 * @param[in]  b Second element.
 * @param[in]  ctx Ring context.
 */
void poly_ring_mul(bigpoly* r, const bigpoly* a, const bigpoly* b,
                   poly_ring* ctx);

/**
 * @brief \f$r = a^{-1} \bmod q\f$ in \f$Z_m[x]/(q)\f$.
 *
 * Uses the extended GCD: \f$x\cdota + y\cdotq = g\f$. If \f$g\f$ is a nonzero
 * constant then
 * \f$a\f$ is a unit and its inverse is \f$x \cdot g^{-1} \bmod q\f$.
 *
 * Let \f$d =\f$ a->deg and \f$d_q =\f$ q->deg.
 *
 * Complexity:
 *   - Time: \f$O(d^2 \cdot n^2)\f$ for the xgcd plus \f$O(d \cdot n)\f$ for the
 * normalisation
 *   - Auxiliary memory: \f$O(d)\f$ bignums
 *   - Output memory: \f$O(d_q)\f$ bignums
 *
 * @param[out] r Receives the inverse if a is a unit.
 * @param[in]  a Element to invert.
 * @param[in]  ctx Ring context.
 *
 * @return true  If \f$a\f$ is a unit; \f$r\f$ is set.
 * @return false If \f$a\f$ and \f$q\f$ are not coprime; \f$r\f$ is left
 * unchanged.
 */
bool poly_ring_inv(bigpoly* r, const bigpoly* a, poly_ring* ctx);

/**
 * @brief \f$r = (a \cdot b^{-1}) \bmod q\f$ in \f$Z_m[x]/(q)\f$.
 *
 * Let \f$d =\f$ max(a->deg, b->deg) and \f$d_q =\f$ q->deg.
 *
 * Complexity:
 *   - Time: see poly_ring_inv() plus one poly_ring_mul()
 *   - Auxiliary memory: \f$O(d)\f$ bignums
 *   - Output memory: \f$O(d_q)\f$ bignums
 *
 * @param[out] r Receives the quotient if \f$b\f$ is a unit.
 * @param[in]  a Numerator.
 * @param[in]  b Denominator.
 * @param[in]  ctx Ring context.
 *
 * @return true  If \f$b\f$ is a unit; \f$r\f$ is set.
 * @return false If \f$b\f$ is not a unit; \f$r\f$ is left unchanged.
 */
bool poly_ring_div(bigpoly* r, const bigpoly* a, const bigpoly* b,
                   poly_ring* ctx);

/**
 * @brief \f$r = a^e \bmod q\f$ in \f$Z_m[x]/(q)\f$.
 *
 * Left-to-right binary exponentiation with ring multiplication.
 *
 * Let \f$d_q =\f$ q->deg.
 *
 * Complexity:
 *   - Time: \f$O(e_{bits} \cdot (d_q^2 \cdot n^2))\f$ where \f$e_{bits} =\f$
 * bit length of \f$e\f$
 *   - Auxiliary memory: \f$O(d_q)\f$ bignums
 *   - Output memory: \f$O(d_q)\f$ bignums
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
 *   - Time: \f$O(d)\f$ where \f$d =\f$ a->deg
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(1)\f$
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
 *   - Time: \f$O(d \cdot n)\f$ where \f$d =\f$ max(a->deg, b->deg)
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(1)\f$
 *
 * @param[in] a First element.
 * @param[in] b Second element.
 *
 * @return true If a and b are equal ring elements.
 */
bool poly_ring_equal(const bigpoly* a, const bigpoly* b);

#endif
