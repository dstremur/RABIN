#ifndef RZFIELD_H
#define RZFIELD_H

#include "rabin_errors.h"
#include "rpol.h"
#include "rz.h"

typedef struct {
  rz_t p;          // The prime modulus
  rz_t r_sq;       // R^2 mod p (used for entering Montgomery form)
  uint64_t p_inv;  // -p^-1 mod 2^64 (Montgomery reduction constant)
  size_t limbs;    // Number of machine words p takes up
} fp_ctx_t;

typedef struct {
  rz_t a;          // Curve parameter 'a'
  rz_t b;          // Curve parameter 'b'
  fp_ctx_t field;  // The underlying finite field
} ec_curve_t;

// TODO: Not yet implemented (placeholder declaration).
rabin_err_t ec_point_add(rz_t* res, const rz_t* p, const rz_t* q,
                         const ec_curve_t* curve);

/*
 * field_ctx_t: arithmetic in the ring Z_m = Z / m Z.
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
  rz_t* q;          /* modulus m (owned heap copy)            */
  bool mont;        /* true if the Montgomery path is active  */
  rz_mont_ctx mctx; /* valid only if mont                     */
} field_ctx_t;

/*
 * poly_ring_t: the ring of polynomials over Z_m modulo q, i.e. Z_m[x]/(q).
 *
 * Coefficients are field elements of fctx (internal representation). The
 * modulus polynomial q is expected to have a unit leading coefficient
 * (in particular it is usually taken monic) so that reduction is always
 * possible.
 */
typedef struct {
  rpol_t* q;         /* modulus polynomial (owned heap copy)    */
  field_ctx_t* fctx; /* coefficient ring Z_m (borrowed, not owned) */
} poly_ring_t;

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
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR,
 * RABIN_ERR_INVALID_ARG, or RABIN_ERR_OUT_OF_MEMORY.
 */
rabin_err_t field_ctx_init(field_ctx_t* ctx, const rz_t* m);

/**
 * @brief Free a \f$Z_m\f$ context.
 *
 * Complexity:
 *   - Time: \f$O(n)\f$
 *   - Auxiliary memory: \f$O(1)\f$
 *   - Output memory: \f$O(1)\f$
 *
 * @param[in,out] ctx Context to free.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR.
 */
rabin_err_t field_ctx_clear(field_ctx_t* ctx);

/**
 * @brief Convert a plain value into the internal representation.
 *
 * Montgomery path: \f$r = a \cdot R \bmod m\f$. Plain path: \f$r = a \bmod
 * m\f$ in \f$[0, m)\f$.
 *
 * Complexity:
 *   - Time: \f$O(n^2)\f$ (Montgomery) or \f$O(n^2)\f$ (one division, plain)
 *   - Auxiliary memory: \f$O(n)\f$ limbs
 *   - Output memory: \f$O(n)\f$ limbs
 *
 * @param[out] r Result in the internal representation.
 * @param[in]  a Plain value.
 * @param[in]  ctx Initialized context.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR /
 * RABIN_ERR_OUT_OF_MEMORY.
 */
rabin_err_t field_in(rz_t* r, const rz_t* a, field_ctx_t* ctx);

/**
 * @brief Convert an internal value back to the plain form \f$[0, m)\f$.
 *
 * Montgomery path: \f$r = a \cdot R^{-1} \bmod m\f$. Plain path: \f$r = a\f$.
 *
 * Complexity:
 *   - Time: \f$O(n^2)\f$ (Montgomery) or \f$O(n)\f$ (plain copy)
 *   - Auxiliary memory: \f$O(n)\f$ limbs
 *   - Output memory: \f$O(n)\f$ limbs
 *
 * @param[out] r Result in the plain form.
 * @param[in]  a Value in the internal representation.
 * @param[in]  ctx Initialized context.
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR /
 * RABIN_ERR_OUT_OF_MEMORY.
 */
rabin_err_t field_out(rz_t* r, const rz_t* a, field_ctx_t* ctx);

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
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR /
 * RABIN_ERR_OUT_OF_MEMORY.
 */
rabin_err_t field_set_u64(rz_t* r, u64 v, field_ctx_t* ctx);

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
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR /
 * RABIN_ERR_OUT_OF_MEMORY.
 */
rabin_err_t field_add(rz_t* r, const rz_t* a, const rz_t* b, field_ctx_t* ctx);

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
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR /
 * RABIN_ERR_OUT_OF_MEMORY.
 */
rabin_err_t field_sub(rz_t* r, const rz_t* a, const rz_t* b, field_ctx_t* ctx);

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
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR /
 * RABIN_ERR_OUT_OF_MEMORY.
 */
rabin_err_t field_neg(rz_t* r, const rz_t* a, field_ctx_t* ctx);

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
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR /
 * RABIN_ERR_OUT_OF_MEMORY.
 */
rabin_err_t field_mul(rz_t* r, const rz_t* a, const rz_t* b, field_ctx_t* ctx);

/**
 * @brief \f$r = a^{-1} \bmod m\f$.
 *
 * The inverse exists exactly when \f$a\f$ is a unit of \f$Z_m\f$ (\f$\gcd(a, m)
 * = 1\f$); for a prime modulus that is every nonzero element. The Montgomery
 * path (odd \f$m\f$) converts the operand to the plain domain and uses the
 * fast binary extended GCD (rz_mod_inverse()); the plain path (even
 * \f$m\f$) uses the classical extended Euclidean algorithm
 * (u64_mod_inverse_euclid()), because the binary GCD is only correct for
 * odd moduli.
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
 * @return RABIN_SUCCESS if \f$a\f$ is a unit and \f$r\f$ is set, or
 * RABIN_ERR_INVALID_ARG if \f$a\f$ is zero or not coprime to \f$m\f$
 * (\f$r\f$ left unchanged).
 */
rabin_err_t field_inv(rz_t* r, const rz_t* a, field_ctx_t* ctx);

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
 * @return RABIN_SUCCESS if \f$b\f$ is a unit and \f$r\f$ is set, or
 * RABIN_ERR_INVALID_ARG if \f$b\f$ is zero or not coprime to \f$m\f$
 * (\f$r\f$ left unchanged).
 */
rabin_err_t field_div(rz_t* r, const rz_t* a, const rz_t* b, field_ctx_t* ctx);

/**
 * @brief \f$r = a^e \bmod m\f$.
 *
 * Montgomery path: binary exponentiation with REDC (the base is already
 * in the Montgomery domain). Plain path: rz_mod_exp(), which handles even
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
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR.
 */
rabin_err_t field_pow(rz_t* r, const rz_t* a, const rz_t* e, field_ctx_t* ctx);

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
bool field_is_zero(const rz_t* a);

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
bool field_equal(const rz_t* a, const rz_t* b);

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
 * @return RABIN_SUCCESS on success, RABIN_ERR_DIV_BY_ZERO if \f$b\f$ is the
 * zero polynomial, or RABIN_ERR_INVALID_ARG if the leading coefficient of
 * \f$b\f$ is not a unit. On both errors \f$r = a\f$ and \f$q = 0\f$ are
 * written (when the corresponding outputs are non-NULL).
 */
rabin_err_t rpol_divmod(rpol_t* q, rpol_t* r, const rpol_t* a, const rpol_t* b,
                        field_ctx_t* fctx);

/**
 * @brief Extended polynomial GCD in \f$Z_m[x]\f$.
 *
 * Let \f$d =\f$ max(a->deg, b->deg).
 *
 * Iterated rpol_divmod() tracking the Bezout coefficients:
 *
 *   \f$(r_0, s_0, t_0) = (a, 1, 0)\f$
 *   \f$(r_1, s_1, t_1) = (b, 0, 1)\f$
 *   while \f$r_1 \neq 0\f$:
 *     \f$(qq, r_2) = divmod(r_0, r_1)\f$
 *     \f$(r_0, r_1) = (r_1, r_2)\f$
 *     \f$(s_0, s_1) = (s_1, s_0 - qq \cdot s_1)\f$
 *     \f$(t_0, t_1) = (t_1, t_0 - qq \cdot t_1)\f$
 *
 * On success \f$g = r_0\f$ is the GCD and \f$x = s_0\f$, \f$y = t_0\f$ satisfy
 * \f$x \cdot a + y \cdot b = g\f$.
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
 * @return RABIN_SUCCESS on success, RABIN_ERR_DIV_BY_ZERO or
 * RABIN_ERR_INVALID_ARG if a division step fails, or
 * RABIN_ERR_OUT_OF_MEMORY.
 */
rabin_err_t rpol_xgcd(rpol_t* g, rpol_t* x, rpol_t* y, const rpol_t* a,
                      const rpol_t* b, field_ctx_t* fctx);

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
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR /
 * RABIN_ERR_OUT_OF_MEMORY.
 */
rabin_err_t poly_ring_init(poly_ring_t* ctx, const rpol_t* q,
                           field_ctx_t* fctx);

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
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR.
 */
rabin_err_t poly_ring_clear(poly_ring_t* ctx);

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
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR /
 * RABIN_ERR_OUT_OF_MEMORY.
 */
rabin_err_t poly_ring_add(rpol_t* r, const rpol_t* a, const rpol_t* b,
                          poly_ring_t* ctx);

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
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR /
 * RABIN_ERR_OUT_OF_MEMORY.
 */
rabin_err_t poly_ring_sub(rpol_t* r, const rpol_t* a, const rpol_t* b,
                          poly_ring_t* ctx);

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
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR /
 * RABIN_ERR_OUT_OF_MEMORY.
 */
rabin_err_t poly_ring_neg(rpol_t* r, const rpol_t* a, poly_ring_t* ctx);

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
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR /
 * RABIN_ERR_OUT_OF_MEMORY.
 */
rabin_err_t poly_ring_mul(rpol_t* r, const rpol_t* a, const rpol_t* b,
                          poly_ring_t* ctx);

/**
 * @brief \f$r = a^{-1} \bmod q\f$ in \f$Z_m[x]/(q)\f$.
 *
 * Uses the extended GCD: \f$x \cdot a + y \cdot q = g\f$. If \f$g\f$ is a
 * nonzero constant then \f$a\f$ is a unit and its inverse is \f$x \cdot
 * g^{-1} \bmod q\f$.
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
 * @return RABIN_SUCCESS if \f$a\f$ is a unit and \f$r\f$ is set, or
 * RABIN_ERR_INVALID_ARG if \f$a\f$ and \f$q\f$ are not coprime (\f$r\f$ left
 * unchanged).
 */
rabin_err_t poly_ring_inv(rpol_t* r, const rpol_t* a, poly_ring_t* ctx);

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
 * @return RABIN_SUCCESS if \f$b\f$ is a unit and \f$r\f$ is set, or
 * RABIN_ERR_INVALID_ARG if \f$b\f$ is not a unit (\f$r\f$ left unchanged).
 */
rabin_err_t poly_ring_div(rpol_t* r, const rpol_t* a, const rpol_t* b,
                          poly_ring_t* ctx);

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
 * @return RABIN_SUCCESS on success, or RABIN_ERR_NULL_PTR /
 * RABIN_ERR_OUT_OF_MEMORY.
 */
rabin_err_t poly_ring_pow(rpol_t* r, const rpol_t* a, const rz_t* e,
                          poly_ring_t* ctx);

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
bool poly_ring_is_zero(const rpol_t* a);

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
bool poly_ring_equal(const rpol_t* a, const rpol_t* b);

#endif
