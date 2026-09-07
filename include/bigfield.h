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

/* ---- context management -------------------------------------------------*/

void field_ctx_init(field_ctx* ctx, const bignum* m);
void field_ctx_free(field_ctx* ctx);

/* ---- Z_m element operations (internal representation) -------------------*/

/* r = a converted into the internal representation (a mod m in [0, m)). */
void field_in(bignum* r, const bignum* a, field_ctx* ctx);
/* r = a converted back to the plain form [0, m). */
void field_out(bignum* r, const bignum* a, field_ctx* ctx);
/* r = v (a 64-bit value) in the internal representation. */
void field_set_u64(bignum* r, u64 v, field_ctx* ctx);
/* r = (a + b) mod m. */
void field_add(bignum* r, const bignum* a, const bignum* b, field_ctx* ctx);
/* r = (a - b) mod m. */
void field_sub(bignum* r, const bignum* a, const bignum* b, field_ctx* ctx);
/* r = (-a) mod m. */
void field_neg(bignum* r, const bignum* a, field_ctx* ctx);
/* r = (a * b) mod m. */
void field_mul(bignum* r, const bignum* a, const bignum* b, field_ctx* ctx);
/* r = a^(-1) mod m; returns false if a is not a unit. */
bool field_inv(bignum* r, const bignum* a, field_ctx* ctx);
/* r = (a * b^(-1)) mod m; returns false if b is not a unit. */
bool field_div(bignum* r, const bignum* a, const bignum* b, field_ctx* ctx);
/* r = a^e mod m (e a plain nonnegative bignum exponent). */
void field_pow(bignum* r, const bignum* a, const bignum* e, field_ctx* ctx);
/* true if a == 0. */
bool field_is_zero(const bignum* a);
/* true if a == b. */
bool field_equal(const bignum* a, const bignum* b);

/* ---- polynomial arithmetic over the coefficient ring --------------------*/

/*
 * q = a / b, r = a mod b in Z_m[x]. Requires the leading coefficient of b
 * to be a unit of Z_m (always true for a monic b). Returns false (and sets
 * r = a, q = 0) if b is zero or its leading coefficient is not invertible.
 * q and/or r may be NULL to discard that output.
 */
bool bigpoly_divmod(bigpoly* q, bigpoly* r, const bigpoly* a, const bigpoly* b,
                    field_ctx* fctx);

/*
 * Extended GCD in Z_m[x]: finds g, x, y with x*a + y*b = g where g is the
 * (monic-normalised) GCD. Returns false if a division step fails (a
 * non-unit leading coefficient). g, x, y may be NULL to discard.
 */
bool bigpoly_xgcd(bigpoly* g, bigpoly* x, bigpoly* y, const bigpoly* a,
                  const bigpoly* b, field_ctx* fctx);

/* ---- ring Z_m[x]/(q) operations -----------------------------------------*/

void poly_ring_init(poly_ring* ctx, const bigpoly* q, field_ctx* fctx);
void poly_ring_free(poly_ring* ctx);
/* r = (a + b) mod q. */
void poly_ring_add(bigpoly* r, const bigpoly* a, const bigpoly* b,
                   poly_ring* ctx);
/* r = (a - b) mod q. */
void poly_ring_sub(bigpoly* r, const bigpoly* a, const bigpoly* b,
                   poly_ring* ctx);
/* r = (-a) mod q. */
void poly_ring_neg(bigpoly* r, const bigpoly* a, poly_ring* ctx);
/* r = (a * b) mod q. */
void poly_ring_mul(bigpoly* r, const bigpoly* a, const bigpoly* b,
                   poly_ring* ctx);
/* r = a^(-1) mod q; returns false if a is not a unit of the ring. */
bool poly_ring_inv(bigpoly* r, const bigpoly* a, poly_ring* ctx);
/* r = (a * b^(-1)) mod q; returns false if b is not a unit. */
bool poly_ring_div(bigpoly* r, const bigpoly* a, const bigpoly* b,
                   poly_ring* ctx);
/* r = a^e mod q (e a plain nonnegative bignum exponent). */
void poly_ring_pow(bigpoly* r, const bigpoly* a, const bignum* e,
                   poly_ring* ctx);
/* true if a is the zero element. */
bool poly_ring_is_zero(const bigpoly* a);
/* true if a and b are equal ring elements. */
bool poly_ring_equal(const bigpoly* a, const bigpoly* b);

#endif
