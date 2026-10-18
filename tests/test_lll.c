/*
 * Extended LLL / MLLL tests.
 *
 * Goes beyond test_lll()/test_mlll() in tests/test_matrix.c by adding:
 *   - deterministic known-answer cases,
 *   - property-based checks over many shapes, entry sizes and delta values,
 *   - an independent lattice-equality check (HNF of the reduced basis must
 *     equal the HNF of the original basis),
 *   - edge cases (negatives, multi-limb entries, zero / duplicate columns,
 *     delta = 1, idempotency).
 *
 * The LLL *conditions* (size reduction |mu| <= 1/2 and the Lovasz
 * inequality) are verified with exact rational Gram-Schmidt, mirroring
 * test_check_lll() in tests/test_matrix.c.
 */
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "../include/rmat.h"
#include "../include/rpol.h"
#include "../include/rq.h"
#include "../include/rz.h"
#include "../include/rzrand.h"

/* ------------------------------------------------------------------ */
/* Small helpers                                                      */
/* ------------------------------------------------------------------ */

static double get_time(void)
{
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return ts.tv_sec + ts.tv_nsec * 1e-9;
}

static void set164(rmat_t* M, u64 r, u64 c, i64 v)
{
  rz_t t;
  rz_init(&t);
  rz_set_i64(&t, v);
  assert(rmat_set(M, &t, r, c) == RABIN_SUCCESS);
  rz_clear(&t);
}

/* ------------------------------------------------------------------ */
/* Exact rational vector helpers (plain arrays of rationals)          */
/* ------------------------------------------------------------------ */

static rq_t* lllv_new(u64 m)
{
  rq_t* v = calloc(m, sizeof(rq_t));
  assert(v != NULL);
  for (u64 i = 0; i < m; i++) {
    assert(rq_init(&v[i]) == RABIN_SUCCESS);
  }
  return v;
}

static void lllv_free(rq_t* v, u64 m)
{
  for (u64 i = 0; i < m; i++) {
    rq_clear(&v[i]);
  }
  free(v);
}

static void lllv_set_col(rq_t* v, const rmat_t* M, u64 c)
{
  for (u64 i = 0; i < M->rows; i++) {
    assert(rq_set(&v[i], RMAT_GET(M, i, c)) == RABIN_SUCCESS);
  }
}

static void lllv_dot(rq_t* r, const rq_t* a, const rq_t* b, u64 m)
{
  rq_t t;
  assert(rq_init(&t) == RABIN_SUCCESS);
  assert(rq_set_u64(r, 0) == RABIN_SUCCESS);
  for (u64 i = 0; i < m; i++) {
    assert(rq_mul(&t, &a[i], &b[i]) == RABIN_SUCCESS);
    assert(rq_add(r, r, &t) == RABIN_SUCCESS);
  }
  rq_clear(&t);
}

/* a - b >= 0 ? (rq values are normalized after every op) */
static bool lllv_cmp_ge(const rq_t* a, const rq_t* b)
{
  rq_t t;
  assert(rq_init(&t) == RABIN_SUCCESS);
  assert(rq_sub(&t, a, b) == RABIN_SUCCESS);
  bool r = !t.num.is_neg;
  rq_clear(&t);
  return r;
}

/*
 * Exact LLL check: the columns of B (m x p) form an LLL-reduced basis for
 * delta = dnum/dden. Runs exact rational Gram-Schmidt and verifies
 *   |mu[k][j]| <= 1/2                            (all j < k)
 *   ||b*_k||^2 >= (delta - mu[k][k-1]^2) * ||b*_{k-1}||^2   (if check_lovasz)
 *
 * The size-reduction condition is exact arithmetic and always required.
 * The Lovasz inequality is checked exactly when check_lovasz is true;
 * rmat_lll()/rmat_mlll() evaluate it with exact integer arithmetic for
 * every double delta, so passing true here asserts the guarantee the
 * algorithm provides (see strict_lovasz for how targets are chosen).
 */
static bool check_lll_conds(const rmat_t* B, u64 p, bool check_lovasz, i64 dnum,
                            i64 dden)
{
  u64 m = B->rows;
  if (p <= 1) return true;

  rq_t** bs = malloc(p * sizeof(rq_t*));
  rq_t** bk = malloc(p * sizeof(rq_t*));
  assert(bs != NULL && bk != NULL);
  for (u64 k = 0; k < p; k++) {
    bs[k] = lllv_new(m);
    bk[k] = lllv_new(m);
  }
  rq_t* mu = lllv_new(p * p);  // flattened: mu[k * p + j]

  bool ok = true;
  rq_t num, two, oneq;
  assert(rq_init_multi(&num, &two, &oneq, NULL) == RABIN_SUCCESS);
  assert(rq_set_u64(&two, 2) == RABIN_SUCCESS);
  assert(rq_set_u64(&oneq, 1) == RABIN_SUCCESS);

  for (u64 k = 0; k < p && ok; k++) {
    lllv_set_col(bk[k], B, k);
    for (u64 j = k; j-- > 0;) {
      rq_t n2, d2;
      assert(rq_init(&n2) == RABIN_SUCCESS);
      assert(rq_init(&d2) == RABIN_SUCCESS);
      lllv_dot(&n2, bk[k], bs[j], m);
      lllv_dot(&d2, bs[j], bs[j], m);
      assert(rq_div(&mu[k * p + j], &n2, &d2) == RABIN_SUCCESS);
      for (u64 i = 0; i < m; i++) {
        assert(rq_mul(&n2, &mu[k * p + j], &bs[j][i]) == RABIN_SUCCESS);
        assert(rq_sub(&bk[k][i], &bk[k][i], &n2) == RABIN_SUCCESS);
      }
      rq_clear(&n2);
      rq_clear(&d2);
    }
    for (u64 i = 0; i < m; i++) {
      assert(rq_copy(&bs[k][i], &bk[k][i]) == RABIN_SUCCESS);
    }

    // size reduction: |2 * mu| <= 1
    for (u64 j = 0; j < k; j++) {
      assert(rq_mul(&num, &two, &mu[k * p + j]) == RABIN_SUCCESS);
      if (!lllv_cmp_ge(&oneq, &num)) ok = false;
    }
  }

  // Lovasz condition (exact, only when required by the caller)
  if (check_lovasz) {
    rq_t s1, s2, dq, dq2, rhs, dn, dd;
    rz_t dnumz, ddenz;
    assert(rq_init_multi(&s1, &s2, &dq, &dq2, &rhs, &dn, &dd, NULL) ==
           RABIN_SUCCESS);
    assert(rz_init_multi(&dnumz, &ddenz, NULL) == RABIN_SUCCESS);
    assert(rz_set_i64(&dnumz, dnum) == RABIN_SUCCESS);
    assert(rz_set_i64(&ddenz, dden) == RABIN_SUCCESS);
    assert(rq_set(&dn, &dnumz) == RABIN_SUCCESS);
    assert(rq_set(&dd, &ddenz) == RABIN_SUCCESS);
    assert(rq_div(&dq, &dn, &dd) == RABIN_SUCCESS);
    for (u64 k = 1; k < p && ok; k++) {
      lllv_dot(&s2, bs[k], bs[k], m);
      lllv_dot(&s1, bs[k - 1], bs[k - 1], m);
      assert(rq_sqr(&dq2, &mu[k * p + (k - 1)]) == RABIN_SUCCESS);
      assert(rq_sub(&dq, &dq, &dq2) == RABIN_SUCCESS);
      assert(rq_mul(&rhs, &dq, &s1) == RABIN_SUCCESS);
      if (!lllv_cmp_ge(&s2, &rhs)) ok = false;
    }
    rz_clear_multi(&dnumz, &ddenz, NULL);
    assert(rq_clear_multi(&s1, &s2, &dq, &dq2, &rhs, &dn, &dd, NULL) ==
           RABIN_SUCCESS);
  }

  assert(rq_clear_multi(&num, &two, &oneq, NULL) == RABIN_SUCCESS);
  for (u64 i = 0; i < p * p; i++) {
    rq_clear(&mu[i]);
  }
  free(mu);
  for (u64 k = 0; k < p; k++) {
    lllv_free(bk[k], m);
    lllv_free(bs[k], m);
  }
  free(bk);
  free(bs);
  return ok;
}

/*
 * Is the integer vector b (length m) an integer combination of the columns
 * of W (m x p, full column rank)? Solves W x = b over Q via the normal
 * equations x = (W^T W)^{-1} W^T b and checks integrality with the adjugate.
 */
static bool in_col_Z(const rmat_t* W, const rmat_t* b, u64 p)
{
  u64 m = W->rows;
  assert(b->rows == m && b->cols == 1);
  if (p == 0) return false;

  rmat_t WT = {0}, G = {0}, N = {0}, Adj = {0}, Anum = {0};
  rpol_t cp = {0};
  rz_t detG, q, r;

  assert(rmat_init(&WT, p, m) == RABIN_SUCCESS);
  for (u64 i = 0; i < p; i++) {
    for (u64 j = 0; j < m; j++) {
      assert(rmat_set(&WT, RMAT_GET(W, j, i), i, j) == RABIN_SUCCESS);
    }
  }
  assert(rmat_init(&G, p, p) == RABIN_SUCCESS);
  assert(rmat_mul(&G, &WT, W) == RABIN_SUCCESS);
  assert(rmat_init(&N, p, 1) == RABIN_SUCCESS);
  for (u64 i = 0; i < p; i++) {
    rz_t acc;
    rz_init(&acc);
    rz_set_u64(&acc, 0);
    for (u64 j = 0; j < m; j++) {
      rz_t t;
      rz_init(&t);
      rz_mul(&t, RMAT_GET(W, j, i), RMAT_GET(b, j, 0));
      rz_add(&acc, &acc, &t);
      rz_clear(&t);
    }
    rmat_set(&N, &acc, i, 0);
    rz_clear(&acc);
  }
  assert(rmat_init(&Adj, p, p) == RABIN_SUCCESS);
  assert(rmat_init(&Anum, p, 1) == RABIN_SUCCESS);
  assert(rmat_charpoly_adj(&cp, &Adj, &G) == RABIN_SUCCESS);
  assert(rmat_mul(&Anum, &Adj, &N) == RABIN_SUCCESS);
  assert(rmat_det(&detG, &G) == RABIN_SUCCESS);
  rz_init_multi(&q, &r, NULL);

  bool ok = true;
  for (u64 i = 0; i < p && ok; i++) {
    assert(rz_divmod(&q, &r, RMAT_GET(&Anum, i, 0), &detG) == RABIN_SUCCESS);
    ok = rz_is_zero(&r);
  }

  rz_clear_multi(&q, &r, NULL);
  rz_clear(&detG);
  rpol_clear(&cp);
  rmat_clear(&Anum);
  rmat_clear(&Adj);
  rmat_clear(&N);
  rmat_clear(&G);
  rmat_clear(&WT);
  return ok;
}

/* ------------------------------------------------------------------ */
/* Invariant checkers                                                 */
/* ------------------------------------------------------------------ */

/*
 * rmat_lll()/rmat_mlll() evaluate the Lovasz condition with exact integer
 * arithmetic for every double delta (delta is handled as its exact binary
 * rational value), so the exact inequality d_{k+1} d_{k-1} < delta d_k^2 -
 * lambda^2 must always hold. It is checked here at the rational target
 * dnum/dden; non-binary targets like 9/10 are only used when the double
 * delta rounds above them (e.g. 0.9 > 9/10), which makes the exact target
 * a weaker, implied condition.
 */
static bool strict_lovasz(i64 dnum, i64 dden, u64 bits)
{
  (void)dnum;
  (void)dden;
  (void)bits;
  return true;
}

/* |det M| == 1 (M must be square) */
static bool is_unimodular(const rmat_t* M)
{
  if (M->rows == 0 || M->cols == 0) return true;
  if (M->rows != M->cols) return false;
  rz_t det;
  rz_init(&det);
  assert(rmat_det(&det, M) == RABIN_SUCCESS);
  det.is_neg = false;
  bool ok = rz_is_eq_i64(&det, 1);
  rz_clear(&det);
  return ok;
}

/* Same lattice? HNF is unique per lattice, so compare HNFs. */
static bool same_lattice(const rmat_t* A, const rmat_t* B)
{
  rmat_t HA = {0}, HB = {0};
  assert(rmat_hermite(&HA, A) == RABIN_SUCCESS);
  assert(rmat_hermite(&HB, B) == RABIN_SUCCESS);
  bool ok = (HA.cols == HB.cols) && rmat_equal(&HA, &HB);
  rmat_clear(&HA);
  rmat_clear(&HB);
  return ok;
}

/*
 * Full-rank LLL invariants for A (m x n, full column rank) -> (B, H):
 *   - B == A * H
 *   - H is unimodular
 *   - lattice(B) == lattice(A)
 *   - columns of B are LLL-reduced for delta = dnum/dden
 */
static bool lll_check(const rmat_t* A, const rmat_t* B, const rmat_t* H, u64 n,
                      bool check_lovasz, i64 dnum, i64 dden)
{
  u64 m = A->rows;
  rmat_t T = {0};
  assert(rmat_init(&T, m, n) == RABIN_SUCCESS);
  assert(rmat_mul(&T, A, H) == RABIN_SUCCESS);
  bool ok = rmat_equal(&T, B);
  rmat_clear(&T);
  if (!ok) return false;
  if (!is_unimodular(H)) return false;
  if (!same_lattice(A, B)) return false;
  return check_lll_conds(B, n, check_lovasz, dnum, dden);
}

/*
 * MLLL invariants for A0 (m x n) -> (B, H), rank p:
 *   - B == A0 * H,  H unimodular
 *   - first n - p columns of B are zero
 *   - lattice(B) == lattice(A0)
 *   - last p columns of B are LLL-reduced for delta = dnum/dden
 *   - last p columns of B are a Z-basis of col_Z(A0)
 */
static bool mlll_check(const rmat_t* A0, const rmat_t* B, const rmat_t* H,
                       u64 p, bool check_lovasz, i64 dnum, i64 dden)
{
  u64 m = A0->rows;
  u64 n = A0->cols;
  if (p > n) return false;

  rmat_t T = {0};
  assert(rmat_init(&T, m, n) == RABIN_SUCCESS);
  assert(rmat_mul(&T, A0, H) == RABIN_SUCCESS);
  bool ok = rmat_equal(&T, B);
  rmat_clear(&T);
  if (!ok) return false;
  if (!is_unimodular(H)) return false;
  if (!same_lattice(A0, B)) return false;

  // first n - p columns of B are zero
  for (u64 x = 0; x < n - p; x++) {
    for (u64 y = 0; y < m; y++) {
      ok = ok && rz_is_zero(RMAT_GET(B, y, x));
    }
  }
  if (!ok) return false;

  if (p > 0) {
    rmat_t Bp = {0};
    assert(rmat_init(&Bp, m, p) == RABIN_SUCCESS);
    for (u64 x = 0; x < p; x++) {
      for (u64 y = 0; y < m; y++) {
        assert(rmat_set(&Bp, RMAT_GET(B, y, n - p + x), y, x) == RABIN_SUCCESS);
      }
    }
    ok = check_lll_conds(&Bp, p, check_lovasz, dnum, dden);
    // basis block is a Z-basis of col_Z(A0): compare against HNF(A0)
    if (ok) {
      rmat_t W = {0};
      assert(rmat_hermite(&W, A0) == RABIN_SUCCESS);
      ok = (W.cols == p);
      if (ok) {
        for (u64 x = 0; x < p && ok; x++) {
          rmat_t col = {0};
          assert(rmat_init(&col, m, 1) == RABIN_SUCCESS);
          for (u64 y = 0; y < m; y++) {
            assert(rmat_set(&col, RMAT_GET(&Bp, y, x), y, 0) == RABIN_SUCCESS);
          }
          ok = in_col_Z(&W, &col, p);
          rmat_clear(&col);
        }
      }
      rmat_clear(&W);
    }
    rmat_clear(&Bp);
  }
  return ok;
}

/* ------------------------------------------------------------------ */
/* Random matrix builders                                             */
/* ------------------------------------------------------------------ */

/*
 * Fill A (m x n, m >= n) with guaranteed full column rank n: the top n x n
 * block is lower triangular with small non-zero (signed) diagonal.
 */
static void build_full_rank(rmat_t* A, u64 m, u64 n, u64 bits, rz_t* val)
{
  for (u64 i = 0; i < m; i++) {
    for (u64 j = 0; j < n; j++) {
      if (i < j) {
        rz_set_u64(val, 0);
      } else if (i == j) {
        assert(rz_gen_random(val, 3) == RABIN_SUCCESS);
        if (rz_is_zero(val)) rz_set_u64(val, 1);
        if (rand() % 2 == 0) val->is_neg = true;
      } else {
        assert(rz_gen_random(val, bits) == RABIN_SUCCESS);
      }
      assert(rmat_set(A, val, i, j) == RABIN_SUCCESS);
    }
  }
}

/*
 * Build A0 = [V | V * C] (m x n) with exact rank p (m >= p):
 *   V is m x p lower triangular with small non-zero (signed) diagonal,
 *   C is p x (n-p) random.  All columns of A0 lie in col(V); the first p
 *   columns are V (rank p), so rank(A0) = p.
 */
static void build_rank_p(rmat_t* A0, u64 m, u64 n, u64 p, u64 bits, rz_t* val)
{
  rmat_t V = {0}, C = {0}, T = {0};
  assert(rmat_init(&V, m, p) == RABIN_SUCCESS);
  assert(rmat_init(&C, p, n - p) == RABIN_SUCCESS);
  for (u64 i = 0; i < m; i++) {
    for (u64 j = 0; j < p; j++) {
      if (i < j) {
        rz_set_u64(val, 0);
      } else if (i == j) {
        assert(rz_gen_random(val, 3) == RABIN_SUCCESS);
        if (rz_is_zero(val)) rz_set_u64(val, 1);
        if (rand() % 2 == 0) val->is_neg = true;
      } else {
        assert(rz_gen_random(val, bits) == RABIN_SUCCESS);
      }
      assert(rmat_set(&V, val, i, j) == RABIN_SUCCESS);
    }
  }
  for (u64 i = 0; i < p; i++) {
    for (u64 j = 0; j < n - p; j++) {
      assert(rz_gen_random(val, bits) == RABIN_SUCCESS);
      assert(rmat_set(&C, val, i, j) == RABIN_SUCCESS);
    }
  }
  // A0[:, 0..p-1] = V,  A0[:, p..n-1] = V * C
  for (u64 i = 0; i < m; i++) {
    for (u64 j = 0; j < p; j++) {
      assert(rmat_set(A0, RMAT_GET(&V, i, j), i, j) == RABIN_SUCCESS);
    }
  }
  assert(rmat_init(&T, m, n - p) == RABIN_SUCCESS);
  assert(rmat_mul(&T, &V, &C) == RABIN_SUCCESS);
  for (u64 i = 0; i < m; i++) {
    for (u64 j = 0; j < n - p; j++) {
      assert(rmat_set(A0, RMAT_GET(&T, i, j), i, p + j) == RABIN_SUCCESS);
    }
  }
  rmat_clear(&V);
  rmat_clear(&C);
  rmat_clear(&T);
}

/* ------------------------------------------------------------------ */
/* Tests                                                              */
/* ------------------------------------------------------------------ */

// Deterministic known-answer cases (exact reduced basis and transform).
static void test_lll_known(void)
{
  // 1. identity 2x2 is already reduced: B = I, H = I
  {
    rmat_t A = {0}, B = {0}, H = {0}, I = {0};
    rmat_init(&A, 2, 2);
    set164(&A, 0, 0, 1);
    set164(&A, 0, 1, 0);
    set164(&A, 1, 0, 0);
    set164(&A, 1, 1, 1);
    assert(rmat_init(&B, 2, 2) == RABIN_SUCCESS);
    assert(rmat_copy(&B, &A) == RABIN_SUCCESS);
    assert(rmat_lll(&B, 2, 0.75, &H) == RABIN_SUCCESS);
    assert(rmat_init(&I, 2, 2) == RABIN_SUCCESS);
    assert(rmat_id(&I, 2) == RABIN_SUCCESS);
    assert(rmat_equal(&B, &I) && rmat_equal(&H, &I));
    rmat_clear(&A);
    rmat_clear(&B);
    rmat_clear(&H);
    rmat_clear(&I);
  }

  // 2. columns [1,0] and [100,1] -> reduced to [1,0] and [0,1],
  //    H = [[1,-100],[0,1]]
  {
    rmat_t A = {0}, B = {0}, H = {0};
    rmat_init(&A, 2, 2);
    set164(&A, 0, 0, 1);
    set164(&A, 0, 1, 100);
    set164(&A, 1, 0, 0);
    set164(&A, 1, 1, 1);
    assert(rmat_init(&B, 2, 2) == RABIN_SUCCESS);
    assert(rmat_copy(&B, &A) == RABIN_SUCCESS);
    assert(rmat_lll(&B, 2, 0.75, &H) == RABIN_SUCCESS);
    assert(rz_is_eq_i64(RMAT_GET(&B, 0, 0), 1));
    assert(rz_is_eq_i64(RMAT_GET(&B, 0, 1), 0));
    assert(rz_is_eq_i64(RMAT_GET(&B, 1, 0), 0));
    assert(rz_is_eq_i64(RMAT_GET(&B, 1, 1), 1));
    assert(rz_is_eq_i64(RMAT_GET(&H, 0, 0), 1));
    assert(rz_is_eq_i64(RMAT_GET(&H, 0, 1), -100));
    assert(rz_is_eq_i64(RMAT_GET(&H, 1, 0), 0));
    assert(rz_is_eq_i64(RMAT_GET(&H, 1, 1), 1));
    rmat_clear(&A);
    rmat_clear(&B);
    rmat_clear(&H);
  }

  // 3. diagonal [[2,0],[0,4]] is already reduced (mu = 0, Lovasz holds)
  {
    rmat_t A = {0}, B = {0}, H = {0}, I = {0};
    rmat_init(&A, 2, 2);
    set164(&A, 0, 0, 2);
    set164(&A, 0, 1, 0);
    set164(&A, 1, 0, 0);
    set164(&A, 1, 1, 4);
    assert(rmat_init(&B, 2, 2) == RABIN_SUCCESS);
    assert(rmat_copy(&B, &A) == RABIN_SUCCESS);
    assert(rmat_lll(&B, 2, 0.75, &H) == RABIN_SUCCESS);
    assert(rmat_init(&I, 2, 2) == RABIN_SUCCESS);
    assert(rmat_id(&I, 2) == RABIN_SUCCESS);
    assert(rmat_equal(&B, &A) && rmat_equal(&H, &I));
    rmat_clear(&A);
    rmat_clear(&B);
    rmat_clear(&H);
    rmat_clear(&I);
  }
  printf("[PASS] lll known answers\n");
}

// Property-based: random full-rank matrices across shapes, bits and delta.
static void test_lll_property(void)
{
  struct {
    u64 m, n;
  } shapes[] = {{2, 2}, {3, 3}, {4, 3}, {5, 3}, {6, 4}, {8, 5}, {4, 4}};
  u64 nshapes = sizeof(shapes) / sizeof(shapes[0]);
  struct {
    double delta;
    i64 dnum, dden;
    // 0.9 as a double is slightly above 9/10, so checking the exact 9/10
    // Lovasz target is a weaker but implied condition; 0.875 = 7/8 is exact
  } deltas[] = {
      {0.75, 3, 4}, {0.5, 1, 2}, {1.0, 1, 1}, {0.9, 9, 10}, {0.875, 7, 8},
  };
  u64 ndeltas = sizeof(deltas) / sizeof(deltas[0]);
  u64 bits_list[] = {8, 16, 64, 128};
  u64 nbits = sizeof(bits_list) / sizeof(bits_list[0]);

  rz_t val;
  rz_init(&val);
  for (u64 s = 0; s < nshapes; s++) {
    for (u64 b = 0; b < nbits; b++) {
      for (u64 d = 0; d < ndeltas; d++) {
        u64 m = shapes[s].m, n = shapes[s].n;
        for (int trial = 0; trial < 2; trial++) {
          rmat_t A = {0}, B = {0}, H = {0};
          assert(rmat_init(&A, m, n) == RABIN_SUCCESS);
          assert(rmat_init(&B, m, n) == RABIN_SUCCESS);
          build_full_rank(&A, m, n, bits_list[b], &val);
          assert(rmat_copy(&B, &A) == RABIN_SUCCESS);
          assert(rmat_lll(&B, n, deltas[d].delta, &H) == RABIN_SUCCESS);
          assert(lll_check(
              &A, &B, &H, n,
              strict_lovasz(deltas[d].dnum, deltas[d].dden, bits_list[b]),
              deltas[d].dnum, deltas[d].dden));
          rmat_clear(&A);
          rmat_clear(&B);
          rmat_clear(&H);
        }
      }
    }
  }
  rz_clear(&val);
  printf("[PASS] lll property (%llu shapes x %llu bits x %llu delta)\n",
         (unsigned long long)nshapes, (unsigned long long)nbits,
         (unsigned long long)ndeltas);
}

// Large (multi-limb) entries with a moderately wide basis.
static void test_lll_large(void)
{
  u64 bits_list[] = {512, 1024};
  u64 nbits = sizeof(bits_list) / sizeof(bits_list[0]);
  struct {
    u64 m, n;
  } shapes[] = {{4, 3}, {6, 4}, {5, 4}};
  u64 nshapes = sizeof(shapes) / sizeof(shapes[0]);
  rz_t val;
  rz_init(&val);
  for (u64 s = 0; s < nshapes; s++) {
    for (u64 b = 0; b < nbits; b++) {
      u64 m = shapes[s].m, n = shapes[s].n;
      rmat_t A = {0}, B = {0}, H = {0};
      assert(rmat_init(&A, m, n) == RABIN_SUCCESS);
      assert(rmat_init(&B, m, n) == RABIN_SUCCESS);
      build_full_rank(&A, m, n, bits_list[b], &val);
      assert(rmat_copy(&B, &A) == RABIN_SUCCESS);
      assert(rmat_lll(&B, n, 0.75, &H) == RABIN_SUCCESS);
      assert(lll_check(&A, &B, &H, n, true, 3, 4));
      rmat_clear(&A);
      rmat_clear(&B);
      rmat_clear(&H);
    }
  }
  rz_clear(&val);
  printf("[PASS] lll large multi-limb entries\n");
}

// Negative entries (full-rank lll) and zero / duplicate columns (mlll).
static void test_lll_edge(void)
{
  // 1. full-rank with many negative entries
  {
    rmat_t A = {0}, B = {0}, H = {0};
    rz_t val;
    rz_init(&val);
    assert(rmat_init(&A, 5, 4) == RABIN_SUCCESS);
    assert(rmat_init(&B, 5, 4) == RABIN_SUCCESS);
    build_full_rank(&A, 5, 4, 32, &val);
    // force the top-left 2x2 block to be negative to stress signs
    set164(&A, 0, 0, -3);
    set164(&A, 1, 0, -11);
    set164(&A, 1, 1, -7);
    assert(rmat_copy(&B, &A) == RABIN_SUCCESS);
    assert(rmat_lll(&B, 4, 0.75, &H) == RABIN_SUCCESS);
    assert(lll_check(&A, &B, &H, 4, true, 3, 4));
    rmat_clear(&A);
    rmat_clear(&B);
    rmat_clear(&H);
    rz_clear(&val);
  }

  // 2. duplicate columns: [1,2,1,2] -> rank 2 (mlll)
  {
    rmat_t A0 = {0}, B = {0}, H = {0};
    rmat_init(&A0, 1, 4);
    set164(&A0, 0, 0, 1);
    set164(&A0, 0, 1, 2);
    set164(&A0, 0, 2, 1);
    set164(&A0, 0, 3, 2);
    assert(rmat_init(&B, 1, 4) == RABIN_SUCCESS);
    assert(rmat_copy(&B, &A0) == RABIN_SUCCESS);
    u64 rank = 99;
    assert(rmat_mlll(&B, 4, 0.75, &H, &rank) == RABIN_SUCCESS);
    assert(rank == 1);
    assert(mlll_check(&A0, &B, &H, 1, true, 3, 4));
    rmat_clear(&A0);
    rmat_clear(&B);
    rmat_clear(&H);
  }

  // 3. zero columns mixed with independent ones (mlll): rank 2
  {
    rmat_t A0 = {0}, B = {0}, H = {0};
    rmat_init(&A0, 2, 4);
    set164(&A0, 0, 0, 3);
    set164(&A0, 1, 0, 1);
    set164(&A0, 0, 1, 0);  // zero column
    set164(&A0, 1, 1, 0);
    set164(&A0, 0, 2, 5);
    set164(&A0, 1, 2, 2);
    set164(&A0, 0, 3, 0);  // zero column
    set164(&A0, 1, 3, 0);
    assert(rmat_init(&B, 2, 4) == RABIN_SUCCESS);
    assert(rmat_copy(&B, &A0) == RABIN_SUCCESS);
    u64 rank = 99;
    assert(rmat_mlll(&B, 4, 0.75, &H, &rank) == RABIN_SUCCESS);
    assert(rank == 2);
    assert(mlll_check(&A0, &B, &H, 2, true, 3, 4));
    rmat_clear(&A0);
    rmat_clear(&B);
    rmat_clear(&H);
  }

  // 4. all-zero matrix -> rank 0, H = I
  {
    rmat_t B = {0}, H = {0}, I4 = {0};
    assert(rmat_init(&B, 2, 4) == RABIN_SUCCESS);
    u64 rank = 99;
    assert(rmat_mlll(&B, 4, 0.75, &H, &rank) == RABIN_SUCCESS);
    assert(rank == 0);
    assert(rmat_init(&I4, 4, 4) == RABIN_SUCCESS);
    assert(rmat_id(&I4, 4) == RABIN_SUCCESS);
    assert(rmat_equal(&H, &I4));
    rmat_clear(&B);
    rmat_clear(&H);
    rmat_clear(&I4);
  }
  printf("[PASS] lll/mlll edge cases\n");
}

// Property-based: random rank-deficient matrices across shapes and delta.
static void test_mlll_property(void)
{
  struct {
    u64 m, n, p;
  } shapes[] = {{2, 2, 2}, {2, 3, 1}, {2, 3, 2}, {3, 4, 2}, {3, 5, 3},
                {4, 6, 2}, {4, 6, 4}, {5, 6, 3}, {3, 3, 1}, {6, 6, 4}};
  u64 nshapes = sizeof(shapes) / sizeof(shapes[0]);
  struct {
    double delta;
    i64 dnum, dden;
  } deltas[] = {
      {0.75, 3, 4},
      {0.5, 1, 2},
      {0.9, 9, 10},
      {0.875, 7, 8},
  };
  u64 ndeltas = sizeof(deltas) / sizeof(deltas[0]);
  u64 bits_list[] = {8, 16, 32};
  u64 nbits = sizeof(bits_list) / sizeof(bits_list[0]);

  rz_t val;
  rz_init(&val);
  for (u64 s = 0; s < nshapes; s++) {
    u64 m = shapes[s].m, n = shapes[s].n, p = shapes[s].p;
    for (u64 b = 0; b < nbits; b++) {
      for (u64 d = 0; d < ndeltas; d++) {
        for (int trial = 0; trial < 2; trial++) {
          rmat_t A0 = {0}, B = {0}, H = {0};
          assert(rmat_init(&A0, m, n) == RABIN_SUCCESS);
          assert(rmat_init(&B, m, n) == RABIN_SUCCESS);
          build_rank_p(&A0, m, n, p, bits_list[b], &val);
          assert(rmat_copy(&B, &A0) == RABIN_SUCCESS);
          u64 rank = 99;
          assert(rmat_mlll(&B, n, deltas[d].delta, &H, &rank) == RABIN_SUCCESS);
          assert(rank == p);
          assert(mlll_check(
              &A0, &B, &H, p,
              strict_lovasz(deltas[d].dnum, deltas[d].dden, bits_list[b]),
              deltas[d].dnum, deltas[d].dden));
          rmat_clear(&A0);
          rmat_clear(&B);
          rmat_clear(&H);
        }
      }
    }
  }
  rz_clear(&val);
  printf("[PASS] mlll property (%llu shapes x %llu bits x %llu delta)\n",
         (unsigned long long)nshapes, (unsigned long long)nbits,
         (unsigned long long)ndeltas);
}

// Idempotency: reducing an already-reduced basis keeps it a valid reduced
// basis of the same lattice (and the composed transform stays unimodular).
static void test_lll_idempotent(void)
{
  struct {
    u64 m, n;
  } shapes[] = {{4, 3}, {5, 4}, {6, 5}};
  u64 nshapes = sizeof(shapes) / sizeof(shapes[0]);
  rz_t val;
  rz_init(&val);
  for (u64 s = 0; s < nshapes; s++) {
    u64 m = shapes[s].m, n = shapes[s].n;
    rmat_t A = {0}, B1 = {0}, H1 = {0}, B2 = {0}, H2 = {0}, T = {0}, C = {0};
    assert(rmat_init(&A, m, n) == RABIN_SUCCESS);
    assert(rmat_init(&B1, m, n) == RABIN_SUCCESS);
    build_full_rank(&A, m, n, 32, &val);
    assert(rmat_copy(&B1, &A) == RABIN_SUCCESS);
    assert(rmat_lll(&B1, n, 0.75, &H1) == RABIN_SUCCESS);
    // second reduction on the result
    assert(rmat_init(&B2, m, n) == RABIN_SUCCESS);
    assert(rmat_copy(&B2, &B1) == RABIN_SUCCESS);
    assert(rmat_lll(&B2, n, 0.75, &H2) == RABIN_SUCCESS);
    assert(lll_check(&B1, &B2, &H2, n, true, 3, 4));
    // composed transform H1 * H2 must also map A to B2 and be unimodular
    assert(rmat_init(&C, n, n) == RABIN_SUCCESS);
    assert(rmat_mul(&C, &H1, &H2) == RABIN_SUCCESS);
    assert(is_unimodular(&C));
    assert(rmat_init(&T, m, n) == RABIN_SUCCESS);
    assert(rmat_mul(&T, &A, &C) == RABIN_SUCCESS);
    assert(rmat_equal(&T, &B2));
    rmat_clear(&T);
    rmat_clear(&C);
    rmat_clear(&A);
    rmat_clear(&B1);
    rmat_clear(&B2);
    rmat_clear(&H1);
    rmat_clear(&H2);
  }
  rz_clear(&val);
  printf("[PASS] lll idempotent\n");
}

// Wall-clock timing of rmat_lll on full-rank square matrices of increasing
// size (fixed entry size), to show how the reduction cost scales with the
// number of basis vectors. Mirrors the run_*_benchmark() style in
// tests/test_matrix.c. The default range stops at 64: the cost grows
// super-cubically, so 128 x 128 already takes ~40s and would make this test
// slow. Bump the last entry to explore larger sizes if needed.
static void test_lll_timing(void)
{
  u64 sizes[] = {4, 8, 16, 32, 64, 128};
  u64 nsizes = sizeof(sizes) / sizeof(sizes[0]);
  u64 bits = 32;
  rz_t val;
  rz_init(&val);
  for (u64 s = 0; s < nsizes; s++) {
    u64 n = sizes[s];
    rmat_t A = {0}, B = {0}, H = {0};
    assert(rmat_init(&A, n, n) == RABIN_SUCCESS);
    assert(rmat_init(&B, n, n) == RABIN_SUCCESS);
    build_full_rank(&A, n, n, bits, &val);
    assert(rmat_copy(&B, &A) == RABIN_SUCCESS);
    printf("Benchmarking %zu x %zu (%zu-bit entries) LLL... ", (size_t)n,
           (size_t)n, (size_t)bits);
    fflush(stdout);
    double start = get_time();
    assert(rmat_lll(&B, n, 0.75, &H) == RABIN_SUCCESS);
    double end = get_time();
    printf("Time: %f seconds\n", end - start);
    rmat_clear(&A);
    rmat_clear(&B);
    rmat_clear(&H);
  }
  rz_clear(&val);
}

void test_lll_3x3()
{
  rmat_t A = {0}, B = {0}, H = {0};
  rmat_init(&A, 3, 3);
  rmat_init(&B, 3, 3);
  rmat_init(&H, 3, 3);
  set164(&A, 0, 0, 2);
  set164(&A, 1, 0, 3);
  set164(&A, 2, 0, 1);
  set164(&A, 0, 1, 1);
  set164(&A, 1, 1, 2);
  set164(&A, 2, 1, -1);
  set164(&A, 0, 2, -2);
  set164(&A, 1, 2, -2);
  set164(&A, 2, 2, 2);

  rmat_print(&A);
  rmat_copy(&B, &A);
  rmat_lll(&B, 3, 0.75, &H);

  rmat_print(&B);

  rmat_print(&H);

  rmat_clear(&A);
  rmat_clear(&B);
  rmat_clear(&H);
}

// A failed reduction must leave the caller's matrix untouched.
static void test_lll_failure_untouched(void)
{
  // dependent columns -> rmat_lll fails with B unchanged
  {
    rmat_t A = {0}, B = {0}, H = {0};
    assert(rmat_init(&A, 2, 3) == RABIN_SUCCESS);
    set164(&A, 0, 0, 1);
    set164(&A, 1, 0, 0);
    set164(&A, 0, 1, 0);
    set164(&A, 1, 1, 1);
    set164(&A, 0, 2, 1);  // col 2 = col 0 + col 1
    set164(&A, 1, 2, 1);
    assert(rmat_init(&B, 2, 3) == RABIN_SUCCESS);
    assert(rmat_copy(&B, &A) == RABIN_SUCCESS);
    assert(rmat_lll(&B, 3, 0.75, &H) == RABIN_ERR_INVALID_ARG);
    assert(rmat_equal(&B, &A));
    rmat_clear(&A);
    rmat_clear(&B);
    rmat_clear(&H);
  }

  // wide matrix (B->rows > n) -> rmat_mlll fails early, B unchanged
  {
    rmat_t A = {0}, B = {0}, H = {0};
    assert(rmat_init(&A, 3, 2) == RABIN_SUCCESS);
    set164(&A, 0, 0, 1);
    set164(&A, 1, 0, 0);
    set164(&A, 2, 0, 5);
    set164(&A, 0, 1, 0);
    set164(&A, 1, 1, 1);
    set164(&A, 2, 1, 7);
    assert(rmat_init(&B, 3, 2) == RABIN_SUCCESS);
    assert(rmat_copy(&B, &A) == RABIN_SUCCESS);
    u64 rank = 99;
    assert(rmat_mlll(&B, 2, 0.75, &H, &rank) == RABIN_ERR_INVALID_ARG);
    assert(rmat_equal(&B, &A));
    rmat_clear(&A);
    rmat_clear(&B);
    rmat_clear(&H);
  }
  printf("[PASS] failed reduction leaves B untouched\n");
}

int main(void)
{
  // srand(0x5eedULL);
  test_lll_known();
  test_lll_property();
  test_lll_large();
  test_lll_edge();
  test_lll_failure_untouched();
  test_mlll_property();
  test_lll_idempotent();
  test_lll_timing();
  test_lll_3x3();

  printf("All extended LLL tests passed\n");
  return 0;
}
