#include <stdio.h>
#include <stdlib.h>

#include "../include/rmat.h"
#include "../include/rmatq.h"
#include "../include/rq.h"
#include "../include/rvec.h"

static int failures = 0;

#define CHECK(cond)                                          \
  do {                                                       \
    if (!(cond)) {                                           \
      printf("FAIL %s:%d: %s\n", __FILE__, __LINE__, #cond); \
      failures++;                                            \
    }                                                        \
  } while (0)

#define MAXN 8

static void set_int(rmatq_t* A, u64 r, u64 c, i64 v)
{
  rz_t z;
  rz_init(&z);
  rz_set_i64(&z, v);
  rmat_set(&A->M, &z, r, c);
  rz_clear(&z);
}

static void set_rat(rmatq_t* A, u64 r, u64 c, i64 num, i64 den)
{
  rq_t a;
  rq_init(&a);
  rz_set_i64(&a.num, num);
  rz_set_i64(&a.den, den);
  rq_normalize(&a);
  rmatq_set(A, &a, r, c);
  rq_clear(&a);
}

static void set_i(rq_t* a, i64 num, i64 den)
{
  rz_set_i64(&a->num, num);
  rz_set_i64(&a->den, den);
  rq_normalize(a);
}

static int rq_eq(const rq_t* a, i64 num, i64 den)
{
  rq_t w;
  rq_init(&w);
  set_i(&w, num, den);
  int r = (rq_cmp(a, &w) == 0);
  rq_clear(&w);
  return r;
}

static void mat_from_ints(rmatq_t* A, const i64* v, u64 n)
{
  for (u64 i = 0; i < n; i++) {
    for (u64 j = 0; j < n; j++) {
      set_int(A, i, j, v[i * n + j]);
    }
  }
}

static void print_vec(const rvec_t* x, u64 n)
{
  printf("[");
  for (u64 i = 0; i < n; i++) {
    if (i) printf(", ");
    rz_print(&x->data[i]);
  }
  printf("]");
}

static void print_rat_mat(const rmatq_t* A, u64 n)
{
  printf("[");
  for (u64 i = 0; i < n; i++) {
    for (u64 j = 0; j < n; j++) {
      printf(j ? ", " : " ");
      rq_t a;
      rq_init(&a);
      rmatq_get(&a, A, i, j);
      rz_print(&a.num);
      if (!rz_is_eq_i64(&a.den, 1)) {
        printf("/");
        rz_print(&a.den);
      }
      rq_clear(&a);
    }
    printf(i + 1 < n ? ";" : " ]");
  }
}

// quadratic form value Q(x) = x^T A x as a rational
static void quad_form(rq_t* out, const rmatq_t* A, const rz_t* x, u64 n)
{
  rq_set_u64(out, 0);
  for (u64 i = 0; i < n; i++) {
    if (rz_is_zero(&x[i])) continue;
    for (u64 j = 0; j < n; j++) {
      if (rz_is_zero(&x[j])) continue;
      rq_t a, t, xj;
      rq_init_multi(&a, &t, &xj, NULL);
      rmatq_get(&a, A, i, j);
      rq_set(&xj, &x[j]);
      rq_mul(&t, &a, &xj);
      rq_set(&xj, &x[i]);
      rq_mul(&t, &t, &xj);
      rq_add(out, out, &t);
      rq_clear_multi(&a, &t, &xj, NULL);
    }
  }
  rq_normalize(out);
}

// the reported norm must equal Q(x) and x must be non-zero
static void check_consistent(const rmatq_t* A, const rvec_t* x,
                             const rq_t* norm, u64 n)
{
  rq_t q, pos;
  rq_init_multi(&q, &pos, NULL);
  quad_form(&q, A, x->data, n);
  if (rq_cmp(&q, norm) != 0) {
    printf("FAIL: Q(x) = ");
    rq_print(&q);
    printf(" != reported norm ");
    rq_print(norm);
    printf("\n");
    failures++;
  }
  int nz = 0;
  for (u64 i = 0; i < n; i++) nz |= !rz_is_zero(&x->data[i]);
  if (!nz) {
    printf("FAIL: zero vector returned\n");
    failures++;
  }
  set_i(&pos, 0, 1);
  CHECK(rq_cmp(norm, &pos) > 0);
  rq_clear_multi(&q, &pos, NULL);
}

// run Fincke-Pohst on an integer matrix and compare against num/den
static void check_expected(const char* name, const i64* v, u64 n, i64 num,
                           i64 den)
{
  rmatq_t A = {0};
  rvec_t x = {0};
  rq_t norm = {0};
  rmatq_init(&A, n, n);
  mat_from_ints(&A, v, n);
  rvec_init(&x, n);
  rq_init(&norm);

  printf("== %s ==\n", name);
  printf("A   : n=%llu, ", (unsigned long long)n);
  print_rat_mat(&A, n);
  printf("\n");

  rabin_err_t err = rmatq_fincke_pohst(&x, &norm, &A);
  if (err != RABIN_SUCCESS) {
    printf("FAIL: err=%d\n", (int)err);
    failures++;
  } else {
    rq_t want;
    rq_init_multi(&want, NULL);
    set_i(&want, num, den);
    printf("x   : ");
    print_vec(&x, n);
    printf("\nnorm: ");
    rq_print(&norm);
    printf("want: ");
    rq_print(&want);
    if (rq_cmp(&norm, &want) != 0) {
      printf("FAIL: norm mismatch\n");
      failures++;
    }
    check_consistent(&A, &x, &norm, n);
    rq_clear(&want);
  }
  rvec_clear(&x);
  rq_clear(&norm);
  rmatq_clear(&A);
}

static void brute_rec(rq_t* best, int* found, const rmatq_t* A, u64 n, u64 i,
                      i64 B, rz_t* x)
{
  if (i == n) {
    int nz = 0;
    for (u64 k = 0; k < n; k++) nz |= !rz_is_zero(&x[k]);
    if (!nz) return;
    rq_t q;
    rq_init(&q);
    quad_form(&q, A, x, n);
    if (!*found || rq_cmp(&q, best) < 0) {
      rq_copy(best, &q);
      *found = 1;
    }
    rq_clear(&q);
    return;
  }
  for (i64 t = -B; t <= B; t++) {
    rz_set_i64(&x[i], t);
    brute_rec(best, found, A, n, i + 1, B, x);
  }
}

// compare Fincke-Pohst against brute force over the box [-B, B]^n
static void check_brute(const rmatq_t* A, u64 n, i64 B)
{
  rvec_t x = {0};
  rq_t norm = {0};
  rvec_init(&x, n);
  rq_init(&norm);

  rabin_err_t err = rmatq_fincke_pohst(&x, &norm, A);
  if (err != RABIN_SUCCESS) {
    printf("FAIL: brute err=%d\n", (int)err);
    failures++;
    goto cleanup;
  }

  printf("A     : n=%llu, ", (unsigned long long)n);
  print_rat_mat(A, n);
  printf("\nx     : ");
  rvec_print(&x);
  printf("\nfp    : ");
  rq_print(&norm);

  rz_t xs[MAXN];
  for (u64 i = 0; i < n; i++) {
    rz_init(&xs[i]);
    rz_set_u64(&xs[i], 0);
  }
  rq_t best;
  rq_init(&best);
  int found = 0;
  brute_rec(&best, &found, A, n, 0, B, xs);
  printf("brute : ");
  rq_print(&best);
  printf("\n");
  CHECK(found);
  if (rq_cmp(&norm, &best) != 0) {
    printf("FAIL: fp norm = ");
    rq_print(&norm);
    printf(", brute = ");
    rq_print(&best);
    printf("\n");
    failures++;
  }
  check_consistent(A, &x, &norm, n);

  for (u64 i = 0; i < n; i++) rz_clear(&xs[i]);
  rq_clear(&best);
cleanup:
  rvec_clear(&x);
  rq_clear(&norm);
}

static void test_known(void)
{
  const i64 A2[4] = {2, -1, -1, 2};
  check_expected("A2 = [[2,-1],[-1,2]], lambda^2 = 2", A2, 2, 2, 1);

  const i64 I3[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
  check_expected("I3 = identity, lambda^2 = 1", I3, 3, 1, 1);

  const i64 one[1] = {5};
  check_expected("[5], lambda^2 = 5", one, 1, 5, 1);

  const i64 D4[16] = {2, -1, 0, 0, -1, 2, -1, 0, 0, -1, 2, -1, 0, 0, -1, 2};
  check_expected("D4 Cartan, lambda^2 = 2", D4, 4, 2, 1);

  // clang-format off
  const i64 E8[64] = {
     2, -1,  0,  0,  0,  0,  0,  0,
    -1,  2, -1,  0,  0,  0,  0,  0,
     0, -1,  2, -1,  0,  0,  0, -1,
     0,  0, -1,  2, -1,  0,  0,  0,
     0,  0,  0, -1,  2, -1,  0,  0,
     0,  0,  0,  0, -1,  2, -1,  0,
     0,  0,  0,  0,  0, -1,  2,  0,
     0,  0, -1,  0,  0,  0,  0,  2,
  };
  // clang-format on
  check_expected("E8 Cartan, lambda^2 = 2", E8, 8, 2, 1);

  // rational Gram matrix: (1/2) * [[3, 1], [1, 3]], lambda = 3/2
  {
    rmatq_t A = {0};
    rvec_t x = {0};
    rq_t norm = {0};
    rmatq_init(&A, 2, 2);
    set_rat(&A, 0, 0, 3, 2);
    set_rat(&A, 1, 1, 3, 2);
    set_rat(&A, 0, 1, 1, 2);
    set_rat(&A, 1, 0, 1, 2);
    rvec_init(&x, 2);
    rq_init(&norm);
    printf("== rational A = (1/2)[[3,1],[1,3]], lambda^2 = 3/2 ==\n");
    printf("A   : n=2, ");
    print_rat_mat(&A, 2);
    printf("\n");
    CHECK(rmatq_fincke_pohst(&x, &norm, &A) == RABIN_SUCCESS);
    rq_t want;
    rq_init_multi(&want, NULL);
    set_i(&want, 3, 2);
    printf("x   : ");
    print_vec(&x, 2);
    printf("\nnorm: ");
    rq_print(&norm);
    printf("want: ");
    rq_print(&want);
    CHECK(rq_cmp(&norm, &want) == 0);
    check_consistent(&A, &x, &norm, 2);
    rq_clear_multi(&want, NULL);
    rvec_clear(&x);
    rq_clear(&norm);
    rmatq_clear(&A);
  }
}

// random positive definite A = M^T M + I with small entries,
// cross-checked against brute force over a guaranteed-safe box
static void test_random(void)
{
  srand(42);
  for (int rep = 0; rep < 25; rep++) {
    u64 n = 2 + (u64)(rand() % 3);  // 2..4
    i64 m[4][4], a[4 * 4];
    for (u64 i = 0; i < 4; i++) {
      for (u64 j = 0; j < 4; j++) m[i][j] = (rand() % 10) - 2;
    }
    for (u64 i = 0; i < n; i++) {
      for (u64 j = 0; j < n; j++) {
        i64 s = (i == j) ? 1 : 0;
        for (u64 k = 0; k < n; k++) s += m[k][i] * m[k][j];
        a[i * n + j] = s;
      }
    }
    rmatq_t A = {0};
    rmatq_init(&A, n, n);
    mat_from_ints(&A, a, n);
    printf("== random rep %d (brute box |x| <= 5) ==\n", rep);
    check_brute(&A, n, 5);
    rmatq_clear(&A);
  }
}

// rational form (1/d) * M, integer symmetric M; run Fincke-Pohst and
// compare the norm against num/den
static void check_rat_expected(const char* name, const i64* m, i64 d, i64 num,
                               i64 den)
{
  rmatq_t A = {0};
  rvec_t x = {0};
  rq_t norm = {0};
  rmatq_init(&A, 2, 2);
  for (u64 i = 0; i < 2; i++) {
    for (u64 j = 0; j < 2; j++) set_rat(&A, i, j, m[i * 2 + j], d);
  }
  rvec_init(&x, 2);
  rq_init(&norm);

  printf("== %s ==\n", name);
  printf("A   : n=2, ");
  print_rat_mat(&A, 2);
  printf("\n");

  rabin_err_t err = rmatq_fincke_pohst(&x, &norm, &A);
  if (err != RABIN_SUCCESS) {
    printf("FAIL: err=%d\n", (int)err);
    failures++;
  } else {
    rq_t want;
    rq_init_multi(&want, NULL);
    set_i(&want, num, den);
    printf("x   : ");
    print_vec(&x, 2);
    printf("\nnorm: ");
    rq_print(&norm);
    printf("want: ");
    rq_print(&want);
    if (rq_cmp(&norm, &want) != 0) {
      printf("FAIL: norm mismatch\n");
      failures++;
    }
    check_consistent(&A, &x, &norm, 2);
    rq_clear(&want);
  }
  rvec_clear(&x);
  rq_clear(&norm);
  rmatq_clear(&A);
}

// Regression: the enumeration bounds floor(sqrt(T/q_ii) -/+ shift) need
// exact rounding of sqrt of non-square rationals plus a fractional shift.
// With the old (cand+c)^2 <= r accept test the shift floor was wrong and
// several of these forms returned a norm strictly above the true minimum.
static void test_rounding(void)
{
  // clang-format off
  const i64 R0[4] = {24,  8,  8, 25}; check_rat_expected("(1/17)[[24,8],[8,25]]", R0, 17, 24, 17);
  const i64 R1[4] = {11, 16, 16, 25}; check_rat_expected("(1/2)[[11,16],[16,25]]", R1, 2, 2, 1);
  const i64 R2[4] = {26,-16,-16, 11}; check_rat_expected("(1/6)[[26,-16],[-16,11]]", R2, 6, 5, 6);
  const i64 R3[4] = {19,  6,  6, 20}; check_rat_expected("(1/16)[[19,6],[6,20]]", R3, 16, 19, 16);
  const i64 R4[4] = {16,-20,-20, 26}; check_rat_expected("(1/11)[[16,-20],[-20,26]]", R4, 11, 2, 11);
  const i64 R5[4] = {11,  3,  3,  1}; check_rat_expected("(1/7)[[11,3],[3,1]]", R5, 7, 1, 7);
  const i64 R6[4] = {22, 16, 16, 12}; check_rat_expected("(1/4)[[22,16],[16,12]]", R6, 4, 1, 2);
  const i64 R7[4] = {18,  6,  6, 19}; check_rat_expected("(1/4)[[18,6],[6,19]]", R7, 4, 9, 2);
  // clang-format on
}

static void test_errors(void)
{
  const i64 A2[4] = {2, -1, -1, 2};
  rmatq_t A = {0};
  rvec_t x = {0};
  rq_t norm = {0};
  rmatq_init(&A, 2, 2);
  mat_from_ints(&A, A2, 2);
  rvec_init(&x, 2);
  rq_init(&norm);

  CHECK(rmatq_fincke_pohst(NULL, &norm, &A) == RABIN_ERR_NULL_PTR);
  CHECK(rmatq_fincke_pohst(&x, NULL, &A) == RABIN_ERR_NULL_PTR);
  CHECK(rmatq_fincke_pohst(&x, &norm, NULL) == RABIN_ERR_NULL_PTR);
  CHECK(rmatq_fincke_pohst(NULL, NULL, NULL) == RABIN_ERR_NULL_PTR);

  rmatq_t sq = {0};
  rmatq_init(&sq, 2, 3);
  CHECK(rmatq_fincke_pohst(&x, &norm, &sq) == RABIN_ERR_MATRIX_DIM);
  rmatq_clear(&sq);

  rvec_t bad = {0};
  rvec_init(&bad, 3);
  CHECK(rmatq_fincke_pohst(&bad, &norm, &A) == RABIN_ERR_MATRIX_DIM);
  rvec_clear(&bad);

  // singular input is rejected
  rmatq_t sing = {0};
  rmatq_init(&sing, 2, 2);
  set_int(&sing, 0, 0, 1);
  set_int(&sing, 0, 1, 1);
  set_int(&sing, 1, 0, 1);
  set_int(&sing, 1, 1, 1);
  CHECK(rmatq_fincke_pohst(&x, &norm, &sing) != RABIN_SUCCESS);
  rmatq_clear(&sing);

  // not positive definite: [[1, 2], [2, 1]] has one negative eigenvalue
  rmatq_t npd = {0};
  rmatq_init(&npd, 2, 2);
  set_int(&npd, 0, 0, 1);
  set_int(&npd, 0, 1, 2);
  set_int(&npd, 1, 0, 2);
  set_int(&npd, 1, 1, 1);
  CHECK(rmatq_fincke_pohst(&x, &norm, &npd) != RABIN_SUCCESS);
  rmatq_clear(&npd);

  rvec_clear(&x);
  rq_clear(&norm);
  rmatq_clear(&A);
}

static void test_helpers(void)
{
  // rmat_transpose
  rmat_t M = {0}, T = {0};
  rmat_init(&M, 2, 3);
  rmat_init(&T, 3, 2);
  for (u64 r = 0; r < 2; r++) {
    for (u64 c = 0; c < 3; c++) {
      rz_t z;
      rz_init(&z);
      rz_set_i64(&z, (i64)(r * 10 + c));
      rmat_set(&M, &z, r, c);
      rz_clear(&z);
    }
  }
  CHECK(rmat_transpose(&T, &M) == RABIN_SUCCESS);
  CHECK(T.rows == 3 && T.cols == 2);
  // T[c][r] == M[r][c]; verify a couple of entries
  {
    rz_t a, b;
    rz_init_multi(&a, &b, NULL);
    rz_copy(&a, RMAT_GET(&M, 0, 2));
    rz_copy(&b, RMAT_GET(&T, 2, 0));
    CHECK(rz_cmp(&a, &b) == 0);
    rz_copy(&a, RMAT_GET(&M, 1, 0));
    rz_copy(&b, RMAT_GET(&T, 0, 1));
    CHECK(rz_cmp(&a, &b) == 0);
    rz_clear_multi(&a, &b, NULL);
  }
  rmat_clear(&M);
  rmat_clear(&T);

  // rmat_inv_unimodular: inverse of [[2, 1], [1, 1]] is [[1, -1], [-1, 2]]
  rmat_t A = {0}, B = {0};
  rmat_init(&A, 2, 2);
  rmat_init(&B, 2, 2);
  i64 av[4] = {2, 1, 1, 1};
  for (u64 r = 0; r < 2; r++) {
    for (u64 c = 0; c < 2; c++) {
      rz_t z;
      rz_init(&z);
      rz_set_i64(&z, av[r * 2 + c]);
      rmat_set(&A, &z, r, c);
      rz_clear(&z);
    }
  }
  CHECK(rmat_inv_unimodular(&B, &A) == RABIN_SUCCESS);
  {
    rz_t z;
    rz_init(&z);
    rz_set_i64(&z, 1);
    CHECK(rz_cmp(RMAT_GET(&B, 0, 0), &z) == 0);
    rz_set_i64(&z, -1);
    CHECK(rz_cmp(RMAT_GET(&B, 0, 1), &z) == 0);
    CHECK(rz_cmp(RMAT_GET(&B, 1, 0), &z) == 0);
    rz_set_i64(&z, 2);
    CHECK(rz_cmp(RMAT_GET(&B, 1, 1), &z) == 0);
    rz_clear(&z);
  }
  rmat_clear(&A);
  rmat_clear(&B);

  // |det| != 1 is rejected
  rmat_init(&A, 2, 2);
  rmat_init(&B, 2, 2);
  i64 dv[4] = {2, 0, 0, 1};
  for (u64 r = 0; r < 2; r++) {
    for (u64 c = 0; c < 2; c++) {
      rz_t z;
      rz_init(&z);
      rz_set_i64(&z, dv[r * 2 + c]);
      rmat_set(&A, &z, r, c);
      rz_clear(&z);
    }
  }
  CHECK(rmat_inv_unimodular(&B, &A) == RABIN_ERR_INVALID_ARG);
  rmat_clear(&A);

  // rmatq_inv on a rational 2x2: inverse of [[1/2, 0], [0, 2]] is [[2, 0],
  // [0, 1/2]]
  rmatq_t Q = {0}, Qi = {0};
  rmatq_init(&Q, 2, 2);
  rmatq_init(&Qi, 2, 2);
  set_rat(&Q, 0, 0, 1, 2);
  set_rat(&Q, 1, 1, 2, 1);
  CHECK(rmatq_inv(&Qi, &Q) == RABIN_SUCCESS);
  {
    rq_t a;
    rq_init_multi(&a, NULL);
    rmatq_get(&a, &Qi, 0, 0);
    CHECK(rq_eq(&a, 2, 1));
    rmatq_get(&a, &Qi, 1, 1);
    CHECK(rq_eq(&a, 1, 2));
    rq_clear(&a);
  }
  rmatq_clear(&Q);
  rmatq_clear(&Qi);

  // rmatq_mul: [[1,2],[3,4]] * [[5,6],[7,8]] = [[19,22],[43,50]]
  rmatq_t X = {0}, Y = {0}, Z = {0};
  rmatq_init(&X, 2, 2);
  rmatq_init(&Y, 2, 2);
  rmatq_init(&Z, 2, 2);
  i64 xv[4] = {1, 2, 3, 4};
  i64 yv[4] = {5, 6, 7, 8};
  for (u64 r = 0; r < 2; r++) {
    for (u64 c = 0; c < 2; c++) {
      set_int(&X, r, c, xv[r * 2 + c]);
      set_int(&Y, r, c, yv[r * 2 + c]);
    }
  }
  CHECK(rmatq_mul(&Z, &X, &Y) == RABIN_SUCCESS);
  {
    rq_t a;
    rq_init_multi(&a, NULL);
    rmatq_get(&a, &Z, 0, 0);
    CHECK(rq_eq(&a, 19, 1));
    rmatq_get(&a, &Z, 0, 1);
    CHECK(rq_eq(&a, 22, 1));
    rmatq_get(&a, &Z, 1, 0);
    CHECK(rq_eq(&a, 43, 1));
    rmatq_get(&a, &Z, 1, 1);
    CHECK(rq_eq(&a, 50, 1));
    rq_clear(&a);
  }
  rmatq_clear(&X);
  rmatq_clear(&Y);
  rmatq_clear(&Z);
}

int main(void)
{
  rz_init_constants();

  test_known();
  test_rounding();
  test_random();
  test_errors();
  test_helpers();

  if (failures) {
    printf("test_fincke_pohst: %d failure(s)\n", failures);
    return EXIT_FAILURE;
  }
  printf("test_fincke_pohst: OK\n");
  return EXIT_SUCCESS;
}
