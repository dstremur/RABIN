#include <stdio.h>

#include "../include/rmat.h"
#include "../include/rmatq.h"
#include "../include/rq.h"

static int failures = 0;

#define CHECK(cond)                                          \
  do {                                                       \
    if (!(cond)) {                                           \
      printf("FAIL %s:%d: %s\n", __FILE__, __LINE__, #cond); \
      failures++;                                            \
    }                                                        \
  } while (0)

static void set_int(rmatq_t* A, u64 r, u64 c, i64 v)
{
  rz_t z;
  rz_init(&z);
  rz_set_i64(&z, v);
  rmat_set(&A->M, &z, r, c);
  rz_clear(&z);
}

static void set_i(rq_t* a, i64 num, i64 den)
{
  rz_set_i64(&a->num, num);
  rz_set_i64(&a->den, den);
  rq_normalize(a);
}

static void check_eq_entry(rmatq_t* A, u64 r, u64 c, i64 num, i64 den)
{
  rq_t a, want;
  rq_init_multi(&a, &want, NULL);
  rmatq_get(&a, A, r, c);
  set_i(&want, num, den);
  if (rq_cmp(&a, &want) != 0) {
    printf("FAIL: entry (%lu, %lu) = ", r, c);
    rq_print(&a);
    printf(" != ");
    rq_print(&want);
    printf("\n");
    failures++;
  }
  rq_clear_multi(&a, &want, NULL);
}

static void check_cholesky(rmatq_t* C, rmatq_t* A)
{
  u64 n = A->M.rows;

  for (u64 i = 0; i < n; i++) {
    for (u64 j = 0; j < i; j++) {
      rq_t a;
      rq_init(&a);
      rmatq_get(&a, C, i, j);
      if (!rq_is_zero(&a)) {
        printf("FAIL: C[%lu][%lu] is not zero below the diagonal\n", i, j);
        failures++;
      }
      rq_clear(&a);
    }
  }

  rq_t s, t, b, a;
  rq_init_multi(&s, &t, &b, &a, NULL);

  for (u64 i = 0; i < n; i++) {
    for (u64 j = 0; j < n; j++) {
      rq_set_u64(&s, 0);
      for (u64 k = 0; k < n; k++) {
        rmatq_get(&t, C, k, i);
        if (rz_is_zero(&t.num)) continue;
        rmatq_get(&b, C, k, j);
        rq_mul(&t, &t, &b);
        rq_add(&s, &s, &t);
      }
      rmatq_get(&a, A, i, j);
      if (rq_is_zero(&s) != rq_is_zero(&a) ||
          (!rq_is_zero(&s) && rq_cmp(&s, &a) != 0)) {
        printf("FAIL: (C^t C)[%lu][%lu] = ", i, j);
        rq_print(&s);
        printf(" != ");
        rq_print(&a);
        printf("\n");
        failures++;
      }
    }
  }

  rq_clear_multi(&s, &t, &b, &a, NULL);
}

int main()
{
  // 1. Exact integer Cholesky of a 3x3 positive definite form.
  // A = [[4,2,0],[2,5,2],[0,2,5]] = R^t R with R = [[2,1,0],[0,2,1],[0,0,2]].
  rmatq_t A = {0}, C = {0};
  rmatq_init(&A, 3, 3);
  rmatq_init(&C, 3, 3);

  set_int(&A, 0, 0, 4);
  set_int(&A, 0, 1, 2);
  set_int(&A, 1, 0, 2);
  set_int(&A, 1, 1, 5);
  set_int(&A, 1, 2, 2);
  set_int(&A, 2, 1, 2);
  set_int(&A, 2, 2, 5);

  CHECK(rmatq_cholesky(&C, &A) == RABIN_SUCCESS);
  CHECK(rz_cmp(&C.den, &A.den) == 0);
  check_eq_entry(&C, 0, 0, 2, 1);
  check_eq_entry(&C, 0, 1, 1, 1);
  check_eq_entry(&C, 1, 1, 2, 1);
  check_eq_entry(&C, 1, 2, 1, 1);
  check_eq_entry(&C, 2, 2, 2, 1);
  check_cholesky(&C, &A);

  rmatq_clear(&A);
  rmatq_clear(&C);

  // 2. Pivot 2 is not a rational square -> no Cholesky over Q.
  rmatq_init(&A, 2, 2);
  rmatq_init(&C, 2, 2);
  set_int(&A, 0, 0, 2);
  set_int(&A, 0, 1, 1);
  set_int(&A, 1, 0, 1);
  set_int(&A, 1, 1, 2);
  CHECK(rmatq_cholesky(&C, &A) == RABIN_ERR_INVALID_ARG);
  CHECK(rz_is_zero(RMAT_GET(&C.M, 0, 0)));
  rmatq_clear(&A);
  rmatq_clear(&C);

  // 3. Zero pivot -> division by zero.
  rmatq_init(&A, 2, 2);
  rmatq_init(&C, 2, 2);
  set_int(&A, 0, 1, 1);
  set_int(&A, 1, 0, 1);
  set_int(&A, 1, 1, 2);
  CHECK(rmatq_cholesky(&C, &A) == RABIN_ERR_DIV_BY_ZERO);
  rmatq_clear(&A);
  rmatq_clear(&C);

  // 4. 1x1 decomposition sqrt(4/25) = 2/5.
  rmatq_init(&A, 1, 1);
  rmatq_init(&C, 1, 1);
  rz_set_u64(&A.den, 25);
  set_int(&A, 0, 0, 4);
  CHECK(rmatq_cholesky(&C, &A) == RABIN_SUCCESS);
  check_eq_entry(&C, 0, 0, 2, 5);
  rmatq_clear(&A);
  rmatq_clear(&C);

  // 5. Element accessors: get/set round-trip with denominator growth.
  rmatq_init(&A, 2, 2);
  rz_set_u64(&A.den, 6);
  set_int(&A, 0, 0, 2);
  set_int(&A, 0, 1, 3);
  set_int(&A, 1, 0, 4);
  set_int(&A, 1, 1, 5);

  check_eq_entry(&A, 0, 0, 1, 3);

  rq_t q;
  rq_init(&q);
  set_i(&q, 1, 4);
  CHECK(rmatq_set(&A, &q, 0, 1) == RABIN_SUCCESS);
  check_eq_entry(&A, 0, 1, 1, 4);
  check_eq_entry(&A, 1, 1, 5, 6);

  CHECK(rmatq_set_zero(&A, 1, 0) == RABIN_SUCCESS);
  check_eq_entry(&A, 1, 0, 0, 1);

  rq_clear(&q);

  // Bounds checks.
  CHECK(rmatq_get(&q, &A, 2, 0) == RABIN_ERR_INVALID_ARG);
  CHECK(rmatq_set(&A, &q, 0, 2) == RABIN_ERR_INVALID_ARG);
  CHECK(rmatq_set_zero(&A, 3, 3) == RABIN_ERR_INVALID_ARG);

  rmatq_clear(&A);

  if (failures == 0) {
    printf("All cholesky tests passed.\n");
    return 0;
  }

  printf("%d failure(s).\n", failures);
  return 1;
}
