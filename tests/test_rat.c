
#include "../include/rmat.h"
#include "../include/rmatq.h"
#include "../include/rq.h"

int main()
{
  rq_t r, q, p, k, l;
  rq_init(&r);
  rq_init(&q);
  rq_init(&p);
  rq_init_multi(&k, &l, NULL);

  rz_set_i64(&r.den, 8905493208);
  rz_set_i64(&r.num, 3);

  rz_set_i64(&q.den, 23534563456);
  rz_set_i64(&q.num, 3);

  rq_print(&r);

  rq_mul(&p, &r, &q);
  rq_mul(&r, &p, &q);

  rq_print(&p);
  rq_print(&r);

  // r = 1 / 3
  rz_set_i64(&r.den, 4);
  rz_set_i64(&r.num, -1);

  // q = 1 / 4
  rz_set_i64(&q.den, 4);
  rz_set_i64(&q.num, 1);

  rq_div(&p, &r, &q);

  rq_print(&p);

  rq_add(&p, &r, &q);

  rq_print(&p);

  /*
  for (u64 x = 0; x < 100; x++) {
    rz_gen_random(&r.den, 100);
    rz_gen_random(&r.num, 100);
    rz_gen_random(&q.den, 100);
    rz_gen_random(&q.num, 100);

    rq_div(&k, &r, &q);

        rq_print(&k);
    rq_div(&l, &q, &r);

    rq_mul(&r, &k, &l);

  //  rq_print(&r);
  }
*/
  rq_clear(&r);
  rq_clear(&p);
  rq_clear(&q);
  rq_clear_multi(&k, &l, NULL);

  rmatq_t A = {0}, B = {0}, R_add = {0}, R_sub = {0};

  // 1. Initialize 5x5 rational matrices
  rmatq_init(&A, 5, 5);
  rmatq_init(&B, 5, 5);
  rmatq_init(&R_add, 5, 5);
  rmatq_init(&R_sub, 5, 5);

  // Set non-unit denominators (e.g., A denominator = 2, B denominator = 3)
  rz_set_u64(&A.den, 2);
  rz_set_u64(&B.den, 3);

  // Temporary rz_t helper for populating values
  rz_t val;
  rz_init(&val);

  // 2. Fill both matrices with sample values
  for (u64 i = 0; i < 5; i++) {
    for (u64 j = 0; j < 5; j++) {
      // Fill A numerator entries with (i + j + 1)
      rz_set_u64(&val, i + j + 1);
      rmat_set(&A.M, &val, i, j);

      // Fill B numerator entries with (i * j + 2)
      rz_set_u64(&val, i * j + 2);
      rmat_set(&B.M, &val, i, j);
    }
  }
  rz_clear(&val);

  // 3. Perform addition and subtraction
  rmatq_add(&R_add, &A, &B);
  rmatq_normalize(&R_add);  // Reduce fractions to lowest terms

  rmatq_sub(&R_sub, &A, &B);
  rmatq_normalize(&R_sub);

  printf("Matrix A:\n");
  rmatq_print(&A);

  printf("Matrix B:\n");
  rmatq_print(&B);
  // 4. Output results
  printf("Matrix A + B:\n");
  rmatq_print(&R_add);

  printf("\nMatrix A - B:\n");
  rmatq_print(&R_sub);

  // 5. Free allocated memory
  rmatq_clear(&A);
  rmatq_clear(&B);
  rmatq_clear(&R_add);
  rmatq_clear(&R_sub);
}
