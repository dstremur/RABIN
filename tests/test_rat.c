
#include "../include/bigmatrix.h"
#include "../include/bigmatrixrat.h"
#include "../include/bigrat.h"

int main()
{
  bigrat r, q, p, k, l;
  br_init(&r);
  br_init(&q);
  br_init(&p);
  br_init_multi(&k, &l, NULL);

  bn_set_i64(&r.den, 8905493208);
  bn_set_i64(&r.num, 3);

  bn_set_i64(&q.den, 23534563456);
  bn_set_i64(&q.num, 3);

  br_print(&r);

  br_mul(&p, &r, &q);
  br_mul(&r, &p, &q);

  br_print(&p);
  br_print(&r);

  // r = 1 / 3
  bn_set_i64(&r.den, 4);
  bn_set_i64(&r.num, -1);

  // q = 1 / 4
  bn_set_i64(&q.den, 4);
  bn_set_i64(&q.num, 1);

  br_div(&p, &r, &q);

  br_print(&p);

  br_add(&p, &r, &q);

  br_print(&p);

  /*
  for (u64 x = 0; x < 100; x++) {
    bn_gen_random(&r.den, 100);
    bn_gen_random(&r.num, 100);
    bn_gen_random(&q.den, 100);
    bn_gen_random(&q.num, 100);

    br_div(&k, &r, &q);

        br_print(&k);
    br_div(&l, &q, &r);

    br_mul(&r, &k, &l);

  //  br_print(&r);
  }
*/
  br_free(&r);
  br_free(&p);
  br_free(&q);
  br_free_multi(&k, &l, NULL);

  bigmatrixrat A, B, R_add, R_sub;

  // 1. Initialize 5x5 rational matrices
  bigmatrixrat_init(&A, 5, 5);
  bigmatrixrat_init(&B, 5, 5);
  bigmatrixrat_init(&R_add, 5, 5);
  bigmatrixrat_init(&R_sub, 5, 5);

  // Set non-unit denominators (e.g., A denominator = 2, B denominator = 3)
  bn_set_u64(&A.den, 2);
  bn_set_u64(&B.den, 3);

  // Temporary bignum helper for populating values
  bignum val;
  bn_init(&val);

  // 2. Fill both matrices with sample values
  for (u64 i = 0; i < 5; i++) {
    for (u64 j = 0; j < 5; j++) {
      // Fill A numerator entries with (i + j + 1)
      bn_set_u64(&val, i + j + 1);
      bigmatrix_set(&A.M, &val, i, j);

      // Fill B numerator entries with (i * j + 2)
      bn_set_u64(&val, i * j + 2);
      bigmatrix_set(&B.M, &val, i, j);
    }
  }
  bn_free(&val);

  // 3. Perform addition and subtraction
  bigmatrixrat_add(&R_add, &A, &B);
  bigmatrixrat_normalize(&R_add);  // Reduce fractions to lowest terms

  bigmatrixrat_sub(&R_sub, &A, &B);
  bigmatrixrat_normalize(&R_sub);

  printf("Matrix A:\n");
  bigmatrixrat_print(&A);

  printf("Matrix B:\n");
  bigmatrixrat_print(&B);
  // 4. Output results
  printf("Matrix A + B:\n");
  bigmatrixrat_print(&R_add);

  printf("\nMatrix A - B:\n");
  bigmatrixrat_print(&R_sub);

  // 5. Free allocated memory
  bigmatrixrat_free(&A);
  bigmatrixrat_free(&B);
  bigmatrixrat_free(&R_add);
  bigmatrixrat_free(&R_sub);
}
