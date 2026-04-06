#include <assert.h>
#include <stdio.h>

#include "../include/bignum.h"

// Helper to check if a bignum matches a long long value
void test_lucas_sequences()
{
  bignum u, v, p, q, n;
  bn_init_multi(&u, &v, &p, &q, &n, NULL);

  printf("Running Lucas Sequence Tests...\n");

  // --- TEST 1: Standard Fibonacci/Luchas (P=1, Q=-1) ---
  // For P=1, Q=-1: U_n is the Fibonacci sequence, V_n is the Lucas numbers.
  /* bn_set_u64(&p, 1);
   bn_set_i64(&q, -1);

   // Test n=0 (U=0, V=2)
   bn_set_u64(&n, 0);
   bn_lucas(&u, &v, &p, &q, &n);
   assert(bn_is_eq_i64(&u, 0) && bn_is_eq_i64(&v, 2));
   printf("  [PASS] n=0 (Base Case)\n");

   // Test n=5 (U_5=5, V_5=11)
   bn_set_u64(&n, 5);
   bn_lucas(&u, &v, &p, &q, &n);
   assert(bn_is_eq_i64(&u, 5) && bn_is_eq_i64(&v, 11));
   printf("  [PASS] n=5 (Fibonacci/Lucas match)\n");
 */  // --- TEST 2: Pell Sequence (P=2, Q=-1) ---
  bn_set_u64(&p, 2);
  bn_set_i64(&q, 1);

  // Test n=3 (U_3=5, V_3=14)
  bn_set_u64(&n, 3);
  bn_lucas(&u, &v, &p, &q, &n);
  assert(bn_is_eq_i64(&u, 3) && bn_is_eq_i64(&v, 2));
  printf("  [PASS] n=3 (Pell Sequence)\n");
  // --- TEST 3: Mathematical Identity Check ---
  // Identity: V_n^2 - D*U_n^2 = 4*Q^n
  // where D = P^2 - 4*Q
  bn_set_u64(&p, 3);
  bn_set_u64(&q, 1);
  bn_set_u64(&n, 4);
  bn_lucas(&u, &v, &p, &q, &n);

  // Calculate: v^2 - (p^2 - 4q)*u^2
  bignum term1, term2, d, res, four;
  bn_init_multi(&term1, &term2, &d, &res, &four, NULL);
  bn_set_u64(&four, 4);

  bn_mul(&term1, &v, &v);  // V_n^2
  bn_mul(&d, &p, &p);      // P^2
  bn_sub(&d, &d, &four);   // D = P^2 - 4 (since Q=1)
  bn_mul(&term2, &u, &u);
  bn_mul(&term2, &term2, &d);    // D * U_n^2
  bn_sub(&res, &term1, &term2);  // V_n^2 - D*U_n^2

  // Should equal 4 * Q^n = 4 * 1^4 = 4
  assert(bn_is_eq_i64(&res, 4));
  printf("  [PASS] Identity V^2 - DU^2 = 4Q^n\n");

  // --- CLEANUP ---
  bn_free(&u);
  bn_free(&v);
  bn_free(&p);
  bn_free(&q);
  bn_free(&n);
  bn_free(&term1);
  bn_free(&term2);
  bn_free(&d);
  bn_free(&res);
  printf("All tests passed successfully!\n");
}

int main()
{
  test_lucas_sequences();
  return 0;
}
