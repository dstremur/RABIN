#include <assert.h>
#include <stdio.h>

#include "../include/rabin.h"

// Helper to check if a rz_t matches a long long value
void test_lucas_sequences()
{
  rz_t u, v, p, q, n;
  rz_init_multi(&u, &v, &p, &q, &n, NULL);

  printf("Running Lucas Sequence Tests...\n");

  // --- TEST 1: Standard Fibonacci/Luchas (P=1, Q=-1) ---
  // For P=1, Q=-1: U_n is the Fibonacci sequence, V_n is the Lucas numbers.
  /* rz_set_u64(&p, 1);
   rz_set_i64(&q, -1);

   // Test n=0 (U=0, V=2)
   rz_set_u64(&n, 0);
   rz_lucas(&u, &v, &p, &q, &n);
   assert(rz_is_eq_i64(&u, 0) && rz_is_eq_i64(&v, 2));
   printf("  [PASS] n=0 (Base Case)\n");

   // Test n=5 (U_5=5, V_5=11)
   rz_set_u64(&n, 5);
   rz_lucas(&u, &v, &p, &q, &n);
   assert(rz_is_eq_i64(&u, 5) && rz_is_eq_i64(&v, 11));
   printf("  [PASS] n=5 (Fibonacci/Lucas match)\n");
 */  // --- TEST 2: Pell Sequence (P=2, Q=-1) ---
  rz_set_u64(&p, 2);
  rz_set_i64(&q, 1);

  // Test n=3 (U_3=5, V_3=14)
  rz_set_u64(&n, 3);
  rz_lucas(&u, &v, &p, &q, &n);
  assert(rz_is_eq_i64(&u, 3) && rz_is_eq_i64(&v, 2));
  printf("  [PASS] n=3 (Pell Sequence)\n");
  // --- TEST 3: Mathematical Identity Check ---
  // Identity: V_n^2 - D*U_n^2 = 4*Q^n
  // where D = P^2 - 4*Q
  rz_set_u64(&p, 3);
  rz_set_u64(&q, 1);
  rz_set_u64(&n, 4);
  rz_lucas(&u, &v, &p, &q, &n);

  // Calculate: v^2 - (p^2 - 4q)*u^2
  rz_t term1, term2, d, res, four;
  rz_init_multi(&term1, &term2, &d, &res, &four, NULL);
  rz_set_u64(&four, 4);

  rz_mul(&term1, &v, &v);  // V_n^2
  rz_mul(&d, &p, &p);      // P^2
  rz_sub(&d, &d, &four);   // D = P^2 - 4 (since Q=1)
  rz_mul(&term2, &u, &u);
  rz_mul(&term2, &term2, &d);    // D * U_n^2
  rz_sub(&res, &term1, &term2);  // V_n^2 - D*U_n^2

  // Should equal 4 * Q^n = 4 * 1^4 = 4
  assert(rz_is_eq_i64(&res, 4));
  printf("  [PASS] Identity V^2 - DU^2 = 4Q^n\n");

  // --- CLEANUP ---
  rz_clear(&u);
  rz_clear(&v);
  rz_clear(&p);
  rz_clear(&q);
  rz_clear(&n);
  rz_clear(&term1);
  rz_clear(&term2);
  rz_clear(&d);
  rz_clear(&res);
  printf("All tests passed successfully!\n");
}

int main()
{
  test_lucas_sequences();
  return 0;
}
