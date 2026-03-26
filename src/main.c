// main.c
#include <assert.h>
#include <fcntl.h>
#include <stdint.h>
#include <stdio.h>
#include <unistd.h>

#include "../include/bignum.h"
#include "ctype.h"
#include "string.h"

typedef struct {
  const char* n_hex;
  const char* base_hex;
  bool expected;
  const char* description;
} test_case;

void bn_from_hex(bignum* r, const char* hex) {
  bn_set_u64(r, 0);
  size_t len = strlen(hex);
  bignum sixteen, val;
  bn_init(&sixteen);
  bn_init(&val);
  bn_set_u64(&sixteen, 16);

  for (size_t i = 0; i < len; i++) {
    char c = tolower(hex[i]);
    uint64_t digit = 0;
    if (c >= '0' && c <= '9')
      digit = c - '0';
    else if (c >= 'a' && c <= 'f')
      digit = c - 'a' + 10;
    else
      continue;

    bn_mul(r, r, &sixteen);
    bn_set_u64(&val, digit);
    bn_add(r, r, &val);
  }
  bn_free(&sixteen);
  bn_free(&val);
}

void test_jacobi_hang() {
  printf("--- Testing Jacobi Zero Trap ---\n");
  bignum a, n;
  bn_init(&a);
  bn_init(&n);

  // If a is a multiple of n, bn_mod in jacobi will produce 0
  bn_set_u64(&a, 5);
  bn_set_u64(&n, 10);  // Candidate is a multiple of D=5

  printf("Testing jacobi(5, 10)... (Should not hang)\n");
  int j = bn_jacobi(&a, &n);
  printf("[PASS] Jacobi result: %d\n", j);

  bn_free(&a);
  bn_free(&n);
}

void run_test(test_case tc) {
  bignum n;
  bn_init(&n);

  // Note: You'll need a bn_from_hex function.
  // If you don't have one, use bn_set_u64 for small tests
  // or see the helper below.
  bn_from_hex(&n, tc.n_hex);

  u64 base = strtoull(tc.base_hex, NULL, 16);

  bool result = bn_millerRabin(&n, base);

  printf("[%s] Testing n=%s, base=%s... ",
         (result == tc.expected ? "PASS" : "FAIL"), tc.n_hex, tc.base_hex);
  printf("(Expected %s, Got %s) - %s\n", tc.expected ? "PRIME" : "COMPOSITE",
         result ? "PRIME" : "COMPOSITE", tc.description);

  bn_free(&n);
}

void test_addition() {
  printf("--- Testing Addition ---\n");
  bignum a, b, res;
  bn_init(&a);
  bn_init(&b);
  bn_init(&res);

  // Test 1: Simple Addition
  bn_set_u64(&a, 500);
  bn_set_u64(&b, 1000);
  bn_add(&res, &a, &b);
  assert(res.size == 1 && res.limbs[0] == 1500);
  printf("[PASS] Simple Addition\n");

  // Test 2: The Carry Test (UINT64_MAX + 1)
  // 0xFFFFFFFFFFFFFFFF + 1 should equal 0x00000000000000010000000000000000
  bn_set_u64(&a, UINT64_MAX);
  bn_set_u64(&b, 1);
  bn_add(&res, &a, &b);

  assert(res.size == 2);      // Must have expanded to 2 limbs
  assert(res.limbs[0] == 0);  // Bottom limb wrapped around to 0
  assert(res.limbs[1] == 1);  // Top limb caught the carry
  printf("[PASS] Cross-Limb Carry Addition\n");

  bn_free(&a);
  bn_free(&b);
  bn_free(&res);
}

void test_multiplication() {
  printf("--- Testing Multiplication ---\n");
  bignum a, b, res;
  bn_init(&a);
  bn_init(&b);
  bn_init(&res);

  // Test 1: UINT64_MAX * UINT64_MAX
  // (2^64 - 1) * (2^64 - 1) = 2^128 - 2*(2^64) + 1
  // In hex: 0xFFFFFFFFFFFFFFFE0000000000000001
  bn_set_u64(&a, UINT64_MAX - 1);
  bn_set_u64(&b, UINT64_MAX);
  bn_mul(&res, &a, &b);
  bn_mul(&res, &res, &b);

  bn_print_hex("Max64 * Max64 Result", &res);

  assert(res.size == 2);
  assert(res.limbs[0] == 1);               // Lower 64 bits
  assert(res.limbs[1] == UINT64_MAX - 1);  // Upper 64 bits
  printf("[PASS] 64-bit x 64-bit Multiplication\n");

  bn_free(&a);
  bn_free(&b);
  bn_free(&res);
}

void bn_smallTest() {
  test_case tests[] = {
      // Small Primes
      {"3", "2", true, "Small prime 3"},
      {"d", "2", true, "Small prime 13"},
      {"1f", "2", true, "Small prime 31"},

      // Small Composites
      {"9", "2", false, "Small composite 9"},
      {"f", "2", false, "Small composite 15"},
      {"5b", "2", false, "Small composite 91 (7*13)"},

      // Carmichael Numbers (Fools Fermat but not Miller-Rabin)
      {"231", "3", false, "Carmichael 561"},
      {"451", "5", false, "Carmichael 1105"},

      // Large Primes (Example: 2^127 - 1, Mersenne Prime)
      {"7fffffffffffffffffffffffffffffff", "3", true, "Mersenne Prime M127"},
      // Large Composites (2^127 - 1) + 1 (which is even)
      {"80000000000000000000000000000000", "2", false, "Large even number"},

      // A specific composite that passes base 2 (Strong Pseudoprime)
      {"801", "2", false, "Strong Pseudoprime base 2 (2049)"},
      {"801", "3", false, "Strong Pseudoprime base 2 caught by base 3"}};

  int num_tests = sizeof(tests) / sizeof(test_case);
  for (int i = 0; i < num_tests; i++) {
    run_test(tests[i]);
  }
}

int main(void) {
  bignum a, b, c;
  bn_init(&a);
  bn_init(&b);
  bn_init(&c);
  bn_init_val(&a,
              "3402823468347696734638673489672893467938246798234672983467432986"
              "72398467289346829346734689366920938463463374607431768211455");
  bn_init_val(&b, "18446744073709551616");
  bn_mod(&c, &a, &b);
  bn_print(&c);

  return 0;
}

