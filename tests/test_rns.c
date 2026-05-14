#include <stdio.h>

#include "../include/bigrns.h"
#include "../include/primes.h"
int main()
{
  // 1. Setup Context
  ctx_rns ctx;
  u64 primes[] = {0xFFFFFFFFFFFFFFC5ULL,
                  0xFFFFFFFFFFFFFF4DULL};  // Example 64-bit primes
  rns_context_init(&ctx, RNS_PRIMES, 20);

  // 2. Prepare Numbers
  bignum n;
  bn_init(&n);
  bn_init_val(&n, "100112324252335235235252350");

  rns_num r;
  r.residues = NULL;  // bignum_to_rns will malloc if NULL

  // 3. Convert and Convert Back
  bignum_to_rns(&r, &n, &ctx);

  bignum reconstructed;
  bn_init(&reconstructed);

  rns_to_bignum(&reconstructed, &r, &ctx);

  // Now: reconstructed == my_big

  printf("original :");
  bn_println(&n);

  printf("rns :");
  bn_println(&reconstructed);

  bn_free_multi(&n, &reconstructed, NULL);

  printf("\n--- Debugging RNS Layer ---\n");

  // TEST 1: Round Trip
  bignum test_val, recon_val;
  bn_init_multi(&test_val, &recon_val, NULL);
  bn_init_val(&test_val, "123456789012345678901234567890");

  r.residues = NULL;
  bignum_to_rns(&r, &test_val, &ctx);
  rns_to_bignum(&recon_val, &r, &ctx);

  if (bn_cmp(&test_val, &recon_val) == 0) {
    printf("[PASS] CRT Round Trip\n");
  } else {
    printf("[FAIL] CRT Round Trip! reconstructed: ");
    bn_println(&recon_val);
  }

  // TEST 2: Sign Flip
  bignum m_val, sign_test;
  bn_init_multi(&m_val, &sign_test, NULL);
  bn_copy(&m_val, &ctx.prod);

  // Create a value that is (M - 5) -> should become -5
  bignum five;
  bn_init(&five);
  bn_set_u64(&five, 5);
  bn_sub(&sign_test, &m_val, &five);

  bignum m_half;
  bn_init(&m_half);
  bn_copy(&m_half, &m_val);
  bn_rshift1(&m_half);

  if (bn_cmp(&sign_test, &m_half) > 0) {
    bn_sub(&sign_test, &sign_test, &m_val);
  }

  printf("Sign test (M-5) result: ");
  bn_println(&sign_test);
  // THIS SHOULD PR&INT -5. If it prints a giant number, your bn_sub is the bug.

  bn_free_multi(&test_val, &recon_val, &m_val, &sign_test, &five, &m_half,
                NULL);
  free(r.residues);
}
