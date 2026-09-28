#include <stdio.h>

#include "../include/rns.h"
#include "../include/rns_primes.h"
int main()
{
  // 1. Setup Context
  rns_ctx_t ctx = {0};
  if (rns_ctx_init(&ctx, RNS_PRIMES, 20) != RABIN_SUCCESS) return 1;

  // 2. Prepare Numbers
  rz_t n;
  rz_init(&n);
  rz_init_val(&n, "100112324252335235235252350");

  rns_num_t r;
  r.residues = NULL;  // rns_import will malloc if NULL

  // 3. Convert and Convert Back
  rns_import(&r, &n, &ctx);

  rz_t reconstructed;
  rz_init(&reconstructed);

  // NOTE: rns_to_bignum_garner never existed; use rns_export.
  rns_export(&reconstructed, &r, &ctx);

  // Now: reconstructed == my_big

  printf("original :");
  rz_println(&n);

  printf("rns :");
  rz_println(&reconstructed);

  rz_clear_multi(&n, &reconstructed, NULL);

  printf("\n--- Debugging RNS Layer ---\n");

  // TEST 1: Round Trip
  rz_t test_val, recon_val;
  rz_init_multi(&test_val, &recon_val, NULL);
  rz_init_val(&test_val, "123456789012345678901234567890");

  r.residues = NULL;
  rns_import(&r, &test_val, &ctx);
  rns_export(&recon_val, &r, &ctx);

  if (rz_cmp(&test_val, &recon_val) == 0) {
    printf("[PASS] CRT Round Trip\n");
  } else {
    printf("[FAIL] CRT Round Trip! reconstructed: ");
    rz_println(&recon_val);
  }

  // TEST 2: Sign Flip
  rz_t m_val, sign_test;
  rz_init_multi(&m_val, &sign_test, NULL);
  rz_copy(&m_val, &ctx.prod);

  // Create a value that is (M - 5) -> should become -5
  rz_t five;
  rz_init(&five);
  rz_set_u64(&five, 5);
  rz_sub(&sign_test, &m_val, &five);

  rz_t m_half;
  rz_init(&m_half);
  rz_copy(&m_half, &m_val);
  rz_rshift1(&m_half);

  if (rz_cmp(&sign_test, &m_half) > 0) {
    rz_sub(&sign_test, &sign_test, &m_val);
  }

  printf("Sign test (M-5) result: ");
  rz_println(&sign_test);
  // THIS SHOULD PR&INT -5. If it prints a giant number, your rz_sub is the bug.

  rz_clear_multi(&test_val, &recon_val, &m_val, &sign_test, &five, &m_half,
                 NULL);
  free(r.residues);
  rns_ctx_clear(&ctx);
}
