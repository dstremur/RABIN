#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include "../include/bignum.h"
// Adjust these include paths to match your project structure
#include "../include/bigntt.h"
#include "../include/bigpoly.h"
#include "../include/u64.h"

// Note: Replace 'bn_println' below with whatever printing function
// your custom bignum library exposes (e.g., bn_print, bn_dump, etc.)

void bigntt_ctx_print_debug(const ntt_ctx* ctx)
{
  if (!ctx) {
    printf("Context is NULL\n");
    return;
  }

  printf("=========================================\n");
  printf("         NTT CONTEXT DEBUG LOG           \n");
  printf("=========================================\n");
  printf("Transform Length (n) : %llu\n", ctx->n);
  printf("Modulus (q)          : ");
  bn_println(&ctx->q);
  printf("N Inverse (n_inv)    : ");
  bn_println(&ctx->n_inv);

  printf("\n--- OMEGA POWERS (w^i mod q) ---\n");
  printf("omega[0]       (Expected: 1)   : ");
  bn_println(&ctx->omega_powers[0]);
  printf("omega[1]       (Base Root)     : ");
  bn_println(&ctx->omega_powers[1]);
  printf("omega[n/2]     (Expected: q-1) : ");
  bn_println(&ctx->omega_powers[ctx->n / 2]);
  printf("omega[n-1]                     : ");
  bn_println(&ctx->omega_powers[ctx->n - 1]);

  printf("\n--- OMEGA INVERSE POWERS (w^-i mod q) ---\n");
  printf("omega_inv[0]   (Expected: 1)   : ");
  bn_println(&ctx->omega_inv_powers[0]);
  printf("omega_inv[1]   (Base Inv Root) : ");
  bn_println(&ctx->omega_inv_powers[1]);
  printf("omega_inv[n/2] (Expected: q-1) : ");
  bn_println(&ctx->omega_inv_powers[ctx->n / 2]);

  if (ctx->psi_powers != NULL) {
    printf("\n--- PSI POWERS (psi^i mod q) ---\n");
    printf("psi[0]         (Expected: 1)   : ");
    bn_println(&ctx->psi_powers[0]);
    printf("psi[1]         (Base Psi Root) : ");
    bn_println(&ctx->psi_powers[1]);
    printf("psi[n/2]       (Square rt of -1): ");
    bn_println(&ctx->psi_powers[ctx->n / 2]);
  }

  printf("\n--- BIT REVERSAL INDICES (First 8 Samples) ---\n");
  u64 max_samples = (ctx->n < 8) ? ctx->n : 8;
  for (u64 i = 0; i < max_samples; i++) {
    printf("  Index [%2llu] ---> Bit-Reversed Index [%2llu]\n", i,
           ctx->bit_rev_indices[i]);
  }
  printf("=========================================\n");
}

int main(void)
{
  bn_init_constants();
  ntt_ctx ctx;

  // k = 3 implies n = 8. c = 1.
  // This will cleanly pick the Proth Prime q = 17 (1 * 2^4 + 1)
  u64 k = 3;
  u64 c = 2738;
  u64 n = 1ULL << k;

  printf("--- Step 1: Initializing NTT Context for n = %llu ---\n",
         (unsigned long long)n);
  if (!bigntt_ctx_init_golden(&ctx, k)) {
    printf("ERROR: Failed to initialize NTT context.\n");
    return 1;
  }

  bigntt_ctx_print_debug(&ctx);

  printf("Context generated successfully!\n");
  printf("Modulus q = ");
  bn_println(&ctx.q);
  printf("Principal root omega = ");
  bn_println(&ctx.omega_powers[1]);  // omega^1

  // --- Step 2: Initialize Test Polynomials ---
  bigpoly a, a_hat, a_res;
  bigpoly_init(&a);
  bigpoly_init(&a_hat);
  bigpoly_init(&a_res);
  bigpoly_alloc(&a, n);
  bigpoly_alloc(&a_hat, n);
  bigpoly_alloc(&a_res, n);

  // Fill input polynomial 'a' with simple sequential values: a[i] = i
  for (u64 i = 0; i < n; i++) {
    bn_set_u64(&a.coeff[i], i);
  }

  printf("\n--- Step 3: Original Input Polynomial (a) ---\n");
  for (u64 i = 0; i < n; i++) {
    printf("a[%llu] = ", (unsigned long long)i);
    bn_println(&a.coeff[i]);
  }

  a.deg = n - 1;

  // --- Step 4: Execute Forward NTT ---
  printf("\n--- Step 4: Running Forward NTT ---\n");
  bigntt_cyclic_forward(&a_hat, &a, &ctx);

  printf("NTT Frequency Domain Representation (a_hat):\n");
  for (u64 i = 0; i < n; i++) {
    printf("a_hat[%llu] = ", (unsigned long long)i);
    bn_println(&a_hat.coeff[i]);
  }

  // --- Step 5: Execute Inverse NTT ---
  printf("\n--- Step 5: Running Inverse NTT ---\n");
  bigntt_cyclic_inverse(&a_res, &a_hat, &ctx);

  printf("Recovered Time Domain Representation (a_res):\n");
  bool round_trip_success = true;
  for (u64 i = 0; i < n; i++) {
    printf("a_res[%llu] = ", (unsigned long long)i);
    bn_println(&a_res.coeff[i]);
  }

  bigpoly_free(&a);
  bigpoly_free(&a_hat);
  bigpoly_free(&a_res);
  bigntt_ctx_free(&ctx);
  bn_free_constants();

  return 0;
}
