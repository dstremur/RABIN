#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include "../include/rabin.h"
// Adjust these include paths to match your project structure
#include "../include/rntt.h"
#include "../include/rpol.h"
#include "../include/u64.h"

// Note: Replace 'rz_println' below with whatever printing function
// your custom rz_t library exposes (e.g., rz_print, rz_dump, etc.)

void rntt_ctx_print_debug(const rntt_ctx_t* ctx)
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
  rz_println(&ctx->q);
  printf("N Inverse (n_inv)    : ");
  rz_println(&ctx->n_inv);

  printf("\n--- OMEGA POWERS (w^i mod q) ---\n");
  printf("omega[0]       (Expected: 1)   : ");
  rz_println(&ctx->omega_powers[0]);
  printf("omega[1]       (Base Root)     : ");
  rz_println(&ctx->omega_powers[1]);
  printf("omega[n/2]     (Expected: q-1) : ");
  rz_println(&ctx->omega_powers[ctx->n / 2]);
  printf("omega[n-1]                     : ");
  rz_println(&ctx->omega_powers[ctx->n - 1]);

  printf("\n--- OMEGA INVERSE POWERS (w^-i mod q) ---\n");
  printf("omega_inv[0]   (Expected: 1)   : ");
  rz_println(&ctx->omega_inv_powers[0]);
  printf("omega_inv[1]   (Base Inv Root) : ");
  rz_println(&ctx->omega_inv_powers[1]);
  printf("omega_inv[n/2] (Expected: q-1) : ");
  rz_println(&ctx->omega_inv_powers[ctx->n / 2]);

  if (ctx->psi_powers != NULL) {
    printf("\n--- PSI POWERS (psi^i mod q) ---\n");
    printf("psi[0]         (Expected: 1)   : ");
    rz_println(&ctx->psi_powers[0]);
    printf("psi[1]         (Base Psi Root) : ");
    rz_println(&ctx->psi_powers[1]);
    printf("psi[n/2]       (Square rt of -1): ");
    rz_println(&ctx->psi_powers[ctx->n / 2]);
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
  rz_init_constants();
  rntt_ctx_t ctx = {0};

  // k = 3 implies n = 8. c = 1.
  // This will cleanly pick the Proth Prime q = 17 (1 * 2^4 + 1)
  u64 k = 3;
  u64 c = 2738;
  u64 n = 1ULL << k;

  printf("--- Step 1: Initializing NTT Context for n = %llu ---\n",
         (unsigned long long)n);
  if (rntt_ctx_init_golden(&ctx, k) != RABIN_SUCCESS) {
    printf("ERROR: Failed to initialize NTT context.\n");
    return 1;
  }

  rntt_ctx_print_debug(&ctx);

  printf("Context generated successfully!\n");
  printf("Modulus q = ");
  rz_println(&ctx.q);
  printf("Principal root omega = ");
  rz_println(&ctx.omega_powers[1]);  // omega^1

  // --- Step 2: Initialize Test Polynomials ---
  rpol_t a = {0}, a_hat = {0}, a_res = {0};
  rpol_init(&a);
  rpol_init(&a_hat);
  rpol_init(&a_res);
  rpol_alloc(&a, n);
  rpol_alloc(&a_hat, n);
  rpol_alloc(&a_res, n);

  // Fill input polynomial 'a' with simple sequential values: a[i] = i
  for (u64 i = 0; i < n; i++) {
    rz_set_u64(&a.coeff[i], i);
  }

  printf("\n--- Step 3: Original Input Polynomial (a) ---\n");
  for (u64 i = 0; i < n; i++) {
    printf("a[%llu] = ", (unsigned long long)i);
    rz_println(&a.coeff[i]);
  }

  a.deg = n - 1;

  // --- Step 4: Execute Forward NTT ---
  printf("\n--- Step 4: Running Forward NTT ---\n");
  rntt_cyclic_forward(&a_hat, &a, &ctx);

  printf("NTT Frequency Domain Representation (a_hat):\n");
  for (u64 i = 0; i < n; i++) {
    printf("a_hat[%llu] = ", (unsigned long long)i);
    rz_println(&a_hat.coeff[i]);
  }

  // --- Step 5: Execute Inverse NTT ---
  printf("\n--- Step 5: Running Inverse NTT ---\n");
  rntt_cyclic_inverse(&a_res, &a_hat, &ctx);

  printf("Recovered Time Domain Representation (a_res):\n");
  bool round_trip_success = true;
  for (u64 i = 0; i < n; i++) {
    printf("a_res[%llu] = ", (unsigned long long)i);
    rz_println(&a_res.coeff[i]);
  }

  rpol_clear(&a);
  rpol_clear(&a_hat);
  rpol_clear(&a_res);
  rntt_ctx_clear(&ctx);
  rz_clear_constants();

  return 0;
}
