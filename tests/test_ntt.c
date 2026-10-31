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

  /* ---- q2 = 69*2^55+1 (second RNS prime) u64 NTT: asserted checks ---- */
  const u64 q2 = 2485986994308513793ULL; /* 69 * 2^55 + 1 */
  int p2_fails = 0;

  /* context: root of unity orders, Montgomery tables, cached identity */
  for (u64 k = 1; k <= 12; k++) {
    u64_ntt_ctx_t* c2 = u64_ntt_ctx_proth2_cached(k);
    if (c2 == NULL || u64_ntt_ctx_proth2_cached(k) != c2) {
      p2_fails++;
      continue;
    }
    u64 nk = 1ULL << k;
    if (c2->n != nk || c2->q != q2 || c2->bit_rev_indices == NULL) p2_fails++;
    /* omega table entry is mont(omega): canonicalize, then check orders */
    u64 om = u64_mont_out(c2->omega_powers[1], &c2->mctx);
    if (u64_mod_pow(om, nk, q2) != 1 || u64_mod_pow(om, nk >> 1, q2) != q2 - 1)
      p2_fails++;
    /* n * n_inv = 1 mod q2 (n_inv stored in Montgomery form) */
    if (u64_mod_mul(nk % q2, u64_mont_out(c2->n_inv, &c2->mctx), q2) != 1)
      p2_fails++;
    /* omega * omega_inv = 1 */
    if (u64_mod_mul(om, u64_mont_out(c2->omega_inv_powers[1], &c2->mctx), q2) !=
        1)
      p2_fails++;
    /* Montgomery residues stay below q2 */
    for (u64 i = 0; i < nk; i++)
      if (c2->omega_powers[i] >= q2 || c2->omega_inv_powers[i] >= q2)
        p2_fails++;
  }

  /* same-prime sanity: q2 is prime, q2-1 = 69 * 2^55 has 2-adicity 55 */
  {
    u64 m = q2 - 1;
    u64 tz = 0;
    while ((m & 1ULL) == 0) {
      m >>= 1;
      tz++;
    }
    if (tz != 55 || m != 69) p2_fails++;
    if (q2 >= (1ULL << 62)) p2_fails++; /* butterfly needs 4q < 2^64 */
    /* q1^{-1} mod q2 used by the CRT must match the kernel constant */
    u64 q1 = 180143985094819841ULL;
    if (u64_mod_mul(q1 % q2, 2291769260378161154ULL, q2) != 1) p2_fails++;
  }

  /* round trip forward (canonical in) -> inverse (canonical out) */
  for (u64 k = 1; k <= 12; k++) {
    u64_ntt_ctx_t* c2 = u64_ntt_ctx_proth2_cached(k);
    u64 nk = 1ULL << k;
    u64* src = malloc(sizeof(u64) * nk);
    u64* mid = malloc(sizeof(u64) * nk);
    u64* dst = malloc(sizeof(u64) * nk);
    if (src == NULL || mid == NULL || dst == NULL) {
      p2_fails++;
    } else {
      for (u64 i = 0; i < nk; i++)
        src[i] = (i == 0) ? q2 - 1 : ((((u64)rand() << 32) ^ (u64)rand()) % q2);
      u64_ntt_cyclic_forward(mid, src, c2);
      u64_ntt_cyclic_inverse_montgomery_in(dst, mid, c2);
      for (u64 i = 0; i < nk; i++)
        if (dst[i] != src[i] || dst[i] >= q2) p2_fails++;
    }
    free(src);
    free(mid);
    free(dst);
  }

  /* cyclic convolution vs naive O(n^2) over the field (n = 16) */
  {
    u64 k = 4, nk = 1ULL << k;
    u64_ntt_ctx_t* c2 = u64_ntt_ctx_proth2_cached(k);
    u64 a_arr[16], b_arr[16], a_hat[16], b_hat[16], c_arr[16], ref[16];
    for (u64 i = 0; i < nk; i++) {
      a_arr[i] = ((u64)rand() << 32) ^ (u64)rand();
      b_arr[i] = ((u64)rand() << 32) ^ (u64)rand();
      a_arr[i] = (i == 0) ? q2 - 1 : a_arr[i] % q2;
      b_arr[i] = (i == nk - 1) ? q2 - 1 : b_arr[i] % q2;
      ref[i] = 0;
    }
    for (u64 i = 0; i < nk; i++)
      for (u64 j = 0; j < nk; j++) {
        u64 idx = i + j;
        if (idx >= nk) idx -= nk;
        ref[idx] =
            u64_mod_add(ref[idx], u64_mod_mul(a_arr[i], b_arr[j], q2), q2);
      }
    u64_ntt_cyclic_forward(a_hat, a_arr, c2);
    u64_ntt_cyclic_forward(b_hat, b_arr, c2);
    for (u64 i = 0; i < nk; i++)
      a_hat[i] = u64_mont_mul(a_hat[i], b_hat[i], &c2->mctx);
    u64_ntt_cyclic_inverse_montgomery_in(c_arr, a_hat, c2);
    for (u64 i = 0; i < nk; i++)
      if (c_arr[i] != ref[i]) p2_fails++;
  }

  /* invalid k rejected, k = 54 accepted */
  {
    u64_ntt_ctx_t tmp;
    if (u64_ntt_ctx_init_proth2(&tmp, 55) != RABIN_ERR_INVALID_ARG) p2_fails++;
    if (u64_ntt_ctx_init_proth2(NULL, 4) != RABIN_ERR_NULL_PTR) p2_fails++;
    if (u64_ntt_ctx_proth2_cached(55) != NULL) p2_fails++;
    if (u64_ntt_ctx_init_proth2(&tmp, 13) != RABIN_SUCCESS)
      p2_fails++;
    else
      u64_ntt_ctx_clear(&tmp);
  }

  if (p2_fails != 0) {
    printf("PROTH2 (q2 = 69*2^55+1) NTT: %d FAILURES\n", p2_fails);
    return 1;
  }
  printf("PROTH2 (q2 = 69*2^55+1) NTT: all checks passed\n");

  /* ---- Goldilocks u64 NTT: asserted checks ---- */
  const u64 gl_p = 0xFFFFFFFF00000001ULL;
  int gl_fails = 0;

  /* fold reduction + multiplication vs u128 modulo oracle */
  static const u64 bvals[] = {0,
                              1,
                              2,
                              0xFFFFFFFFULL,
                              0x100000000ULL,
                              0x100000001ULL,
                              1ULL << 32,
                              1ULL << 63,
                              0xFFFFFFFEFFFFFFFFULL,
                              0xFFFFFFFF00000002ULL % 0xFFFFFFFF00000001ULL,
                              0x8000000000000001ULL};
  for (size_t i = 0; i < sizeof bvals / sizeof bvals[0]; i++) {
    for (size_t j = 0; j < sizeof bvals / sizeof bvals[0]; j++) {
      u64 x = bvals[i], y = bvals[j];
      if (x >= gl_p || y >= gl_p) continue;
      u128 prod = (u128)x * y;
      if (u64_mul_goldilocks(x, y) != (u64)(prod % gl_p)) gl_fails++;
      if (u64_goldilock_red(prod) != (u64)(prod % gl_p)) gl_fails++;
    }
  }
  /* products at the fold boundaries: hi word all-ones, lo word all-ones */
  {
    u128 edges[] = {0,
                    ~0ULL,
                    (u128)~0ULL << 64,
                    ~((u128)0),
                    ((u128)gl_p << 64) | 7,
                    (u128)(gl_p - 1) * (gl_p - 1)};
    for (size_t i = 0; i < sizeof edges / sizeof edges[0]; i++)
      if (u64_goldilock_red(edges[i]) != (u64)(edges[i] % (u128)gl_p))
        gl_fails++;
  }

  /* context: root of unity orders, canonical tables, cached identity */
  for (u64 k = 1; k <= 12; k++) {
    u64_ntt_ctx_t* gl_ctx = u64_ntt_ctx_goldilocks_cached(k);
    if (gl_ctx == NULL || u64_ntt_ctx_goldilocks_cached(k) != gl_ctx) {
      gl_fails++;
      continue;
    }
    u64 n = 1ULL << k;
    if (gl_ctx->n != n || gl_ctx->q != gl_p || gl_ctx->bit_rev_indices != NULL)
      gl_fails++;
    u64 om = gl_ctx->omega_powers[1];
    /* omega^n = 1 and omega^(n/2) = p-1 (primitive n-th root) */
    if (u64_mod_pow(om, n, gl_p) != 1 ||
        u64_mod_pow(om, n >> 1, gl_p) != gl_p - 1)
      gl_fails++;
    /* n * n_inv = 1 mod p */
    if (u64_mul_goldilocks(n % gl_p, gl_ctx->n_inv) != 1) gl_fails++;
    /* omega * omega_inv = 1, tables canonical */
    if (u64_mul_goldilocks(om, gl_ctx->omega_inv_powers[1]) != 1) gl_fails++;
    for (u64 i = 0; i < n; i++)
      if (gl_ctx->omega_powers[i] >= gl_p ||
          gl_ctx->omega_inv_powers[i] >= gl_p)
        gl_fails++;
  }

  /* round trip fwd -> inv on random data for several sizes */
  for (u64 k = 1; k <= 12; k++) {
    u64_ntt_ctx_t* gl_ctx = u64_ntt_ctx_goldilocks_cached(k);
    u64 n = 1ULL << k;
    u64* src = malloc(sizeof(u64) * n);
    u64* mid = malloc(sizeof(u64) * n);
    u64* dst = malloc(sizeof(u64) * n);
    if (src == NULL || mid == NULL || dst == NULL) {
      gl_fails++;
    } else {
      for (u64 i = 0; i < n; i++) src[i] = ((u64)rand() << 32) ^ (u64)rand();
      /* keep inputs reduced; also exercise the p-1 edge */
      for (u64 i = 0; i < n; i++) src[i] = (i == 0) ? gl_p - 1 : src[i] % gl_p;
      u64_ntt_goldilocks_forward(mid, src, gl_ctx);
      for (u64 i = 0; i < n; i++)
        if (mid[i] >= gl_p) gl_fails++;
      u64_ntt_goldilocks_inverse(dst, mid, gl_ctx);
      for (u64 i = 0; i < n; i++)
        if (dst[i] != src[i]) gl_fails++;
    }
    free(src);
    free(mid);
    free(dst);
  }

  /* cyclic convolution vs naive O(n^2) over the field (n = 16) */
  {
    u64 k = 4, n = 1ULL << k;
    u64_ntt_ctx_t* gl_ctx = u64_ntt_ctx_goldilocks_cached(k);
    u64 a_arr[16], b_arr[16], a_hat[16], b_hat[16], c_arr[16], ref[16];
    for (u64 i = 0; i < n; i++) {
      a_arr[i] = ((u64)rand() << 32) ^ (u64)rand();
      b_arr[i] = ((u64)rand() << 32) ^ (u64)rand();
      a_arr[i] = (i == 0) ? gl_p - 1 : a_arr[i] % gl_p;
      b_arr[i] = (i == n - 1) ? gl_p - 1 : b_arr[i] % gl_p;
      ref[i] = 0;
    }
    for (u64 i = 0; i < n; i++)
      for (u64 j = 0; j < n; j++) {
        u64 idx = i + j;
        if (idx >= n) idx -= n;
        ref[idx] =
            u64_mod_add(ref[idx], u64_mul_goldilocks(a_arr[i], b_arr[j]), gl_p);
      }
    u64_ntt_goldilocks_forward(a_hat, a_arr, gl_ctx);
    u64_ntt_goldilocks_forward(b_hat, b_arr, gl_ctx);
    for (u64 i = 0; i < n; i++)
      a_hat[i] = u64_mul_goldilocks(a_hat[i], b_hat[i]);
    u64_ntt_goldilocks_inverse(c_arr, a_hat, gl_ctx);
    for (u64 i = 0; i < n; i++)
      if (c_arr[i] != ref[i]) gl_fails++;
  }

  /* invalid k rejected, k = 31 accepted */
  {
    u64_ntt_ctx_t tmp;
    if (u64_ntt_ctx_init_goldilocks(&tmp, 32) != RABIN_ERR_INVALID_ARG)
      gl_fails++;
    if (u64_ntt_ctx_init_goldilocks(NULL, 4) != RABIN_ERR_NULL_PTR) gl_fails++;
    if (u64_ntt_ctx_goldilocks_cached(32) != NULL) gl_fails++;
  }

  if (gl_fails != 0) {
    printf("GOLDILOCKS NTT: %d FAILURES\n", gl_fails);
    return 1;
  }
  printf("GOLDILOCKS NTT: all checks passed\n");

  rz_clear_constants();

  return 0;
}
