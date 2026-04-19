
#include <stdio.h>

#include "../../include/bignum.h"

void bn_ntt_init(bn_ntt* cfg, u64 n, const char* prime)
{
  cfg->n = n;
  bn_init(&cfg->q);
  bn_init_val(&cfg->q, prime);
  bn_mont_ctx_init(&cfg->ctx, &cfg->q);

  cfg->psi_pow = malloc(n * sizeof(bignum));
  cfg->inv_psi_pow = malloc(n * sizeof(bignum));

  bignum psi, psi_inv, temp, temp2, g, phi;
  bn_init_multi(&psi, &psi_inv, &temp, &temp2, &g, &phi, NULL);

  // 1. Calculate phi = (q-1) / 2n
  bn_set_u64(&temp, 1);
  bn_sub(&phi, &cfg->q, &temp);  // phi = q - 1
  bn_set_u64(&temp, 2 * n);
  bn_div(&phi, &phi, &temp);  // phi = (q-1)/2n

  // 2. Find psi (a primitive 2n-th root of unity)
  u64 g_val = 2;
  while (true) {
    bn_set_u64(&g, g_val++);
    bn_mod_exp(&psi, &g, &phi, &cfg->q);

    // For negacyclic, psi^n must be -1 mod q, NOT 1 mod q.
    bn_set_u64(&temp, n);
    bn_mod_exp(&temp, &psi, &temp, &cfg->q);
    if (!bn_is_eq_i64(&temp, 1)) break;  // Found a valid 2n-th root
  }

  // 3. Calculate psi_inv = psi^(q-2) mod q
  bn_set_u64(&temp, 2);
  bn_sub(&temp, &cfg->q, &temp);  // temp = q - 2
  bn_mod_exp(&psi_inv, &psi, &temp, &cfg->q);

  // 4. Precompute powers
  for (u64 i = 0; i < n; i++) {
    bn_init(&cfg->psi_pow[i]);
    bn_init(&cfg->inv_psi_pow[i]);

    bn_set_u64(&temp, i);
    bn_mod_exp(&cfg->psi_pow[i], &psi, &temp, &cfg->q);
    bn_mont_in(&cfg->psi_pow[i], &cfg->psi_pow[i], &cfg->ctx);

    bn_mod_exp(&cfg->inv_psi_pow[i], &psi_inv, &temp, &cfg->q);
    bn_mont_in(&cfg->inv_psi_pow[i], &cfg->inv_psi_pow[i], &cfg->ctx);
  }

  // 5. Calculate n_inv = n^(q-2) mod q
  bn_init(&cfg->n_inv);
  bn_set_u64(&temp, n);
  bn_set_u64(&temp2, 2);
  bn_sub(&temp2, &cfg->q, &temp2);  // temp2 = q - 2
  bn_mod_exp(&cfg->n_inv, &temp, &temp2, &cfg->q);
  bn_mont_in(&cfg->n_inv, &cfg->n_inv, &cfg->ctx);

  bn_free_multi(&g, &psi, &psi_inv, &temp, &temp2, &phi, NULL);
}
// write in bit-reversed order
void ntt_bit_reverse(bignum* a, u64 n)
{
  for (u64 i = 1, j = 0; i < n; i++) {
    u64 bit = n >> 1;
    for (; j & bit; bit >>= 1) j ^= bit;
    j ^= bit;
    if (i < j) {
      // Simple swap of bignum pointers or values
      bignum temp;
      temp = a[i];
      a[i] = a[j];
      a[j] = temp;
    }
  }
}

void ntt_forward(bignum* a, bn_ntt* cfg)
{
  u64 n = cfg->n;
  bignum t1, t2;
  bn_init_multi(&t1, &t2, NULL);

  // 1. Pre-weighting (Negacyclic Step)
  // Allows to use a standard Cooley-Tukey algo
  for (u64 i = 0; i < n; i++) {
    bn_mont_in(&a[i], &a[i], &cfg->ctx);
    bn_mont_mul(&a[i], &a[i], &cfg->psi_pow[i], &cfg->ctx);
  }

  // 2. Bit-reverse
  ntt_bit_reverse(a, n);

  // 3. Cooley-Tukey NTT
  for (u64 len = 2; len <= n; len <<= 1) {
    u64 half = len >> 1;
    u64 step = n / len;
    for (u64 i = 0; i < n; i += len) {
      for (u64 j = 0; j < half; j++) {
        bignum* u = &a[i + j];
        bignum* v = &a[i + j + half];
        // Omega power: psi^(2 * j * step) Twiddle factor
        bignum* w = &cfg->psi_pow[2 * j * step];

        // Butterfly: T = V * W
        bn_mont_mul(&t1, v, w, &cfg->ctx);

        // V = U - T
        if (bn_cmp(u, &t1) < 0) bn_add(u, u, &cfg->q);  // Handle borrow
        bn_sub(v, u, &t1);

        // U = U + T
        bn_add(u, u, &t1);
        if (bn_cmp(u, &cfg->q) >= 0) bn_sub(u, u, &cfg->q);
      }
    }
  }
  bn_free_multi(&t1, &t2, NULL);
}

void ntt_inverse(bignum* a, bn_ntt* cfg)
{
  u64 n = cfg->n;
  bignum t1, t2;
  bn_init_multi(&t1, &t2, NULL);

  ntt_bit_reverse(a, n);

  for (u64 len = 2; len <= n; len <<= 1) {
    u64 half = len >> 1;
    u64 step = n / len;
    for (u64 i = 0; i < n; i += len) {
      for (u64 j = 0; j < half; j++) {
        bignum* u = &a[i + j];
        bignum* v = &a[i + j + half];
        // Inverse Omega power: psi^(-2 * j * step)
        bignum* w = &cfg->inv_psi_pow[2 * j * step];

        bn_mont_mul(&t1, v, w, &cfg->ctx);

        if (bn_cmp(u, &t1) < 0) bn_add(u, u, &cfg->q);
        bn_sub(v, u, &t1);

        bn_add(u, u, &t1);
        if (bn_cmp(u, &cfg->q) >= 0) bn_sub(u, u, &cfg->q);
      }
    }
  }
  bn_free_multi(&t1, &t2, NULL);
}

void bn_poly_mul_negacyclic(bignum* r, bignum* a, bignum* b, bn_ntt* cfg)
{
  u64 n = cfg->n;

  ntt_forward(a, cfg);
  ntt_forward(b, cfg);

  // Point-wise multiplication in Montgomery space
  for (u64 i = 0; i < n; i++) {
    bn_mont_mul(&r[i], &a[i], &b[i], &cfg->ctx);
  }

  ntt_inverse(r, cfg);

  // Post-processing: Scaling by 1/n and psi^-i
  for (u64 i = 0; i < n; i++) {
    bn_mont_mul(&r[i], &r[i], &cfg->n_inv, &cfg->ctx);
    bn_mont_mul(&r[i], &r[i], &cfg->inv_psi_pow[i], &cfg->ctx);
    bn_mont_out(&r[i], &r[i], &cfg->ctx);  // Return to normal domain
  }
}

void print_poly(const char* name, bignum* p, int n)
{
  printf("%s = [ ", name);
  for (int i = 0; i < n; i++) {
    char* s = bn_to_string(&p[i]);
    printf("%s%s", s, (i == n - 1) ? "" : ", ");
    free(s);
  }
  printf(" ]\n");
}

void nntest()
{
  u64 n = 4;
  const char* q_str = "17";

  // 1. Initialize NTT Configuration
  bn_ntt cfg;
  bn_ntt_init(&cfg, n, q_str);

  // 2. Initialize Input Polynomials
  bignum a[4], b[4], r[4];
  u64 val_a[] = {1, 2, 3, 4};
  u64 val_b[] = {5, 6, 7, 8};

  for (u64 i = 0; i < n; i++) {
    bn_init(&a[i]);
    bn_init(&b[i]);
    bn_init(&r[i]);
    bn_set_u64(&a[i], val_a[i]);
    bn_set_u64(&b[i], val_b[i]);
  }

  printf("--- NTT Negacyclic Convolution Test ---\n");
  printf("Modulus q = 17, n = 4 (Ring: x^4 + 1)\n");
  print_poly("Input A", a, n);
  print_poly("Input B", b, n);

  // 3. Perform Negacyclic Convolution
  // Note: a and b are modified in-place during the transform
  bn_poly_mul_negacyclic(r, a, b, &cfg);

  // 4. Print Result
  print_poly("Result R", r, n);

  // 5. Verification
  u64 expected[] = {12, 15, 2, 9};
  int passed = 1;
  for (u64 i = 0; i < n; i++) {
    if (!bn_is_eq_i64(&r[i], expected[i])) {
      passed = 0;
      printf("Error at index %llu: Expected %lu\n", (unsigned long long)i,
             expected[i]);
    }
  }

  if (passed) {
    printf("\nRESULT: SUCCESS (Matches mathematical expectation)\n");
  } else {
    printf("\nRESULT: FAILED\n");
  }

  // Cleanup
  for (u64 i = 0; i < n; i++) {
    bn_free(&a[i]);
    bn_free(&b[i]);
    bn_free(&r[i]);
  }
}
