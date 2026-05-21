// main.c
#include <assert.h>
#include <fcntl.h>
#include <stdint.h>
#include <stdio.h>
#include <unistd.h>

#include "../include/bignum.h"
#include "../include/bigntt.h"
#include "../include/bigmatrix.h" 
#include "ctype.h"
#include "string.h"

typedef struct {
  const char* n_hex;
  const char* base_hex;
  bool expected;
  const char* description;
} test_case;

void bigntt_ctx_print_debug(const ntt_ctx *ctx)
{
    if (!ctx) {
        printf("Context is NULL\n");
        return;
    }

    printf("=========================================\n");
    printf("         NTT CONTEXT DEBUG LOG           \n");
    printf("=========================================\n");
    printf("Transform Length (n) : %llu\n", ctx->n);
    printf("Modulus (q)          : "); bn_println(&ctx->q);
    printf("N Inverse (n_inv)    : "); bn_println(&ctx->n_inv);
    
    printf("\n--- OMEGA POWERS (w^i mod q) ---\n");
    printf("omega[0]       (Expected: 1)   : "); bn_println(&ctx->omega_powers[0]);
    printf("omega[1]       (Base Root)     : "); bn_println(&ctx->omega_powers[1]);
    printf("omega[n/2]     (Expected: q-1) : "); bn_println(&ctx->omega_powers[ctx->n / 2]);
    printf("omega[n-1]                     : "); bn_println(&ctx->omega_powers[ctx->n - 1]);

    printf("\n--- OMEGA INVERSE POWERS (w^-i mod q) ---\n");
    printf("omega_inv[0]   (Expected: 1)   : "); bn_println(&ctx->omega_inv_powers[0]);
    printf("omega_inv[1]   (Base Inv Root) : "); bn_println(&ctx->omega_inv_powers[1]);
    printf("omega_inv[n/2] (Expected: q-1) : "); bn_println(&ctx->omega_inv_powers[ctx->n / 2]);

    if (ctx->psi_powers != NULL) {
        printf("\n--- PSI POWERS (psi^i mod q) ---\n");
        printf("psi[0]         (Expected: 1)   : "); bn_println(&ctx->psi_powers[0]);
        printf("psi[1]         (Base Psi Root) : "); bn_println(&ctx->psi_powers[1]);
        printf("psi[n/2]       (Square rt of -1): "); bn_println(&ctx->psi_powers[ctx->n / 2]);
    }

    printf("\n--- BIT REVERSAL INDICES (First 8 Samples) ---\n");
    u64 max_samples = (ctx->n < 8) ? ctx->n : 8;
    for (u64 i = 0; i < max_samples; i++) {
        printf("  Index [%2llu] ---> Bit-Reversed Index [%2llu]\n", i, ctx->bit_rev_indices[i]);
    }
    printf("=========================================\n");
}

int main()
{
  bn_init_constants();

  //gen_rns_primes(500);
  bignum a, b, c, psi, omega;
  bn_init_multi(&a, &b, &c, &psi, &omega, NULL);
  char buf1[1024];
  char buf2[1024];

  printf("Enter first big number (a): ");
  if (!fgets(buf1, sizeof(buf1), stdin)) return 1;
  buf1[strcspn(buf1, "\n")] = 0;  // Remove newline

  printf("Enter second big number (b, p): ");
  if (!fgets(buf2, sizeof(buf2), stdin)) return 1;
  buf2[strcspn(buf2, "\n")] = 0;  // Remove newline

  bn_init_val(&a, buf1);
  bn_init_val(&b, buf2);
  printf("Inputs: \n");
  printf("a: ");
  bn_println(&a);
  printf("b: ");
  bn_println(&b);
  printf("ADD: ");
  bn_add(&c, &a, &b);
  bn_println(&c);
  bn_sub(&c, &a, &b);
  printf("SUB: ");
  bn_println(&c);
  bn_isqrt(&c, &a);
  printf("ISQRT (a): ");
  bn_println(&c);
  bn_ln(&c, &a);
  printf("LN (a): "); 
  bn_println(&c); 
  bn_mul(&c, &a, &b);
  printf("MUL: ");
  bn_println(&c);
  bn_div(&c, &a, &b);
  printf("DIV: ");
  bn_println(&c);
  bn_mod(&c, &a, &b);
  printf("MOD: ");
  bn_println(&c);
  printf("RSHIFT: ");
  bn_rshift(&c, &a, 64);
  bn_println(&c);
  printf("LSHIFT: ");
  bn_lshift(&c, &a, 32);
  bn_println(&c);
  printf("Pow: ");
  bn_pow(&c, &a, &b);
  bn_println(&c);

  bigpoly_test();

  bn_gen_proth_ntt(&c, &psi, &c, &c, 20, 10000);
  printf("proth prime: ");
  bn_println(&psi);
  printf("proth generator: ");
  bn_println(&c);

  printf("factorize a: ");
  bigvector v;
  bigvector_init_dynamic(&v);

  bn_factorize(&v, &a);

  bigvector_println(&v);
  bigvector_free(&v);

  bignum n1;
  bn_init(&n1);

  bn_init_val(&n1, "18446744069414584321");
  bn_init_val(&omega, "1803076106186727246");
  bn_init_val(&psi, "11353340290879379826");
  ntt_ctx ctx;

  //bigntt_ctx_init(&ctx, 512, &n1, &omega, &psi);
  bigntt_ctx_init_simple(&ctx, 25, 10001);
  bigntt_ctx_print_debug(&ctx);

  bn_free(&n1);

  bigntt_ctx_free(&ctx);

  if (bn_pollard_rho(&c, &a)) {
    printf("Factor found: ");
    bn_println(&c);
  } else {
    printf("No non-trivial factor found (Number might be prime).\n");
  }

  if (bn_pollard_p_minus_one(&c, &a)){
	  printf("Factor found p - 1: ");
	  bn_println(&c);
  } else {
	  printf("p - 1 fail \n" );
  }



  bignum n, d;
  bn_init_multi(&n, &d, NULL);
  /*
  u64 n_1 = 10;
    for (u64 i = 0; i < n_1; i++) {
      bn_set_u64(&n, i);
      bn_lucas(&c, &d, &a, &b, &n);
      printf("Lucas\n");
      bn_println(&c);
      bn_println(&d);
    }
  */
  if (bn_bpsw(&a)) {
    printf("a is Prime \n");
  } else {
    printf("a is not prime to base b \n");
  }

  i64 j = bn_jacobi(&a, &b);

  printf("JACOBI: %lli\n", (long long)j);

  bignum r;
  bn_init(&r);

  tonelli_shanks(&r, &a, &b);

  bn_print(&r);
  printf("^2 = ");
  bn_print(&a);
  printf(" mod ");
  bn_println(&b);
  char* s = bn_to_string(&b);
  printf("Result: %s\n", s);
  free(s);

  bignum my_prime;
  bn_init(&my_prime);
  printf("Generating a 2048-bit prime... (this may take a minute)\n");

  if (bn_gen_prime(&my_prime, 2048)) {
    printf("\nSuccess! Your 2048-bit prime is:\n");
    printf(
        "----------------------------------------------------------------\n");
    bn_print(&my_prime);
    printf(
        "----------------------------------------------------------------\n");
  } else {
    fprintf(stderr, "Failed to generate prime or read from /dev/urandom\n");
  }

//  printf("Generating 512-bit provable prime...\n");
 //// bn_provable_prime(&my_prime, 512);

 // printf("Result: ");
 // bn_println(&my_prime);

  // avxtest();

  // 
  //
  bigpoly_test();
  //
  bigpoly p;
  bigpoly_init(&p);

  bn_free(&my_prime);
  bn_free_multi(&a, &b, &c, &n, &d, &r, &omega, &psi, NULL);
  bn_free_constants();
  return 0;
}
