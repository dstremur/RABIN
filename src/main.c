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

int main()
{
  bignum a, b, c;
  bn_init(&a);
  bn_init(&b);
  bn_init(&c);
  char buf1[1024];
  char buf2[1024];

  printf("Enter first big number (a): ");
  if (!fgets(buf1, sizeof(buf1), stdin)) return 1;
  buf1[strcspn(buf1, "\n")] = 0;  // Remove newline

  printf("Enter second big number (b): ");
  if (!fgets(buf2, sizeof(buf2), stdin)) return 1;
  buf2[strcspn(buf2, "\n")] = 0;  // Remove newline

  bn_init_val(&a, buf1);
  bn_init_val(&b, buf2);
  printf("Inputs: \n");
  printf("a: ");
  bn_print(&a);
  printf("b: ");
  bn_print(&b);
  printf("ADD: \n");
  bn_add(&c, &a, &b);
  bn_print(&c);
  bn_init(&c);
  bn_sub(&c, &a, &b);
  printf("SUB: \n");
  bn_print(&c);
  bn_init(&c);
  bn_mul(&c, &a, &b);
  printf("MUL: \n");
  bn_print(&c);
  bn_init(&c);
  bn_div(&c, &a, &b);
  printf("DIV: \n");
  bn_print(&c);
  printf("RSHIFT: \n");
  bn_rshift(&c, &a, 64);
  bn_print(&c);
  printf("LSHIFT: \n");
  bn_lshift(&c, &a, 32);
  bn_print(&c);
  printf("Pow: \n");
  bn_pow(&c, &a, &b);
  bn_print(&c);

  i64 j = bn_jacobi(&a, &b);

  printf("JACOBI: %lli\n", (long long)j);

  bignum p, n, r;
  bn_init_multi(&p, &n, &r);

  bn_set_u64(&p, 13);
  bn_set_u64(&n, 13434634);

  tonelli_shanks(&r, &p, &n);

  printf("root = ");
  bn_print(&r);
  printf("\n");

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

  bn_free(&my_prime);
  return 0;
}
