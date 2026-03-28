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

int main() {
  bignum a, b, c;
  bn_init(&a);
  bn_init(&b);
  bn_init(&c);
  bn_init_val(&a, "36346346346839476893763473469837689734896739846794353486798347689347683476347698347693749683458476384");
  bn_init_val(&b, "4000000");
  bn_print(&a);
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

  return 0;
}

