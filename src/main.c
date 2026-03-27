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
  bn_init_val(&a, "34633745745745784686458568548456456486");
  bn_init_val(&b, "476486548548548785769569569799569");
  bn_div(&c, &a, &b);
  bn_print(&c);

  return 0;
}

