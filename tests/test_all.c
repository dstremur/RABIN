// tests all available functions

#include <stdio.h>

#include "../include/biglogic.h"
#include "../include/bignum.h"
#include "../include/bigvector.h"

int main()
{
  bignum a, b, res;
  bn_init_multi(&a, &b, &res, NULL);

  bn_init_val(&a,
              "2142409230825209382098293582908583529035932854353453453463463464"
              "56456456456456464256324636346346346");
  bn_init_val(&b, "2344");

  printf("A: ");
  bn_println(&a);
  printf("B: ");
  bn_println(&b);

  bn_add(&res, &a, &b);
  printf("ADD: ");
  bn_println(&res);

  bn_sub(&res, &a, &b);
  printf("SUB: ");
  bn_println(&res);

  bn_mul(&res, &a, &b);
  printf("MUL: ");
  bn_println(&res);

  bn_div(&res, &a, &b);
  printf("DIV: ");
  bn_println(&res);

  bn_and(&res, &a, &b);
  printf("a & b: ");
  bn_println(&res);

  bn_or(&res, &a, &b);
  printf("a or b: ");
  bn_println(&res);

  bn_xor(&res, &a, &b);
  printf("a xor b: ");
  bn_println(&res);

  bn_mod(&res, &a, &b);
  printf("MOD: ");
  bn_println(&res);

  // bn_pow(&res, &a, &b);
  printf("POW: ");
  bn_println(&res);

  bn_mod_exp(&res, &a, &b, &b);
  printf("mod b POW: ");
  bn_println(&res);

  bn_isqrt(&res, &a);
  printf("sqrt of a: ");
  bn_println(&res);

  printf("Prime BPSW 2048 bits: ");
  bn_gen_prime(&res, 2048);
  bn_println(&res);

  printf("Provable prime 100 bits: ");
  bn_provable_prime(&res, 512);
  bn_println(&res);

  bn_pollard_rho(&res, &a);
  printf("Factor (pollard_rho) of a: ");
  bn_println(&res);

  printf("JACOBI (b / a): %i \n", bn_jacobi(&b, &a));

  tonelli_shanks(&res, &b, &a);

  bigvector v;
  bigvector_init(&v, 5);

  bigvector_set(&v, &a, 0);
  bigvector_set(&v, &b, 1);
  bigvector_set(&v, &a, 2);
  bigvector_set(&v, &b, 3);
  bigvector_set(&v, &b, 4);

  bigvector_println(&v);

  bigvector_dot(&res, &v, &v);

  bn_println(&res);

  bigvector_norm(&res, &v);

  bn_println(&res);

  bigvector_free(&v);
  bn_free_multi(&a, &b, &res);
}
