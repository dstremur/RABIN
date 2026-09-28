// tests all available functions

#include <stdio.h>

#include "../include/rabin.h"
#include "../include/rvec.h"
#include "../include/rzlogic.h"

int main()
{
  rz_t a, b, res;
  rz_init_multi(&a, &b, &res, NULL);

  rz_init_val(&a,
              "2142409230825209382098293582908583529035932854353453453463463464"
              "56456456456456464256324636346346346");
  rz_init_val(&b, "2344");

  printf("A: ");
  rz_println(&a);
  printf("B: ");
  rz_println(&b);

  rz_add(&res, &a, &b);
  printf("ADD: ");
  rz_println(&res);

  rz_sub(&res, &a, &b);
  printf("SUB: ");
  rz_println(&res);

  rz_mul(&res, &a, &b);
  printf("MUL: ");
  rz_println(&res);

  rz_div(&res, &a, &b);
  printf("DIV: ");
  rz_println(&res);

  rz_and(&res, &a, &b);
  printf("a & b: ");
  rz_println(&res);

  rz_or(&res, &a, &b);
  printf("a or b: ");
  rz_println(&res);

  rz_xor(&res, &a, &b);
  printf("a xor b: ");
  rz_println(&res);

  rz_mod(&res, &a, &b);
  printf("MOD: ");
  rz_println(&res);

  // rz_pow(&res, &a, &b);
  printf("POW: ");
  rz_println(&res);

  rz_mod_exp(&res, &a, &b, &b);
  printf("mod b POW: ");
  rz_println(&res);

  rz_isqrt(&res, &a);
  printf("sqrt of a: ");
  rz_println(&res);

  printf("Prime BPSW 2048 bits: ");
  rz_gen_prime(&res, 2048);
  rz_println(&res);

  printf("Provable prime 100 bits: ");
  rz_provable_prime(&res, 512);
  rz_println(&res);

  rz_pollard_rho(&res, &a);
  printf("Factor (pollard_rho) of a: ");
  rz_println(&res);

  printf("JACOBI (b / a): %i \n", rz_jacobi(&b, &a));

  rz_tonelli_shanks(&res, &b, &a);

  rvec_t v = {0};
  rvec_init(&v, 5);

  rvec_set(&v, &a, 0);
  rvec_set(&v, &b, 1);
  rvec_set(&v, &a, 2);
  rvec_set(&v, &b, 3);
  rvec_set(&v, &b, 4);

  rvec_println(&v);

  rvec_dot(&res, &v, &v);

  rz_println(&res);

  rvec_norm(&res, &v);

  rz_println(&res);

  rvec_clear(&v);
  rz_clear_multi(&a, &b, &res);
}
