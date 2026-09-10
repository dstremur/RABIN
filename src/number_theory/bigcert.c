#include "../../include/bigcert.h"

#include "stdio.h"

void print_pocklington_cert(pocklington_cert* cert)
{
  if (!cert) return;

  printf("Prime N: ");
  bn_print(&cert->N);
  printf("\n");

  if (cert->size == 0) {
    printf("-> Base case (verified via bpsw)\n");
    return;
  }

  for (u64 i = 0; i < cert->size; i++) {
    printf("-> Factor q: ");
    bn_print(&cert->data[i].q);
    printf(" | Base alpha: ");
    bn_print(&cert->data[i].alpha_q);
    printf("\n");

    // Recursively print the certificate for q
    print_pocklington_cert(cert->data[i].q_cert);
  }
}

void pocklington_cert_free(pocklington_cert* cert)
{
  if (!cert) return;

  bn_free(&cert->N);

  if (cert->data) {
    for (u64 i = 0; i < cert->size; i++) {
      bn_free(&cert->data[i].q);
      bn_free(&cert->data[i].alpha_q);
      pocklington_cert_free(cert->data[i].q_cert);
    }
    free(cert->data);
  }
  free(cert);
}

// bool pocklington_cert_verify(const pocklington_cert* cert, u64
// base_case_bound)
// {
//   // 1. a^(N - 1) = 1 mod N

//   // 2. p | N - 1 & p > sqrt(N) - 1

//   // 3. gcd(a^((N - 1) / p), N) = 1
// }

// bool pocklington_cert_verify2(const pocklington_cert* cert, u64
// base_case_bound)
// {
//   // 1. BASE CASE: If N is small enough, verify without a certificate
//   if (bignum_cmp_u64(&cert->N, base_case_bound) <= 0) {
//     return verify_base_case(&cert->N);  // Implement trial division here
//   }

//   bignum F, N_minus_1;
//   bignum_init_set_u64(&F, 1);
//   bignum_init(&N_minus_1);
//   bignum_sub_u64(&N_minus_1, &cert->N, 1);

//   // Verify each prime factor in the certificate
//   for (u64 i = 0; i < cert->size; i++) {
//     pocklington_cert_elem* elem = &cert->data[i];

//     // 2. RECURSIVE VERIFICATION: Prove 'q' is actually prime
//     if (elem->q_cert != NULL) {
//       // Verify the child certificate
//       if (!pocklington_cert_verify(elem->q_cert, base_case_bound)) return
//       false;
//       // Ensure the child certificate is actually for 'q'
//       if (bignum_cmp(&elem->q_cert->N, &elem->q) != 0) return false;
//     } else {
//       // If leaf node, it must be <= base_case_bound
//       if (bignum_cmp_u64(&elem->q, base_case_bound) > 0) return false;
//       if (!verify_base_case(&elem->q)) return false;
//     }

//     // 3. ACCUMULATE F: Find highest power of q dividing N-1
//     bignum q_power;
//     get_highest_power_dividing(&q_power, &N_minus_1, &elem->q);
//     bignum_mul(&F, &F, &q_power);  // F = F * q^e

//     // 4. POCKLINGTON CONDITIONS
//     bignum fermat_res, exp, gcd_res, gcd_input;

//     // A. Fermat condition: alpha_q^(N-1) == 1 mod N
//     bignum_mod_exp(&fermat_res, &elem->alpha_q, &N_minus_1, &cert->N);
//     if (bignum_cmp_u64(&fermat_res, 1) != 0) return false;

//     // B. GCD condition: gcd(alpha_q^((N-1)/q) - 1, N) == 1
//     bignum_div(&exp, &N_minus_1, &elem->q);  // exp = (N-1) / q
//     bignum_mod_exp(&gcd_input, &elem->alpha_q, &exp, &cert->N);
//     bignum_sub_u64(&gcd_input, &gcd_input, 1);  // alpha_q^((N-1)/q) - 1

//     bignum_gcd(&gcd_res, &gcd_input, &cert->N);
//     if (bignum_cmp_u64(&gcd_res, 1) != 0) return false;
//   }

//   // 5. FINAL BOUND CHECK: F > sqrt(N) - 1
//   // Mathematically equivalent to: (F + 1)^2 > N
//   bignum F_plus_1, F_squared;
//   bignum_add_u64(&F_plus_1, &F, 1);
//   bignum_mul(&F_squared, &F_plus_1, &F_plus_1);

//   if (bignum_cmp(&F_squared, &cert->N) <= 0) {
//     return false;  // Factored part isn't large enough
//   }

//   return true;  // All conditions met, N is prime!
// }