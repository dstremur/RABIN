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

void free_pocklington_cert(pocklington_cert* cert)
{
  if (!cert) return;

  bn_free(&cert->N);

  if (cert->data) {
    for (u64 i = 0; i < cert->size; i++) {
      bn_free(&cert->data[i].q);
      bn_free(&cert->data[i].alpha_q);
      free_pocklington_cert(cert->data[i].q_cert);
    }
    free(cert->data);
  }
  free(cert);
}
