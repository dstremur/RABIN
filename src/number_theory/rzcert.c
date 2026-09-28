#include "../../include/rzcert.h"

#include "stdio.h"

rabin_err_t print_pocklington_cert(rzcert_t* cert)
{
  if (cert == NULL) return RABIN_ERR_NULL_PTR;

  printf("Prime N: ");
  rz_print(&cert->N);
  printf("\n");

  if (cert->size == 0) {
    printf("-> Base case (verified via bpsw)\n");
    return RABIN_SUCCESS;
  }

  rabin_err_t err = RABIN_SUCCESS;
  for (u64 i = 0; i < cert->size && err == RABIN_SUCCESS; i++) {
    printf("-> Factor q: ");
    rz_print(&cert->data[i].q);
    printf(" | Base alpha: ");
    rz_print(&cert->data[i].alpha_q);
    printf("\n");

    // Recursively print the certificate for q
    err = print_pocklington_cert(cert->data[i].q_cert);
  }

  return err;
}

rabin_err_t rzcert_clear(rzcert_t* cert)
{
  if (cert == NULL) return RABIN_ERR_NULL_PTR;

  rz_clear(&cert->N);

  if (cert->data) {
    for (u64 i = 0; i < cert->size; i++) {
      rz_clear(&cert->data[i].q);
      rz_clear(&cert->data[i].alpha_q);
      rzcert_clear(cert->data[i].q_cert);
    }
    free(cert->data);
  }
  free(cert);
  return RABIN_SUCCESS;
}
