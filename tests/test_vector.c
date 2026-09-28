#include "../include/rvec.h"

int main()
{
  rz_init_constants();

  rvec_t a = {0};

  rvec_init_dynamic(&a);

  rvec_println(&a);

  rvec_append(&a, &RZ_ONE);
  rvec_append(&a, &RZ_ONE);
  rvec_append(&a, &RZ_ONE);
  rvec_append(&a, &RZ_ONE);
  rvec_append(&a, &RZ_ONE);
  rvec_append(&a, &RZ_ONE);
  rvec_append(&a, &RZ_ONE);
  rvec_append(&a, &RZ_ONE);

  rvec_println(&a);

  rvec_neg(&a, &a);
  rvec_println(&a);

  rvec_clear(&a);

  rz_clear_constants();
}
