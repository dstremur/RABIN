#include "../include/bigvector.h"

int main()
{
  bn_init_constants();

  bigvector a;

  bigvector_init_dynamic(&a);

  bigvector_println(&a);

  bigvector_append(&a, &BN_ONE);
  bigvector_append(&a, &BN_ONE);
  bigvector_append(&a, &BN_ONE);
  bigvector_append(&a, &BN_ONE);
  bigvector_append(&a, &BN_ONE);
  bigvector_append(&a, &BN_ONE);
  bigvector_append(&a, &BN_ONE);
  bigvector_append(&a, &BN_ONE);

  bigvector_println(&a);

  bigvector_neg(&a, &a);
  bigvector_println(&a);

  bigvector_free(&a);

  bn_free_constants();
}
