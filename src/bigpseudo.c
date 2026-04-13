#include <fcntl.h>
#include <stdarg.h>
#include <stddef.h>
#include <unistd.h>

#include "../include/bignum.h"

// returns a k-bit strong pseudoprime to base 2
bool bn_gen_strps(bignum* p, u64 k)
{
  bignum two;
  bn_init(&two);
  bn_set_u64(&two, 2);
  int fd = open("/dev/urandom", O_RDONLY);
  if (fd < 0) return false;

  do {
    bn_gen_random_with_fd(p, k, fd);

  } while (bn_bpsw(p) || !bn_rabin_mont(p, &two));

  close(fd);
  bn_free(&two);
  return true;
}
