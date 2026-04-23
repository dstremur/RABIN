#include <ctype.h>
#include <fcntl.h>
#include <math.h>
#include <stdarg.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "../../include/bignum.h"

// generates a random bits long odd number
bool bn_gen_random(bignum* r, u64 bits)
{
  int fd = open("/dev/urandom", O_RDONLY);
  if (fd < 0) return false;

  bn_gen_random_with_fd(r, bits, fd);
  close(fd);

  return true;
}

// generate random bits long odd number with given dev/urandom
bool bn_gen_random_with_fd(bignum* r, u64 bits, int fd)
{
  int limbs_needed = (bits + 63) / 64;

  if (!bn_alloc(r, limbs_needed)) return false;
  r->size = limbs_needed;

  // Read random bytes directly using the open file descriptor
  if (read(fd, r->limbs, limbs_needed * sizeof(u64)) !=
      (ssize_t)(limbs_needed * sizeof(u64))) {
    return false;
  }

  // Mask the top limb to fit the exact bit length
  int top_bits = bits % 64;
  if (top_bits != 0) {
    u64 mask = ((u64)1 << top_bits) - 1;
    r->limbs[r->size - 1] &= mask;
  }

  // Ensure it's exactly 'bits' long by setting the MSB
  r->limbs[r->size - 1] |= ((u64)1 << ((bits - 1) % 64));

  // Ensure it's odd by setting the LSB
  r->limbs[0] |= 1;

  return true;
}

// Helper to generate a random number in range [low, high]
void bn_gen_random_range(bignum* r, const bignum* low, const bignum* high)
{
  bignum range;
  bn_init(&range);
  bn_sub(&range, high, low);

  if (bn_is_zero(&range)) {
    bn_copy(r, low);
    bn_free(&range);
    return;
  }

  int bits = bn_bit_length(&range);

  do {
    bn_gen_random(r, bits);
  } while (bn_cmp(r, &range) > 0);

  bn_add(r, r, low);
  bn_free(&range);
}
