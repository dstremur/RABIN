#include "../../include/bigntt.h"

#include "../../include/u64.h"

static inline bool is_pow_2(u64 n) { return n && !(n & (n - 1)); }

static inline u64 ilog2_u64(u64 x)
{
  u64 r = 0;
  while (x >>= 1) r++;
  return r;
}

static u64 reverse_bits(u64 x, u64 bits)
{
  u64 r = 0;

  for (u64 i = 0; i < bits; i++) {
    r <<= 1;
    r |= (x & 1);
    x >>= 1;
  }

  return r;
}

bool bigntt_find_prime(bignum* q, u64 n, u64 bits)
{
  // n must be a power of 2
  if (!is_pow_2(n)) {
    return false;
  }

  // q needs to be of the form
  // q = k * (2n) + 1

  u64 stride = 2 * n;

  bignum cand;
  bn_init_multi(&cand, NULL);

  // cand = k * 2^(bits - 1) + 1
  // k = ceil(2^(bits-1) / 2n)
  bn_lshift(&cand, &BN_ONE, bits - 1);

  u64 rem = bn_mod_u64(&cand, stride);

  if (rem != 1) {
    u64 delta = (stride + 1 - rem) % stride;
    bn_add_u64(&cand, &cand, delta);
  }
  while (true) {
    if (bn_bpsw(&cand)) {
      bn_copy(q, &cand);
      bn_free(&cand);
      return true;
    }

    bn_add_u64(&cand, &cand, stride);
  }

  bn_free(&cand);
  return false;
}

bool bigntt_find_generator(bignum* g, const bignum* q)
{
  bignum q_min_1, exp, tmp;
  bn_init_multi(&q_min_1, &exp, &tmp, NULL);

  bn_sub(&q_min_1, q, &BN_ONE);

  bool found = false;

  // Loop through potential generator bases
  for (u64 a = 2; a < 10000; a++) {
    bn_set_u64(g, a);

    // Test 1: g^((q-1)/2) != 1 mod q (The most basic check, handles the factor
    // 2)
    bn_divmod_u64(&exp, &q_min_1, 2);
    bn_mod_exp(&tmp, g, &exp, q);
    if (bn_cmp(&tmp, &BN_ONE) == 0) {
      continue;  // Not a generator
    }

    // Test 2: If q-1 has other obvious small prime factors, you'd check them
    // here. However, for the purpose of finding a 2n-th root of unity, if a
    // base passes the check for the 2-power component, it will yield the
    // correct primitive roots in bigntt_compute_roots.

    found = true;
    break;
  }

  bn_free_multi(&q_min_1, &exp, &tmp, NULL);
  return found;
}

bool bigntt_compute_roots(bignum* omega, bignum* psi, u64 n, const bignum* q)
{
  if (!is_pow_2(n)) {
    return false;
  }

  bignum q_min_1, k, exp, g;
  bn_init_multi(&q_min_1, &k, &exp, &g, NULL);

  bn_sub(&q_min_1, q, &BN_ONE);
  bool res = false;

  // Ensure 2n evenly divides q-1
  // (Fixing the logic check: if remainder is NOT 0, it's an error)
  if (bn_mod_u64(&q_min_1, 2 * n) != 0) {
    goto cleanup;
  }

  // k = (q-1) / 2n
  if (!bn_divmod_u64(&k, &q_min_1, 2 * n)) {
    goto cleanup;
  }

  // Step 1: Find a generator g mod q
  if (!bigntt_find_generator(&g, q)) {
    goto cleanup;
  }

  // Step 2: psi = g^k mod q  (where k = (q-1)/2n)
  // By definition of a generator, this guarantees psi is a primitive 2n-th
  // root.
  bn_mod_exp(psi, &g, &k, q);

  // Step 3: omega = psi^2 mod q (primitive n-th root)
  bn_set_u64(&exp, 2);
  bn_mod_exp(omega, psi, &exp, q);

  res = true;

cleanup:
  bn_free_multi(&q_min_1, &k, &exp, &g, NULL);
  return res;
}
