#include "../../include/u64.h"

void ntt_u64_cyclic_forward(u64* a_hat, const u64* a, ntt_ctx_u64* ctx)
{
  u64 n = ctx->n;
  u64 q = ctx->q;

  // 1. Bit-reversal permutation AND entering Montgomery space

  for (u64 i = 0; i < n; i++) {
    u64 rev = ctx->bit_rev_indices[i];
    // Convert to Montgomery form right as we load the data
    a_hat[rev] = mont_in(a[i], &ctx->mctx);
  }

  // 2. Cooley-Tukey Butterfly
  for (u64 len = 2; len <= n; len <<= 1) {
    u64 half = len >> 1;
    u64 step = n / len;

    for (u64 i = 0; i < n; i += len) {
      for (u64 j = 0; j < half; j++) {
        u64 twiddle =
            ctx->omega_powers[j * step];  // Already in Montgomery form

        u64 even = i + j;
        u64 odd = i + j + half;

        // t = a_hat[odd] * twiddle (Montgomery multiplication)
        u64 t = mont_mul(a_hat[odd], twiddle, &ctx->mctx);
        u64 u = a_hat[even];

        // Montgomery form preserves addition and subtraction natively
        a_hat[even] = mod_add(u, t, q);
        a_hat[odd] = mod_sub(u, t, q);
      }
    }
  }
}

bool ntt_ctx_u64_init_golden(ntt_ctx_u64* ctx, u64 k)
{
  if (k > 54) return false;

  // p = 5 * 2^55 + 1
  u64 p = 180143985094819841ULL;
  u64 g = 3;

  // Calculate the exponent for psi: c * 2^(55 - (k+1))
  // We use 5ULL to ensure it shifts as a 64-bit integer
  u64 exp = 5ULL << (55 - (k + 1));

  // psi = g^exp mod p
  u64 psi = mod_pow(g, exp, p);

  // omega = psi^2 mod p
  u64 omega = mod_mul(psi, psi, p);

  // Initialize the Montgomery NTT context
  return ntt_ctx_u64_init(ctx, p, k, omega, psi);
}

bool ntt_ctx_u64_init(ntt_ctx_u64* ctx, u64 p, u64 k, u64 omega, u64 psi)
{
  u64 n = (1ULL << k);
  ctx->n = n;
  ctx->k = k;
  ctx->q = p;

  mont_init(&ctx->mctx, p);

  ctx->omega_powers = malloc(sizeof(u64) * n);
  ctx->omega_inv_powers = malloc(sizeof(u64) * n);
  ctx->bit_rev_indices = malloc(sizeof(u64) * n);

  // Calculate n^-1 mod q and put it in Montgomery form
  u64 n_inv_standard = mod_inverse_euclid(n, p);
  ctx->n_inv = mont_in(n_inv_standard, &ctx->mctx);

  // Precompute powers and move to Montgomery form immediately
  u64 current_omega = 1;
  u64 omega_inv = mod_inverse_euclid(omega, p);
  u64 current_omega_inv = 1;

  for (u64 i = 0; i < n; i++) {
    ctx->omega_powers[i] = mont_in(current_omega, &ctx->mctx);
    ctx->omega_inv_powers[i] = mont_in(current_omega_inv, &ctx->mctx);

    current_omega =
        mod_mul(current_omega, omega, p);  // Standard mod_mul for precalc
    current_omega_inv = mod_mul(current_omega_inv, omega_inv, p);
  }

  // Precompute Bit-Reversal
  u64 bits = 0;
  while (((u64)1 << bits) < n) bits++;
  for (u64 i = 0; i < n; i++) {
    u64 rev = 0, temp = i;
    for (u64 j = 0; j < bits; j++) {
      rev = (rev << 1) | (temp & 1);
      temp >>= 1;
    }
    ctx->bit_rev_indices[i] = rev;
  }

  return true;
}

void ntt_ctx_u64_free(ntt_ctx_u64* ctx)
{
  if (ctx->omega_powers) free(ctx->omega_powers);
  if (ctx->omega_inv_powers) free(ctx->omega_inv_powers);
  if (ctx->bit_rev_indices) free(ctx->bit_rev_indices);
}

// Inverse NTT using Montgomery Reduction
void ntt_u64_cyclic_inverse(u64* a_hat, const u64* a, ntt_ctx_u64* ctx)
{
  u64 n = ctx->n;
  u64 q = ctx->q;

  // 1. Bit-reversal permutation AND entering Montgomery space
  for (u64 i = 0; i < n; i++) {
    u64 rev = ctx->bit_rev_indices[i];
    a_hat[rev] = mont_in(a[i], &ctx->mctx);
  }

  // 2. Cooley-Tukey Butterfly
  for (u64 len = 2; len <= n; len <<= 1) {
    u64 half = len >> 1;
    u64 step = n / len;

    for (u64 i = 0; i < n; i += len) {
      for (u64 j = 0; j < half; j++) {
        u64 twiddle =
            ctx->omega_inv_powers[j * step];  // Already in Montgomery form

        u64 even = i + j;
        u64 odd = i + j + half;

        u64 t = mont_mul(a_hat[odd], twiddle, &ctx->mctx);
        u64 u = a_hat[even];

        a_hat[even] = mod_add(u, t, q);
        a_hat[odd] = mod_sub(u, t, q);
      }
    }
  }

  // 3. Scale by n^-1 mod q AND exit Montgomery space
  for (u64 i = 0; i < n; i++) {
    // Multiply by n_inv (which is in Montgomery form)
    a_hat[i] = mont_mul(a_hat[i], ctx->n_inv, &ctx->mctx);

    // Exit Montgomery space to get the final standard integer
    a_hat[i] = mont_out(a_hat[i], &ctx->mctx);
  }
}

void ntt_u64_cyclic_inverse_montgomery_in(u64* a_hat, const u64* a,
                                          ntt_ctx_u64* ctx)
{
  u64 n = ctx->n;
  u64 q = ctx->q;

  // 1. Bit-reversal permutation AND entering Montgomery space
  for (u64 i = 0; i < n; i++) {
    u64 rev = ctx->bit_rev_indices[i];
    a_hat[rev] = a[i];
  }

  // 2. Cooley-Tukey Butterfly
  for (u64 len = 2; len <= n; len <<= 1) {
    u64 half = len >> 1;
    u64 step = n / len;

    for (u64 i = 0; i < n; i += len) {
      for (u64 j = 0; j < half; j++) {
        u64 twiddle =
            ctx->omega_inv_powers[j * step];  // Already in Montgomery form

        u64 even = i + j;
        u64 odd = i + j + half;

        u64 t = mont_mul(a_hat[odd], twiddle, &ctx->mctx);
        u64 u = a_hat[even];

        a_hat[even] = mod_add(u, t, q);
        a_hat[odd] = mod_sub(u, t, q);
      }
    }
  }

  // 3. Scale by n^-1 mod q AND exit Montgomery space
  for (u64 i = 0; i < n; i++) {
    // Multiply by n_inv (which is in Montgomery form)
    a_hat[i] = mont_mul(a_hat[i], ctx->n_inv, &ctx->mctx);

    // Exit Montgomery space to get the final standard integer
    a_hat[i] = mont_out(a_hat[i], &ctx->mctx);
  }
}