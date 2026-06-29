#include "../../include/u64.h"

void mont_init(mont_ctx* ctx, u64 p)
{
  ctx->p = p;

  ctx->p_inv = mod_inverse_u64(p);

  // Calculate 2^128 mod p
  unsigned __int128 r2 = ((unsigned __int128)1 << 64) % p;
  r2 = (r2 * r2) % p;
  ctx->r2_mod_p = (u64)r2;
}

u64 mont_redc(unsigned __int128 T, mont_ctx* ctx)
{
  u64 p = ctx->p;

  // T mod R
  u64 T_modR = (u64)T;

  u64 m = T_modR * ctx->p_inv;

  unsigned __int128 t_wide = (T + (unsigned __int128)m * p);
  u64 t = (u64)(t_wide >> 64);

  return (t >= p) ? (t - p) : t;
}

u64 mont_mul(u64 a, u64 b, mont_ctx* ctx)
{
  unsigned __int128 R = (unsigned __int128)a * b;

  return mont_redc(R, ctx);
}

u64 mont_in(u64 a, mont_ctx* ctx) { return mont_mul(a, ctx->r2_mod_p, ctx); }

u64 mont_out(u64 a_hat, mont_ctx* ctx) { return mont_mul(a_hat, 1, ctx); }

u64 mont_inverse(u64 a_mont, const mont_ctx* ctx)
{
  u64 a = mont_out(a_mont, ctx);

  u64 inv = mod_inverse_euclid(a, ctx->p);

  return mont_in(inv, ctx);
}