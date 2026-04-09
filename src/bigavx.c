#include <immintrin.h>
#include <stdio.h>

#include "../include/bignum.h"

void mac_52bit_avx512(__m512i* lo, __m512i* hi, __m512i A, u64 b)
{
  // broadcast to all 8 lanes
  __m512i B = _mm512_set1_epi64(b);

  *lo = _mm512_madd52lo_epu64(*lo, A, B);

  *hi = _mm512_madd52hi_epu64(*hi, A, B);
}

void normalize(__m512i* lo, __m512i* hi)
{
  __m512i mask = _mm512_set1_epi64((1ULL << 52) - 1);

  // Extract overflow from low part
  __m512i overflow = _mm512_srli_epi64(*lo, 52);

  // Keep only 52 bits in the low part
  *lo = _mm512_and_si512(*lo, mask);

  // Add overflow to the high part (the carry-ins for the next limb set)
  *hi = _mm512_add_epi64(*hi, overflow);
}

void avxtest()
{
  // Example: Multiply [10, 20, 30, 40, 50, 60, 70, 80] by 5
  uint64_t val_a[8] = {10, 20, 30, 40, 50, 60, 70, 80000000000};
  uint64_t b = 5;

  __m512i acc_lo = _mm512_setzero_si512();
  __m512i acc_hi = _mm512_setzero_si512();
  __m512i vec_a = _mm512_loadu_si512((__m512i*)val_a);

  mac_52bit_avx512(&acc_lo, &acc_hi, vec_a, b);
  normalize(&acc_lo, &acc_hi);

  // Store results back to print
  uint64_t res_lo[8], res_hi[8];
  _mm512_storeu_si512((__m512i*)res_lo, acc_lo);
  _mm512_storeu_si512((__m512i*)res_hi, acc_hi);

  for (int i = 0; i < 8; i++) {
    printf("Limb %d: Result = %lu, Carry for next = %lu\n", i, res_lo[i],
           res_hi[i]);
  }
}
