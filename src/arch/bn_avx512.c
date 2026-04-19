#include <immintrin.h>
#include <stdlib.h>
#include <string.h>

#include "../include/bignum.h"
/*
void convert_to_52(bn52* dest, const bignum* src)
{
  if (src->size == 0) {
    dest->size = 0;
    dest->limbs = NULL;
    return;
  }

  u64 num_limbs52 = (src->size * 64 + 51) / 52;
  dest->limbs = (u64*)calloc(num_limbs52 + 8, sizeof(u64));
  dest->size = num_limbs52;

  u64 bit_acc = 0;
  int bits_in_acc = 0;
  u64 src_idx = 0;
  u64 dst_idx = 0;

  while (src_idx < src->size || bits_in_acc > 0) {
    // Fill the accumulator if we have space and source data
    if (bits_in_acc < 52 && src_idx < src->size) {
      u64 word = src->limbs[src_idx++];
      bit_acc |= (word << bits_in_acc);

      // Extract 52 bits
      dest->limbs[dst_idx++] = bit_acc & 0xFFFFFFFFFFFFF;

      // The remaining bits (64 - (52 - old_bits))
      int consumed = 52 - bits_in_acc;
      bit_acc = word >> consumed;
      bits_in_acc = 64 - consumed;
    } else {
      // Flush remaining bits in accumulator
      dest->limbs[dst_idx++] = bit_acc & 0xFFFFFFFFFFFFF;
      bit_acc >>= 52;
      bits_in_acc = (bits_in_acc > 52) ? bits_in_acc - 52 : 0;
    }

    if (dst_idx >= num_limbs52) break;
  }
}

static void bn52_alloc(bn52* n, uint64_t size)
{
  n->size = size;
  // We allocate 8 extra limbs so _mm512_loadu/storeu
  // never access unmapped memory at the end of the loop.
  n->limbs = (uint64_t*)calloc(size + 8, sizeof(uint64_t));
}

void bn52_free(bn52* n)
{
  if (n->limbs) free(n->limbs);
  n->limbs = NULL;
  n->size = 0;
}

// multiply using avx512 IFMA, only works for 52 bit limbs
void bn_mul_add_avx512(bn52* r, bn52* a, bn52* b)
{
  u64 max = a->size + b->size;
  bn52_alloc(r, max);

  // Main loop
  for (u64 i = 0; i < a->size; i++) {
    if (a->limbs[i] == 0) continue;

    // --- NEW: PERIODIC CARRY RIPPLE ---
    // Prevents 64-bit overflow when summing many 52-bit products
    if (i > 0 && i % 2048 == 0) {
      u64 carry = 0;
      for (u64 k = 0; k < r->size; k++) {
        u64 val = r->limbs[k] + carry;
        r->limbs[k] = val & 0xFFFFFFFFFFFFF;
        carry = val >> 52;
      }
    }

    // set curr limb of a
    __m512i vA = _mm512_set1_epi64(a->limbs[i]);
    u64 j = 0;

    for (; j + 7 < b->size; j += 8) {
      u64 idx = i + j;

      // load 8 limbs of b
      __m512i vB = _mm512_loadu_si512((__m512i*)&b->limbs[j]);

      // Low 52 bits

      __m512i vR_lo = _mm512_loadu_si512((__m512i*)&r->limbs[idx]);
      vR_lo = _mm512_madd52lo_epu64(vR_lo, vA, vB);
      _mm512_storeu_si512((__m512i*)&r->limbs[idx], vR_lo);

      // high 52 bits
      __m512i vR_hi = _mm512_loadu_si512((__m512i*)&r->limbs[idx + 1]);
      vR_hi = _mm512_madd52hi_epu64(vR_hi, vA, vB);
      _mm512_storeu_si512((__m512i*)&r->limbs[idx + 1], vR_hi);
    }
    // Scalar Fallback for tail of B
    for (; j < b->size; j++) {
      uint64_t idx = i + j;
      unsigned __int128 prod = (unsigned __int128)a->limbs[i] * b->limbs[j];

      r->limbs[idx] += (uint64_t)(prod & 0xFFFFFFFFFFFFF);
      r->limbs[idx + 1] += (uint64_t)(prod >> 52);
    }
  }

  // 3. Normalization (The Carry Pass)
  // Convert the accumulated sums back into proper base 2^52
  uint64_t carry = 0;
  for (uint64_t k = 0; k < r->size; k++) {
    uint64_t val = r->limbs[k] + carry;
    r->limbs[k] = val & 0xFFFFFFFFFFFFF;
    carry = val >> 52;
  }

  // If there is an ultimate carry, we use the padding space
  if (carry) {
    r->limbs[r->size] = carry;
    r->size++;
  }
}

void convert_back_to_64(bignum* r, const bn52* src)
{
  if (src->size == 0) {
    r->size = 0;
    return;
  }

  u64 total_bits = src->size * 52;
  u64 num_limbs64 = (total_bits + 63) / 64;

  // Use your library's allocator
  bn_alloc(r, num_limbs64);
  memset(r->limbs, 0, r->capacity * sizeof(u64));
  r->size = num_limbs64;

  for (u64 i = 0; i < src->size; i++) {
    u64 val = src->limbs[i];
    u64 current_bit = i * 52;

    u64 limb_idx = current_bit / 64;
    u64 bit_offset = current_bit % 64;

    r->limbs[limb_idx] |= (val << bit_offset);

    // If the 52 bits cross a 64-bit boundary
    if (bit_offset > 12 && (limb_idx + 1) < r->size) {
      r->limbs[limb_idx + 1] |= (val >> (64 - bit_offset));
    }
  }
  bn_trim(r);
}

void bn_mul_512(bignum* r, bignum* a, bignum* b)
{
  if (a->size == 0 || b->size == 0) {
    r->size = 0;
    return;
  }

  if (a->size >= 256 && b->size >= 256) {
    bn52 wa, wb;
    convert_to_52(&wa, a);
    convert_to_52(&wb, b);

    bn52 wr;
    bn_mul_add_avx512(&wr, &wa, &wb);

    convert_back_to_64(r, &wr);

    free(wa.limbs);
    free(wb.limbs);
    free(wr.limbs);

  } else {
    bn_mul(r, a, b);
  }

  r->is_neg = a->is_neg ^ b->is_neg;
}

*/
