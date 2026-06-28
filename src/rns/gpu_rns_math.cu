#include <stdint.h>
#include <cuda_runtime.h>

// __device__ means this runs on the GPU and is called by the GPU
__device__ uint64_t gpu_mont_mul(uint64_t a, uint64_t b, uint64_t p, uint64_t p_inv) {
    // 1. Calculate 128-bit product T = a * b
    uint64_t T_lo = a * b;
    // CUDA intrinsic for high 64 bits of 64x64 multiply
    uint64_t T_hi = __umul64hi(a, b); 

    // 2. Calculate m = T_lo * p_inv mod 2^64
    uint64_t m = T_lo * p_inv;

    // 3. Calculate m * p (128-bit)
    uint64_t mp_lo = m * p;
    uint64_t mp_hi = __umul64hi(m, p);

    // 4. Add T + m * p (128-bit addition)
    // We only care about the high 64 bits of the result because we shift by 2^64
    uint64_t res_hi = T_hi + mp_hi + (T_lo + mp_lo < T_lo ? 1 : 0);

    // 5. Final reduction
    if (res_hi >= p) {
        res_hi -= p;
    }
    return res_hi;
}

// BLOCK_SIZE is our tile size. 16 or 32 is standard.
#define BLOCK_SIZE 16

__global__ void batched_rns_matmul_kernel(
    uint64_t* C, const uint64_t* A, const uint64_t* B, 
    uint64_t n, const uint64_t* primes, const uint64_t* p_invs) 
{
    // Block row and column
    int row = blockIdx.y * blockDim.y + threadIdx.y;
    int col = blockIdx.x * blockDim.x + threadIdx.x;
    
    // The Z dimension tells us WHICH prime (and which matrix batch) we are working on
    int batch_idx = blockIdx.z; 

    uint64_t p = primes[batch_idx];
    uint64_t p_inv = p_invs[batch_idx];

    // Offsets for this specific batch
    uint64_t matrix_offset = batch_idx * n * n;
    const uint64_t* A_batch = A + matrix_offset;
    const uint64_t* B_batch = B + matrix_offset;
    uint64_t* C_batch = C + matrix_offset;

    uint64_t sum = 0;

    if (row < n && col < n) {
        // Simple global memory implementation
        // (For max performance, you'd load chunks of A and B into __shared__ memory here)
        for (int k = 0; k < n; ++k) {
            uint64_t a_val = A_batch[row * n + k];
            uint64_t b_val = B_batch[k * n + col];
            
            uint64_t prod = gpu_mont_mul(a_val, b_val, p, p_inv);
            
            // mod_add
            sum += prod;
            if (sum >= p || sum < prod) {
                sum -= p;
            }
        }
        C_batch[row * n + col] = sum;
    }
}
