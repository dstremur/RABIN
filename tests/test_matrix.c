#include <assert.h>
#include <stdio.h>
#include <time.h>
#include "../include/bigmatrix.h"

// Helper to print a matrix for manual debugging
void bigmatrix_print22(bigmatrix* M, const char* name)
{
  printf("Matrix %s (%llu x %llu):\n", name, M->r_size, M->c_size);
  for (u64 i = 0; i < M->r_size; i++) {
    for (u64 j = 0; j < M->c_size; j++) {
      // Replace 'bn_print' with your library's actual print function
      bn_print(&M->data[i * M->c_size + j]);
      printf("\t");
    }
    printf("\n");
  }
  printf("\n");
}

double get_time() {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

void run_det_benchmark(u64 size, u64 bits) {
    bigmatrix M;
    bignum det, val;
    
    bn_init(&det);
    bn_init(&val);
    bigmatrix_init(&M, size, size);

    // Populate with random data to prevent "easy" zeros
    for (u64 i = 0; i < size; i++) {
        for (u64 j = 0; j < size; j++) {
            bn_gen_random(&val, bits);
            bigmatrix_set(&M, &val, i, j);
        }
    }

    printf("Benchmarking %llu x %llu (%llu-bit entries)... ", size, size, bits);
    fflush(stdout);

    double start = get_time();
    bigmatrix_det(&det, &M);
    double end = get_time();

    printf("Time: %f seconds\n", end - start);

    bigmatrix_free(&M);
    bn_free(&det);
    bn_free(&val);
}

int main()
{
  printf("--- Starting BigMatrix Test Suite ---\n");

  // 1. Test Initialization and Free
  bigmatrix A, B, R;
  bigmatrix_init(&A, 2, 2);
  bigmatrix_init(&B, 2, 2);
  bigmatrix_init(&R, 2, 2);

  assert(A.data != NULL);
  assert(A.r_size == 2 && A.c_size == 2);
  printf("[PASS] Initialization\n");

  // 2. Setup Data: A = [[1, 2], [3, 4]]
  bignum val;
  bn_init(&val);

  bn_set_u64(&val, 1);
  bigmatrix_set(&A, &val, 0, 0);
  bn_set_u64(&val, 2);
  bigmatrix_set(&A, &val, 0, 1);
  bn_set_u64(&val, 3);
  bigmatrix_set(&A, &val, 1, 0);
  bn_set_u64(&val, 4);
  bigmatrix_set(&A, &val, 1, 1);
  // 3. Test bigmatrix_get
  bignum check;
  bn_init(&check);
  bigmatrix_get(&check, &A, 0, 1);  // Should be 2
  // Replace with your library's comparison function, e.g., bn_cmp_u64
  assert(bn_is_eq_i64(&check, 2));
  printf("[PASS] Get/Set\n");

  // 4. Test bigmatrix_copy
  bigmatrix_copy(&B, &A);
  bigmatrix_get(&check, &B, 1, 1);
  assert(bn_is_eq_i64(&check, 4));
  printf("[PASS] Copy\n");

  // 5. Test bigmatrix_add: R = A + B (Since B=A, R should be [[2, 4], [6, 8]])
  bigmatrix_add(&R, &A, &B);
  bigmatrix_get(&check, &R, 1, 0);
  assert(bn_is_eq_i64(&check, 6));
  printf("[PASS] Addition\n");

  // 6. Test bigmatrix_mul: R = A * Identity
  // Setup B as Identity: [[1, 0], [0, 1]]
  bn_set_u64(&val, 1);
  bigmatrix_set(&B, &val, 0, 0);
  bn_set_u64(&val, 0);
  bigmatrix_set(&B, &val, 0, 1);
  bn_set_u64(&val, 0);
  bigmatrix_set(&B, &val, 1, 0);
  bn_set_u64(&val, 1);
  bigmatrix_set(&B, &val, 1, 1);

  bigmatrix_mul(&R, &A, &B);  // Result should be A
  bigmatrix_get(&check, &R, 1, 1);
  assert(bn_is_eq_i64(&check, 4));
  printf("[PASS] Multiplication (Identity)\n");

  // 7. Test bigmatrix_mul: Rectangular Multiplication
  // A (2x3) * B (3x2) = R (2x2)
  bigmatrix RectA, RectB, RectR;
  bigmatrix_init(&RectA, 2, 3);
  bigmatrix_init(&RectB, 3, 2);
  bigmatrix_init(&RectR, 2, 2);

  // Populate with 1s
  bn_set_u64(&val, 1);
  for (u64 i = 0; i < 6; i++) {
    bn_copy(&RectA.data[i], &val);
    bn_copy(&RectB.data[i], &val);
  }

  bigmatrix_mul(&RectR, &RectA, &RectB);
  // Each cell in a 2x3 * 3x2 of all 1s should be 3
  bigmatrix_get(&check, &RectR, 0, 0);
  assert(bn_is_eq_i64(&check, 3));
  printf("[PASS] Multiplication (Rectangular)\n");
 
  bigmatrix_print(&R); 

  bigmatrix_det(&val, &A); 
  bn_println(&val); 


  // 8. Cleanup
  bn_free(&val);
  bn_free(&check);
  bigmatrix_free(&A);
  bigmatrix_free(&B);
  bigmatrix_free(&R);
  bigmatrix_free(&RectA);
  bigmatrix_free(&RectB);
  bigmatrix_free(&RectR);

  printf("[PASS] Cleanup / Free\n");
  printf("--- All Tests Passed! ---\n");


  u64 n = 100; 

  bigmatrix X;
  bignum t;
  bn_init(&t); 
  bigmatrix_init(&X, n, n); 
  for (u64 i = 0; i < n; i++) {
	  for (u64 j = 0; j < n; j++) {
		//bn_set_i64(&t, rand()); 
		bn_gen_random(&t, 32); 
		bigmatrix_set(&X, &t, i, j);
	  }
  }

  bigmatrix_print22(&X, "343");

  
	bigmatrix_det(&t, &X);

	bn_println(&t);

  u64 sizes[] = {2, 4, 8, 16, 32, 64, 128, 256, 300, 512};
    int num_tests = sizeof(sizes) / sizeof(sizes[0]);

    for (int i = 0; i < num_tests; i++) {
        run_det_benchmark(sizes[i], 32); // 32-bit random entries
    }

	bn_free(&t);
	bigmatrix_free(&X);

  return 0;
}
