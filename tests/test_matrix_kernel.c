/*
 * Runtime benchmark for the Cohen 2.7 kernel routines: rmat_kernel(),
 * rmat_kernel_lll() and rmat_kernel_image(). For each (shape, bits) pair
 * the HNF-route kernel and the LLL-route kernel are first checked to be
 * equal (same lattice, compared via HNF), then every function is timed on
 * random matrices of that shape until a minimum wall-clock budget is
 * spent; the reported value is the mean seconds per call.
 */
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "../include/rmat.h"
#include "../include/rz.h"
#include "../include/rzrand.h"

#define BENCH_TARGET_SEC 0.1

static double now_sec(void)
{
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec;
}

typedef struct {
  u64 rows;
  u64 cols;
  u64 bits;
} bench_case_t;

static rmat_t A;

static void build_case(u64 rows, u64 cols, u64 bits)
{
  if (rmat_init(&A, rows, cols) != RABIN_SUCCESS) exit(1);
  rz_t e;
  if (rz_init(&e) != RABIN_SUCCESS) exit(1);
  for (u64 r = 0; r < rows; r++) {
    for (u64 c = 0; c < cols; c++) {
      if (rz_gen_random(&e, bits) != RABIN_SUCCESS) exit(1);
      if (rand() & 1) rz_neg(&e, &e);
      if (rmat_set(&A, &e, r, c) != RABIN_SUCCESS) exit(1);
    }
  }
  rz_clear(&e);
}

static bool same_lattice(const rmat_t* X, const rmat_t* Y)
{
  if (X->cols == 0 || Y->cols == 0) return X->cols == Y->cols;
  rmat_t HX = {0}, HY = {0};
  if (rmat_hermite(&HX, X) != RABIN_SUCCESS ||
      rmat_hermite(&HY, Y) != RABIN_SUCCESS)
    exit(1);
  const bool ok = (HX.cols == HY.cols) && rmat_equal(&HX, &HY);
  rmat_clear(&HX);
  rmat_clear(&HY);
  return ok;
}

/* mode: 0 = rmat_kernel, 1 = rmat_kernel_lll, 2 = rmat_kernel_image */
static double bench(int mode, const char* name)
{
  u64 calls = 0;
  double total = 0.0;
  do {
    rmat_t K = {0};
    u64 p = 0;
    const double t0 = now_sec();
    rabin_err_t err = mode == 0   ? rmat_kernel(&A, &K)
                      : mode == 1 ? rmat_kernel_lll(&A, &K)
                                  : rmat_kernel_image(&A, &K, &p);
    total += now_sec() - t0;
    rmat_clear(&K);
    if (err != RABIN_SUCCESS) {
      fprintf(stderr, "[FAIL] %s returned %d\n", name, (int)err);
      exit(1);
    }
    calls++;
  } while (total < BENCH_TARGET_SEC);
  return total / (double)calls;
}

int main(void)
{
  srand(0xBEEF);

  static const bench_case_t cases[] = {
      {4, 6, 32},   {6, 8, 32},  {8, 10, 32},  {12, 15, 32}, {4, 6, 64},
      {6, 8, 64},   {8, 10, 64}, {12, 15, 64}, {4, 6, 128},  {6, 8, 128},
      {8, 10, 128}, {4, 6, 256}, {6, 8, 256},
  };
  const u64 nc = sizeof(cases) / sizeof(cases[0]);

  printf("+---------------------+----------+------+-------------------+\n");
  printf("| %-19s | %-8s | %4s | %17s |\n", "function", "shape", "bits",
         "avg time (s)");
  printf("+---------------------+----------+------+-------------------+\n");
  for (u64 i = 0; i < nc; i++) {
    char shape[16];
    snprintf(shape, sizeof(shape), "%lux%lu", (u64)cases[i].rows,
             (u64)cases[i].cols);
    build_case(cases[i].rows, cases[i].cols, cases[i].bits);

    // the HNF-route and the LLL-route kernel must agree
    rmat_t K_hnf = {0}, K_lll = {0};
    if (rmat_kernel(&A, &K_hnf) != RABIN_SUCCESS ||
        rmat_kernel_lll(&A, &K_lll) != RABIN_SUCCESS) {
      fprintf(stderr, "[FAIL] kernel call failed for %s/%lu bits\n", shape,
              (u64)cases[i].bits);
      return 1;
    }
    if (!same_lattice(&K_hnf, &K_lll)) {
      fprintf(stderr, "[FAIL] kernel routes disagree for %s/%lu bits\n", shape,
              (u64)cases[i].bits);
      return 1;
    }
    rmat_clear(&K_hnf);
    rmat_clear(&K_lll);

    for (int m = 0; m < 3; m++) {
      const char* name = m == 0   ? "rmat_kernel"
                         : m == 1 ? "rmat_kernel_lll"
                                  : "rmat_kernel_image";
      const double t = bench(m, name);
      printf("| %-19s | %-8s | %4lu | %17.9f |\n", name, shape,
             (u64)cases[i].bits, t);
    }
    rmat_clear(&A);
  }
  printf("+---------------------+----------+------+-------------------+\n");
  printf("Benchmark finished\n");
  return 0;
}
