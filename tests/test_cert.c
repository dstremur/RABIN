/*
 * test_cert.c
 *
 * Pocklington certificate tests (Phase 3):
 *
 *   1. Hand-built certificate chains for small primes must verify.
 *   2. Forged certificates must be rejected: tampered N, q, alpha,
 *      truncated tree, out-of-bound leaf, composite leaf factor,
 *      composite root, child/factor mismatch.
 *   3. Tree helpers: is_base_case, get_depth, get_node_count,
 *      get_byte_size.
 *   4. Round-trip: rz_gen_provable_arithmetic() output must pass
 *      rzcert_verify(); a tampered generated certificate must fail.
 *
 * Copyright (C) 2026 Diego Strebel
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <stdio.h>
#include <stdlib.h>

#include "../include/rabin.h"
#include "../include/rzcert.h"
#include "test_util.h"

/* Build a leaf certificate (size 0) for the value n. */
static rzcert_t* cert_leaf(u64 n)
{
  rzcert_t* c = calloc(1, sizeof(rzcert_t));
  if (c == NULL) {
    fprintf(stderr, "OOM in cert_leaf\n");
    exit(EXIT_FAILURE);
  }
  rz_init(&c->N);
  if (rz_set_u64(&c->N, n) != RABIN_SUCCESS) {
    fprintf(stderr, "rz_set_u64 failed\n");
    exit(EXIT_FAILURE);
  }
  return c;
}

/* Build an interior certificate node N with one element (q, alpha) whose
   child proves q. */
static rzcert_t* cert_node(u64 n, u64 q, u64 alpha, rzcert_t* child)
{
  rzcert_t* c = cert_leaf(n);
  rz_t bq, ba;
  rz_init_multi(&bq, &ba, NULL);
  rz_set_u64(&bq, q);
  rz_set_u64(&ba, alpha);
  if (!rzcert_add_elem(c, &bq, &ba, child)) {
    fprintf(stderr, "rzcert_add_elem failed\n");
    exit(EXIT_FAILURE);
  }
  rz_clear_multi(&bq, &ba, NULL);
  return c;
}

/* Append a second element without freeing (used for multi-element checks).
   q and alpha are u64 for convenience. */
static void cert_push(rzcert_t* c, u64 q, u64 alpha, rzcert_t* child)
{
  rz_t bq, ba;
  rz_init_multi(&bq, &ba, NULL);
  rz_set_u64(&bq, q);
  rz_set_u64(&ba, alpha);
  if (!rzcert_add_elem(c, &bq, &ba, child)) {
    fprintf(stderr, "rzcert_add_elem failed\n");
    exit(EXIT_FAILURE);
  }
  rz_clear_multi(&bq, &ba, NULL);
}

/* =====================================================================
   1. Hand-built valid chains
   ===================================================================== */

TEST(hand_built_chains_verify)
{
  /* N = 7:  N-1 = 6, q = 3, alpha = 3; leaf 3.
     (3+1)^2 = 16 > 7; 3^6 = 729 = 104*7 + 1; 3^2 - 1 = 8, gcd(8,7) = 1. */
  rzcert_t* c7 = cert_node(7, 3, 3, cert_leaf(3));
  TEST_ASSERT(rzcert_verify(c7, 8), "chain 7<-3 verifies (bound 8)");
  TEST_ASSERT(!rzcert_verify(c7, 1), "chain 7<-3 fails with bound 1");
  rzcert_clear(c7);

  /* N = 13: N-1 = 12, q = 3, alpha = 2; leaf 3.
     16 > 13; 2^12 = 1 mod 13; 2^4 - 1 = 15, gcd(15,13) = 1. */
  rzcert_t* c13 = cert_node(13, 3, 2, cert_leaf(3));
  TEST_ASSERT(rzcert_verify(c13, 8), "chain 13<-3 verifies");
  rzcert_clear(c13);

  /* Chain 47 <- 23 <- 11 <- 5 <- 2, all with checked witnesses:
     47: q=23, alpha=5   (24^2 = 576 > 47; 5^46 = 1; gcd(5^2-1,47) = 1)
     23: q=11, alpha=2   (12^2 = 144 > 23; 2^22 = 1; gcd(2^2-1,23) = 1)
     11: q=5,  alpha=2   (6^2  =  36 > 11; 2^10 = 1; gcd(2^2-1,11) = 1)
      5: q=2,  alpha=2   (3^2  =   9 >  5; 2^4  = 1; gcd(2^2-1,5)  = 1) */
  rzcert_t* c47 = cert_node(
      47, 23, 5,
      cert_node(23, 11, 2,
                cert_node(11, 5, 2, cert_node(5, 2, 2, cert_leaf(2)))));
  TEST_ASSERT(rzcert_verify(c47, 4), "chain 47<-23<-11<-5<-2 verifies");
  TEST_ASSERT(!rzcert_verify(c47, 1), "chain fails when leaf 2 exceeds bound");
  rzcert_clear(c47);
}

/* =====================================================================
   2. Forged certificates must fail
   ===================================================================== */

TEST(forged_certificates_rejected)
{
  /* Tampered N: 7 -> 9 breaks q | N-1 (3 does not divide 8). */
  rzcert_t* c = cert_node(7, 3, 3, cert_leaf(3));
  TEST_ASSERT(rz_add_u64(&c->N, &c->N, 2) == RABIN_SUCCESS, "bump N");
  TEST_ASSERT(!rzcert_verify(c, 64), "tampered N rejected");
  rzcert_clear(c);

  /* Tampered q: 3 -> 5 breaks q | N-1. */
  c = cert_node(7, 3, 3, cert_leaf(3));
  TEST_ASSERT(rz_add_u64(&c->data[0].q, &c->data[0].q, 2) == RABIN_SUCCESS,
              "bump q");
  TEST_ASSERT(!rzcert_verify(c, 64), "tampered q rejected");
  rzcert_clear(c);

  /* Tampered alpha: alpha = N is 0 mod N, fails Fermat. */
  c = cert_node(13, 3, 2, cert_leaf(3));
  TEST_ASSERT(rz_copy(&c->data[0].alpha_q, &c->N) == RABIN_SUCCESS,
              "alpha = N");
  TEST_ASSERT(!rzcert_verify(c, 64), "alpha = 0 mod N rejected");
  rzcert_clear(c);

  /* Truncated tree: drop the child for q = 23 (5 bits) with bound 4. */
  rzcert_t* leaf2 = cert_leaf(2);
  c = cert_node(
      47, 23, 5,
      cert_node(23, 11, 2, cert_node(11, 5, 2, cert_node(5, 2, 2, leaf2))));
  rzcert_t* dropped = c->data[0].q_cert;
  c->data[0].q_cert = NULL;
  TEST_ASSERT(!rzcert_verify(c, 4), "truncated tree rejected (q above bound)");
  /* With a bound that trusts 23, BPSW(23) accepts it: documented behavior. */
  TEST_ASSERT(rzcert_verify(c, 64), "truncated tree OK once q is in bound");
  c->data[0].q_cert = dropped; /* re-attach so the tree frees fully */
  rzcert_clear(c);

  /* Child proves the wrong value: child N = 5 while q = 3. */
  c = cert_node(7, 3, 3, cert_leaf(5));
  TEST_ASSERT(!rzcert_verify(c, 64), "child/factor mismatch rejected");
  rzcert_clear(c);

  /* Composite leaf factor: N = 19, fake q = 9 divides 18 and (10)^2 > 19,
     and alpha = 2 satisfies both congruence conditions, but BPSW(9) is
     false, so the certificate must fail on the leaf alone. */
  c = cert_node(19, 9, 2, cert_leaf(9));
  TEST_ASSERT(!rzcert_verify(c, 64), "composite leaf factor rejected");
  rzcert_clear(c);

  /* Composite root: 341 = 11 * 31, q = 17 divides 340 but (18)^2 = 324
     is not greater than 341, so the Pocklington size condition fails. */
  c = cert_node(341, 17, 3, cert_leaf(17));
  TEST_ASSERT(!rzcert_verify(c, 64), "composite root with small q rejected");
  rzcert_clear(c);

  /* Degenerate roots: N = 0, negative N, empty element array. */
  c = cert_leaf(0);
  TEST_ASSERT(!rzcert_verify(c, 64), "N = 0 rejected");
  rzcert_clear(c);

  c = cert_leaf(3);
  c->N.is_neg = true;
  TEST_ASSERT(!rzcert_verify(c, 64), "negative N rejected");
  rzcert_clear(c);

  /* size > 0 but data == NULL is malformed, not a base case. */
  c = cert_leaf(7);
  c->size = 1;
  c->data = NULL;
  TEST_ASSERT(!rzcert_verify(c, 64), "size>0 with data==NULL rejected");
  c->size = 0;
  rzcert_clear(c);

  /* Even interior N can never satisfy the conditions. */
  c = cert_node(8, 7, 3, cert_leaf(7));
  TEST_ASSERT(!rzcert_verify(c, 64), "even interior N rejected");
  rzcert_clear(c);
}

/* =====================================================================
   3. Tree helpers
   ===================================================================== */

TEST(cert_helpers)
{
  rzcert_t* leaf = cert_leaf(2);
  TEST_ASSERT(rzcert_is_base_case(leaf), "leaf is base case");
  TEST_ASSERT(rzcert_get_depth(leaf) == 1, "leaf depth 1");
  TEST_ASSERT(rzcert_get_node_count(leaf) == 1, "leaf node count 1");
  TEST_ASSERT(rzcert_get_byte_size(leaf) >= sizeof(rzcert_t),
              "leaf byte size sane");

  rzcert_t* c47 = cert_node(
      47, 23, 5,
      cert_node(23, 11, 2, cert_node(11, 5, 2, cert_node(5, 2, 2, leaf))));
  TEST_ASSERT(!rzcert_is_base_case(c47), "interior is not base case");
  TEST_ASSERT(rzcert_get_depth(c47) == 5, "chain depth 5");
  TEST_ASSERT(rzcert_get_node_count(c47) == 5, "chain node count 5");
  TEST_ASSERT(rzcert_get_byte_size(c47) > 5 * rzcert_get_byte_size(leaf),
              "byte size accumulates over the tree");

  /* Multi-element growth path of rzcert_add_elem. */
  cert_push(c47, 23, 5, NULL);
  TEST_ASSERT(c47->size == 2, "add_elem grows data array");
  TEST_ASSERT(rzcert_get_node_count(c47) == 5,
              "node count unchanged for leaf factor");
  rzcert_clear(c47);

  TEST_ASSERT(!rzcert_is_base_case(NULL), "NULL is not base case");
  TEST_ASSERT(rzcert_get_depth(NULL) == 0, "NULL depth 0");
  TEST_ASSERT(rzcert_get_node_count(NULL) == 0, "NULL node count 0");
  TEST_ASSERT(rzcert_get_byte_size(NULL) == 0, "NULL byte size 0");
}

/* =====================================================================
   4. Generation round-trip
   ===================================================================== */

TEST(generated_certificates_verify)
{
  const int bits_list[] = {256, 384, 512};
  const u64 bound = 64;

  for (size_t k = 0; k < sizeof(bits_list) / sizeof(bits_list[0]); k++) {
    rz_t p;
    rz_init(&p);
    rzcert_t* cert = NULL;

    rabin_err_t err = rz_gen_provable_arithmetic(&p, (u64)bits_list[k], &cert);
    TEST_ASSERT(err == RABIN_SUCCESS, "generation succeeds");
    TEST_ASSERT(cert != NULL, "certificate emitted");

    TEST_ASSERT(rz_cmp(&cert->N, &p) == 0, "certificate root matches p");
    TEST_ASSERT(rzcert_verify(cert, bound), "generated cert verifies");
    TEST_ASSERT(rzcert_get_depth(cert) >= 2, "generated cert has depth");
    TEST_ASSERT(rzcert_get_node_count(cert) == rzcert_get_depth(cert),
                "generated chain is linear");
    TEST_ASSERT(rzcert_get_byte_size(cert) > 0, "byte size nonzero");

    /* Forged variant: bumping the root factor breaks divisibility. */
    TEST_ASSERT(
        rz_add_u64(&cert->data[0].q, &cert->data[0].q, 2) == RABIN_SUCCESS,
        "tamper q");
    TEST_ASSERT(!rzcert_verify(cert, bound), "tampered generated cert fails");

    rzcert_clear(cert);
    rz_clear(&p);
  }
}

int main(void)
{
  rz_init_constants();

  RUN_TEST(hand_built_chains_verify);
  RUN_TEST(forged_certificates_rejected);
  RUN_TEST(cert_helpers);
  RUN_TEST(generated_certificates_verify);

  rz_clear_constants();
  return TEST_REPORT();
}
