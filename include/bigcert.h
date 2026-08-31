#ifndef BIGCERT_H
#define BIGCERT_H

#include "bignum.h"

typedef struct pocklington_cert pocklington_cert;
typedef struct pocklington_cert_elem pocklington_cert_elem;

typedef struct pocklington_cert_elem {
  bignum q;        // prime
  bignum alpha_q;  // base for pocklington theorem
  pocklington_cert*
      q_cert;  // certificate for q's primality, if NULL then q is leaf
};

// todo add enum for bpsw, trial div, ECPP
struct pocklington_cert {
  bignum N;  // prime number
  pocklington_cert_elem*
      data;  // pointer to start node of pocklington certificates
  u64 size;
  u64 capacity;
};

// base_case bound = 2^64

void gen_provable_primes_arithmetic(bignum* p, u64 n,
                                    pocklington_cert** cert_out);

void pocklington_cert_init(pocklington_cert* cert, const bignum* N);
void pocklington_cert_free(pocklington_cert* cert);
void pocklington_cert_copy(pocklington_cert* cert_a,
                           const pocklington_cert* cert_b);

void print_pocklington_cert(pocklington_cert* cert);

bool pocklington_cert_add_elem(pocklington_cert* cert, const bignum* q,
                               const bignum* alpha_q, pocklington_cert* q_cert);

bool pocklington_cert_verify(const pocklington_cert* cert, u64 base_case_bound);
bool pocklington_cert_verify_step(const pocklington_cert* cert);

bool pocklington_cert_is_base_case(const pocklington_cert* cert);
u64 pocklington_cert_get_depth(const pocklington_cert* cert);
u64 pocklington_cert_get_node_count(const pocklington_cert* cert);
size_t pocklington_cert_get_byte_size(const pocklington_cert* cert);

#endif