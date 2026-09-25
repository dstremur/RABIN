# RABIN

**RABIN** (Rapid Arithmetic for Big Integers and Number Theory) is a C library designed for arbitrary-precision arithmetic and number-theoretic algorithms.

---

> [!WARNING]
> This library is not safe for cryptography. It doesn't even attempt to be constant-time. If you seriously think about using this for crypto, you probably aren't smart enough to write your own crypto code.

---

## Features

- Arbitrary-precision integers and rationals (experimental)
- Higher-level structs for arbitrary-precision matrices over $\mathbb{Z}$ and $\mathbb{Q}$, polynomials, and vectors
- Various number theory algorithms and primality tests
- Implementation of a provable prime generation algorithm including Pocklington certificate

---

## Plan

- [ ] Finish implementation of rationals
- [ ] General hardening of the library
- [ ] Implement all algorithms in [1]
- [ ] Optimize u64 NTT, add negacyclic variant, and implement Schönhage-Strassen
- [ ] Rename functions to `rbn_` ($\mathbb{Z}$), `rbq_` ($\mathbb{Q}$)
- [ ] Implement Frobenius primality test
- [ ] Write CUDA/HIP kernels for highly parallel operations
- ...

---

## Getting Started

### Prerequisites

- A C compiler, preferably `gcc`
- `make`
- `gmp` for running the tests

### Building from Source

```bash
https://github.com/dstremur/RABIN.git
cd rabin
make
