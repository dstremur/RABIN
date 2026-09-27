# Project Overview
This repository contains RABIN, a low-level C bignum library for arbitrary precision integers, rationals, matrices and number theory functions. It is designed for high performance, and zero external dependencies.
It is also a learning project that aims to implement as many algorithms as possible.

# Build & Test Commands
- **Build Library:** `make clean && make`
- **Run Unit Tests:** `make test_` there are several available, make test_complete runs all, but takes ~15 min
- **Format Code:** `./format.sh`

# Naming Conventions
- **Types:** Use `_t` suffix.
  - Integers: `rz_t`
  - Rationals: `rq_t`
  - Integer Matrices: `rmat_t`
  - Rational Matrices: `rmatq_t`
- **Functions:** Lowercase `snake_case` with domain prefixes (e.g., `rz_add`, `rmatq_invert`).
- **Macros/Constants:** UPPERCASE with prefixes (e.g., `RABIN_SUCCESS`, `RZ_MAX`).

# Code Style & Architecture
- **Language Standard:** C11.
- **Formatting:** Use the formatter. Braces on new lines for functions (Allman style), same-line for control statements.
- **Structure:** `include/` for public API headers, `src/` for internal implementation files, `tests/` for unit tests. Do not expose internal helper functions in public headers.
- `src` has `arch` for asm routines, `core` for core functions, etc. 

# C Coding & Safety Guidelines
- **Memory Allocation:** Always check memory allocation results. Never leak memory on early returns or error states.
- **Pointer Safety:** Validate all public API pointer inputs against NULL. Return appropriate error codes if invalid.
- **Integer Safety:** Explicitly check for arithmetic overflows before performing allocations or buffer manipulations.
- **Side Effects:** Functions must not print to `stdout` or `stderr`. Output status via return codes.

# Doxygen Documentation Standards
- **Placement:** All Doxygen comments must reside exclusively in public header files (`include/*.h`), not in `.c` files.
- **Structure:** Use the standard `/** ... */` block syntax. Keep descriptions punchy and avoid redundant memory lifecycle explanations if they follow standard library rules.
- **Standard Tags:**
  - `@brief`: Single-line summary of the function, struct, or enum.
  - `@param[in]` / `@param[out]`: Document all parameters explicitly.
  - `@return`: Must specify the `rabin_err_t` return code (e.g., "RABIN_SUCCESS on success, or RABIN_ERR_OUT_OF_MEMORY").
  - `@pre` / `@warning`: Group multiple preconditions into a single `@pre` block to avoid vertical bloat.
- **Math & LaTeX:** Use `\f$ ... \f$` for inline equations. Ensure proper spacing (e.g., use `\cdot a`, never `\cdota`).
- **Literature References:** Use `@par Algorithm Reference:` to cite academic papers, algorithms, or cryptography standards at the bottom of the description.
- **Cross-Referencing:** Always include trailing parentheses when mentioning functions (e.g., `rz_add()`) to trigger auto-linking. Use `@see` at the end of the block for related functions.

# Agent Workflow Requirements
- Run tests after completing any modification to ensure zero regressions. DONT RUN TEST_LOGIC, IF A TEST TIMES OUT; STOP AND ASK FOR HELP
- Keep modifications scoped strictly to the requested feature or bug fix.
- Do not introduce external library dependencies (e.g., OpenSSL, GMP) under any circumstances.
