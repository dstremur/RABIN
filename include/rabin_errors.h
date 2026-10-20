#ifndef RABIN_ERRORS_H
#define RABIN_ERRORS_H

/**
 * @brief Error codes returned by all public RABIN functions.
 *
 * Every public mathematical, vector, polynomial and matrix operation
 * returns a \c rabin_err_t. \c RABIN_SUCCESS signals a completed
 * operation; all other values indicate the rejected precondition or the
 * resource failure that aborted it.
 */
typedef enum {
  RABIN_SUCCESS = 0,             ///< Operation completed.
  RABIN_ERR_NULL_PTR = -1,       ///< A required pointer argument was NULL.
  RABIN_ERR_OUT_OF_MEMORY = -2,  ///< An internal allocation or growth failed.
  RABIN_ERR_DIV_BY_ZERO = -3,    ///< Division or inversion by zero.
  RABIN_ERR_OVERFLOW = -4,       ///< An integer or dimension overflow.
  RABIN_ERR_INVALID_ARG = -5,    ///< An argument is out of domain.
  RABIN_ERR_MATRIX_DIM = -6      ///< Matrix dimensions do not match.
} rabin_err_t;

/**
 * @brief Run a fallible call, forwarding its status to a local cleanup path.
 *
 * Evaluates \p expr once, and on failure records the code in a local
 * `rabin_err_t err` and jumps to the enclosing function's `cleanup` label.
 * `err` is left untouched when the call succeeds, so a status set earlier
 * survives. The enclosing function must declare `err` and provide a
 * `cleanup` label that releases every resource and returns `err`.
 */
#define TRY(expr)                \
  do {                           \
    rabin_err_t _err = (expr);   \
    if (_err != RABIN_SUCCESS) { \
      err = _err;                \
      goto cleanup;              \
    }                            \
  } while (0)

#endif
