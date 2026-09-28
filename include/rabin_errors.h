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

#endif
