#ifndef COMB_CUDA_H
#define COMB_CUDA_H

/**
 * \file comb_cuda.h
 *
 * Additive CUDA interface for Comb.  This header does not change the main CPU
 * ABI in comb_external.h.  Every input pointer in comb_problem_t is
 * non-owning and only needs to remain valid until comb_cuda_context_create
 * returns; the context owns its normalized host and device copies.
 */

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#include "comb_external.h"

#define COMB_CUDA_ABI_VERSION 1U
#define COMB_CUDA_RANK_NONE UINT64_MAX

/** Host/runtime outcome of a CUDA API call (not the scientific decision). */
typedef enum {
  COMB_CUDA_STATUS_SUCCESS = 0,
  COMB_CUDA_STATUS_INVALID_ARGUMENT = 1,
  COMB_CUDA_STATUS_INVALID_PROBLEM = 2,
  COMB_CUDA_STATUS_PROBLEM_TOO_LARGE = 3,
  COMB_CUDA_STATUS_CUDA_UNAVAILABLE = 4,
  COMB_CUDA_STATUS_CUDA_ERROR = 5,
  COMB_CUDA_STATUS_OUT_OF_MEMORY = 6,
  COMB_CUDA_STATUS_NOT_IMPLEMENTED = 7,
  COMB_CUDA_STATUS_INTERNAL_ERROR = 8
} comb_cuda_status_t;

/** Witness-selection policy. */
typedef enum {
  /** Use the policy stored in comb_cuda_options_t. */
  COMB_CUDA_WITNESS_DEFAULT = 0,
  /** Complete the requested range and return its smallest bad global rank. */
  COMB_CUDA_WITNESS_DETERMINISTIC = 1,
  /** Permit early termination after any independently valid bad witness. */
  COMB_CUDA_WITNESS_FAST_ANY = 2
} comb_cuda_witness_mode_t;

/**
 * Versioned, non-owning description of one verification problem.
 *
 * Lookup-table lengths are explicit because the main CPU entry point cannot
 * validate the allocation behind a raw pointer.  Lengths are in logical bits,
 * not bytes or comb_bit_word_t elements.
 */
typedef struct {
  uint32_t struct_size;
  uint32_t abi_version;

  uint32_t n_all;
  uint32_t n;
  uint32_t n_x;
  uint32_t n_z;
  uint32_t k;

  const comb_pauli_x_t *x_stabilisers;
  const comb_bit_word_t *x_lookup_table;
  uint64_t x_lookup_bit_count;
  const comb_pauli_z_t *z_stabilisers;
  const comb_bit_word_t *z_lookup_table;
  uint64_t z_lookup_bit_count;
  const comb_dcp_t *dcp_stabilisers;

  uint32_t fault_location_count;
  const uint32_t *fault_location_partitions;
  uint32_t s;

  uint32_t output_error_count;
  const comb_dcp_t *output_errors;
  uint32_t x_measurement_count;
  const comb_pauli_x_t *x_measurement_paulis;
  uint32_t z_measurement_count;
  const comb_pauli_z_t *z_measurement_paulis;
} comb_problem_t;

/** Context creation and default execution policy. */
typedef struct {
  uint32_t struct_size;
  uint32_t abi_version;
  /** CUDA ordinal, or -1 to use the runtime's current device. */
  int32_t cuda_device;
  /** comb_cuda_witness_mode_t; DEFAULT is normalized to DETERMINISTIC. */
  uint32_t witness_mode;
  /** Maximum ranks assigned to one bounded device epoch; 0 selects a default. */
  uint64_t epoch_size;
  /** Upper bound for reusable per-context device scratch; 0 selects a default. */
  uint64_t scratch_limit_bytes;
  uint64_t reserved[4];
} comb_cuda_options_t;

/** Verification result.  Timings are milliseconds measured on the host. */
typedef struct {
  uint32_t struct_size;
  uint32_t abi_version;
  comb_code_t verification_code;
  uint32_t reserved0;
  uint64_t total_combination_count;
  uint64_t checked_combination_count;
  uint64_t minimum_bad_rank;
  double setup_milliseconds;
  double verification_milliseconds;
} comb_cuda_result_t;

/**
 * Optional decoded bad witness supplied by the caller.
 *
 * fault_indices is non-owning output storage.  Its capacity must be at least
 * the problem's s when a witness is requested.
 */
typedef struct {
  uint32_t struct_size;
  uint32_t abi_version;
  uint64_t rank;
  uint32_t fault_count;
  uint32_t fault_capacity;
  uint32_t *fault_indices;
} comb_cuda_witness_t;

typedef struct comb_cuda_context comb_cuda_context_t;

#define COMB_PROBLEM_INITIALIZER                                          \
  {                                                                      \
    sizeof(comb_problem_t), COMB_CUDA_ABI_VERSION, 0, 0, 0, 0, 0, NULL, \
        NULL, 0, NULL, NULL, 0, NULL, 0, NULL, 0, 0, NULL, 0, NULL, 0,   \
        NULL                                                             \
  }
#define COMB_CUDA_OPTIONS_INITIALIZER                                      \
  {                                                                       \
    sizeof(comb_cuda_options_t), COMB_CUDA_ABI_VERSION, -1,               \
        COMB_CUDA_WITNESS_DETERMINISTIC, 0, 0, {0, 0, 0, 0}               \
  }
#define COMB_CUDA_RESULT_INITIALIZER                                       \
  {                                                                       \
    sizeof(comb_cuda_result_t), COMB_CUDA_ABI_VERSION, OTHER_ERROR, 0, 0, \
        0, COMB_CUDA_RANK_NONE, 0.0, 0.0                                  \
  }
#define COMB_CUDA_WITNESS_INITIALIZER                                   \
  {                                                                    \
    sizeof(comb_cuda_witness_t), COMB_CUDA_ABI_VERSION,                \
        COMB_CUDA_RANK_NONE, 0, 0, NULL                                \
  }

comb_cuda_status_t comb_cuda_context_create(
    const comb_problem_t *problem, const comb_cuda_options_t *options,
    comb_cuda_context_t **context);

comb_cuda_status_t comb_cuda_verify_range(
    comb_cuda_context_t *context, uint64_t rank_begin, uint64_t rank_end,
    comb_cuda_witness_mode_t mode, comb_cuda_result_t *result,
    comb_cuda_witness_t *witness);

comb_cuda_status_t comb_cuda_verify(comb_cuda_context_t *context,
                                    comb_cuda_witness_mode_t mode,
                                    comb_cuda_result_t *result,
                                    comb_cuda_witness_t *witness);

void comb_cuda_context_destroy(comb_cuda_context_t *context);

const char *comb_cuda_status_string(comb_cuda_status_t status);

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // COMB_CUDA_H
