#include "comb_cuda.h"
#include "comb_cuda_internal.h"

#include <limits.h>

/*
 * CPU-only builds retain the complete CUDA ABI and fail explicitly.  CUDA
 * targets must compile this translation unit with COMB_ENABLE_CUDA_BACKEND or
 * omit it, then link comb_cuda.cu instead.
 */
#ifndef COMB_ENABLE_CUDA_BACKEND

comb_cuda_status_t comb_cuda_context_create(
    const comb_problem_t *problem, const comb_cuda_options_t *options,
    comb_cuda_context_t **context) {
  void *blob = NULL;
  uint64_t blob_size = 0;
  comb_cuda_status_t status;
  if (context == NULL) return COMB_CUDA_STATUS_INVALID_ARGUMENT;
  *context = NULL;
  status = comb_cuda_flat_problem_pack(problem, &blob, &blob_size);
  if (status != COMB_CUDA_STATUS_SUCCESS) return status;
  status = comb_cuda_flat_problem_validate_science(blob, blob_size);
  comb_cuda_flat_problem_free(blob);
  if (status != COMB_CUDA_STATUS_SUCCESS) return status;
  if (options != NULL &&
      (options->struct_size < sizeof(*options) ||
       options->abi_version != COMB_CUDA_ABI_VERSION ||
       options->cuda_device < -1 ||
       (options->witness_mode != COMB_CUDA_WITNESS_DEFAULT &&
        options->witness_mode != COMB_CUDA_WITNESS_DETERMINISTIC &&
        options->witness_mode != COMB_CUDA_WITNESS_FAST_ANY))) {
    return COMB_CUDA_STATUS_INVALID_ARGUMENT;
  }
  return COMB_CUDA_STATUS_CUDA_UNAVAILABLE;
}

comb_cuda_status_t comb_cuda_device_pci_bus_id(int32_t device, char *buffer,
                                               size_t buffer_size) {
  if (device < 0 || buffer == NULL || buffer_size < 16U ||
      buffer_size > (size_t)INT_MAX) {
    return COMB_CUDA_STATUS_INVALID_ARGUMENT;
  }
  buffer[0] = '\0';
  return COMB_CUDA_STATUS_CUDA_UNAVAILABLE;
}

comb_cuda_status_t comb_cuda_verify_range(
    comb_cuda_context_t *context, uint64_t rank_begin, uint64_t rank_end,
    comb_cuda_witness_mode_t mode, comb_cuda_result_t *result,
    comb_cuda_witness_t *witness) {
  (void)rank_begin;
  (void)rank_end;
  (void)mode;
  (void)result;
  (void)witness;
  if (context == NULL) return COMB_CUDA_STATUS_INVALID_ARGUMENT;
  return COMB_CUDA_STATUS_CUDA_UNAVAILABLE;
}

comb_cuda_status_t comb_cuda_verify(comb_cuda_context_t *context,
                                    comb_cuda_witness_mode_t mode,
                                    comb_cuda_result_t *result,
                                    comb_cuda_witness_t *witness) {
  (void)mode;
  (void)result;
  (void)witness;
  if (context == NULL) return COMB_CUDA_STATUS_INVALID_ARGUMENT;
  return COMB_CUDA_STATUS_CUDA_UNAVAILABLE;
}

void comb_cuda_context_destroy(comb_cuda_context_t *context) { (void)context; }

const char *comb_cuda_status_string(comb_cuda_status_t status) {
  switch (status) {
    case COMB_CUDA_STATUS_SUCCESS:
      return "success";
    case COMB_CUDA_STATUS_INVALID_ARGUMENT:
      return "invalid argument";
    case COMB_CUDA_STATUS_INVALID_PROBLEM:
      return "invalid problem";
    case COMB_CUDA_STATUS_PROBLEM_TOO_LARGE:
      return "problem too large";
    case COMB_CUDA_STATUS_CUDA_UNAVAILABLE:
      return "CUDA unavailable";
    case COMB_CUDA_STATUS_CUDA_ERROR:
      return "CUDA runtime error";
    case COMB_CUDA_STATUS_OUT_OF_MEMORY:
      return "out of memory";
    case COMB_CUDA_STATUS_NOT_IMPLEMENTED:
      return "CUDA scientific pipeline not implemented";
    case COMB_CUDA_STATUS_INTERNAL_ERROR:
      return "internal error";
    default:
      return "unknown CUDA status";
  }
}

#endif  // COMB_ENABLE_CUDA_BACKEND
