/* Test-only CUDA API adapter. Runs the same allocation-free candidate code on
 * the host to test distributed scientific composition and launcher behavior.
 * This file must never be linked into the production CUDA library or launcher.
 */
#include "comb_cuda_pipeline.h"
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

struct comb_cuda_context {
  uint8_t *blob;
  uint8_t *prepared;
  void *scratch;
};

comb_cuda_status_t comb_cuda_device_pci_bus_id(int32_t device, char *buffer,
                                               size_t size) {
  int rank = -1;
  if (device < 0 || buffer == NULL || size < 32)
    return COMB_CUDA_STATUS_INVALID_ARGUMENT;
  if (MPI_Comm_rank(MPI_COMM_WORLD, &rank) != MPI_SUCCESS)
    return COMB_CUDA_STATUS_INTERNAL_ERROR;
  snprintf(buffer, size, "host-test-rank-%d", rank);
  return COMB_CUDA_STATUS_SUCCESS;
}

void comb_cuda_context_destroy(comb_cuda_context_t *context) {
  if (context == NULL) return;
  free(context->blob);
  free(context->prepared);
  free(context->scratch);
  free(context);
}

comb_cuda_status_t comb_cuda_context_create_from_flat_blob(
    const void *blob, uint64_t blob_size, const comb_cuda_options_t *options,
    comb_cuda_context_t **context) {
  comb_cuda_status_t status;
  uint64_t prepared_size = 0;
  uint64_t scratch_bytes = 0;
  comb_cuda_context_t *created;
  (void)options;
  if (context == NULL) return COMB_CUDA_STATUS_INVALID_ARGUMENT;
  *context = NULL;
  status = comb_cuda_flat_problem_validate_science(blob, blob_size);
  if (status != COMB_CUDA_STATUS_SUCCESS) return status;
  created = calloc(1, sizeof(*created));
  if (created == NULL) return COMB_CUDA_STATUS_OUT_OF_MEMORY;
  created->blob = malloc((size_t)blob_size);
  if (created->blob == NULL) {
    comb_cuda_context_destroy(created);
    return COMB_CUDA_STATUS_OUT_OF_MEMORY;
  }
  memcpy(created->blob, blob, (size_t)blob_size);
  status = comb_cuda_prepare(blob, blob_size, (void **)&created->prepared,
                             &prepared_size, &scratch_bytes);
  if (status == COMB_CUDA_STATUS_SUCCESS) {
    created->scratch = calloc(1, (size_t)scratch_bytes);
    if (created->scratch == NULL) status = COMB_CUDA_STATUS_OUT_OF_MEMORY;
  }
  if (status != COMB_CUDA_STATUS_SUCCESS) {
    comb_cuda_context_destroy(created);
    return status;
  }
  *context = created;
  return COMB_CUDA_STATUS_SUCCESS;
}

comb_cuda_status_t comb_cuda_verify_range(
    comb_cuda_context_t *context, uint64_t begin, uint64_t end,
    comb_cuda_witness_mode_t mode, comb_cuda_result_t *result,
    comb_cuda_witness_t *witness) {
  const comb_cuda_flat_problem_t *flat;
  uint32_t *locations = NULL;
  uint32_t *faults = NULL;
  uint64_t minimum = COMB_CUDA_RANK_NONE;
  (void)mode;
  if (context == NULL || result == NULL) return COMB_CUDA_STATUS_INVALID_ARGUMENT;
  flat = (const comb_cuda_flat_problem_t *)context->blob;
  if (begin > end || end > flat->total_combination_count)
    return COMB_CUDA_STATUS_INVALID_ARGUMENT;
  for (uint64_t rank = begin; rank < end; ++rank) {
    int bad = 0;
    const int status = comb_gpu_candidate(context->blob, context->prepared,
                                          rank, context->scratch, &bad);
    if (status != COMB_CUDA_STATUS_SUCCESS) return (comb_cuda_status_t)status;
    if (bad && minimum == COMB_CUDA_RANK_NONE) minimum = rank;
  }
  result->verification_code = minimum == COMB_CUDA_RANK_NONE
                                  ? ALL_FAULTS_GOOD : BAD_FAULT_FOUND;
  result->total_combination_count = flat->total_combination_count;
  result->checked_combination_count = end - begin;
  result->minimum_bad_rank = minimum;
  if (witness != NULL) {
    witness->rank = minimum;
    witness->fault_count = minimum == COMB_CUDA_RANK_NONE ? 0 : flat->s;
    if (minimum != COMB_CUDA_RANK_NONE) {
      locations = calloc(flat->s, sizeof(*locations));
      faults = calloc(flat->s, sizeof(*faults));
      if (locations == NULL || faults == NULL) {
        free(locations);
        free(faults);
        return COMB_CUDA_STATUS_OUT_OF_MEMORY;
      }
      if (!comb_gpu_unrank(context->blob, minimum, locations, faults)) {
        free(locations);
        free(faults);
        return COMB_CUDA_STATUS_INTERNAL_ERROR;
      }
      memcpy(witness->fault_indices, faults, flat->s * sizeof(*faults));
      free(locations);
      free(faults);
    }
  }
  return COMB_CUDA_STATUS_SUCCESS;
}

const char *comb_cuda_status_string(comb_cuda_status_t status) {
  switch (status) {
    case COMB_CUDA_STATUS_SUCCESS: return "success";
    case COMB_CUDA_STATUS_INVALID_ARGUMENT: return "invalid argument";
    case COMB_CUDA_STATUS_INVALID_PROBLEM: return "invalid problem";
    case COMB_CUDA_STATUS_OUT_OF_MEMORY: return "out of memory";
    default: return "host test backend error";
  }
}
