#include <stdlib.h>
#include <string.h>

#include "comb_cuda_pipeline.h"

static int reserve_bytes(uint64_t *size, uint64_t count, uint64_t element,
                         uint64_t *offset) {
  uint64_t aligned;
  if (*size > UINT64_MAX - 7) return 0;
  aligned = (*size + 7) & ~UINT64_C(7);
  if (count > (UINT64_MAX - aligned) / element) return 0;
  *offset = aligned;
  *size = aligned + count * element;
  return *size <= SIZE_MAX;
}

static void extract_columns(const uint64_t *rows, uint32_t count, uint32_t n,
                            uint32_t words, uint32_t *pivots,
                            uint32_t *free_columns) {
  uint32_t pivot = 0, free_count = 0;
  for (uint32_t row = 0; row < count; ++row) {
    uint32_t column = 0;
    while (!comb_gpu_bit(rows + (uint64_t)row * words, column)) ++column;
    pivots[row] = column;
  }
  for (uint32_t column = 0; column < n; ++column) {
    if (pivot < count && pivots[pivot] == column)
      ++pivot;
    else
      free_columns[free_count++] = column;
  }
}

comb_cuda_status_t comb_cuda_prepare(const void *blob, uint64_t blob_size,
                                     void **prepared, uint64_t *prepared_size,
                                     uint64_t *scratch_bytes) {
  const uint8_t *bytes = (const uint8_t *)blob;
  const comb_cuda_flat_problem_t *p = (const comb_cuda_flat_problem_t *)blob;
  comb_gpu_prepared_t q = {0};
  uint64_t size = sizeof(q), scratch = 0, ignored;
  uint8_t *storage;
  uint32_t *z_pivots;
  uint64_t *z_sequence, *x_sequence;
  comb_cuda_status_t status;
  if (!prepared || !prepared_size || !scratch_bytes)
    return COMB_CUDA_STATUS_INVALID_ARGUMENT;
  *prepared = NULL;
  *prepared_size = 0;
  *scratch_bytes = 0;
  status = comb_cuda_flat_problem_validate_science(blob, blob_size);
  if (status != COMB_CUDA_STATUS_SUCCESS) return status;
  if (p->n_all == UINT32_MAX ||
      (uint64_t)p->z_measurement_count + p->n_z + p->x_measurement_count >
          UINT32_MAX)
    return COMB_CUDA_STATUS_PROBLEM_TOO_LARGE;
  q.x_free_count = p->n - p->n_x;
  q.z_free_count = p->n - p->n_z;
  q.z_rows = p->z_measurement_count + p->n_z;
  q.x_rows = q.z_rows + p->x_measurement_count;
  if (!reserve_bytes(&size, p->n_x, 4, &q.x_pivots_offset) ||
      !reserve_bytes(&size, p->n_z, 4, &q.z_pivots_offset) ||
      !reserve_bytes(&size, q.x_free_count, 4, &q.x_free_offset) ||
      !reserve_bytes(&size, q.z_free_count, 4, &q.z_free_offset) ||
      !reserve_bytes(&size, (uint64_t)q.z_rows * p->words_all, 8,
                     &q.z_sequence_offset) ||
      !reserve_bytes(&size, (uint64_t)q.x_rows * p->words_all, 8,
                     &q.x_sequence_offset) ||
      !reserve_bytes(&scratch, p->s, 8, &ignored) ||
      !reserve_bytes(&scratch, 4 * ((uint64_t)p->n_all + 2) * p->words_all, 8,
                     &ignored) ||
      !reserve_bytes(&scratch, 2 * ((uint64_t)p->n_all + 1) * p->words_all, 8,
                     &ignored) ||
      !reserve_bytes(&scratch, 2 * ((uint64_t)p->n_all + 1), 1, &ignored) ||
      !reserve_bytes(&scratch, p->words_all, 8, &ignored) ||
      !reserve_bytes(&scratch, p->n_all, 8, &ignored))
    return COMB_CUDA_STATUS_PROBLEM_TOO_LARGE;
  q.scratch_bytes = scratch;
  storage = (uint8_t *)calloc(1, (size_t)size);
  if (!storage) return COMB_CUDA_STATUS_OUT_OF_MEMORY;
  memcpy(storage, &q, sizeof(q));
  extract_columns((const uint64_t *)(bytes + p->x_stabilisers_offset), p->n_x,
                  p->n, p->words_code,
                  (uint32_t *)(storage + q.x_pivots_offset),
                  (uint32_t *)(storage + q.x_free_offset));
  z_pivots = (uint32_t *)(storage + q.z_pivots_offset);
  extract_columns((const uint64_t *)(bytes + p->z_stabilisers_offset), p->n_z,
                  p->n, p->words_code, z_pivots,
                  (uint32_t *)(storage + q.z_free_offset));
  z_sequence = (uint64_t *)(storage + q.z_sequence_offset);
  x_sequence = (uint64_t *)(storage + q.x_sequence_offset);
  for (uint32_t row = 0; row < p->z_measurement_count; ++row) {
    const uint64_t *source =
        (const uint64_t *)(bytes + p->z_measurements_offset) +
        (uint64_t)row * p->words_all;
    memcpy(z_sequence + (uint64_t)row * p->words_all, source,
           (size_t)p->words_all * 8);
    memcpy(x_sequence + (uint64_t)row * p->words_all, source,
           (size_t)p->words_all * 8);
  }
  for (uint32_t row = 0; row < p->n_z; ++row) {
    const uint32_t dest = p->z_measurement_count + row;
    memcpy(z_sequence + (uint64_t)dest * p->words_all,
           (const uint64_t *)(bytes + p->z_stabilisers_offset) +
               (uint64_t)row * p->words_code,
           (size_t)p->words_code * 8);
    comb_gpu_set(x_sequence + (uint64_t)dest * p->words_all, z_pivots[row], 1);
  }
  for (uint32_t row = 0; row < p->x_measurement_count; ++row)
    memcpy(x_sequence + ((uint64_t)q.z_rows + row) * p->words_all,
           (const uint64_t *)(bytes + p->x_measurements_offset) +
               (uint64_t)row * p->words_all,
           (size_t)p->words_all * 8);
  *prepared = storage;
  *prepared_size = size;
  *scratch_bytes = scratch;
  return COMB_CUDA_STATUS_SUCCESS;
}
