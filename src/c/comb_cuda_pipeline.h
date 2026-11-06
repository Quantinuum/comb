#ifndef COMB_CUDA_PIPELINE_H
#define COMB_CUDA_PIPELINE_H

/* Candidate mathematics shared by the actual CUDA kernel and host differential
 * tests. No allocation, CUDA runtime calls, or CPU fallback occurs here. */
#include "comb_cuda_internal.h"
#include "comb_cuda_math.h"

typedef struct {
  uint64_t x_pivots_offset, z_pivots_offset;
  uint64_t x_free_offset, z_free_offset;
  uint64_t z_sequence_offset, x_sequence_offset;
  uint64_t scratch_bytes;
  uint32_t x_free_count, z_free_count, z_rows, x_rows;
} comb_gpu_prepared_t;

#ifdef __cplusplus
extern "C" {
#endif
comb_cuda_status_t comb_cuda_prepare(const void *blob, uint64_t blob_size,
                                     void **prepared, uint64_t *prepared_size,
                                     uint64_t *scratch_bytes);
#ifdef __cplusplus
}
#endif

#if defined(__CUDACC__) || defined(__CUDA__)
#define COMB_PIPELINE_HD static __host__ __device__ inline
#else
#define COMB_PIPELINE_HD static inline
#endif

COMB_PIPELINE_HD int comb_gpu_unrank(const uint8_t *blob, uint64_t rank,
                                     uint32_t *locations, uint32_t *faults) {
  const comb_cuda_flat_problem_t *p = (const comb_cuda_flat_problem_t *)blob;
  const uint32_t *partitions =
      (const uint32_t *)(blob + p->fault_partitions_offset);
  const uint64_t *prefixes =
      (const uint64_t *)(blob + p->fault_prefixes_offset);
  const uint64_t *dp = (const uint64_t *)(blob + p->rank_dp_offset);
  uint32_t m = p->fault_location_count, r = p->s;
  uint64_t probe = rank, block_begin = 0, local;
  if (rank >= p->total_combination_count) return 0;
  while (r != 0) {
    uint64_t without_last;
    if (m == 0) return 0;
    without_last = dp[(uint64_t)(m - 1) * p->rank_dp_columns + r];
    if (probe < without_last) {
      --m;
      continue;
    }
    locations[r - 1] = m - 1;
    probe = (probe - without_last) / partitions[m - 1];
    --m;
    --r;
  }
  for (uint32_t i = 0; i < p->s; ++i) {
    const uint32_t location = locations[i];
    block_begin = dp[(uint64_t)location * p->rank_dp_columns + i + 1] +
                  (uint64_t)partitions[location] * block_begin;
  }
  local = rank - block_begin;
  for (uint32_t i = 0; i < p->s; ++i) {
    const uint32_t location = locations[i];
    faults[i] = (uint32_t)(prefixes[location] + local % partitions[location]);
    local /= partitions[location];
  }
  return local == 0;
}

COMB_PIPELINE_HD comb_gpu_dcp_t comb_gpu_dcp_view(uint64_t *storage, uint32_t n,
                                                  uint32_t w) {
  comb_gpu_dcp_t result;
  result.qubits = n;
  result.words = w;
  result.f = storage;
  result.hx = storage + (uint64_t)n * w;
  result.hz = result.hx + w;
  return result;
}

COMB_PIPELINE_HD comb_gpu_dcp_t comb_gpu_flat_dcp(const uint8_t *blob,
                                                  uint32_t index, int logical) {
  const comb_cuda_flat_problem_t *p = (const comb_cuda_flat_problem_t *)blob;
  comb_gpu_dcp_t result;
  result.qubits = p->n_all;
  result.words = p->words_all;
  result.f = (uint64_t *)(blob + (logical ? p->logical_dcp_f_offset
                                          : p->output_dcp_f_offset)) +
             (uint64_t)index * p->n_all * p->words_all;
  result.hx = (uint64_t *)(blob + (logical ? p->logical_dcp_hx_offset
                                           : p->output_dcp_hx_offset)) +
              (uint64_t)index * p->words_all;
  result.hz = (uint64_t *)(blob + (logical ? p->logical_dcp_hz_offset
                                           : p->output_dcp_hz_offset)) +
              (uint64_t)index * p->words_all;
  return result;
}

COMB_PIPELINE_HD void comb_gpu_copy_dcp(comb_gpu_dcp_t *to,
                                        const comb_gpu_dcp_t *from) {
  for (uint64_t i = 0; i < (uint64_t)to->qubits * to->words; ++i)
    to->f[i] = from->f[i];
  for (uint32_t i = 0; i < to->words; ++i) {
    to->hx[i] = from->hx[i];
    to->hz[i] = from->hz[i];
  }
}

COMB_PIPELINE_HD int comb_gpu_variant_good(const uint8_t *blob,
                                           const uint8_t *prepared,
                                           const comb_gpu_dcp_t *error,
                                           comb_gpu_tableau_t *tableau,
                                           uint64_t *x_error,
                                           uint64_t *generators, int *good) {
  const comb_cuda_flat_problem_t *p = (const comb_cuda_flat_problem_t *)blob;
  const comb_gpu_prepared_t *q = (const comb_gpu_prepared_t *)prepared;
  const uint32_t *x_pivots = (const uint32_t *)(prepared + q->x_pivots_offset);
  const uint32_t *x_free = (const uint32_t *)(prepared + q->x_free_offset);
  const uint32_t *z_free = (const uint32_t *)(prepared + q->z_free_offset);
  const uint64_t *z_sequence =
      (const uint64_t *)(prepared + q->z_sequence_offset);
  const uint64_t *x_sequence =
      (const uint64_t *)(prepared + q->x_sequence_offset);
  const uint64_t *x_stabs = (const uint64_t *)(blob + p->x_stabilisers_offset);
  const uint64_t *x_lookup = (const uint64_t *)(blob + p->x_lookup_offset);
  const uint64_t *z_lookup = (const uint64_t *)(blob + p->z_lookup_offset);
  uint64_t reduced = 0;
  unsigned rejected = 0, coset_good = 0;
  int status;
  *good = 0;
  for (uint32_t w = 0; w < p->words_all; ++w) x_error[w] = error->hx[w];
  for (uint32_t row = 0; row < p->z_measurement_count; ++row) {
    if (comb_gpu_dot(x_error, z_sequence + (uint64_t)row * p->words_all,
                     p->words_all)) {
      *good = 1;
      return COMB_CUDA_STATUS_SUCCESS;
    }
  }
  for (uint32_t row = 0; row < p->n_x; ++row) {
    if (comb_gpu_bit(x_error, x_pivots[row])) {
      for (uint32_t w = 0; w < p->words_code; ++w)
        x_error[w] ^= x_stabs[(uint64_t)row * p->words_code + w];
    }
  }
  for (uint32_t i = 0; i < q->x_free_count; ++i)
    reduced |= (uint64_t)comb_gpu_bit(x_error, x_free[i]) << i;
  if (!comb_gpu_bit(x_lookup, reduced)) return COMB_CUDA_STATUS_SUCCESS;

  comb_gpu_tableau_from_dcp(tableau, error);
  status = comb_gpu_postselect(tableau, z_sequence, q->z_rows, 1, &rejected);
  if (status != 0) return COMB_CUDA_STATUS_INVALID_PROBLEM;
  if (rejected) {
    *good = 1;
    return COMB_CUDA_STATUS_SUCCESS;
  }
  status = comb_gpu_postselect(tableau, x_sequence, q->x_rows, 0, &rejected);
  if (status != 0) return COMB_CUDA_STATUS_INVALID_PROBLEM;
  if (rejected) {
    *good = 1;
    return COMB_CUDA_STATUS_SUCCESS;
  }
  status = comb_gpu_coset_good(tableau, z_free, q->z_free_count, z_lookup,
                               p->z_lookup_bit_count, generators, &coset_good);
  if (status == 3) return COMB_CUDA_STATUS_PROBLEM_TOO_LARGE;
  if (status != 0) return COMB_CUDA_STATUS_INVALID_PROBLEM;
  *good = (int)coset_good;
  return COMB_CUDA_STATUS_SUCCESS;
}

/* One independent candidate owns one disjoint reusable scratch slice. The
 * initial implementation parallelizes candidates, retaining deterministic
 * serial pivot selection within a candidate. */
COMB_PIPELINE_HD int comb_gpu_candidate(const uint8_t *blob,
                                        const uint8_t *prepared, uint64_t rank,
                                        void *scratch, int *bad) {
  const comb_cuda_flat_problem_t *p = (const comb_cuda_flat_problem_t *)blob;
  uint32_t *locations = (uint32_t *)scratch;
  uint32_t *faults = locations + p->s;
  uint64_t *cursor = (uint64_t *)scratch + p->s;
  const uint64_t dcp_words = ((uint64_t)p->n_all + 2) * p->words_all;
  comb_gpu_dcp_t dcps[4], combined, error, temporary, source;
  comb_gpu_tableau_t tableau;
  uint64_t *x_error, *generators;
  if (!comb_gpu_unrank(blob, rank, locations, faults))
    return COMB_CUDA_STATUS_INTERNAL_ERROR;
  for (uint32_t i = 0; i < 4; ++i) {
    dcps[i] = comb_gpu_dcp_view(cursor, p->n_all, p->words_all);
    cursor += dcp_words;
  }
  tableau.qubits = p->n_all;
  tableau.words = p->words_all;
  tableau.rows = p->n_all;
  tableau.x = cursor;
  cursor += ((uint64_t)p->n_all + 1) * p->words_all;
  tableau.z = cursor;
  cursor += ((uint64_t)p->n_all + 1) * p->words_all;
  tableau.delta = (uint8_t *)cursor;
  tableau.epsilon = tableau.delta + (uint64_t)p->n_all + 1;
  cursor += (2 * ((uint64_t)p->n_all + 1) + 7) / 8;
  x_error = cursor;
  generators = cursor + p->words_all;

  combined = dcps[0];
  temporary = dcps[1];
  source = comb_gpu_flat_dcp(blob, faults[0], 0);
  comb_gpu_copy_dcp(&combined, &source);
  for (uint32_t i = 1; i < p->s; ++i) {
    comb_gpu_dcp_t swap;
    source = comb_gpu_flat_dcp(blob, faults[i], 0);
    comb_gpu_dcp_multiply(&temporary, &combined, &source);
    swap = combined;
    combined = temporary;
    temporary = swap;
  }
  *bad = 1;
  for (uint64_t mask = 0; mask < (UINT64_C(1) << p->k); ++mask) {
    int good = 0, status;
    error = dcps[2];
    temporary = dcps[3];
    comb_gpu_copy_dcp(&error, &combined);
    for (uint32_t i = 0; i < p->k; ++i) {
      if ((mask >> i) & 1U) {
        comb_gpu_dcp_t swap;
        source = comb_gpu_flat_dcp(blob, i, 1);
        comb_gpu_dcp_multiply(&temporary, &source, &error);
        swap = error;
        error = temporary;
        temporary = swap;
      }
    }
    status = comb_gpu_variant_good(blob, prepared, &error, &tableau, x_error,
                                   generators, &good);
    if (status != COMB_CUDA_STATUS_SUCCESS) return status;
    if (good) {
      *bad = 0;
      break;
    }
  }
  return COMB_CUDA_STATUS_SUCCESS;
}

#undef COMB_PIPELINE_HD
#endif
