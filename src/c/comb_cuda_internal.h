#ifndef COMB_CUDA_INTERNAL_H
#define COMB_CUDA_INTERNAL_H

#include <stddef.h>
#include <stdint.h>

#include "comb_cuda.h"

/*
 * Pointer-free descriptor copied at byte zero of the immutable device blob.
 * Every nonzero data member named *_offset is a byte offset from the blob
 * base.  Rows use dense uint64_t words with bit i in word i / 64; the final
 * word of every row is masked to its declared logical width.
 */

#define COMB_CUDA_FLAT_MAGIC UINT64_C(0x434f4d42464c4154) /* "COMBFLAT" */
#define COMB_CUDA_FLAT_VERSION 1U

typedef struct {
  uint64_t magic;
  uint32_t version;
  uint32_t header_size;
  uint64_t blob_size;

  uint32_t n_all;
  uint32_t n;
  uint32_t n_x;
  uint32_t n_z;
  uint32_t k;
  uint32_t words_all;
  uint32_t words_code;
  uint32_t fault_location_count;
  uint32_t s;
  uint32_t rank_dp_rows;
  uint32_t rank_dp_columns;
  uint32_t output_error_count;
  uint32_t x_measurement_count;
  uint32_t z_measurement_count;
  uint32_t reserved0;

  uint64_t total_combination_count;
  uint64_t x_lookup_bit_count;
  uint64_t z_lookup_bit_count;

  uint64_t fault_partitions_offset;
  uint64_t fault_prefixes_offset;
  /** Checked W[m][r] table, row-major, for weighted colex rank/unrank. */
  uint64_t rank_dp_offset;
  uint64_t x_stabilisers_offset;
  uint64_t z_stabilisers_offset;
  uint64_t x_lookup_offset;
  uint64_t z_lookup_offset;

  /* Structure-of-arrays DCP storage. */
  uint64_t logical_dcp_f_offset;
  uint64_t logical_dcp_hx_offset;
  uint64_t logical_dcp_hz_offset;
  uint64_t output_dcp_f_offset;
  uint64_t output_dcp_hx_offset;
  uint64_t output_dcp_hz_offset;

  uint64_t x_measurements_offset;
  uint64_t z_measurements_offset;
  uint64_t reserved[8];
} comb_cuda_flat_problem_t;

typedef struct {
  uint64_t next_rank;
  uint64_t rank_end;
  uint64_t scheduled_count;
  uint64_t minimum_bad_rank;
  uint32_t stop_requested;
  uint32_t pipeline_incomplete;
  uint32_t runtime_error;
  uint32_t reserved;
} comb_cuda_scheduler_state_t;

#if defined(__cplusplus)
static_assert(sizeof(comb_cuda_flat_problem_t) == 296,
              "FlatProblem v1 layout changed");
static_assert(offsetof(comb_cuda_flat_problem_t, total_combination_count) == 88,
              "FlatProblem v1 offset changed");
static_assert(sizeof(comb_cuda_scheduler_state_t) == 48,
              "CUDA scheduler layout changed");
#else
_Static_assert(sizeof(comb_cuda_flat_problem_t) == 296,
               "FlatProblem v1 layout changed");
_Static_assert(offsetof(comb_cuda_flat_problem_t, total_combination_count) ==
                   88,
               "FlatProblem v1 offset changed");
_Static_assert(sizeof(comb_cuda_scheduler_state_t) == 48,
               "CUDA scheduler layout changed");
#endif

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Internal MPI-facing flat-blob API.  Rank zero may pack once, broadcast the
 * returned bytes, and have each rank validate and upload the same immutable
 * representation.  The allocation returned by pack is owned by the caller.
 */
comb_cuda_status_t comb_cuda_flat_problem_pack(const comb_problem_t *problem,
                                               void **blob,
                                               uint64_t *blob_size);
comb_cuda_status_t comb_cuda_flat_problem_validate(const void *blob,
                                                   uint64_t blob_size);
/* Validate the scientific contract required before CPU preprocessing or GPU
 * upload: CSS dimensions/commutation, RREF stabilisers, symmetric DCPs, and
 * supported measurement dimensions. */
comb_cuda_status_t comb_cuda_flat_problem_validate_science(
    const void *blob, uint64_t blob_size);
void comb_cuda_flat_problem_free(void *blob);

comb_cuda_status_t comb_cuda_context_create_from_flat_blob(
    const void *blob, uint64_t blob_size, const comb_cuda_options_t *options,
    comb_cuda_context_t **context);

/*
 * Return the selected runtime device's physical PCI bus identifier.  The MPI
 * layer uses this after applying CUDA_VISIBLE_DEVICES so isolated ranks that
 * all see ordinal zero can still prove that they own different GPUs.
 */
comb_cuda_status_t comb_cuda_device_pci_bus_id(int32_t device, char *buffer,
                                               size_t buffer_size);

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // COMB_CUDA_INTERNAL_H
