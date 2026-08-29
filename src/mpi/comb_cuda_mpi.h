#ifndef COMB_CUDA_MPI_H
#define COMB_CUDA_MPI_H

/**
 * \file comb_cuda_mpi.h
 *
 * MPI orchestration for the candidate-parallel CUDA verifier.  The CUDA core
 * remains MPI-independent: this layer broadcasts one immutable flat problem,
 * creates one local CUDA context per rank, and schedules disjoint rank ranges.
 */

#include <mpi.h>
#include <stddef.h>
#include <stdint.h>

#include "comb_cuda.h"

#ifdef __cplusplus
extern "C" {
#endif

#define COMB_CUDA_MPI_ABI_VERSION 1U
#define COMB_CUDA_MPI_CHECKPOINT_VERSION 1U

/** Orchestration-layer outcome.  The scientific decision is in the result. */
typedef enum {
  COMB_CUDA_MPI_STATUS_SUCCESS = 0,
  COMB_CUDA_MPI_STATUS_INVALID_ARGUMENT = 1,
  COMB_CUDA_MPI_STATUS_MPI_NOT_INITIALIZED = 2,
  COMB_CUDA_MPI_STATUS_MPI_ERROR = 3,
  COMB_CUDA_MPI_STATUS_DEVICE_MAPPING_ERROR = 4,
  COMB_CUDA_MPI_STATUS_CUDA_ERROR = 5,
  COMB_CUDA_MPI_STATUS_CHECKPOINT_IO_ERROR = 6,
  COMB_CUDA_MPI_STATUS_CHECKPOINT_MISMATCH = 7,
  COMB_CUDA_MPI_STATUS_OUT_OF_MEMORY = 8,
  COMB_CUDA_MPI_STATUS_INTERNAL_ERROR = 9
} comb_cuda_mpi_status_t;

/**
 * Root-authoritative distributed execution options.
 *
 * Only rank zero's value is read.  In particular, checkpoint_path is never
 * dereferenced by another rank.  cuda_options.cuda_device must be -1 because
 * this layer assigns devices from each rank's shared-memory-local rank.
 */
typedef struct {
  uint32_t struct_size;
  uint32_t abi_version;
  comb_cuda_options_t cuda_options;
  /** Globally contiguous checkpoint and deterministic ordering unit. */
  uint64_t epoch_size;
  /** Fine cyclic scheduling unit within an epoch. */
  uint64_t chunk_size;
  /** comb_cuda_witness_mode_t; DEFAULT becomes DETERMINISTIC. */
  uint32_t witness_mode;
  /** Per-rank chunks processed between fast-mode nonblocking reductions. */
  uint32_t fast_poll_chunks;
  /** Read checkpoint_path before starting when nonzero. */
  uint32_t resume_from_checkpoint;
  uint32_t reserved0;
  /**
   * Optional checkpoint path.  Rank zero atomically rewrites it after every
   * completely good epoch.  A resumed checkpoint must match the flat blob.
   */
  const char *checkpoint_path;
  uint64_t reserved[4];
} comb_cuda_mpi_options_t;

/** Global result plus this process's MPI/GPU placement metadata. */
typedef struct {
  uint32_t struct_size;
  uint32_t abi_version;
  comb_code_t verification_code;
  comb_cuda_status_t cuda_status;
  comb_cuda_mpi_status_t orchestration_status;
  int32_t mpi_error_code;

  uint64_t total_combination_count;
  uint64_t checked_combination_count;
  uint64_t minimum_bad_rank;
  uint64_t start_rank;
  uint64_t next_epoch;
  uint64_t problem_hash;

  int32_t world_rank;
  int32_t world_size;
  int32_t local_rank;
  int32_t local_size;
  int32_t cuda_device;
  uint32_t reserved0;
  double elapsed_milliseconds;
} comb_cuda_mpi_result_t;

#define COMB_CUDA_MPI_OPTIONS_INITIALIZER                                 \
  {                                                                       \
    sizeof(comb_cuda_mpi_options_t), COMB_CUDA_MPI_ABI_VERSION,           \
        COMB_CUDA_OPTIONS_INITIALIZER, 0, 0,                               \
        COMB_CUDA_WITNESS_DETERMINISTIC, 1, 0, 0, NULL, {0, 0, 0, 0}     \
  }

#define COMB_CUDA_MPI_RESULT_INITIALIZER                                  \
  {                                                                       \
    sizeof(comb_cuda_mpi_result_t), COMB_CUDA_MPI_ABI_VERSION,            \
        OTHER_ERROR, COMB_CUDA_STATUS_SUCCESS,                            \
        COMB_CUDA_MPI_STATUS_SUCCESS, MPI_SUCCESS, 0, 0,                  \
        COMB_CUDA_RANK_NONE, 0, 0, 0, -1, 0, -1, 0, -1, 0, 0.0          \
  }

/**
 * Verify one problem collectively on communicator.
 *
 * MPI must already be initialized and every rank in communicator must call
 * this function.  problem, options, and witness are read or written only on
 * communicator rank zero; problem is non-owning and may be NULL elsewhere.
 * result may be supplied on every rank and receives identical global counts
 * together with rank-local placement metadata.
 */
comb_cuda_mpi_status_t comb_cuda_mpi_verify_comm(
    MPI_Comm communicator, const comb_problem_t *problem,
    const comb_cuda_mpi_options_t *options, comb_cuda_mpi_result_t *result,
    comb_cuda_witness_t *witness);

/**
 * Convenience wrapper around MPI_COMM_WORLD.
 *
 * If MPI is not initialized, this function initializes it with
 * MPI_THREAD_FUNNELED and finalizes it before returning.  Applications that
 * load a problem only on rank zero should normally own MPI initialization,
 * perform that loading, and call comb_cuda_mpi_verify_comm directly.
 */
comb_cuda_mpi_status_t comb_cuda_mpi_run(
    int *argc, char ***argv, const comb_problem_t *problem,
    const comb_cuda_mpi_options_t *options, comb_cuda_mpi_result_t *result,
    comb_cuda_witness_t *witness);

const char *comb_cuda_mpi_status_string(comb_cuda_mpi_status_t status);

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // COMB_CUDA_MPI_H
