#define _POSIX_C_SOURCE 200809L

#include "comb_cuda_internal.h"
#include "comb_cuda_mpi.h"

#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#define CHECK(condition)                                                       \
  do {                                                                         \
    if (!(condition)) {                                                        \
      fprintf(stderr, "check failed at %s:%d: %s\n", __FILE__, __LINE__,     \
              #condition);                                                     \
      exit(EXIT_FAILURE);                                                      \
    }                                                                          \
  } while (0)

enum { MOCK_MAX_BAD_RANKS = 4 };

typedef struct {
  uint64_t bad_ranks[MOCK_MAX_BAD_RANKS];
  size_t bad_rank_count;
  unsigned identity_calls;
  unsigned context_creations;
  unsigned context_destructions;
  unsigned verify_calls;
  int32_t identity_device;
  int32_t context_device;
} mock_cuda_state_t;

static mock_cuda_state_t mock_cuda;

struct comb_cuda_context {
  uint64_t total_combination_count;
  uint32_t s;
  int32_t device;
};

static void mock_cuda_reset(const uint64_t *bad_ranks, size_t count) {
  CHECK(count <= MOCK_MAX_BAD_RANKS);
  memset(&mock_cuda, 0, sizeof(mock_cuda));
  mock_cuda.identity_device = -1;
  mock_cuda.context_device = -1;
  if (count != 0U) {
    CHECK(bad_ranks != NULL);
    memcpy(mock_cuda.bad_ranks, bad_ranks, count * sizeof(*bad_ranks));
  }
  mock_cuda.bad_rank_count = count;
}

static uint64_t mock_cuda_minimum_bad_rank(uint64_t rank_begin,
                                           uint64_t rank_end) {
  uint64_t minimum = COMB_CUDA_RANK_NONE;
  size_t index;
  for (index = 0; index < mock_cuda.bad_rank_count; ++index) {
    const uint64_t candidate = mock_cuda.bad_ranks[index];
    if (candidate >= rank_begin && candidate < rank_end && candidate < minimum) {
      minimum = candidate;
    }
  }
  return minimum;
}

comb_cuda_status_t comb_cuda_device_pci_bus_id(int32_t device, char *buffer,
                                               size_t buffer_size) {
  static const char identity[] = "0000:00:00.0";
  ++mock_cuda.identity_calls;
  mock_cuda.identity_device = device;
  if (device != 0 || buffer == NULL || buffer_size < sizeof(identity)) {
    return COMB_CUDA_STATUS_INVALID_ARGUMENT;
  }
  memcpy(buffer, identity, sizeof(identity));
  return COMB_CUDA_STATUS_SUCCESS;
}

comb_cuda_status_t comb_cuda_context_create_from_flat_blob(
    const void *blob, uint64_t blob_size, const comb_cuda_options_t *options,
    comb_cuda_context_t **context) {
  const comb_cuda_flat_problem_t *flat;
  comb_cuda_context_t *created;
  if (context == NULL) return COMB_CUDA_STATUS_INVALID_ARGUMENT;
  *context = NULL;
  if (blob == NULL || options == NULL ||
      options->struct_size < sizeof(*options) ||
      options->abi_version != COMB_CUDA_ABI_VERSION ||
      comb_cuda_flat_problem_validate(blob, blob_size) !=
          COMB_CUDA_STATUS_SUCCESS) {
    return COMB_CUDA_STATUS_INVALID_ARGUMENT;
  }
  flat = (const comb_cuda_flat_problem_t *)blob;
  created = (comb_cuda_context_t *)calloc(1, sizeof(*created));
  if (created == NULL) return COMB_CUDA_STATUS_OUT_OF_MEMORY;
  created->total_combination_count = flat->total_combination_count;
  created->s = flat->s;
  created->device = options->cuda_device;
  ++mock_cuda.context_creations;
  mock_cuda.context_device = options->cuda_device;
  *context = created;
  return COMB_CUDA_STATUS_SUCCESS;
}

comb_cuda_status_t comb_cuda_verify_range(
    comb_cuda_context_t *context, uint64_t rank_begin, uint64_t rank_end,
    comb_cuda_witness_mode_t mode, comb_cuda_result_t *result,
    comb_cuda_witness_t *witness) {
  uint64_t bad_rank;
  if (context == NULL || result == NULL ||
      result->struct_size < sizeof(*result) ||
      result->abi_version != COMB_CUDA_ABI_VERSION || rank_begin > rank_end ||
      rank_end > context->total_combination_count ||
      (mode != COMB_CUDA_WITNESS_DETERMINISTIC &&
       mode != COMB_CUDA_WITNESS_FAST_ANY) ||
      witness == NULL || witness->struct_size < sizeof(*witness) ||
      witness->abi_version != COMB_CUDA_ABI_VERSION ||
      witness->fault_capacity < context->s ||
      (context->s != 0U && witness->fault_indices == NULL)) {
    return COMB_CUDA_STATUS_INVALID_ARGUMENT;
  }

  ++mock_cuda.verify_calls;
  result->verification_code = ALL_FAULTS_GOOD;
  result->total_combination_count = context->total_combination_count;
  result->checked_combination_count = rank_end - rank_begin;
  result->minimum_bad_rank = COMB_CUDA_RANK_NONE;
  result->setup_milliseconds = 0.0;
  result->verification_milliseconds = 0.0;
  witness->rank = COMB_CUDA_RANK_NONE;
  witness->fault_count = 0;

  bad_rank = mock_cuda_minimum_bad_rank(rank_begin, rank_end);
  if (bad_rank == COMB_CUDA_RANK_NONE) return COMB_CUDA_STATUS_SUCCESS;

  result->verification_code = BAD_FAULT_FOUND;
  result->minimum_bad_rank = bad_rank;
  if (mode == COMB_CUDA_WITNESS_FAST_ANY) {
    result->checked_combination_count = bad_rank - rank_begin + 1U;
  }
  witness->rank = bad_rank;
  witness->fault_count = context->s;
  if (context->s != 0U) {
    uint32_t index;
    CHECK(bad_rank <= UINT32_MAX);
    for (index = 0; index < context->s; ++index) {
      witness->fault_indices[index] = (uint32_t)bad_rank;
    }
  }
  return COMB_CUDA_STATUS_SUCCESS;
}

void comb_cuda_context_destroy(comb_cuda_context_t *context) {
  if (context != NULL) {
    ++mock_cuda.context_destructions;
    free(context);
  }
}

typedef struct {
  comb_problem_t problem;
  uint32_t partitions[3];
  comb_dcp_t logicals[1];
  comb_dcp_t errors[7];
  comb_bit_word_t logical_f[1];
  comb_bit_word_t logical_hx[1];
  comb_bit_word_t logical_hz[1];
  comb_bit_word_t error_f[7];
  comb_bit_word_t error_hx[7];
  comb_bit_word_t error_hz[7];
  comb_bit_word_t x_lookup[1];
  comb_bit_word_t z_lookup[1];
} fixture_t;

static void fixture_init(fixture_t *fixture) {
  uint32_t index;
  const comb_problem_t empty_problem = COMB_PROBLEM_INITIALIZER;
  CHECK(fixture != NULL);
  memset(fixture, 0, sizeof(*fixture));
  fixture->problem = empty_problem;
  fixture->partitions[0] = 3U;
  fixture->partitions[1] = 2U;
  fixture->partitions[2] = 2U;
  fixture->logicals[0].F.m = 1U;
  fixture->logicals[0].F.n = 1U;
  fixture->logicals[0].F.data = fixture->logical_f;
  fixture->logicals[0].h_x.n = 1U;
  fixture->logicals[0].h_x.data = fixture->logical_hx;
  fixture->logicals[0].h_z.n = 1U;
  fixture->logicals[0].h_z.data = fixture->logical_hz;
  for (index = 0; index < 7U; ++index) {
    fixture->errors[index].F.m = 1U;
    fixture->errors[index].F.n = 1U;
    fixture->errors[index].F.data = &fixture->error_f[index];
    fixture->errors[index].h_x.n = 1U;
    fixture->errors[index].h_x.data = &fixture->error_hx[index];
    fixture->errors[index].h_z.n = 1U;
    fixture->errors[index].h_z.data = &fixture->error_hz[index];
  }
  fixture->problem.n_all = 1U;
  fixture->problem.n = 1U;
  fixture->problem.k = 1U;
  fixture->problem.dcp_stabilisers = fixture->logicals;
  fixture->problem.x_lookup_table = fixture->x_lookup;
  fixture->problem.x_lookup_bit_count = 2U;
  fixture->problem.z_lookup_table = fixture->z_lookup;
  fixture->problem.z_lookup_bit_count = 2U;
  fixture->problem.fault_location_count = 3U;
  fixture->problem.fault_location_partitions = fixture->partitions;
  fixture->problem.s = 1U;
  fixture->problem.output_error_count = 7U;
  fixture->problem.output_errors = fixture->errors;
}

static comb_cuda_witness_t fresh_witness(uint32_t *fault) {
  comb_cuda_witness_t witness = COMB_CUDA_WITNESS_INITIALIZER;
  witness.fault_capacity = 1U;
  witness.fault_indices = fault;
  return witness;
}

static void read_checkpoint(const char *path, uint64_t *problem_hash,
                            uint64_t *epoch_size, uint64_t *next_epoch) {
  FILE *stream = fopen(path, "r");
  unsigned version = 0U;
  int matched;
  CHECK(stream != NULL);
  matched = fscanf(stream,
                   "COMB_CUDA_MPI_CHECKPOINT %u\nproblem_hash %" SCNx64
                   "\nepoch_size %" SCNu64 "\nnext_epoch %" SCNu64 "\n",
                   &version, problem_hash, epoch_size, next_epoch);
  CHECK(matched == 4);
  CHECK(version == COMB_CUDA_MPI_CHECKPOINT_VERSION);
  CHECK(fclose(stream) == 0);
}

static void check_common_placement(const comb_cuda_mpi_result_t *result) {
  CHECK(result->world_rank == 0);
  CHECK(result->world_size == 1);
  CHECK(result->local_rank == 0);
  CHECK(result->local_size == 1);
  CHECK(result->cuda_device == 0);
  CHECK(mock_cuda.identity_calls == 1U);
  CHECK(mock_cuda.identity_device == 0);
  CHECK(mock_cuda.context_creations == 1U);
  CHECK(mock_cuda.context_destructions == 1U);
  CHECK(mock_cuda.context_device == 0);
}

static void test_deterministic_bad_and_checkpoint(const comb_problem_t *problem,
                                                  const char *checkpoint_path,
                                                  uint64_t *problem_hash) {
  static const uint64_t bad_ranks[] = {5U, 4U};
  comb_cuda_mpi_options_t options = COMB_CUDA_MPI_OPTIONS_INITIALIZER;
  comb_cuda_mpi_result_t result = COMB_CUDA_MPI_RESULT_INITIALIZER;
  uint32_t fault = UINT32_MAX;
  comb_cuda_witness_t witness = fresh_witness(&fault);
  uint64_t checkpoint_hash = 0;
  uint64_t checkpoint_epoch_size = 0;
  uint64_t checkpoint_next_epoch = 0;

  mock_cuda_reset(bad_ranks, sizeof(bad_ranks) / sizeof(bad_ranks[0]));
  options.epoch_size = 3U;
  options.chunk_size = 2U;
  options.witness_mode = COMB_CUDA_WITNESS_DETERMINISTIC;
  options.checkpoint_path = checkpoint_path;
  CHECK(comb_cuda_mpi_verify_comm(MPI_COMM_WORLD, problem, &options, &result,
                                  &witness) ==
        COMB_CUDA_MPI_STATUS_SUCCESS);
  CHECK(result.verification_code == BAD_FAULT_FOUND);
  CHECK(result.cuda_status == COMB_CUDA_STATUS_SUCCESS);
  CHECK(result.orchestration_status == COMB_CUDA_MPI_STATUS_SUCCESS);
  CHECK(result.total_combination_count == 7U);
  CHECK(result.checked_combination_count == 6U);
  CHECK(result.minimum_bad_rank == 4U);
  CHECK(result.start_rank == 0U);
  CHECK(result.next_epoch == 3U);
  CHECK(witness.rank == 4U);
  CHECK(witness.fault_count == 1U);
  CHECK(fault == 4U);
  CHECK(mock_cuda.verify_calls == 4U);
  check_common_placement(&result);

  read_checkpoint(checkpoint_path, &checkpoint_hash, &checkpoint_epoch_size,
                  &checkpoint_next_epoch);
  CHECK(checkpoint_hash == result.problem_hash);
  CHECK(checkpoint_epoch_size == 3U);
  CHECK(checkpoint_next_epoch == 3U);
  *problem_hash = result.problem_hash;
}

static void test_fast_any_witness(const comb_problem_t *problem) {
  static const uint64_t bad_ranks[] = {4U};
  comb_cuda_mpi_options_t options = COMB_CUDA_MPI_OPTIONS_INITIALIZER;
  comb_cuda_mpi_result_t result = COMB_CUDA_MPI_RESULT_INITIALIZER;
  uint32_t fault = UINT32_MAX;
  comb_cuda_witness_t witness = fresh_witness(&fault);

  mock_cuda_reset(bad_ranks, sizeof(bad_ranks) / sizeof(bad_ranks[0]));
  options.epoch_size = 7U;
  options.chunk_size = 2U;
  options.witness_mode = COMB_CUDA_WITNESS_FAST_ANY;
  options.fast_poll_chunks = 1U;
  CHECK(comb_cuda_mpi_verify_comm(MPI_COMM_WORLD, problem, &options, &result,
                                  &witness) ==
        COMB_CUDA_MPI_STATUS_SUCCESS);
  CHECK(result.verification_code == BAD_FAULT_FOUND);
  CHECK(result.checked_combination_count == 6U);
  CHECK(result.minimum_bad_rank == 4U);
  CHECK(result.start_rank == 0U);
  CHECK(result.next_epoch == 0U);
  CHECK(witness.rank == 4U);
  CHECK(witness.fault_count == 1U);
  CHECK(fault == 4U);
  CHECK(mock_cuda.verify_calls == 4U);
  check_common_placement(&result);
}

static void test_all_good_resume(const comb_problem_t *problem,
                                 const char *checkpoint_path,
                                 uint64_t expected_hash) {
  comb_cuda_mpi_options_t options = COMB_CUDA_MPI_OPTIONS_INITIALIZER;
  comb_cuda_mpi_result_t result = COMB_CUDA_MPI_RESULT_INITIALIZER;
  uint32_t fault = 6U;
  comb_cuda_witness_t witness = fresh_witness(&fault);
  uint64_t checkpoint_hash = 0;
  uint64_t checkpoint_epoch_size = 0;
  uint64_t checkpoint_next_epoch = 0;

  mock_cuda_reset(NULL, 0U);
  /* Reusing a caller-owned witness after a bad result must not leak stale
   * witness metadata into an all-good resumed result.  The unused tuple buffer
   * itself remains caller-owned and need not be overwritten. */
  witness.rank = 6U;
  witness.fault_count = 1U;
  options.epoch_size = 3U;
  options.chunk_size = 2U;
  options.witness_mode = COMB_CUDA_WITNESS_DETERMINISTIC;
  options.resume_from_checkpoint = 1U;
  options.checkpoint_path = checkpoint_path;
  CHECK(comb_cuda_mpi_verify_comm(MPI_COMM_WORLD, problem, &options, &result,
                                  &witness) ==
        COMB_CUDA_MPI_STATUS_SUCCESS);
  CHECK(result.verification_code == ALL_FAULTS_GOOD);
  CHECK(result.checked_combination_count == 4U);
  CHECK(result.minimum_bad_rank == COMB_CUDA_RANK_NONE);
  CHECK(result.start_rank == 3U);
  CHECK(result.next_epoch == 7U);
  CHECK(result.problem_hash == expected_hash);
  CHECK(witness.rank == COMB_CUDA_RANK_NONE);
  CHECK(witness.fault_count == 0U);
  CHECK(fault == 6U);
  CHECK(mock_cuda.verify_calls == 3U);
  check_common_placement(&result);

  read_checkpoint(checkpoint_path, &checkpoint_hash, &checkpoint_epoch_size,
                  &checkpoint_next_epoch);
  CHECK(checkpoint_hash == expected_hash);
  CHECK(checkpoint_epoch_size == 3U);
  CHECK(checkpoint_next_epoch == 7U);
}

int main(int argc, char **argv) {
  fixture_t fixture;
  char checkpoint_path[] = "/tmp/comb-cuda-mpi-single-rank-XXXXXX";
  uint64_t problem_hash = 0;
  int descriptor;
  int world_size = 0;

  fixture_init(&fixture);
  CHECK(setenv("CUDA_VISIBLE_DEVICES", "comb-mock-device", 1) == 0);
  CHECK(MPI_Init(&argc, &argv) == MPI_SUCCESS);
  CHECK(MPI_Comm_size(MPI_COMM_WORLD, &world_size) == MPI_SUCCESS);
  CHECK(world_size == 1);

  descriptor = mkstemp(checkpoint_path);
  CHECK(descriptor >= 0);
  CHECK(close(descriptor) == 0);
  CHECK(unlink(checkpoint_path) == 0);

  test_deterministic_bad_and_checkpoint(&fixture.problem, checkpoint_path,
                                        &problem_hash);
  test_fast_any_witness(&fixture.problem);
  test_all_good_resume(&fixture.problem, checkpoint_path, problem_hash);

  CHECK(unlink(checkpoint_path) == 0);
  CHECK(MPI_Finalize() == MPI_SUCCESS);
  puts("single-rank MPI orchestration mock tests passed");
  return EXIT_SUCCESS;
}
