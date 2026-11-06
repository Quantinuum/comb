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
static uint64_t candidate_hits[7];
static int duplicate_device;
static int fail_rank = -1;

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
  int rank = -1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  ++mock_cuda.identity_calls;
  mock_cuda.identity_device = device;
  if (device < 0 || buffer == NULL || buffer_size < 32)
    return COMB_CUDA_STATUS_INVALID_ARGUMENT;
  snprintf(buffer, buffer_size, "0000:%02x:00.0", duplicate_device ? 0 : rank);
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

  {
    int rank = -1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == fail_rank) return COMB_CUDA_STATUS_CUDA_ERROR;
    for (uint64_t i = rank_begin; i < rank_end; ++i) {
      CHECK(i < 7);
      candidate_hits[i]++;
    }
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


static void verify_case(const comb_problem_t *problem,
                        comb_cuda_witness_mode_t mode, int bad, int rank,
                        int world_size) {
  const uint64_t bad_ranks[] = {5, 1};
  uint64_t global_hits[7] = {0};
  uint32_t fault = UINT32_MAX;
  comb_cuda_mpi_options_t options = COMB_CUDA_MPI_OPTIONS_INITIALIZER;
  comb_cuda_mpi_result_t result = COMB_CUDA_MPI_RESULT_INITIALIZER;
  comb_cuda_witness_t witness = COMB_CUDA_WITNESS_INITIALIZER;
  witness.fault_capacity = 1;
  witness.fault_indices = &fault;
  mock_cuda_reset(bad ? bad_ranks : NULL, bad ? 2U : 0U);
  memset(candidate_hits, 0, sizeof(candidate_hits));
  options.epoch_size = 7;
  options.chunk_size = 1;
  options.witness_mode = (uint32_t)mode;
  CHECK(comb_cuda_mpi_verify_comm(
      MPI_COMM_WORLD, rank == 0 ? problem : NULL,
      rank == 0 ? &options : NULL, &result, rank == 0 ? &witness : NULL) ==
        COMB_CUDA_MPI_STATUS_SUCCESS);
  CHECK(result.world_rank == rank && result.world_size == world_size);
  CHECK(result.cuda_device == 0);
  CHECK(result.total_combination_count == 7);
  CHECK(result.verification_code == (bad ? BAD_FAULT_FOUND : ALL_FAULTS_GOOD));
  if (bad) {
    CHECK(result.minimum_bad_rank == 1 || result.minimum_bad_rank == 5);
    if (mode == COMB_CUDA_WITNESS_DETERMINISTIC)
      CHECK(result.minimum_bad_rank == 1);
    if (rank == 0) {
      CHECK(witness.rank == result.minimum_bad_rank);
      CHECK(witness.fault_count == 1 && fault == witness.rank);
    }
  } else {
    CHECK(result.minimum_bad_rank == COMB_CUDA_RANK_NONE);
    if (rank == 0) CHECK(witness.fault_count == 0);
  }
  if (mode == COMB_CUDA_WITNESS_DETERMINISTIC || !bad) {
    CHECK(result.checked_combination_count == 7);
    CHECK(MPI_Allreduce(candidate_hits, global_hits, 7, MPI_UINT64_T,
                        MPI_SUM, MPI_COMM_WORLD) == MPI_SUCCESS);
    for (unsigned i = 0; i < 7; ++i) CHECK(global_hits[i] == 1);
  }
  CHECK(mock_cuda.context_creations == 1 && mock_cuda.context_destructions == 1);
}

int main(int argc, char **argv) {
  fixture_t fixture;
  int rank = -1;
  int world_size = 0;
  fixture_init(&fixture);
  CHECK(setenv("CUDA_VISIBLE_DEVICES", "comb-mock-device", 1) == 0);
  /* An omitted optional result must not dereference NULL in error paths. */
  CHECK(comb_cuda_mpi_verify_comm(MPI_COMM_NULL, NULL, NULL, NULL, NULL) ==
        COMB_CUDA_MPI_STATUS_MPI_NOT_INITIALIZED);
  CHECK(MPI_Init(&argc, &argv) == MPI_SUCCESS);
  CHECK(MPI_Comm_rank(MPI_COMM_WORLD, &rank) == MPI_SUCCESS);
  CHECK(MPI_Comm_size(MPI_COMM_WORLD, &world_size) == MPI_SUCCESS);
  CHECK(world_size >= 2);
  CHECK(comb_cuda_mpi_verify_comm(MPI_COMM_NULL, NULL, NULL, NULL, NULL) ==
        COMB_CUDA_MPI_STATUS_INVALID_ARGUMENT);
  verify_case(&fixture.problem, COMB_CUDA_WITNESS_DETERMINISTIC, 0,
              rank, world_size);
  verify_case(&fixture.problem, COMB_CUDA_WITNESS_DETERMINISTIC, 1,
              rank, world_size);
  verify_case(&fixture.problem, COMB_CUDA_WITNESS_FAST_ANY, 0,
              rank, world_size);
  verify_case(&fixture.problem, COMB_CUDA_WITNESS_FAST_ANY, 1,
              rank, world_size);
  {
    comb_cuda_mpi_options_t options = COMB_CUDA_MPI_OPTIONS_INITIALIZER;
    comb_cuda_mpi_result_t result = COMB_CUDA_MPI_RESULT_INITIALIZER;
    options.chunk_size = 1;
    options.epoch_size = 7;
    mock_cuda_reset(NULL, 0);
    fail_rank = 1;
    CHECK(comb_cuda_mpi_verify_comm(
        MPI_COMM_WORLD, rank == 0 ? &fixture.problem : NULL,
        rank == 0 ? &options : NULL, &result, NULL) ==
          COMB_CUDA_MPI_STATUS_CUDA_ERROR);
    CHECK(result.cuda_status == COMB_CUDA_STATUS_CUDA_ERROR);
    fail_rank = -1;
    duplicate_device = 1;
    CHECK(comb_cuda_mpi_verify_comm(
        MPI_COMM_WORLD, rank == 0 ? &fixture.problem : NULL,
        rank == 0 ? &options : NULL, &result, NULL) ==
          COMB_CUDA_MPI_STATUS_DEVICE_MAPPING_ERROR);
    duplicate_device = 0;
  }
  {
    char checkpoint_path[] = "/tmp/comb-cuda-mpi-multi-rank-XXXXXX";
    comb_cuda_mpi_options_t options = COMB_CUDA_MPI_OPTIONS_INITIALIZER;
    comb_cuda_mpi_result_t result = COMB_CUDA_MPI_RESULT_INITIALIZER;
    if (rank == 0) {
      const int descriptor = mkstemp(checkpoint_path);
      CHECK(descriptor >= 0);
      CHECK(close(descriptor) == 0);
      CHECK(unlink(checkpoint_path) == 0);
    }
    CHECK(MPI_Bcast(checkpoint_path, (int)sizeof(checkpoint_path), MPI_CHAR, 0,
                   MPI_COMM_WORLD) == MPI_SUCCESS);
    options.epoch_size = 3;
    options.chunk_size = 1;
    options.checkpoint_path = checkpoint_path;
    mock_cuda_reset(NULL, 0);
    /* Also exercise the full-node visible-device mapping. */
    CHECK(unsetenv("CUDA_VISIBLE_DEVICES") == 0);
    CHECK(comb_cuda_mpi_verify_comm(
        MPI_COMM_WORLD, rank == 0 ? &fixture.problem : NULL,
        rank == 0 ? &options : NULL, &result, NULL) ==
          COMB_CUDA_MPI_STATUS_SUCCESS);
    CHECK(result.cuda_device == result.local_rank);
    CHECK(result.verification_code == ALL_FAULTS_GOOD);
    CHECK(result.checked_combination_count == 7);
    options.resume_from_checkpoint = 1;
    mock_cuda_reset(NULL, 0);
    CHECK(comb_cuda_mpi_verify_comm(
        MPI_COMM_WORLD, rank == 0 ? &fixture.problem : NULL,
        rank == 0 ? &options : NULL, &result, NULL) ==
          COMB_CUDA_MPI_STATUS_SUCCESS);
    CHECK(result.start_rank == 7 && result.next_epoch == 7);
    CHECK(result.checked_combination_count == 0);
    CHECK(mock_cuda.verify_calls == 0);
    if (rank == 0) CHECK(unlink(checkpoint_path) == 0);
  }
  CHECK(MPI_Finalize() == MPI_SUCCESS);
  if (rank == 0) puts("Multi-rank MPI scheduling, witnesses, and failures passed");
  return EXIT_SUCCESS;
}
