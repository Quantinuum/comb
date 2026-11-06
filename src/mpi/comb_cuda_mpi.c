#define _POSIX_C_SOURCE 200809L

#include "comb_cuda_mpi.h"

#include <fcntl.h>
#include <inttypes.h>
#include <limits.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "comb_cuda_internal.h"

#define COMB_MPI_DEFAULT_EPOCH_SIZE UINT64_C(1048576)
#define COMB_MPI_DEFAULT_CHUNK_SIZE UINT64_C(16384)
#define COMB_MPI_BCAST_CHUNK INT_MAX
#define COMB_MPI_CHECKPOINT_TAG "COMB_CUDA_MPI_CHECKPOINT"

typedef struct {
  uint64_t epoch_size;
  uint64_t chunk_size;
  uint64_t scratch_limit_bytes;
  uint32_t witness_mode;
  uint32_t fast_poll_chunks;
  uint32_t checkpoint_enabled;
  uint32_t resume_from_checkpoint;
} comb_mpi_wire_options_t;

typedef struct {
  uint64_t checked;
  uint64_t minimum_bad_rank;
  comb_cuda_status_t cuda_status;
} comb_mpi_epoch_local_t;

static uint64_t comb_mpi_min_u64(uint64_t a, uint64_t b) {
  return a < b ? a : b;
}

static int comb_mpi_u64_add(uint64_t a, uint64_t b, uint64_t *sum) {
  if (sum == NULL || UINT64_MAX - a < b) return 0;
  *sum = a + b;
  return 1;
}

static int comb_mpi_u64_multiply(uint64_t a, uint64_t b,
                                 uint64_t *product) {
  if (product == NULL || (a != 0 && b > UINT64_MAX / a)) return 0;
  *product = a * b;
  return 1;
}

/* Keep the count at its wire width until this check succeeds.  Comparing a
 * narrowed uint32_t count directly with SIZE_MAX/sizeof(T) is diagnosed as
 * always false on 64-bit GCC, and narrowing first would lose overflow on
 * hosts whose size_t is smaller than the wire count. */
static int comb_mpi_allocation_fits(uint64_t count, size_t element_size) {
  return element_size != 0 && count <= (uint64_t)(SIZE_MAX / element_size);
}

static uint64_t comb_mpi_hash_bytes(const void *data, uint64_t size) {
  const unsigned char *bytes = (const unsigned char *)data;
  uint64_t hash = UINT64_C(14695981039346656037);
  uint64_t i;

  for (i = 0; i < size; ++i) {
    hash ^= (uint64_t)bytes[i];
    hash *= UINT64_C(1099511628211);
  }
  hash ^= size;
  hash *= UINT64_C(1099511628211);
  return hash;
}

static void comb_mpi_result_reset(comb_cuda_mpi_result_t *result) {
  if (result == NULL) return;
  memset(result, 0, sizeof(*result));
  result->struct_size = sizeof(*result);
  result->abi_version = COMB_CUDA_MPI_ABI_VERSION;
  result->verification_code = OTHER_ERROR;
  result->cuda_status = COMB_CUDA_STATUS_SUCCESS;
  result->orchestration_status = COMB_CUDA_MPI_STATUS_SUCCESS;
  result->mpi_error_code = MPI_SUCCESS;
  result->minimum_bad_rank = COMB_CUDA_RANK_NONE;
  result->world_rank = -1;
  result->local_rank = -1;
  result->cuda_device = -1;
}

static int comb_mpi_result_is_valid(const comb_cuda_mpi_result_t *result) {
  return result == NULL ||
         (result->struct_size >= sizeof(*result) &&
          result->abi_version == COMB_CUDA_MPI_ABI_VERSION);
}

static int comb_mpi_witness_is_valid(const comb_cuda_witness_t *witness,
                                     uint32_t s) {
  if (witness == NULL) return 1;
  if (witness->struct_size < sizeof(*witness) ||
      witness->abi_version != COMB_CUDA_ABI_VERSION ||
      witness->fault_capacity < s) {
    return 0;
  }
  return s == 0 || witness->fault_indices != NULL;
}

static int comb_mpi_sync_status(MPI_Comm communicator,
                                comb_cuda_mpi_status_t local,
                                comb_cuda_mpi_status_t *global) {
  int local_value = (int)local;
  int global_value = 0;
  int mpi_error = MPI_Allreduce(&local_value, &global_value, 1, MPI_INT,
                                MPI_MAX, communicator);
  if (mpi_error == MPI_SUCCESS && global != NULL) {
    *global = (comb_cuda_mpi_status_t)global_value;
  }
  return mpi_error;
}

static int comb_mpi_sync_cuda_status(MPI_Comm communicator,
                                     comb_cuda_status_t local,
                                     comb_cuda_status_t *global) {
  int local_value = (int)local;
  int global_value = 0;
  int mpi_error = MPI_Allreduce(&local_value, &global_value, 1, MPI_INT,
                                MPI_MAX, communicator);
  if (mpi_error == MPI_SUCCESS && global != NULL) {
    *global = (comb_cuda_status_t)global_value;
  }
  return mpi_error;
}

static int comb_mpi_bcast_bytes(MPI_Comm communicator, void *buffer,
                                uint64_t size, int root) {
  uint64_t offset = 0;
  while (offset < size) {
    uint64_t remaining = size - offset;
    int count = remaining > (uint64_t)COMB_MPI_BCAST_CHUNK
                    ? COMB_MPI_BCAST_CHUNK
                    : (int)remaining;
    int mpi_error = MPI_Bcast((unsigned char *)buffer + (size_t)offset, count,
                              MPI_BYTE, root, communicator);
    if (mpi_error != MPI_SUCCESS) return mpi_error;
    offset += (uint64_t)count;
  }
  return MPI_SUCCESS;
}

static int comb_mpi_token_count(const char *visible) {
  int count = 0;
  const char *cursor = visible;

  if (visible == NULL) return -1;
  while (*cursor != '\0') {
    const char *begin = cursor;
    const char *end;
    while (*cursor != '\0' && *cursor != ',') ++cursor;
    end = cursor;
    while (begin < end && (*begin == ' ' || *begin == '\t')) ++begin;
    while (end > begin && (end[-1] == ' ' || end[-1] == '\t')) --end;
    if (begin == end || count == INT_MAX) return 0;
    ++count;
    if (*cursor == ',') ++cursor;
  }
  return count;
}

static comb_cuda_mpi_status_t comb_mpi_map_device(
    MPI_Comm local_communicator, int local_rank, int local_size,
    int *cuda_device, int *mpi_error_out) {
  enum { COMB_MPI_PCI_ID_BYTES = 32 };
  const char *visible = getenv("CUDA_VISIBLE_DEVICES");
  int visible_count = comb_mpi_token_count(visible);
  int selected = -1;
  char pci_bus_id[COMB_MPI_PCI_ID_BYTES] = {0};
  char *identities = NULL;
  int mpi_error;
  int i;
  int j;
  comb_cuda_mpi_status_t local_status = COMB_CUDA_MPI_STATUS_SUCCESS;
  comb_cuda_mpi_status_t global_status = COMB_CUDA_MPI_STATUS_SUCCESS;

  if (visible_count == 0) {
    local_status = COMB_CUDA_MPI_STATUS_DEVICE_MAPPING_ERROR;
  } else if (visible_count < 0) {
    selected = local_rank;
  } else if (visible_count == 1) {
    selected = 0;
  } else {
    if (visible_count < local_size || local_rank >= visible_count) {
      local_status = COMB_CUDA_MPI_STATUS_DEVICE_MAPPING_ERROR;
    } else {
      selected = local_rank;
    }
  }

  if (local_status == COMB_CUDA_MPI_STATUS_SUCCESS) {
    comb_cuda_status_t cuda_status = comb_cuda_device_pci_bus_id(
        selected, pci_bus_id, sizeof(pci_bus_id));
    if (cuda_status != COMB_CUDA_STATUS_SUCCESS || pci_bus_id[0] == '\0') {
      local_status = COMB_CUDA_MPI_STATUS_DEVICE_MAPPING_ERROR;
    }
  }

  if (local_size < 0 ||
      !comb_mpi_allocation_fits((uint64_t)local_size, COMB_MPI_PCI_ID_BYTES)) {
    local_status = COMB_CUDA_MPI_STATUS_OUT_OF_MEMORY;
  } else if (local_status == COMB_CUDA_MPI_STATUS_SUCCESS) {
    identities = (char *)calloc((size_t)local_size, COMB_MPI_PCI_ID_BYTES);
    if (identities == NULL)
      local_status = COMB_CUDA_MPI_STATUS_OUT_OF_MEMORY;
  }

  mpi_error = comb_mpi_sync_status(local_communicator, local_status,
                                   &global_status);
  if (mpi_error != MPI_SUCCESS) {
    free(identities);
    if (mpi_error_out != NULL) *mpi_error_out = mpi_error;
    return COMB_CUDA_MPI_STATUS_MPI_ERROR;
  }
  if (global_status != COMB_CUDA_MPI_STATUS_SUCCESS) {
    free(identities);
    return global_status;
  }
  if (identities == NULL) return COMB_CUDA_MPI_STATUS_OUT_OF_MEMORY;

  mpi_error = MPI_Allgather(pci_bus_id, COMB_MPI_PCI_ID_BYTES, MPI_CHAR,
                            identities, COMB_MPI_PCI_ID_BYTES, MPI_CHAR,
                            local_communicator);
  if (mpi_error != MPI_SUCCESS) {
    free(identities);
    if (mpi_error_out != NULL) *mpi_error_out = mpi_error;
    return COMB_CUDA_MPI_STATUS_MPI_ERROR;
  }

  for (i = 0; i < local_size; ++i) {
    for (j = i + 1; j < local_size; ++j) {
      if (memcmp(identities + (size_t)i * COMB_MPI_PCI_ID_BYTES,
                 identities + (size_t)j * COMB_MPI_PCI_ID_BYTES,
                 COMB_MPI_PCI_ID_BYTES) == 0) {
        free(identities);
        return COMB_CUDA_MPI_STATUS_DEVICE_MAPPING_ERROR;
      }
    }
  }

  free(identities);
  *cuda_device = selected;
  return COMB_CUDA_MPI_STATUS_SUCCESS;
}

static int comb_mpi_line_is(const char *line, const char *expected) {
  size_t line_length = strlen(line);
  size_t expected_length = strlen(expected);
  while (line_length > 0 &&
         (line[line_length - 1] == '\n' || line[line_length - 1] == '\r')) {
    --line_length;
  }
  return line_length == expected_length &&
         memcmp(line, expected, expected_length) == 0;
}

static int comb_mpi_parse_hash_line(const char *line, uint64_t *value) {
  char trailing;
  return sscanf(line, "problem_hash %" SCNx64 " %c", value, &trailing) == 1;
}

static int comb_mpi_parse_epoch_line(const char *line, uint64_t *value) {
  char trailing;
  return sscanf(line, "next_epoch %" SCNu64 " %c", value, &trailing) == 1;
}

static int comb_mpi_parse_epoch_size_line(const char *line, uint64_t *value) {
  char trailing;
  return sscanf(line, "epoch_size %" SCNu64 " %c", value, &trailing) == 1;
}

static comb_cuda_mpi_status_t comb_mpi_checkpoint_read(
    const char *path, uint64_t expected_hash, uint64_t epoch_size,
    uint64_t total_count, uint64_t *next_epoch) {
  FILE *stream;
  char line[256];
  char expected_header[64];
  uint64_t stored_hash;
  uint64_t stored_epoch_size;
  uint64_t stored_epoch;
  int header_length;

  stream = fopen(path, "r");
  if (stream == NULL) return COMB_CUDA_MPI_STATUS_CHECKPOINT_IO_ERROR;

  header_length = snprintf(expected_header, sizeof(expected_header), "%s %u",
                           COMB_MPI_CHECKPOINT_TAG,
                           COMB_CUDA_MPI_CHECKPOINT_VERSION);
  if (header_length < 0 || (size_t)header_length >= sizeof(expected_header) ||
      fgets(line, sizeof(line), stream) == NULL ||
      !comb_mpi_line_is(line, expected_header) ||
      fgets(line, sizeof(line), stream) == NULL ||
      !comb_mpi_parse_hash_line(line, &stored_hash) ||
      fgets(line, sizeof(line), stream) == NULL ||
      !comb_mpi_parse_epoch_size_line(line, &stored_epoch_size) ||
      fgets(line, sizeof(line), stream) == NULL ||
      !comb_mpi_parse_epoch_line(line, &stored_epoch)) {
    fclose(stream);
    return COMB_CUDA_MPI_STATUS_CHECKPOINT_MISMATCH;
  }

  while (fgets(line, sizeof(line), stream) != NULL) {
    const char *cursor = line;
    while (*cursor == ' ' || *cursor == '\t' || *cursor == '\n' ||
           *cursor == '\r') {
      ++cursor;
    }
    if (*cursor != '\0') {
      fclose(stream);
      return COMB_CUDA_MPI_STATUS_CHECKPOINT_MISMATCH;
    }
  }
  {
    int read_failed = ferror(stream);
    int close_failed = fclose(stream) != 0;
    if (read_failed || close_failed) {
      return COMB_CUDA_MPI_STATUS_CHECKPOINT_IO_ERROR;
    }
  }

  if (stored_hash != expected_hash || stored_epoch_size != epoch_size ||
      stored_epoch > total_count ||
      (stored_epoch != total_count && stored_epoch % epoch_size != 0)) {
    return COMB_CUDA_MPI_STATUS_CHECKPOINT_MISMATCH;
  }

  *next_epoch = stored_epoch;
  return COMB_CUDA_MPI_STATUS_SUCCESS;
}

static comb_cuda_mpi_status_t comb_mpi_checkpoint_write(
    const char *path, uint64_t problem_hash, uint64_t epoch_size,
    uint64_t next_epoch) {
  static const char suffix[] = ".tmp.XXXXXX";
  size_t path_length = strlen(path);
  size_t temporary_length;
  char *temporary_path;
  int descriptor;
  FILE *stream;
  int failed = 0;

  if (path_length > SIZE_MAX - sizeof(suffix)) {
    return COMB_CUDA_MPI_STATUS_CHECKPOINT_IO_ERROR;
  }
  temporary_length = path_length + sizeof(suffix);
  temporary_path = (char *)malloc(temporary_length);
  if (temporary_path == NULL) return COMB_CUDA_MPI_STATUS_OUT_OF_MEMORY;
  memcpy(temporary_path, path, path_length);
  memcpy(temporary_path + path_length, suffix, sizeof(suffix));

  descriptor = mkstemp(temporary_path);
  if (descriptor < 0) {
    free(temporary_path);
    return COMB_CUDA_MPI_STATUS_CHECKPOINT_IO_ERROR;
  }

  stream = fdopen(descriptor, "w");
  if (stream == NULL) {
    close(descriptor);
    unlink(temporary_path);
    free(temporary_path);
    return COMB_CUDA_MPI_STATUS_CHECKPOINT_IO_ERROR;
  }

  if (fprintf(stream, "%s %u\nproblem_hash %016" PRIx64
                      "\nepoch_size %" PRIu64 "\nnext_epoch %" PRIu64
                      "\n",
              COMB_MPI_CHECKPOINT_TAG, COMB_CUDA_MPI_CHECKPOINT_VERSION,
              problem_hash, epoch_size, next_epoch) < 0 ||
      fflush(stream) != 0 || fsync(fileno(stream)) != 0) {
    failed = 1;
  }
  if (fclose(stream) != 0) failed = 1;

  if (!failed && rename(temporary_path, path) != 0) failed = 1;
  if (!failed) {
    const char *slash = strrchr(path, '/');
    const char *directory = ".";
    char *allocated_directory = NULL;
    int directory_descriptor = -1;
    int directory_flags = O_RDONLY;
#ifdef O_DIRECTORY
    directory_flags |= O_DIRECTORY;
#endif
    if (slash != NULL) {
      size_t directory_length = slash == path ? 1U : (size_t)(slash - path);
      allocated_directory = (char *)malloc(directory_length + 1U);
      if (allocated_directory == NULL) {
        failed = 1;
      } else {
        memcpy(allocated_directory, path, directory_length);
        allocated_directory[directory_length] = '\0';
        directory = allocated_directory;
      }
    }
    if (!failed) {
      directory_descriptor = open(directory, directory_flags);
      if (directory_descriptor < 0 || fsync(directory_descriptor) != 0) {
        failed = 1;
      }
      if (directory_descriptor >= 0 && close(directory_descriptor) != 0) {
        failed = 1;
      }
    }
    free(allocated_directory);
  }
  if (failed) unlink(temporary_path);
  free(temporary_path);
  return failed ? COMB_CUDA_MPI_STATUS_CHECKPOINT_IO_ERROR
                : COMB_CUDA_MPI_STATUS_SUCCESS;
}

static comb_cuda_mpi_status_t comb_mpi_normalize_root_options(
    const comb_cuda_mpi_options_t *options, comb_mpi_wire_options_t *wire,
    const char **checkpoint_path) {
  comb_cuda_mpi_options_t defaults = COMB_CUDA_MPI_OPTIONS_INITIALIZER;
  const comb_cuda_mpi_options_t *selected = options != NULL ? options : &defaults;
  uint32_t mode;

  if (selected->struct_size < sizeof(*selected) ||
      selected->abi_version != COMB_CUDA_MPI_ABI_VERSION ||
      selected->cuda_options.struct_size < sizeof(selected->cuda_options) ||
      selected->cuda_options.abi_version != COMB_CUDA_ABI_VERSION ||
      selected->cuda_options.cuda_device != -1) {
    return COMB_CUDA_MPI_STATUS_INVALID_ARGUMENT;
  }

  mode = selected->witness_mode;
  if (mode == COMB_CUDA_WITNESS_DEFAULT) {
    mode = selected->cuda_options.witness_mode;
  }
  if (mode == COMB_CUDA_WITNESS_DEFAULT) {
    mode = COMB_CUDA_WITNESS_DETERMINISTIC;
  }
  if (mode != COMB_CUDA_WITNESS_DETERMINISTIC &&
      mode != COMB_CUDA_WITNESS_FAST_ANY) {
    return COMB_CUDA_MPI_STATUS_INVALID_ARGUMENT;
  }

  memset(wire, 0, sizeof(*wire));
  wire->epoch_size = selected->epoch_size != 0
                         ? selected->epoch_size
                         : selected->cuda_options.epoch_size;
  if (wire->epoch_size == 0) wire->epoch_size = COMB_MPI_DEFAULT_EPOCH_SIZE;
  wire->chunk_size = selected->chunk_size != 0
                         ? selected->chunk_size
                         : COMB_MPI_DEFAULT_CHUNK_SIZE;
  if (wire->chunk_size > wire->epoch_size) {
    wire->chunk_size = wire->epoch_size;
  }
  wire->scratch_limit_bytes = selected->cuda_options.scratch_limit_bytes;
  wire->witness_mode = mode;
  wire->fast_poll_chunks = selected->fast_poll_chunks != 0
                               ? selected->fast_poll_chunks
                               : 1;
  wire->checkpoint_enabled = selected->checkpoint_path != NULL;
  wire->resume_from_checkpoint = selected->resume_from_checkpoint != 0;

  if (wire->resume_from_checkpoint && !wire->checkpoint_enabled) {
    return COMB_CUDA_MPI_STATUS_INVALID_ARGUMENT;
  }
  if (wire->checkpoint_enabled && selected->checkpoint_path[0] == '\0') {
    return COMB_CUDA_MPI_STATUS_INVALID_ARGUMENT;
  }
  *checkpoint_path = selected->checkpoint_path;
  return COMB_CUDA_MPI_STATUS_SUCCESS;
}

static void comb_mpi_verify_one_range(
    comb_cuda_context_t *context, uint64_t rank_begin, uint64_t rank_end,
    comb_cuda_witness_mode_t mode, uint32_t s, uint32_t *chunk_faults,
    uint32_t *best_faults, comb_mpi_epoch_local_t *local) {
  comb_cuda_result_t range_result = COMB_CUDA_RESULT_INITIALIZER;
  comb_cuda_witness_t range_witness = COMB_CUDA_WITNESS_INITIALIZER;
  comb_cuda_status_t status;
  uint64_t updated_checked;

  if (local->cuda_status != COMB_CUDA_STATUS_SUCCESS) return;

  range_witness.fault_capacity = s;
  range_witness.fault_indices = chunk_faults;
  status = comb_cuda_verify_range(context, rank_begin, rank_end, mode,
                                  &range_result, &range_witness);
  if (status != COMB_CUDA_STATUS_SUCCESS) {
    local->cuda_status = status;
    return;
  }
  if (range_result.checked_combination_count > rank_end - rank_begin ||
      !comb_mpi_u64_add(local->checked,
                        range_result.checked_combination_count,
                        &updated_checked)) {
    local->cuda_status = COMB_CUDA_STATUS_INTERNAL_ERROR;
    return;
  }
  local->checked = updated_checked;

  /*
   * Deterministic distributed ordering relies on every assigned range being
   * complete.  The CUDA backend may still evaluate candidates internally in
   * parallel, but it must report the full range count before MIN reduction.
   */
  if (mode == COMB_CUDA_WITNESS_DETERMINISTIC &&
      range_result.checked_combination_count != rank_end - rank_begin) {
    local->cuda_status = COMB_CUDA_STATUS_INTERNAL_ERROR;
    return;
  }

  if (range_result.verification_code == ALL_FAULTS_GOOD) {
    if (range_result.checked_combination_count != rank_end - rank_begin ||
        range_result.minimum_bad_rank != COMB_CUDA_RANK_NONE) {
      local->cuda_status = COMB_CUDA_STATUS_INTERNAL_ERROR;
    }
    return;
  }
  if (range_result.verification_code != BAD_FAULT_FOUND ||
      range_result.minimum_bad_rank < rank_begin ||
      range_result.minimum_bad_rank >= rank_end ||
      range_witness.rank != range_result.minimum_bad_rank ||
      range_witness.fault_count != s) {
    local->cuda_status = COMB_CUDA_STATUS_INTERNAL_ERROR;
    return;
  }

  if (range_result.minimum_bad_rank < local->minimum_bad_rank) {
    local->minimum_bad_rank = range_result.minimum_bad_rank;
    if (s != 0) memcpy(best_faults, chunk_faults, (size_t)s * sizeof(*best_faults));
  }
}

static int comb_mpi_deterministic_epoch(
    MPI_Comm communicator, int world_rank, int world_size,
    comb_cuda_context_t *context, uint64_t epoch_begin, uint64_t epoch_end,
    uint64_t chunk_size, uint32_t s, uint32_t *chunk_faults,
    uint32_t *best_faults, comb_mpi_epoch_local_t *local,
    uint64_t *global_checked, uint64_t *global_bad,
    comb_cuda_status_t *global_cuda_status) {
  uint64_t chunk_index = 0;
  uint64_t chunk_begin;
  int mpi_error;

  for (chunk_begin = epoch_begin; chunk_begin < epoch_end; ++chunk_index) {
    uint64_t chunk_end = comb_mpi_min_u64(epoch_end, chunk_begin +
                                                        comb_mpi_min_u64(
                                                            chunk_size,
                                                            UINT64_MAX -
                                                                chunk_begin));
    if ((int)(chunk_index % (uint64_t)world_size) == world_rank) {
      comb_mpi_verify_one_range(context, chunk_begin, chunk_end,
                                COMB_CUDA_WITNESS_DETERMINISTIC, s,
                                chunk_faults, best_faults, local);
    }
    chunk_begin = chunk_end;
  }

  mpi_error = comb_mpi_sync_cuda_status(communicator, local->cuda_status,
                                        global_cuda_status);
  if (mpi_error != MPI_SUCCESS) return mpi_error;
  mpi_error = MPI_Allreduce(&local->checked, global_checked, 1, MPI_UINT64_T,
                            MPI_SUM, communicator);
  if (mpi_error != MPI_SUCCESS) return mpi_error;
  return MPI_Allreduce(&local->minimum_bad_rank, global_bad, 1, MPI_UINT64_T,
                       MPI_MIN, communicator);
}

static int comb_mpi_wait_fast_reductions(MPI_Request requests[2],
                                         int *pending) {
  int mpi_error;
  if (!*pending) return MPI_SUCCESS;
  mpi_error = MPI_Waitall(2, requests, MPI_STATUSES_IGNORE);
  *pending = 0;
  return mpi_error;
}

static int comb_mpi_fast_epoch(
    MPI_Comm communicator, int world_rank, int world_size,
    comb_cuda_context_t *context, uint64_t epoch_begin, uint64_t epoch_end,
    uint64_t chunk_size, uint32_t poll_chunks, uint32_t s,
    uint32_t *chunk_faults, uint32_t *best_faults,
    comb_mpi_epoch_local_t *local, uint64_t *global_checked,
    uint64_t *global_bad, comb_cuda_status_t *global_cuda_status) {
  uint64_t epoch_length = epoch_end - epoch_begin;
  uint64_t chunk_count = epoch_length / chunk_size +
                         (epoch_length % chunk_size != 0);
  uint64_t cyclic_rounds = chunk_count / (uint64_t)world_size +
                           (chunk_count % (uint64_t)world_size != 0);
  uint64_t group_count = cyclic_rounds / (uint64_t)poll_chunks +
                         (cyclic_rounds % (uint64_t)poll_chunks != 0);
  uint64_t group;
  uint64_t send_bad = COMB_CUDA_RANK_NONE;
  uint64_t reduced_bad = COMB_CUDA_RANK_NONE;
  int send_cuda = COMB_CUDA_STATUS_SUCCESS;
  int reduced_cuda = COMB_CUDA_STATUS_SUCCESS;
  MPI_Request requests[2] = {MPI_REQUEST_NULL, MPI_REQUEST_NULL};
  int pending = 0;
  int mpi_error = MPI_SUCCESS;

  for (group = 0; group < group_count; ++group) {
    uint64_t group_base = 0;
    uint32_t item;
    if (!comb_mpi_u64_multiply(group, (uint64_t)poll_chunks,
                               &group_base)) {
      local->cuda_status = COMB_CUDA_STATUS_INTERNAL_ERROR;
    }
    for (item = 0; item < poll_chunks; ++item) {
      uint64_t cyclic_round;
      uint64_t chunk_base;
      uint64_t chunk_index;
      uint64_t chunk_begin;
      uint64_t chunk_end;
      if (local->cuda_status != COMB_CUDA_STATUS_SUCCESS) break;
      if (!comb_mpi_u64_add(group_base, (uint64_t)item, &cyclic_round)) {
        local->cuda_status = COMB_CUDA_STATUS_INTERNAL_ERROR;
        break;
      }
      if (cyclic_round >= cyclic_rounds) break;
      if (!comb_mpi_u64_multiply(cyclic_round, (uint64_t)world_size,
                                 &chunk_base) ||
          !comb_mpi_u64_add(chunk_base, (uint64_t)world_rank,
                            &chunk_index)) {
        local->cuda_status = COMB_CUDA_STATUS_INTERNAL_ERROR;
        break;
      }
      if (chunk_index >= chunk_count) continue;
      if (chunk_index > (UINT64_MAX - epoch_begin) / chunk_size) {
        local->cuda_status = COMB_CUDA_STATUS_INTERNAL_ERROR;
        continue;
      }
      chunk_begin = epoch_begin + chunk_index * chunk_size;
      chunk_end = comb_mpi_min_u64(epoch_end, chunk_begin +
                                                  comb_mpi_min_u64(
                                                      chunk_size,
                                                      UINT64_MAX -
                                                          chunk_begin));
      comb_mpi_verify_one_range(context, chunk_begin, chunk_end,
                                COMB_CUDA_WITNESS_FAST_ANY, s, chunk_faults,
                                best_faults, local);
    }

    mpi_error = comb_mpi_wait_fast_reductions(requests, &pending);
    if (mpi_error != MPI_SUCCESS) return mpi_error;
    if (reduced_bad != COMB_CUDA_RANK_NONE ||
        reduced_cuda != COMB_CUDA_STATUS_SUCCESS) {
      break;
    }

    send_bad = local->minimum_bad_rank;
    send_cuda = (int)local->cuda_status;
    mpi_error = MPI_Iallreduce(&send_bad, &reduced_bad, 1, MPI_UINT64_T,
                               MPI_MIN, communicator, &requests[0]);
    if (mpi_error != MPI_SUCCESS) return mpi_error;
    mpi_error = MPI_Iallreduce(&send_cuda, &reduced_cuda, 1, MPI_INT, MPI_MAX,
                               communicator, &requests[1]);
    if (mpi_error != MPI_SUCCESS) {
      MPI_Wait(&requests[0], MPI_STATUS_IGNORE);
      return mpi_error;
    }
    pending = 1;
  }

  mpi_error = comb_mpi_wait_fast_reductions(requests, &pending);
  if (mpi_error != MPI_SUCCESS) return mpi_error;

  mpi_error = comb_mpi_sync_cuda_status(communicator, local->cuda_status,
                                        global_cuda_status);
  if (mpi_error != MPI_SUCCESS) return mpi_error;
  mpi_error = MPI_Allreduce(&local->checked, global_checked, 1, MPI_UINT64_T,
                            MPI_SUM, communicator);
  if (mpi_error != MPI_SUCCESS) return mpi_error;
  return MPI_Allreduce(&local->minimum_bad_rank, global_bad, 1, MPI_UINT64_T,
                       MPI_MIN, communicator);
}

static int comb_mpi_broadcast_witness(MPI_Comm communicator, int world_rank,
                                      uint64_t global_bad, uint32_t s,
                                      uint64_t local_bad,
                                      const uint32_t *best_faults,
                                      uint32_t *winner_faults) {
  int local_owner = local_bad == global_bad ? world_rank : INT_MAX;
  int owner = INT_MAX;
  int mpi_error = MPI_Allreduce(&local_owner, &owner, 1, MPI_INT, MPI_MIN,
                                communicator);
  if (mpi_error != MPI_SUCCESS) return mpi_error;
  if (owner == INT_MAX) return MPI_ERR_OTHER;
  if (world_rank == owner && s != 0) {
    memcpy(winner_faults, best_faults, (size_t)s * sizeof(*winner_faults));
  }
  return MPI_Bcast(winner_faults, (int)s, MPI_UINT32_T, owner, communicator);
}

comb_cuda_mpi_status_t comb_cuda_mpi_verify_comm(
    MPI_Comm communicator, const comb_problem_t *problem,
    const comb_cuda_mpi_options_t *options, comb_cuda_mpi_result_t *result,
    comb_cuda_witness_t *witness) {
  MPI_Comm work_communicator = MPI_COMM_NULL;
  MPI_Comm local_communicator = MPI_COMM_NULL;
  comb_cuda_mpi_status_t local_status = COMB_CUDA_MPI_STATUS_SUCCESS;
  comb_cuda_mpi_status_t global_status = COMB_CUDA_MPI_STATUS_SUCCESS;
  comb_cuda_status_t local_cuda_status = COMB_CUDA_STATUS_SUCCESS;
  comb_cuda_status_t global_cuda_status = COMB_CUDA_STATUS_SUCCESS;
  comb_mpi_wire_options_t wire;
  const char *checkpoint_path = NULL;
  comb_cuda_context_t *context = NULL;
  void *flat_blob = NULL;
  uint64_t flat_blob_size = 0;
  uint64_t problem_hash = 0;
  uint64_t minimum_hash = 0;
  uint64_t maximum_hash = 0;
  uint64_t total_count = 0;
  uint64_t start_rank = 0;
  uint64_t next_epoch = 0;
  uint64_t checked_count = 0;
  uint64_t global_bad = COMB_CUDA_RANK_NONE;
  uint32_t s = 0;
  uint32_t *chunk_faults = NULL;
  uint32_t *best_faults = NULL;
  uint32_t *winner_faults = NULL;
  int mpi_initialized = 0;
  int mpi_finalized = 0;
  int mpi_error = MPI_SUCCESS;
  int world_rank = -1;
  int world_size = 0;
  int local_rank = -1;
  int local_size = 0;
  int cuda_device = -1;
  int root_status_value = 0;
  int root_cuda_status_value = 0;
  double start_time = 0.0;

  memset(&wire, 0, sizeof(wire));

  if (MPI_Initialized(&mpi_initialized) != MPI_SUCCESS || !mpi_initialized ||
      MPI_Finalized(&mpi_finalized) != MPI_SUCCESS || mpi_finalized) {
    if (result != NULL && comb_mpi_result_is_valid(result)) {
      comb_mpi_result_reset(result);
      result->orchestration_status =
          COMB_CUDA_MPI_STATUS_MPI_NOT_INITIALIZED;
    }
    return COMB_CUDA_MPI_STATUS_MPI_NOT_INITIALIZED;
  }
  if (communicator == MPI_COMM_NULL) {
    if (result != NULL && comb_mpi_result_is_valid(result)) {
      comb_mpi_result_reset(result);
      result->orchestration_status = COMB_CUDA_MPI_STATUS_INVALID_ARGUMENT;
    }
    return COMB_CUDA_MPI_STATUS_INVALID_ARGUMENT;
  }

  mpi_error = MPI_Comm_dup(communicator, &work_communicator);
  if (mpi_error != MPI_SUCCESS) {
    if (result != NULL && comb_mpi_result_is_valid(result)) {
      comb_mpi_result_reset(result);
      result->orchestration_status = COMB_CUDA_MPI_STATUS_MPI_ERROR;
      result->mpi_error_code = mpi_error;
    }
    return COMB_CUDA_MPI_STATUS_MPI_ERROR;
  }
  mpi_error = MPI_Comm_set_errhandler(work_communicator, MPI_ERRORS_RETURN);
  if (mpi_error != MPI_SUCCESS) {
    local_status = COMB_CUDA_MPI_STATUS_MPI_ERROR;
    goto cleanup;
  }
  start_time = MPI_Wtime();

  if ((mpi_error = MPI_Comm_rank(work_communicator, &world_rank)) !=
          MPI_SUCCESS ||
      (mpi_error = MPI_Comm_size(work_communicator, &world_size)) !=
          MPI_SUCCESS) {
    local_status = COMB_CUDA_MPI_STATUS_MPI_ERROR;
    goto cleanup;
  }
  {
    int is_intercommunicator = 0;
    mpi_error =
        MPI_Comm_test_inter(work_communicator, &is_intercommunicator);
    if (mpi_error != MPI_SUCCESS) {
      local_status = COMB_CUDA_MPI_STATUS_MPI_ERROR;
      goto cleanup;
    }
    if (is_intercommunicator) {
      local_status = COMB_CUDA_MPI_STATUS_INVALID_ARGUMENT;
      goto cleanup;
    }
  }

  if (!comb_mpi_result_is_valid(result)) {
    local_status = COMB_CUDA_MPI_STATUS_INVALID_ARGUMENT;
  } else {
    comb_mpi_result_reset(result);
  }
  mpi_error = comb_mpi_sync_status(work_communicator, local_status,
                                   &global_status);
  if (mpi_error != MPI_SUCCESS) {
    local_status = COMB_CUDA_MPI_STATUS_MPI_ERROR;
    goto cleanup;
  }
  if (global_status != COMB_CUDA_MPI_STATUS_SUCCESS) {
    local_status = global_status;
    goto cleanup;
  }

  if (world_rank == 0) {
    if (problem == NULL) {
      global_status = COMB_CUDA_MPI_STATUS_INVALID_ARGUMENT;
    } else {
      global_status = comb_mpi_normalize_root_options(
          options, &wire, &checkpoint_path);
    }
    root_status_value = (int)global_status;
  }
  mpi_error = MPI_Bcast(&root_status_value, 1, MPI_INT, 0,
                        work_communicator);
  if (mpi_error != MPI_SUCCESS) {
    local_status = COMB_CUDA_MPI_STATUS_MPI_ERROR;
    goto cleanup;
  }
  global_status = (comb_cuda_mpi_status_t)root_status_value;
  if (global_status != COMB_CUDA_MPI_STATUS_SUCCESS) {
    local_status = global_status;
    goto cleanup;
  }
  mpi_error = MPI_Bcast(&wire, (int)sizeof(wire), MPI_BYTE, 0,
                        work_communicator);
  if (mpi_error != MPI_SUCCESS) {
    local_status = COMB_CUDA_MPI_STATUS_MPI_ERROR;
    goto cleanup;
  }

  mpi_error = MPI_Comm_split_type(work_communicator, MPI_COMM_TYPE_SHARED,
                                  world_rank, MPI_INFO_NULL,
                                  &local_communicator);
  if (mpi_error != MPI_SUCCESS) {
    local_status = COMB_CUDA_MPI_STATUS_MPI_ERROR;
    goto cleanup;
  }
  mpi_error = MPI_Comm_set_errhandler(local_communicator, MPI_ERRORS_RETURN);
  if (mpi_error != MPI_SUCCESS) {
    local_status = COMB_CUDA_MPI_STATUS_MPI_ERROR;
    goto cleanup;
  }
  if ((mpi_error = MPI_Comm_rank(local_communicator, &local_rank)) !=
          MPI_SUCCESS ||
      (mpi_error = MPI_Comm_size(local_communicator, &local_size)) !=
          MPI_SUCCESS) {
    local_status = COMB_CUDA_MPI_STATUS_MPI_ERROR;
    goto cleanup;
  }
  local_status = comb_mpi_map_device(local_communicator, local_rank,
                                     local_size, &cuda_device, &mpi_error);
  if (local_status == COMB_CUDA_MPI_STATUS_MPI_ERROR) goto cleanup;
  mpi_error = comb_mpi_sync_status(work_communicator, local_status,
                                   &global_status);
  if (mpi_error != MPI_SUCCESS) {
    local_status = COMB_CUDA_MPI_STATUS_MPI_ERROR;
    goto cleanup;
  }
  if (global_status != COMB_CUDA_MPI_STATUS_SUCCESS) {
    local_status = global_status;
    goto cleanup;
  }

  if (world_rank == 0) {
    local_cuda_status =
        comb_cuda_flat_problem_pack(problem, &flat_blob, &flat_blob_size);
    if (local_cuda_status == COMB_CUDA_STATUS_SUCCESS) {
      local_cuda_status =
          comb_cuda_flat_problem_validate_science(flat_blob, flat_blob_size);
    }
    root_cuda_status_value = (int)local_cuda_status;
  }
  mpi_error = MPI_Bcast(&root_cuda_status_value, 1, MPI_INT, 0,
                        work_communicator);
  if (mpi_error != MPI_SUCCESS) {
    local_status = COMB_CUDA_MPI_STATUS_MPI_ERROR;
    goto cleanup;
  }
  global_cuda_status = (comb_cuda_status_t)root_cuda_status_value;
  if (global_cuda_status != COMB_CUDA_STATUS_SUCCESS) {
    local_status = COMB_CUDA_MPI_STATUS_CUDA_ERROR;
    goto cleanup;
  }
  mpi_error = MPI_Bcast(&flat_blob_size, 1, MPI_UINT64_T, 0,
                        work_communicator);
  if (mpi_error != MPI_SUCCESS) {
    local_status = COMB_CUDA_MPI_STATUS_MPI_ERROR;
    goto cleanup;
  }
  if (flat_blob_size < sizeof(comb_cuda_flat_problem_t) ||
      !comb_mpi_allocation_fits(flat_blob_size, 1)) {
    local_status = COMB_CUDA_MPI_STATUS_INTERNAL_ERROR;
  } else if (world_rank != 0) {
    flat_blob = malloc((size_t)flat_blob_size);
    if (flat_blob == NULL) local_status = COMB_CUDA_MPI_STATUS_OUT_OF_MEMORY;
  }
  mpi_error = comb_mpi_sync_status(work_communicator, local_status,
                                   &global_status);
  if (mpi_error != MPI_SUCCESS) {
    local_status = COMB_CUDA_MPI_STATUS_MPI_ERROR;
    goto cleanup;
  }
  if (global_status != COMB_CUDA_MPI_STATUS_SUCCESS) {
    local_status = global_status;
    goto cleanup;
  }
  mpi_error = comb_mpi_bcast_bytes(work_communicator, flat_blob,
                                   flat_blob_size, 0);
  if (mpi_error != MPI_SUCCESS) {
    local_status = COMB_CUDA_MPI_STATUS_MPI_ERROR;
    goto cleanup;
  }

  local_cuda_status =
      comb_cuda_flat_problem_validate_science(flat_blob, flat_blob_size);
  mpi_error = comb_mpi_sync_cuda_status(work_communicator, local_cuda_status,
                                        &global_cuda_status);
  if (mpi_error != MPI_SUCCESS) {
    local_status = COMB_CUDA_MPI_STATUS_MPI_ERROR;
    goto cleanup;
  }
  if (global_cuda_status != COMB_CUDA_STATUS_SUCCESS) {
    local_status = COMB_CUDA_MPI_STATUS_CUDA_ERROR;
    goto cleanup;
  }

  problem_hash = comb_mpi_hash_bytes(flat_blob, flat_blob_size);
  mpi_error = MPI_Allreduce(&problem_hash, &minimum_hash, 1, MPI_UINT64_T,
                            MPI_MIN, work_communicator);
  if (mpi_error == MPI_SUCCESS) {
    mpi_error = MPI_Allreduce(&problem_hash, &maximum_hash, 1, MPI_UINT64_T,
                              MPI_MAX, work_communicator);
  }
  if (mpi_error != MPI_SUCCESS) {
    local_status = COMB_CUDA_MPI_STATUS_MPI_ERROR;
    goto cleanup;
  }
  if (minimum_hash != maximum_hash) {
    local_status = COMB_CUDA_MPI_STATUS_INTERNAL_ERROR;
    goto cleanup;
  }

  {
    const comb_cuda_flat_problem_t *flat =
        (const comb_cuda_flat_problem_t *)flat_blob;
    if (flat == NULL) {
      local_status = COMB_CUDA_MPI_STATUS_INTERNAL_ERROR;
      goto cleanup;
    }
    total_count = flat->total_combination_count;
    s = flat->s;
  }
  if (s > (uint32_t)INT_MAX ||
      !comb_mpi_allocation_fits((uint64_t)s, sizeof(*chunk_faults))) {
    global_cuda_status = COMB_CUDA_STATUS_PROBLEM_TOO_LARGE;
    local_status = COMB_CUDA_MPI_STATUS_CUDA_ERROR;
    goto cleanup;
  }

  if (world_rank == 0 && !comb_mpi_witness_is_valid(witness, s)) {
    local_status = COMB_CUDA_MPI_STATUS_INVALID_ARGUMENT;
  } else {
    local_status = COMB_CUDA_MPI_STATUS_SUCCESS;
  }
  mpi_error = comb_mpi_sync_status(work_communicator, local_status,
                                   &global_status);
  if (mpi_error != MPI_SUCCESS) {
    local_status = COMB_CUDA_MPI_STATUS_MPI_ERROR;
    goto cleanup;
  }
  if (global_status != COMB_CUDA_MPI_STATUS_SUCCESS) {
    local_status = global_status;
    goto cleanup;
  }
  if (world_rank == 0 && witness != NULL) {
    witness->rank = COMB_CUDA_RANK_NONE;
    witness->fault_count = 0;
  }
  chunk_faults = (uint32_t *)calloc(s != 0 ? (size_t)s : 1,
                                    sizeof(*chunk_faults));
  best_faults = (uint32_t *)calloc(s != 0 ? (size_t)s : 1,
                                   sizeof(*best_faults));
  winner_faults = (uint32_t *)calloc(s != 0 ? (size_t)s : 1,
                                     sizeof(*winner_faults));
  if (chunk_faults == NULL || best_faults == NULL || winner_faults == NULL) {
    local_status = COMB_CUDA_MPI_STATUS_OUT_OF_MEMORY;
  }
  mpi_error = comb_mpi_sync_status(work_communicator, local_status,
                                   &global_status);
  if (mpi_error != MPI_SUCCESS) {
    local_status = COMB_CUDA_MPI_STATUS_MPI_ERROR;
    goto cleanup;
  }
  if (global_status != COMB_CUDA_MPI_STATUS_SUCCESS) {
    local_status = global_status;
    goto cleanup;
  }
  if (chunk_faults == NULL || best_faults == NULL || winner_faults == NULL) {
    local_status = COMB_CUDA_MPI_STATUS_OUT_OF_MEMORY;
    goto cleanup;
  }

  if (world_rank == 0 && wire.resume_from_checkpoint) {
    global_status = comb_mpi_checkpoint_read(
        checkpoint_path, problem_hash, wire.epoch_size, total_count,
        &start_rank);
    root_status_value = (int)global_status;
  } else if (world_rank == 0) {
    root_status_value = COMB_CUDA_MPI_STATUS_SUCCESS;
  }
  mpi_error = MPI_Bcast(&root_status_value, 1, MPI_INT, 0,
                        work_communicator);
  if (mpi_error != MPI_SUCCESS) {
    local_status = COMB_CUDA_MPI_STATUS_MPI_ERROR;
    goto cleanup;
  }
  global_status = (comb_cuda_mpi_status_t)root_status_value;
  if (global_status != COMB_CUDA_MPI_STATUS_SUCCESS) {
    local_status = global_status;
    goto cleanup;
  }
  mpi_error = MPI_Bcast(&start_rank, 1, MPI_UINT64_T, 0, work_communicator);
  if (mpi_error != MPI_SUCCESS) {
    local_status = COMB_CUDA_MPI_STATUS_MPI_ERROR;
    goto cleanup;
  }
  next_epoch = start_rank;

  {
    comb_cuda_options_t cuda_options = COMB_CUDA_OPTIONS_INITIALIZER;
    cuda_options.cuda_device = cuda_device;
    cuda_options.witness_mode = wire.witness_mode;
    cuda_options.epoch_size = wire.epoch_size;
    cuda_options.scratch_limit_bytes = wire.scratch_limit_bytes;
    local_cuda_status = comb_cuda_context_create_from_flat_blob(
        flat_blob, flat_blob_size, &cuda_options, &context);
  }
  mpi_error = comb_mpi_sync_cuda_status(work_communicator, local_cuda_status,
                                        &global_cuda_status);
  if (mpi_error != MPI_SUCCESS) {
    local_status = COMB_CUDA_MPI_STATUS_MPI_ERROR;
    goto cleanup;
  }
  if (global_cuda_status != COMB_CUDA_STATUS_SUCCESS) {
    local_status = COMB_CUDA_MPI_STATUS_CUDA_ERROR;
    goto cleanup;
  }

  while (next_epoch < total_count) {
    uint64_t epoch_begin = next_epoch;
    uint64_t epoch_end =
        comb_mpi_min_u64(total_count,
                         epoch_begin + comb_mpi_min_u64(
                                           wire.epoch_size,
                                           UINT64_MAX - epoch_begin));
    uint64_t epoch_checked = 0;
    uint64_t epoch_bad = COMB_CUDA_RANK_NONE;
    uint64_t updated_checked;
    comb_mpi_epoch_local_t epoch_local;

    epoch_local.checked = 0;
    epoch_local.minimum_bad_rank = COMB_CUDA_RANK_NONE;
    epoch_local.cuda_status = COMB_CUDA_STATUS_SUCCESS;
    memset(best_faults, 0, (s != 0 ? (size_t)s : 1) * sizeof(*best_faults));

    if (wire.witness_mode == COMB_CUDA_WITNESS_DETERMINISTIC) {
      mpi_error = comb_mpi_deterministic_epoch(
          work_communicator, world_rank, world_size, context, epoch_begin,
          epoch_end, wire.chunk_size, s, chunk_faults, best_faults,
          &epoch_local, &epoch_checked, &epoch_bad, &global_cuda_status);
    } else {
      mpi_error = comb_mpi_fast_epoch(
          work_communicator, world_rank, world_size, context, epoch_begin,
          epoch_end, wire.chunk_size, wire.fast_poll_chunks, s, chunk_faults,
          best_faults, &epoch_local, &epoch_checked, &epoch_bad,
          &global_cuda_status);
    }
    if (mpi_error != MPI_SUCCESS) {
      local_status = COMB_CUDA_MPI_STATUS_MPI_ERROR;
      goto cleanup;
    }
    if (!comb_mpi_u64_add(checked_count, epoch_checked, &updated_checked)) {
      local_status = COMB_CUDA_MPI_STATUS_INTERNAL_ERROR;
      goto cleanup;
    }
    checked_count = updated_checked;
    if (global_cuda_status != COMB_CUDA_STATUS_SUCCESS) {
      local_status = COMB_CUDA_MPI_STATUS_CUDA_ERROR;
      goto cleanup;
    }
    if (epoch_bad != COMB_CUDA_RANK_NONE) {
      global_bad = epoch_bad;
      mpi_error = comb_mpi_broadcast_witness(
          work_communicator, world_rank, global_bad, s,
          epoch_local.minimum_bad_rank, best_faults, winner_faults);
      if (mpi_error != MPI_SUCCESS) {
        local_status = COMB_CUDA_MPI_STATUS_MPI_ERROR;
        goto cleanup;
      }
      break;
    }
    if (epoch_checked != epoch_end - epoch_begin) {
      local_status = COMB_CUDA_MPI_STATUS_INTERNAL_ERROR;
      goto cleanup;
    }

    next_epoch = epoch_end;
    if (world_rank == 0 && wire.checkpoint_enabled) {
      global_status = comb_mpi_checkpoint_write(
          checkpoint_path, problem_hash, wire.epoch_size, next_epoch);
      root_status_value = (int)global_status;
    } else if (world_rank == 0) {
      root_status_value = COMB_CUDA_MPI_STATUS_SUCCESS;
    }
    mpi_error = MPI_Bcast(&root_status_value, 1, MPI_INT, 0,
                          work_communicator);
    if (mpi_error != MPI_SUCCESS) {
      local_status = COMB_CUDA_MPI_STATUS_MPI_ERROR;
      goto cleanup;
    }
    global_status = (comb_cuda_mpi_status_t)root_status_value;
    if (global_status != COMB_CUDA_MPI_STATUS_SUCCESS) {
      local_status = global_status;
      goto cleanup;
    }
  }

  if (world_rank == 0 && witness != NULL && global_bad != COMB_CUDA_RANK_NONE) {
    witness->rank = global_bad;
    witness->fault_count = s;
    if (s != 0) {
      memcpy(witness->fault_indices, winner_faults,
             (size_t)s * sizeof(*winner_faults));
    }
  }
  local_status = COMB_CUDA_MPI_STATUS_SUCCESS;

cleanup:
  if (context != NULL) comb_cuda_context_destroy(context);
  if (flat_blob != NULL) {
    if (world_rank == 0) {
      comb_cuda_flat_problem_free(flat_blob);
    } else {
      free(flat_blob);
    }
  }
  free(chunk_faults);
  free(best_faults);
  free(winner_faults);

  if (result != NULL && comb_mpi_result_is_valid(result)) {
    result->verification_code =
        local_status == COMB_CUDA_MPI_STATUS_SUCCESS
            ? (global_bad == COMB_CUDA_RANK_NONE ? ALL_FAULTS_GOOD
                                                 : BAD_FAULT_FOUND)
            : OTHER_ERROR;
    result->cuda_status = global_cuda_status;
    result->orchestration_status = local_status;
    result->mpi_error_code = mpi_error;
    result->total_combination_count = total_count;
    result->checked_combination_count = checked_count;
    result->minimum_bad_rank = global_bad;
    result->start_rank = start_rank;
    result->next_epoch = next_epoch;
    result->problem_hash = problem_hash;
    result->world_rank = world_rank;
    result->world_size = world_size;
    result->local_rank = local_rank;
    result->local_size = local_size;
    result->cuda_device = cuda_device;
    result->elapsed_milliseconds =
        world_rank >= 0 ? (MPI_Wtime() - start_time) * 1000.0 : 0.0;
  }

  if (local_communicator != MPI_COMM_NULL) MPI_Comm_free(&local_communicator);
  if (work_communicator != MPI_COMM_NULL) MPI_Comm_free(&work_communicator);
  return local_status;
}

comb_cuda_mpi_status_t comb_cuda_mpi_run(
    int *argc, char ***argv, const comb_problem_t *problem,
    const comb_cuda_mpi_options_t *options, comb_cuda_mpi_result_t *result,
    comb_cuda_witness_t *witness) {
  int initialized = 0;
  int finalized = 0;
  int provided = MPI_THREAD_SINGLE;
  int owns_mpi = 0;
  int mpi_error;
  comb_cuda_mpi_status_t status;

  mpi_error = MPI_Initialized(&initialized);
  if (mpi_error != MPI_SUCCESS) return COMB_CUDA_MPI_STATUS_MPI_ERROR;
  mpi_error = MPI_Finalized(&finalized);
  if (mpi_error != MPI_SUCCESS) return COMB_CUDA_MPI_STATUS_MPI_ERROR;
  if (finalized) return COMB_CUDA_MPI_STATUS_MPI_NOT_INITIALIZED;

  if (!initialized) {
    mpi_error = MPI_Init_thread(argc, argv, MPI_THREAD_FUNNELED, &provided);
    if (mpi_error != MPI_SUCCESS) return COMB_CUDA_MPI_STATUS_MPI_ERROR;
    owns_mpi = 1;
    if (provided < MPI_THREAD_FUNNELED) {
      mpi_error = MPI_Finalize();
      if (result != NULL && comb_mpi_result_is_valid(result)) {
        comb_mpi_result_reset(result);
        result->orchestration_status = COMB_CUDA_MPI_STATUS_MPI_ERROR;
        result->mpi_error_code =
            mpi_error == MPI_SUCCESS ? MPI_ERR_OTHER : mpi_error;
      }
      return COMB_CUDA_MPI_STATUS_MPI_ERROR;
    }
  }

  status = comb_cuda_mpi_verify_comm(MPI_COMM_WORLD, problem, options, result,
                                     witness);
  if (owns_mpi) {
    mpi_error = MPI_Finalize();
    if (mpi_error != MPI_SUCCESS && status == COMB_CUDA_MPI_STATUS_SUCCESS) {
      status = COMB_CUDA_MPI_STATUS_MPI_ERROR;
      if (result != NULL && comb_mpi_result_is_valid(result)) {
        result->orchestration_status = status;
        result->mpi_error_code = mpi_error;
        result->verification_code = OTHER_ERROR;
      }
    }
  }
  return status;
}

const char *comb_cuda_mpi_status_string(comb_cuda_mpi_status_t status) {
  switch (status) {
    case COMB_CUDA_MPI_STATUS_SUCCESS:
      return "success";
    case COMB_CUDA_MPI_STATUS_INVALID_ARGUMENT:
      return "invalid argument";
    case COMB_CUDA_MPI_STATUS_MPI_NOT_INITIALIZED:
      return "MPI is not initialized";
    case COMB_CUDA_MPI_STATUS_MPI_ERROR:
      return "MPI error";
    case COMB_CUDA_MPI_STATUS_DEVICE_MAPPING_ERROR:
      return "non-unique or invalid rank-to-GPU mapping";
    case COMB_CUDA_MPI_STATUS_CUDA_ERROR:
      return "CUDA verifier error";
    case COMB_CUDA_MPI_STATUS_CHECKPOINT_IO_ERROR:
      return "checkpoint I/O error";
    case COMB_CUDA_MPI_STATUS_CHECKPOINT_MISMATCH:
      return "checkpoint does not match the problem or epoch schedule";
    case COMB_CUDA_MPI_STATUS_OUT_OF_MEMORY:
      return "out of memory";
    case COMB_CUDA_MPI_STATUS_INTERNAL_ERROR:
      return "internal orchestration error";
    default:
      return "unknown MPI orchestration status";
  }
}
