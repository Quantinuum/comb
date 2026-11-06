#include <cuda_runtime.h>

#include <chrono>
#include <climits>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <new>

#include "comb_cuda.h"
#include "comb_cuda_internal.h"
#include "comb_cuda_pipeline.h"

namespace {

constexpr uint64_t kContextMagic = UINT64_C(0x434f4d4243554441); /* COMBCUDA */
constexpr uint64_t kDefaultEpochSize = UINT64_C(1) << 20;

comb_cuda_status_t normalize_options(const comb_cuda_options_t *options,
                                     comb_cuda_options_t *normalized) {
  *normalized = COMB_CUDA_OPTIONS_INITIALIZER;
  if (options != nullptr) {
    if (options->struct_size < sizeof(comb_cuda_options_t) ||
        options->abi_version != COMB_CUDA_ABI_VERSION) {
      return COMB_CUDA_STATUS_INVALID_ARGUMENT;
    }
    *normalized = *options;
  }
  if (normalized->cuda_device < -1 ||
      (normalized->witness_mode != COMB_CUDA_WITNESS_DEFAULT &&
       normalized->witness_mode != COMB_CUDA_WITNESS_DETERMINISTIC &&
       normalized->witness_mode != COMB_CUDA_WITNESS_FAST_ANY)) {
    return COMB_CUDA_STATUS_INVALID_ARGUMENT;
  }
  if (normalized->witness_mode == COMB_CUDA_WITNESS_DEFAULT) {
    normalized->witness_mode = COMB_CUDA_WITNESS_DETERMINISTIC;
  }
  if (normalized->epoch_size == 0) normalized->epoch_size = kDefaultEpochSize;
  return COMB_CUDA_STATUS_SUCCESS;
}

comb_cuda_status_t status_from_cuda(cudaError_t error) {
  if (error == cudaSuccess) return COMB_CUDA_STATUS_SUCCESS;
  if (error == cudaErrorMemoryAllocation) return COMB_CUDA_STATUS_OUT_OF_MEMORY;
  if (error == cudaErrorNoDevice || error == cudaErrorInsufficientDriver ||
      error == cudaErrorInitializationError) {
    return COMB_CUDA_STATUS_CUDA_UNAVAILABLE;
  }
  return COMB_CUDA_STATUS_CUDA_ERROR;
}

}  // namespace

struct comb_cuda_context {
  uint64_t magic;
  int device;
  uint32_t default_witness_mode;
  uint64_t epoch_size;
  uint64_t scratch_limit_bytes;
  uint64_t total_combination_count;
  uint32_t s;
  double setup_milliseconds;
  uint8_t *host_blob;
  uint64_t blob_size;
  void *device_blob;
  comb_cuda_scheduler_state_t *device_scheduler;
  void *device_scratch;
  uint64_t device_scratch_size;
  void *device_prepared;
  uint64_t scratch_per_worker;
  uint32_t worker_count;
};

/* Saturating rank claim: an atomicAdd could wrap at UINT64_MAX. */
__device__ bool claim_rank(uint64_t *counter_value, uint64_t end,
                           uint64_t *rank) {
  auto *counter = reinterpret_cast<unsigned long long *>(counter_value);
  unsigned long long current = atomicAdd(counter, 0ULL);
  while (current < end) {
    const unsigned long long observed =
        atomicCAS(counter, current, current + 1ULL);
    if (observed == current) {
      *rank = current;
      return true;
    }
    current = observed;
  }
  return false;
}

__global__ void verify_candidates(const uint8_t *blob, const uint8_t *prepared,
                                  uint8_t *scratch, uint64_t scratch_stride,
                                  uint32_t workers, bool fast,
                                  comb_cuda_scheduler_state_t *state) {
  const uint32_t worker = blockIdx.x * blockDim.x + threadIdx.x;
  if (worker >= workers) return;
  void *workspace = scratch + static_cast<uint64_t>(worker) * scratch_stride;
  for (;;) {
    if (atomicAdd(&state->runtime_error, 0U) != 0U) return;
    if (fast && atomicAdd(reinterpret_cast<unsigned long long *>(
                              &state->minimum_bad_rank),
                          0ULL) != COMB_CUDA_RANK_NONE)
      return;
    uint64_t rank;
    if (!claim_rank(&state->next_rank, state->rank_end, &rank)) return;
    int bad = 0;
    const int status =
        comb_gpu_candidate(blob, prepared, rank, workspace, &bad);
    if (status != COMB_CUDA_STATUS_SUCCESS) {
      atomicCAS(&state->runtime_error, 0U, static_cast<unsigned int>(status));
      return;
    }
    if (bad)
      atomicMin(
          reinterpret_cast<unsigned long long *>(&state->minimum_bad_rank),
          static_cast<unsigned long long>(rank));
    atomicAdd(reinterpret_cast<unsigned long long *>(&state->scheduled_count),
              1ULL);
  }
}

/* Allocate a bounded number of independent workspaces once per context. */
comb_cuda_status_t allocate_workspace(comb_cuda_context_t *context,
                                      const void *blob, uint64_t blob_size) {
  void *prepared = nullptr;
  uint64_t prepared_size = 0;
  comb_cuda_status_t status = comb_cuda_prepare(
      blob, blob_size, &prepared, &prepared_size, &context->scratch_per_worker);
  if (status != COMB_CUDA_STATUS_SUCCESS) return status;
  status = status_from_cuda(cudaMalloc(&context->device_prepared,
                                       static_cast<size_t>(prepared_size)));
  if (status == COMB_CUDA_STATUS_SUCCESS)
    status = status_from_cuda(cudaMemcpy(context->device_prepared, prepared,
                                         static_cast<size_t>(prepared_size),
                                         cudaMemcpyHostToDevice));
  std::free(prepared);
  if (status != COMB_CUDA_STATUS_SUCCESS) return status;
  size_t free_bytes = 0, total_bytes = 0;
  status = status_from_cuda(cudaMemGetInfo(&free_bytes, &total_bytes));
  if (status != COMB_CUDA_STATUS_SUCCESS) return status;
  uint64_t budget = context->scratch_limit_bytes;
  if (budget == 0) {
    budget = static_cast<uint64_t>(free_bytes) / 4;
    if (budget > (UINT64_C(256) << 20)) budget = UINT64_C(256) << 20;
  }
  if (budget > static_cast<uint64_t>(free_bytes)) budget = free_bytes;
  uint64_t workers = budget / context->scratch_per_worker;
  if (workers > 8192) workers = 8192;
  if (workers > context->total_combination_count)
    workers = context->total_combination_count;
  if (workers == 0) return COMB_CUDA_STATUS_OUT_OF_MEMORY;
  context->worker_count = static_cast<uint32_t>(workers);
  context->device_scratch_size = workers * context->scratch_per_worker;
  return status_from_cuda(
      cudaMalloc(&context->device_scratch,
                 static_cast<size_t>(context->device_scratch_size)));
}

extern "C" comb_cuda_status_t comb_cuda_context_create_from_flat_blob(
    const void *blob, uint64_t blob_size, const comb_cuda_options_t *options,
    comb_cuda_context_t **context) {
  if (context == nullptr) return COMB_CUDA_STATUS_INVALID_ARGUMENT;
  *context = nullptr;
  const auto started = std::chrono::steady_clock::now();
  comb_cuda_status_t status =
      comb_cuda_flat_problem_validate_science(blob, blob_size);
  if (status != COMB_CUDA_STATUS_SUCCESS) return status;

  comb_cuda_options_t normalized{};
  status = normalize_options(options, &normalized);
  if (status != COMB_CUDA_STATUS_SUCCESS) return status;

  int device_count = 0;
  status = status_from_cuda(cudaGetDeviceCount(&device_count));
  if (status != COMB_CUDA_STATUS_SUCCESS) return status;
  if (device_count <= 0) return COMB_CUDA_STATUS_CUDA_UNAVAILABLE;

  int device = normalized.cuda_device;
  if (device == -1) {
    status = status_from_cuda(cudaGetDevice(&device));
    if (status != COMB_CUDA_STATUS_SUCCESS) return status;
  }
  if (device < 0 || device >= device_count) {
    return COMB_CUDA_STATUS_INVALID_ARGUMENT;
  }
  status = status_from_cuda(cudaSetDevice(device));
  if (status != COMB_CUDA_STATUS_SUCCESS) return status;

  comb_cuda_context_t *created = new (std::nothrow) comb_cuda_context_t{};
  if (created == nullptr) return COMB_CUDA_STATUS_OUT_OF_MEMORY;
  created->magic = kContextMagic;
  created->device = device;
  created->default_witness_mode = normalized.witness_mode;
  created->epoch_size = normalized.epoch_size;
  created->scratch_limit_bytes = normalized.scratch_limit_bytes;
  created->blob_size = blob_size;
  const auto *header = reinterpret_cast<const comb_cuda_flat_problem_t *>(blob);
  created->total_combination_count = header->total_combination_count;
  created->s = header->s;
  created->host_blob =
      static_cast<uint8_t *>(std::malloc(static_cast<size_t>(blob_size)));
  if (created->host_blob == nullptr) {
    delete created;
    return COMB_CUDA_STATUS_OUT_OF_MEMORY;
  }
  std::memcpy(created->host_blob, blob, static_cast<size_t>(blob_size));

  status = status_from_cuda(
      cudaMalloc(&created->device_blob, static_cast<size_t>(blob_size)));
  if (status == COMB_CUDA_STATUS_SUCCESS) {
    status = status_from_cuda(cudaMemcpy(created->device_blob, blob,
                                         static_cast<size_t>(blob_size),
                                         cudaMemcpyHostToDevice));
  }
  if (status == COMB_CUDA_STATUS_SUCCESS) {
    status = status_from_cuda(
        cudaMalloc(reinterpret_cast<void **>(&created->device_scheduler),
                   sizeof(comb_cuda_scheduler_state_t)));
  }
  if (status == COMB_CUDA_STATUS_SUCCESS)
    status = allocate_workspace(created, blob, blob_size);
  if (status != COMB_CUDA_STATUS_SUCCESS) {
    if (created->device_scratch != nullptr) cudaFree(created->device_scratch);
    if (created->device_prepared != nullptr) cudaFree(created->device_prepared);
    if (created->device_scheduler != nullptr)
      cudaFree(created->device_scheduler);
    if (created->device_blob != nullptr) cudaFree(created->device_blob);
    std::free(created->host_blob);
    created->magic = 0;
    delete created;
    return status;
  }

  const auto finished = std::chrono::steady_clock::now();
  created->setup_milliseconds =
      std::chrono::duration<double, std::milli>(finished - started).count();
  *context = created;
  return COMB_CUDA_STATUS_SUCCESS;
}

extern "C" comb_cuda_status_t comb_cuda_device_pci_bus_id(int32_t device,
                                                          char *buffer,
                                                          size_t buffer_size) {
  if (device < 0 || buffer == nullptr || buffer_size < 16U ||
      buffer_size > static_cast<size_t>(INT_MAX)) {
    return COMB_CUDA_STATUS_INVALID_ARGUMENT;
  }
  buffer[0] = '\0';
  return status_from_cuda(cudaDeviceGetPCIBusId(
      buffer, static_cast<int>(buffer_size), static_cast<int>(device)));
}

extern "C" comb_cuda_status_t comb_cuda_context_create(
    const comb_problem_t *problem, const comb_cuda_options_t *options,
    comb_cuda_context_t **context) {
  if (context == nullptr) return COMB_CUDA_STATUS_INVALID_ARGUMENT;
  *context = nullptr;
  const auto started = std::chrono::steady_clock::now();
  void *blob = nullptr;
  uint64_t blob_size = 0;
  comb_cuda_status_t status =
      comb_cuda_flat_problem_pack(problem, &blob, &blob_size);
  if (status != COMB_CUDA_STATUS_SUCCESS) return status;
  status = comb_cuda_context_create_from_flat_blob(blob, blob_size, options,
                                                   context);
  comb_cuda_flat_problem_free(blob);
  if (status == COMB_CUDA_STATUS_SUCCESS) {
    const auto finished = std::chrono::steady_clock::now();
    (*context)->setup_milliseconds =
        std::chrono::duration<double, std::milli>(finished - started).count();
  }
  return status;
}

namespace {

comb_cuda_status_t validate_result(comb_cuda_result_t *result) {
  if (result == nullptr || result->struct_size < sizeof(comb_cuda_result_t) ||
      result->abi_version != COMB_CUDA_ABI_VERSION) {
    return COMB_CUDA_STATUS_INVALID_ARGUMENT;
  }
  return COMB_CUDA_STATUS_SUCCESS;
}

comb_cuda_status_t prepare_witness(comb_cuda_context_t *context,
                                   comb_cuda_witness_t *witness) {
  if (witness == nullptr) return COMB_CUDA_STATUS_SUCCESS;
  if (witness->struct_size < sizeof(comb_cuda_witness_t) ||
      witness->abi_version != COMB_CUDA_ABI_VERSION ||
      witness->fault_capacity < context->s ||
      witness->fault_indices == nullptr) {
    return COMB_CUDA_STATUS_INVALID_ARGUMENT;
  }
  witness->rank = COMB_CUDA_RANK_NONE;
  witness->fault_count = 0;
  return COMB_CUDA_STATUS_SUCCESS;
}

comb_cuda_status_t normalize_mode(const comb_cuda_context_t *context,
                                  comb_cuda_witness_mode_t requested,
                                  comb_cuda_witness_mode_t *normalized) {
  if (requested == COMB_CUDA_WITNESS_DEFAULT) {
    *normalized =
        static_cast<comb_cuda_witness_mode_t>(context->default_witness_mode);
    return COMB_CUDA_STATUS_SUCCESS;
  }
  if (requested != COMB_CUDA_WITNESS_DETERMINISTIC &&
      requested != COMB_CUDA_WITNESS_FAST_ANY) {
    return COMB_CUDA_STATUS_INVALID_ARGUMENT;
  }
  *normalized = requested;
  return COMB_CUDA_STATUS_SUCCESS;
}

void initialize_result(const comb_cuda_context_t *context,
                       comb_cuda_result_t *result) {
  result->verification_code = OTHER_ERROR;
  result->reserved0 = 0;
  result->total_combination_count = context->total_combination_count;
  result->checked_combination_count = 0;
  result->minimum_bad_rank = COMB_CUDA_RANK_NONE;
  result->setup_milliseconds = context->setup_milliseconds;
  result->verification_milliseconds = 0.0;
}

}  // namespace

extern "C" comb_cuda_status_t comb_cuda_verify_range(
    comb_cuda_context_t *context, uint64_t rank_begin, uint64_t rank_end,
    comb_cuda_witness_mode_t mode, comb_cuda_result_t *result,
    comb_cuda_witness_t *witness) {
  if (context == nullptr || context->magic != kContextMagic) {
    return COMB_CUDA_STATUS_INVALID_ARGUMENT;
  }
  comb_cuda_status_t status = validate_result(result);
  if (status != COMB_CUDA_STATUS_SUCCESS) return status;
  initialize_result(context, result);
  status = prepare_witness(context, witness);
  if (status != COMB_CUDA_STATUS_SUCCESS) return status;
  comb_cuda_witness_mode_t normalized_mode = COMB_CUDA_WITNESS_DEFAULT;
  status = normalize_mode(context, mode, &normalized_mode);
  if (status != COMB_CUDA_STATUS_SUCCESS) return status;
  if (rank_begin > rank_end || rank_end > context->total_combination_count) {
    return COMB_CUDA_STATUS_INVALID_ARGUMENT;
  }

  const auto started = std::chrono::steady_clock::now();
  if (rank_begin == rank_end) {
    result->verification_code = ALL_FAULTS_GOOD;
    const auto finished = std::chrono::steady_clock::now();
    result->verification_milliseconds =
        std::chrono::duration<double, std::milli>(finished - started).count();
    return COMB_CUDA_STATUS_SUCCESS;
  }

  status = status_from_cuda(cudaSetDevice(context->device));
  if (status != COMB_CUDA_STATUS_SUCCESS) return status;
  const bool fast = normalized_mode == COMB_CUDA_WITNESS_FAST_ANY;
  uint64_t begin = rank_begin;
  while (begin < rank_end) {
    const uint64_t remaining = rank_end - begin;
    const uint64_t count =
        remaining < context->epoch_size ? remaining : context->epoch_size;
    const uint64_t end = begin + count;
    comb_cuda_scheduler_state_t scheduler{};
    scheduler.next_rank = begin;
    scheduler.rank_end = end;
    scheduler.minimum_bad_rank = COMB_CUDA_RANK_NONE;
    status =
        status_from_cuda(cudaMemcpy(context->device_scheduler, &scheduler,
                                    sizeof(scheduler), cudaMemcpyHostToDevice));
    if (status != COMB_CUDA_STATUS_SUCCESS) return status;
    const uint32_t workers = count < context->worker_count
                                 ? static_cast<uint32_t>(count)
                                 : context->worker_count;
    verify_candidates<<<(workers + 63U) / 64U, 64>>>(
        static_cast<const uint8_t *>(context->device_blob),
        static_cast<const uint8_t *>(context->device_prepared),
        static_cast<uint8_t *>(context->device_scratch),
        context->scratch_per_worker, workers, fast, context->device_scheduler);
    status = status_from_cuda(cudaGetLastError());
    if (status == COMB_CUDA_STATUS_SUCCESS)
      status = status_from_cuda(cudaDeviceSynchronize());
    if (status == COMB_CUDA_STATUS_SUCCESS)
      status = status_from_cuda(
          cudaMemcpy(&scheduler, context->device_scheduler, sizeof(scheduler),
                     cudaMemcpyDeviceToHost));
    if (status != COMB_CUDA_STATUS_SUCCESS) return status;
    if (scheduler.runtime_error != 0)
      return static_cast<comb_cuda_status_t>(scheduler.runtime_error);
    result->checked_combination_count += scheduler.scheduled_count;
    if (scheduler.minimum_bad_rank < result->minimum_bad_rank)
      result->minimum_bad_rank = scheduler.minimum_bad_rank;
    if (!fast && scheduler.scheduled_count != count)
      return COMB_CUDA_STATUS_INTERNAL_ERROR;
    if (fast && result->minimum_bad_rank != COMB_CUDA_RANK_NONE) break;
    begin = end;
  }
  if (witness != nullptr && result->minimum_bad_rank != COMB_CUDA_RANK_NONE) {
    auto *locations = static_cast<uint32_t *>(
        std::malloc(static_cast<size_t>(context->s) * sizeof(uint32_t)));
    if (locations == nullptr) return COMB_CUDA_STATUS_OUT_OF_MEMORY;
    const bool decoded =
        comb_gpu_unrank(context->host_blob, result->minimum_bad_rank, locations,
                        witness->fault_indices) != 0;
    std::free(locations);
    if (!decoded) return COMB_CUDA_STATUS_INTERNAL_ERROR;
    witness->rank = result->minimum_bad_rank;
    witness->fault_count = context->s;
  }
  result->verification_code = result->minimum_bad_rank == COMB_CUDA_RANK_NONE
                                  ? ALL_FAULTS_GOOD
                                  : BAD_FAULT_FOUND;
  const auto finished = std::chrono::steady_clock::now();
  result->verification_milliseconds =
      std::chrono::duration<double, std::milli>(finished - started).count();
  return COMB_CUDA_STATUS_SUCCESS;
}

extern "C" comb_cuda_status_t comb_cuda_verify(comb_cuda_context_t *context,
                                               comb_cuda_witness_mode_t mode,
                                               comb_cuda_result_t *result,
                                               comb_cuda_witness_t *witness) {
  if (context == nullptr || context->magic != kContextMagic) {
    return COMB_CUDA_STATUS_INVALID_ARGUMENT;
  }
  return comb_cuda_verify_range(context, 0, context->total_combination_count,
                                mode, result, witness);
}

extern "C" void comb_cuda_context_destroy(comb_cuda_context_t *context) {
  if (context == nullptr || context->magic != kContextMagic) return;
  cudaSetDevice(context->device);
  if (context->device_scratch != nullptr) cudaFree(context->device_scratch);
  if (context->device_prepared != nullptr) cudaFree(context->device_prepared);
  if (context->device_scheduler != nullptr) cudaFree(context->device_scheduler);
  if (context->device_blob != nullptr) cudaFree(context->device_blob);
  std::free(context->host_blob);
  context->magic = 0;
  delete context;
}

extern "C" const char *comb_cuda_status_string(comb_cuda_status_t status) {
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
      return "operation not implemented";
    case COMB_CUDA_STATUS_INTERNAL_ERROR:
      return "internal error";
    default:
      return "unknown CUDA status";
  }
}
