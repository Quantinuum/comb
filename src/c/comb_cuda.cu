#include "comb_cuda.h"
#include "comb_cuda_internal.h"

#include <cuda_runtime.h>

#include <chrono>
#include <climits>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <new>

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
};

template <typename T>
__device__ __forceinline__ const T *flat_section(const uint8_t *blob,
                                                 uint64_t offset) {
  return reinterpret_cast<const T *>(blob + offset);
}

/* Canonical weighted-colex unranking.  W[m][r] is copied into FlatProblem by
 * the checked pure-C packer.  The first selected tuple element varies fastest
 * when local fault offsets are decoded. */
__device__ bool comb_cuda_unrank_weighted_colex(
    const uint8_t *blob, uint64_t rank, uint32_t *selected_locations,
    uint32_t *fault_indices) {
  const auto *problem =
      reinterpret_cast<const comb_cuda_flat_problem_t *>(blob);
  if (rank >= problem->total_combination_count) return false;
  const uint32_t *partitions =
      flat_section<uint32_t>(blob, problem->fault_partitions_offset);
  const uint64_t *prefixes =
      flat_section<uint64_t>(blob, problem->fault_prefixes_offset);
  const uint64_t *dp = flat_section<uint64_t>(blob, problem->rank_dp_offset);
  const uint32_t columns = problem->rank_dp_columns;

  uint32_t m = problem->fault_location_count;
  uint32_t r = problem->s;
  uint64_t probe = rank;
  while (r != 0) {
    if (m == 0) return false;
    const uint64_t without_last =
        dp[static_cast<uint64_t>(m - 1) * columns + r];
    if (probe < without_last) {
      --m;
      continue;
    }
    const uint32_t location = m - 1;
    const uint32_t partition = partitions[location];
    selected_locations[r - 1] = location;
    probe = (probe - without_last) / partition;
    --m;
    --r;
  }

  uint64_t block_begin = 0;
  for (uint32_t index = 0; index < problem->s; ++index) {
    const uint32_t location = selected_locations[index];
    block_begin = dp[static_cast<uint64_t>(location) * columns + index + 1] +
                  static_cast<uint64_t>(partitions[location]) * block_begin;
  }
  uint64_t local = rank - block_begin;
  for (uint32_t index = 0; index < problem->s; ++index) {
    const uint32_t location = selected_locations[index];
    const uint32_t partition = partitions[location];
    const uint32_t offset = static_cast<uint32_t>(local % partition);
    local /= partition;
    fault_indices[index] = static_cast<uint32_t>(prefixes[location] + offset);
  }
  return local == 0;
}

__device__ __forceinline__ uint32_t comb_cuda_warp_xor(uint32_t value) {
  for (int offset = 16; offset != 0; offset >>= 1) {
    value ^= __shfl_down_sync(0xffffffffU, value, offset);
  }
  return value;
}

__device__ uint32_t comb_cuda_block_and_parity(const uint64_t *left,
                                               const uint64_t *right,
                                               uint32_t words,
                                               uint32_t *warp_scratch) {
  uint32_t parity = 0;
  for (uint32_t word = threadIdx.x; word < words; word += blockDim.x) {
    parity ^= static_cast<uint32_t>(__popcll(left[word] & right[word]) & 1U);
  }
  parity = comb_cuda_warp_xor(parity);
  const uint32_t lane = threadIdx.x & 31U;
  const uint32_t warp = threadIdx.x >> 5U;
  if (lane == 0) warp_scratch[warp] = parity;
  __syncthreads();
  if (warp == 0) {
    const uint32_t warp_count = (blockDim.x + 31U) >> 5U;
    parity = lane < warp_count ? warp_scratch[lane] : 0;
    parity = comb_cuda_warp_xor(parity);
    if (lane == 0) warp_scratch[0] = parity;
  }
  __syncthreads();
  return warp_scratch[0];
}

__device__ void comb_cuda_row_xor(uint64_t *destination,
                                  const uint64_t *source, uint32_t words,
                                  uint64_t final_mask) {
  for (uint32_t word = threadIdx.x; word < words; word += blockDim.x) {
    destination[word] ^= source[word];
  }
  __syncthreads();
  if (threadIdx.x == 0 && words != 0) destination[words - 1] &= final_mask;
  __syncthreads();
}

/* Matches add_row_to: compute x_source . z_destination before either packed
 * destination row mutates, then apply the Dehaene-De Moor phase correction. */
__device__ void comb_cuda_phase_aware_row_add(
    uint64_t *destination_x, uint64_t *destination_z,
    uint32_t *destination_delta, uint32_t *destination_epsilon,
    const uint64_t *source_x, const uint64_t *source_z,
    uint32_t source_delta, uint32_t source_epsilon, uint32_t words,
    uint64_t final_mask, uint32_t *warp_scratch) {
  const uint32_t destination_delta_before = *destination_delta;
  const uint32_t a2_t_u_a1 = comb_cuda_block_and_parity(
      source_x, destination_z, words, warp_scratch);
  comb_cuda_row_xor(destination_x, source_x, words, final_mask);
  comb_cuda_row_xor(destination_z, source_z, words, final_mask);
  if (threadIdx.x == 0) {
    *destination_delta = destination_delta_before ^ source_delta;
    *destination_epsilon =
        *destination_epsilon ^ source_epsilon ^
        (source_delta & destination_delta_before) ^ a2_t_u_a1;
  }
  __syncthreads();
}

__device__ __forceinline__ void comb_cuda_atomic_min_rank(
    uint64_t *minimum_bad_rank, uint64_t rank) {
  atomicMin(reinterpret_cast<unsigned long long *>(minimum_bad_rank),
            static_cast<unsigned long long>(rank));
}

/* Saturating work claim.  atomicAdd would wrap when rank_end is UINT64_MAX,
 * allowing later teams to reschedule rank zero. */
__device__ bool comb_cuda_claim_rank(uint64_t *next_rank, uint64_t rank_end,
                                     uint64_t *rank) {
  auto *counter = reinterpret_cast<unsigned long long *>(next_rank);
  unsigned long long current = atomicAdd(counter, 0ULL);
  while (current < rank_end) {
    const unsigned long long desired = current + 1ULL;
    const unsigned long long observed = atomicCAS(counter, current, desired);
    if (observed == current) {
      *rank = static_cast<uint64_t>(current);
      return true;
    }
    current = observed;
  }
  return false;
}

/* Persistent scheduling skeleton.  It deliberately has no scientific badness
 * predicate yet and is not launched by the public verifier. */
template <bool OneWord>
__global__ void comb_cuda_candidate_scheduler_scaffold(
    const uint8_t *flat_blob, comb_cuda_scheduler_state_t *scheduler) {
  const auto *problem =
      reinterpret_cast<const comb_cuda_flat_problem_t *>(flat_blob);
  __shared__ uint64_t rank;
  __shared__ uint32_t have_rank;
  __shared__ uint32_t warp_scratch[32];
  for (;;) {
    if (threadIdx.x == 0) {
      have_rank = comb_cuda_claim_rank(&scheduler->next_rank,
                                       scheduler->rank_end, &rank)
                      ? 1U
                      : 0U;
    }
    __syncthreads();
    if (have_rank == 0U || scheduler->stop_requested != 0) return;

    uint32_t parity = 0;
    if (OneWord) {
      if (threadIdx.x == 0 && problem->words_all == 1) parity = 0;
    } else if (problem->words_all != 0) {
      const uint64_t *x =
          flat_section<uint64_t>(flat_blob, problem->output_dcp_hx_offset);
      const uint64_t *z =
          flat_section<uint64_t>(flat_blob, problem->output_dcp_hz_offset);
      parity = comb_cuda_block_and_parity(x, z, problem->words_all,
                                          warp_scratch);
    }
    (void)parity;
    if (threadIdx.x == 0) {
      atomicAdd(reinterpret_cast<unsigned long long *>(
                    &scheduler->scheduled_count),
                1ULL);
      atomicExch(reinterpret_cast<unsigned int *>(
                     &scheduler->pipeline_incomplete),
                 1U);
    }
    __syncthreads();
  }
}

/* Kept dormant until the scientific predicate is complete.  Its two branches
 * make the <=64-bit warp path and general multiword block path explicit and
 * compile both scheduler specializations in CUDA builds. */
[[maybe_unused]] comb_cuda_status_t comb_cuda_launch_scheduler_scaffold(
    comb_cuda_context_t *context, uint64_t rank_begin, uint64_t rank_end,
    uint64_t *scheduled_count) {
  comb_cuda_scheduler_state_t scheduler{};
  scheduler.next_rank = rank_begin;
  scheduler.rank_end = rank_end;
  scheduler.minimum_bad_rank = COMB_CUDA_RANK_NONE;
  comb_cuda_status_t status = status_from_cuda(cudaMemcpy(
      context->device_scheduler, &scheduler, sizeof(scheduler),
      cudaMemcpyHostToDevice));
  if (status != COMB_CUDA_STATUS_SUCCESS) return status;

  const auto *problem = reinterpret_cast<const comb_cuda_flat_problem_t *>(
      context->host_blob);
  if (problem->n_all <= 64) {
    comb_cuda_candidate_scheduler_scaffold<true><<<128, 32>>>(
        static_cast<const uint8_t *>(context->device_blob),
        context->device_scheduler);
  } else {
    comb_cuda_candidate_scheduler_scaffold<false><<<128, 128>>>(
        static_cast<const uint8_t *>(context->device_blob),
        context->device_scheduler);
  }
  status = status_from_cuda(cudaGetLastError());
  if (status == COMB_CUDA_STATUS_SUCCESS) {
    status = status_from_cuda(cudaDeviceSynchronize());
  }
  if (status == COMB_CUDA_STATUS_SUCCESS) {
    status = status_from_cuda(cudaMemcpy(&scheduler, context->device_scheduler,
                                         sizeof(scheduler),
                                         cudaMemcpyDeviceToHost));
  }
  if (scheduled_count != nullptr) {
    *scheduled_count = scheduler.scheduled_count;
  }
  return status;
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
  const auto *header =
      reinterpret_cast<const comb_cuda_flat_problem_t *>(blob);
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
    status = status_from_cuda(cudaMalloc(
        reinterpret_cast<void **>(&created->device_scheduler),
        sizeof(comb_cuda_scheduler_state_t)));
  }
  if (status != COMB_CUDA_STATUS_SUCCESS) {
    if (created->device_scheduler != nullptr) cudaFree(created->device_scheduler);
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

extern "C" comb_cuda_status_t comb_cuda_device_pci_bus_id(
    int32_t device, char *buffer, size_t buffer_size) {
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
      witness->fault_capacity < context->s || witness->fault_indices == nullptr) {
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
  (void)normalized_mode;
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

  /* Do not return a scientific decision until FlatProblem includes the CPU's
   * RREF pivots/permutations, C_z/C_x/P_x transforms, PLUQ permutation and row
   * schedule, followed by the complete DCP, postselection, phase-aware GE and
   * universal-coset device pipeline. */
  const auto finished = std::chrono::steady_clock::now();
  result->verification_milliseconds =
      std::chrono::duration<double, std::milli>(finished - started).count();
  return COMB_CUDA_STATUS_NOT_IMPLEMENTED;
}

extern "C" comb_cuda_status_t comb_cuda_verify(
    comb_cuda_context_t *context, comb_cuda_witness_mode_t mode,
    comb_cuda_result_t *result, comb_cuda_witness_t *witness) {
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
      return "CUDA scientific pipeline not implemented";
    case COMB_CUDA_STATUS_INTERNAL_ERROR:
      return "internal error";
    default:
      return "unknown CUDA status";
  }
}
