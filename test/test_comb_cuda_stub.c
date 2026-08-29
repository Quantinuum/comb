#include "comb_external.h"
#include "comb_cuda.h"
#include "comb_cuda_internal.h"

#ifdef NDEBUG
#undef NDEBUG
#endif
#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

int main(void) {
  comb_problem_t problem = COMB_PROBLEM_INITIALIZER;
  comb_cuda_options_t options = COMB_CUDA_OPTIONS_INITIALIZER;
  comb_cuda_result_t result = COMB_CUDA_RESULT_INITIALIZER;
  comb_cuda_witness_t witness = COMB_CUDA_WITNESS_INITIALIZER;
  comb_cuda_context_t *context = (comb_cuda_context_t *)(uintptr_t)1;
  const uint32_t partitions[1] = {1U};
  comb_bit_word_t x_lookup[1] = {UINT64_C(3)};
  comb_bit_word_t z_lookup[1] = {UINT64_C(3)};
  comb_bit_word_t logical_f[1] = {0U};
  comb_bit_word_t logical_hx[1] = {0U};
  comb_bit_word_t logical_hz[1] = {0U};
  comb_bit_word_t output_f[1] = {0U};
  comb_bit_word_t output_hx[1] = {0U};
  comb_bit_word_t output_hz[1] = {0U};
  const comb_dcp_t logical = {
      {1U, 1U, logical_f}, {1U, logical_hx}, {1U, logical_hz}};
  const comb_dcp_t output = {
      {1U, 1U, output_f}, {1U, output_hx}, {1U, output_hz}};
  char pci_bus_id[32] = "not-cleared";

  assert(problem.struct_size == sizeof(problem));
  assert(problem.abi_version == COMB_CUDA_ABI_VERSION);
  assert(options.struct_size == sizeof(options));
  assert(options.cuda_device == -1);
  assert(options.witness_mode == COMB_CUDA_WITNESS_DETERMINISTIC);
  assert(result.verification_code == OTHER_ERROR);
  assert(result.minimum_bad_rank == COMB_CUDA_RANK_NONE);
  assert(witness.rank == COMB_CUDA_RANK_NONE);

  /* The stub preserves host-side validation before reporting no runtime. */
  assert(comb_cuda_context_create(&problem, &options, &context) ==
         COMB_CUDA_STATUS_INVALID_PROBLEM);
  assert(context == NULL);

  problem.n_all = 1U;
  problem.n = 1U;
  problem.k = 1U;
  problem.x_lookup_table = x_lookup;
  problem.x_lookup_bit_count = 2U;
  problem.z_lookup_table = z_lookup;
  problem.z_lookup_bit_count = 2U;
  problem.dcp_stabilisers = &logical;
  problem.fault_location_count = 1U;
  problem.fault_location_partitions = partitions;
  problem.s = 1U;
  problem.output_error_count = 1U;
  problem.output_errors = &output;

  context = (comb_cuda_context_t *)(uintptr_t)1;
  assert(comb_cuda_context_create(&problem, &options, &context) ==
         COMB_CUDA_STATUS_CUDA_UNAVAILABLE);
  assert(context == NULL);
  assert(comb_cuda_context_create(&problem, &options, NULL) ==
         COMB_CUDA_STATUS_INVALID_ARGUMENT);
  assert(comb_cuda_verify_range(NULL, 0U, 0U,
                                COMB_CUDA_WITNESS_DETERMINISTIC, &result,
                                &witness) == COMB_CUDA_STATUS_INVALID_ARGUMENT);
  assert(comb_cuda_verify(NULL, COMB_CUDA_WITNESS_DETERMINISTIC, &result,
                          &witness) == COMB_CUDA_STATUS_INVALID_ARGUMENT);
  assert(comb_cuda_device_pci_bus_id(-1, pci_bus_id, sizeof(pci_bus_id)) ==
         COMB_CUDA_STATUS_INVALID_ARGUMENT);
  assert(comb_cuda_device_pci_bus_id(0, pci_bus_id, sizeof(pci_bus_id)) ==
         COMB_CUDA_STATUS_CUDA_UNAVAILABLE);
  assert(pci_bus_id[0] == '\0');
  assert(strcmp(comb_cuda_status_string(COMB_CUDA_STATUS_CUDA_UNAVAILABLE),
                "CUDA unavailable") == 0);
  assert(strcmp(comb_cuda_status_string((comb_cuda_status_t)999),
                "unknown CUDA status") == 0);

  comb_cuda_context_destroy(NULL);
  puts("comb_cuda_stub: all tests passed");
  return 0;
}
