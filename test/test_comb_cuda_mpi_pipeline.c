#define _POSIX_C_SOURCE 200809L
#include "comb_cuda_mpi.h"
#include "comb_reference.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define CHECK(x) do { if (!(x)) { \
  fprintf(stderr, "%s:%d: %s\n", __FILE__, __LINE__, #x); \
  MPI_Abort(MPI_COMM_WORLD, 1); } } while (0)

int main(int argc, char **argv) {
  int rank = -1;
  unsigned good_cases = 0, bad_cases = 0;
  comb_problem_t problem = COMB_PROBLEM_INITIALIZER;
  uint32_t partitions[] = {2, 1, 3};
  comb_dcp_t errors[6] = {0};
  comb_dcp_t logical = {0};
  comb_vector_t measurement = {0};
  comb_bit_word_t f[6] = {0}, hx[6] = {0}, hz[6] = {0};
  comb_bit_word_t logical_f = 0, logical_hx = 0, logical_hz = 1;
  comb_bit_word_t measure_z = 2;
  comb_bit_word_t xl = 3, zl = 3;
  CHECK(MPI_Init(&argc, &argv) == MPI_SUCCESS);
  CHECK(MPI_Comm_rank(MPI_COMM_WORLD, &rank) == MPI_SUCCESS);
#ifndef COMB_TEST_CUDA_DEVICE
  CHECK(setenv("CUDA_VISIBLE_DEVICES", "host-test-gpu", 1) == 0);
#endif
  problem.n_all = 2;
  problem.n = 1;
  problem.k = 1;
  problem.fault_location_count = 3;
  problem.fault_location_partitions = partitions;
  problem.s = 2;
  problem.output_error_count = 6;
  problem.output_errors = errors;
  problem.dcp_stabilisers = &logical;
  problem.x_lookup_table = &xl;
  problem.z_lookup_table = &zl;
  problem.x_lookup_bit_count = problem.z_lookup_bit_count = 2;
  problem.z_measurement_count = 1;
  problem.z_measurement_paulis = &measurement;
  measurement.n = 2;
  measurement.data = &measure_z;
  logical.F.m = logical.F.n = logical.h_x.n = logical.h_z.n = 2;
  logical.F.data = &logical_f;
  logical.h_x.data = &logical_hx;
  logical.h_z.data = &logical_hz;
  for (unsigned i = 0; i < 6; ++i) {
    errors[i].F.m = errors[i].F.n = 2;
    errors[i].h_x.n = errors[i].h_z.n = 2;
    errors[i].F.data = &f[i];
    errors[i].h_x.data = &hx[i];
    errors[i].h_z.data = &hz[i];
  }
  for (unsigned case_index = 0; case_index < 12; ++case_index) {
    comb_reference_context_t *reference = NULL;
    comb_cuda_result_t expected = COMB_CUDA_RESULT_INITIALIZER;
    comb_cuda_mpi_result_t actual = COMB_CUDA_MPI_RESULT_INITIALIZER;
    comb_cuda_mpi_options_t options = COMB_CUDA_MPI_OPTIONS_INITIALIZER;
    comb_cuda_witness_t expected_witness = COMB_CUDA_WITNESS_INITIALIZER;
    comb_cuda_witness_t actual_witness = COMB_CUDA_WITNESS_INITIALIZER;
    uint32_t expected_faults[2] = {0}, actual_faults[2] = {0};
    expected_witness.fault_capacity = actual_witness.fault_capacity = 2;
    expected_witness.fault_indices = expected_faults;
    actual_witness.fault_indices = actual_faults;
    xl = case_index % 3 == 0 ? 3 : 1;
    zl = case_index % 3 == 1 || case_index == 0 ? 3 : 1;
    for (unsigned i = 0; i < 6; ++i) {
      const unsigned v = case_index * 17 + i * 7;
      f[i] = (v & 1U) | ((v & 2U) ? 6U : 0U) | ((v & 4U) ? 8U : 0U);
      hx[i] = (v >> 2) & 3U;
      hz[i] = (v >> 4) & 3U;
    }
    options.epoch_size = 4;
    options.chunk_size = 1;
    CHECK(comb_reference_context_create(&problem, &reference) ==
          COMB_CUDA_STATUS_SUCCESS);
    CHECK(comb_reference_verify(reference, COMB_CUDA_WITNESS_DETERMINISTIC,
                                &expected, &expected_witness) ==
          COMB_CUDA_STATUS_SUCCESS);
    CHECK(comb_cuda_mpi_verify_comm(
        MPI_COMM_WORLD, rank == 0 ? &problem : NULL,
        rank == 0 ? &options : NULL, &actual,
        rank == 0 ? &actual_witness : NULL) == COMB_CUDA_MPI_STATUS_SUCCESS);
    CHECK(actual.verification_code == expected.verification_code);
    if (expected.verification_code == BAD_FAULT_FOUND) ++bad_cases;
    else ++good_cases;
    CHECK(actual.minimum_bad_rank == expected.minimum_bad_rank);
    CHECK(actual.total_combination_count == expected.total_combination_count);
    {
      uint64_t checked = expected.total_combination_count;
      if (expected.verification_code == BAD_FAULT_FOUND) {
        checked = (expected.minimum_bad_rank / options.epoch_size + 1) * options.epoch_size;
        if (checked > expected.total_combination_count)
          checked = expected.total_combination_count;
      }
      CHECK(actual.checked_combination_count == checked);
    }
    if (rank == 0 && expected.verification_code == BAD_FAULT_FOUND) {
      CHECK(actual_witness.rank == expected_witness.rank);
      CHECK(actual_witness.fault_count == 2);
      CHECK(memcmp(actual_faults, expected_faults, sizeof(actual_faults)) == 0);
    }
    comb_reference_context_destroy(reference);
  }
  CHECK(good_cases != 0 && bad_cases != 0);
  CHECK(MPI_Finalize() == MPI_SUCCESS);
  if (rank == 0) {
#ifdef COMB_TEST_CUDA_DEVICE
    printf("Native CUDA MPI matches scientific reference: %u good and %u bad cases\n",
           good_cases, bad_cases);
#else
    printf("Distributed candidate code matches scientific reference (host execution): %u good and %u bad cases\n",
           good_cases, bad_cases);
#endif
  }
  return 0;
}
