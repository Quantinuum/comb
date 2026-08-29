#ifdef NDEBUG
#undef NDEBUG
#endif
#include <assert.h>
#include <stdint.h>
#include <stdlib.h>

#include "comb_cuda_internal.h"

static comb_bit_word_t *legacy_alloc(uint64_t bit_count) {
  uint64_t word_count = (bit_count + sizeof(comb_bit_word_t) - 1) /
                        sizeof(comb_bit_word_t);
  if (word_count == 0) word_count = 1;
  comb_bit_word_t *data =
      (comb_bit_word_t *)calloc((size_t)word_count, sizeof(*data));
  assert(data != NULL);
  return data;
}

static void legacy_set(comb_bit_word_t *data, uint64_t bit_index) {
  data[bit_index / sizeof(comb_bit_word_t)] |=
      UINT64_C(1) << (bit_index % sizeof(comb_bit_word_t));
}

static comb_vector_t *vectors_create(uint32_t count, uint32_t width) {
  comb_vector_t *vectors =
      (comb_vector_t *)calloc(count == 0 ? 1 : count, sizeof(*vectors));
  assert(vectors != NULL);
  for (uint32_t index = 0; index < count; ++index) {
    vectors[index].n = width;
    vectors[index].data = legacy_alloc(width);
  }
  return vectors;
}

static void vectors_destroy(comb_vector_t *vectors, uint32_t count) {
  for (uint32_t index = 0; index < count; ++index) free(vectors[index].data);
  free(vectors);
}

static comb_dcp_t *dcps_create(uint32_t count, uint32_t width) {
  comb_dcp_t *dcps =
      (comb_dcp_t *)calloc(count == 0 ? 1 : count, sizeof(*dcps));
  assert(dcps != NULL);
  for (uint32_t index = 0; index < count; ++index) {
    dcps[index].F.m = width;
    dcps[index].F.n = width;
    dcps[index].F.data = legacy_alloc((uint64_t)width * width);
    dcps[index].h_x.n = width;
    dcps[index].h_x.data = legacy_alloc(width);
    dcps[index].h_z.n = width;
    dcps[index].h_z.data = legacy_alloc(width);
  }
  return dcps;
}

static void dcps_destroy(comb_dcp_t *dcps, uint32_t count) {
  for (uint32_t index = 0; index < count; ++index) {
    free(dcps[index].F.data);
    free(dcps[index].h_x.data);
    free(dcps[index].h_z.data);
  }
  free(dcps);
}

static void expect_pack_status(const comb_problem_t *problem,
                               comb_cuda_status_t expected) {
  void *blob = NULL;
  uint64_t blob_size = 0;
  assert(comb_cuda_flat_problem_pack(problem, &blob, &blob_size) == expected);
  comb_cuda_flat_problem_free(blob);
}

int main(void) {
  enum {
    N_ALL = 65,
    N = 33,
    N_X = 31,
    N_Z = 32,
    K = 1,
    OUTPUT_COUNT = 6
  };
  const uint32_t partitions[3] = {2, 3, 1};
  comb_vector_t *x_stabilisers = vectors_create(N_X, N);
  comb_vector_t *z_stabilisers = vectors_create(N_Z, N);
  comb_vector_t *x_measurements = vectors_create(1, N_ALL);
  comb_vector_t *z_measurements = vectors_create(1, N_ALL);
  comb_dcp_t *logical = dcps_create(K, N_ALL);
  comb_dcp_t *output = dcps_create(OUTPUT_COUNT, N_ALL);
  comb_bit_word_t *x_lookup = legacy_alloc(4);
  comb_bit_word_t *z_lookup = legacy_alloc(2);

  legacy_set(x_stabilisers[0].data, 32);
  legacy_set(z_stabilisers[0].data, 0);
  legacy_set(x_measurements[0].data, 64);
  legacy_set(z_measurements[0].data, 64);
  legacy_set(x_lookup, 0);
  legacy_set(x_lookup, 3);
  legacy_set(z_lookup, 0);
  legacy_set(z_lookup, 1);

  /* Row 1 begins at legacy bit 65, so this exercises a non-byte-aligned row. */
  legacy_set(output[0].F.data, (uint64_t)N_ALL + 64);
  legacy_set(output[0].h_x.data, 64);
  legacy_set(logical[0].F.data, (uint64_t)2 * N_ALL + 32);

  comb_problem_t problem = COMB_PROBLEM_INITIALIZER;
  problem.n_all = N_ALL;
  problem.n = N;
  problem.n_x = N_X;
  problem.n_z = N_Z;
  problem.k = K;
  problem.x_stabilisers = x_stabilisers;
  problem.x_lookup_table = x_lookup;
  problem.x_lookup_bit_count = 4;
  problem.z_stabilisers = z_stabilisers;
  problem.z_lookup_table = z_lookup;
  problem.z_lookup_bit_count = 2;
  problem.dcp_stabilisers = logical;
  problem.fault_location_count = 3;
  problem.fault_location_partitions = partitions;
  problem.s = 2;
  problem.output_error_count = OUTPUT_COUNT;
  problem.output_errors = output;
  problem.x_measurement_count = 1;
  problem.x_measurement_paulis = x_measurements;
  problem.z_measurement_count = 1;
  problem.z_measurement_paulis = z_measurements;

  void *blob = NULL;
  uint64_t blob_size = 0;
  assert(comb_cuda_flat_problem_pack(&problem, &blob, &blob_size) ==
         COMB_CUDA_STATUS_SUCCESS);
  assert(comb_cuda_flat_problem_validate(blob, blob_size) ==
         COMB_CUDA_STATUS_SUCCESS);
  /* Structural FlatProblem validation is intentionally distinct from the
   * scientific contract; this packing fixture has non-RREF generators and an
   * asymmetric DCP matrix. */
  assert(comb_cuda_flat_problem_validate_science(blob, blob_size) ==
         COMB_CUDA_STATUS_INVALID_PROBLEM);

  comb_cuda_flat_problem_t *header = (comb_cuda_flat_problem_t *)blob;
  assert(header->n == 33);
  assert(header->n_all == 65);
  assert(header->words_code == 1);
  assert(header->words_all == 2);
  assert(header->rank_dp_rows == 4);
  assert(header->rank_dp_columns == 3);
  assert(header->total_combination_count == 11);

  const uint64_t expected_dp[12] = {1, 0, 0, 1, 2, 0,
                                    1, 5, 6, 1, 6, 11};
  uint64_t *rank_dp =
      (uint64_t *)((unsigned char *)blob + header->rank_dp_offset);
  for (uint32_t index = 0; index < 12; ++index) {
    assert(rank_dp[index] == expected_dp[index]);
  }

  const uint64_t *prefixes =
      (const uint64_t *)((const unsigned char *)blob +
                         header->fault_prefixes_offset);
  assert(prefixes[0] == 0 && prefixes[1] == 2 && prefixes[2] == 5 &&
         prefixes[3] == 6);
  const uint64_t *packed_x_stabilisers =
      (const uint64_t *)((const unsigned char *)blob +
                         header->x_stabilisers_offset);
  assert(packed_x_stabilisers[0] == (UINT64_C(1) << 32));
  const uint64_t *packed_x_measurements =
      (const uint64_t *)((const unsigned char *)blob +
                         header->x_measurements_offset);
  assert(packed_x_measurements[0] == 0);
  assert(packed_x_measurements[1] == 1);
  const uint64_t *packed_output_f =
      (const uint64_t *)((const unsigned char *)blob +
                         header->output_dcp_f_offset);
  assert(packed_output_f[2] == 0);
  assert(packed_output_f[3] == 1);
  const uint64_t *packed_output_hx =
      (const uint64_t *)((const unsigned char *)blob +
                         header->output_dcp_hx_offset);
  assert(packed_output_hx[0] == 0 && packed_output_hx[1] == 1);
  const uint64_t *packed_x_lookup =
      (const uint64_t *)((const unsigned char *)blob + header->x_lookup_offset);
  assert(packed_x_lookup[0] == UINT64_C(9));

  const uint64_t saved_magic = header->magic;
  header->magic = 0;
  assert(comb_cuda_flat_problem_validate(blob, blob_size) ==
         COMB_CUDA_STATUS_INVALID_PROBLEM);
  header->magic = saved_magic;

  const uint64_t saved_lookup_offset = header->x_lookup_offset;
  header->x_lookup_offset += 8;
  assert(comb_cuda_flat_problem_validate(blob, blob_size) ==
         COMB_CUDA_STATUS_INVALID_PROBLEM);
  header->x_lookup_offset = saved_lookup_offset;

  rank_dp[11] ^= 1;
  assert(comb_cuda_flat_problem_validate(blob, blob_size) ==
         COMB_CUDA_STATUS_INVALID_PROBLEM);
  rank_dp[11] ^= 1;

  uint64_t *mutable_x_stabilisers =
      (uint64_t *)((unsigned char *)blob + header->x_stabilisers_offset);
  mutable_x_stabilisers[0] |= UINT64_C(1) << 40;
  assert(comb_cuda_flat_problem_validate(blob, blob_size) ==
         COMB_CUDA_STATUS_INVALID_PROBLEM);
  mutable_x_stabilisers[0] &= ~(UINT64_C(1) << 40);
  assert(comb_cuda_flat_problem_validate(blob, blob_size) ==
         COMB_CUDA_STATUS_SUCCESS);

  comb_cuda_flat_problem_free(blob);

  {
    const uint32_t scientific_partitions[1] = {1};
    comb_vector_t *scientific_x_stabilisers = vectors_create(1, 2);
    comb_dcp_t *scientific_logical = dcps_create(1, 2);
    comb_dcp_t *scientific_output = dcps_create(1, 2);
    comb_bit_word_t *scientific_x_lookup = legacy_alloc(2);
    comb_bit_word_t *scientific_z_lookup = legacy_alloc(4);
    comb_problem_t scientific_problem = COMB_PROBLEM_INITIALIZER;
    void *scientific_blob = NULL;
    uint64_t scientific_blob_size = 0;
    comb_cuda_flat_problem_t *scientific_header;
    uint64_t *scientific_output_f;
    uint64_t *scientific_x_rows;

    legacy_set(scientific_x_stabilisers[0].data, 0);
    scientific_problem.n_all = 2;
    scientific_problem.n = 2;
    scientific_problem.n_x = 1;
    scientific_problem.n_z = 0;
    scientific_problem.k = 1;
    scientific_problem.x_stabilisers = scientific_x_stabilisers;
    scientific_problem.x_lookup_table = scientific_x_lookup;
    scientific_problem.x_lookup_bit_count = 2;
    scientific_problem.z_stabilisers = NULL;
    scientific_problem.z_lookup_table = scientific_z_lookup;
    scientific_problem.z_lookup_bit_count = 4;
    scientific_problem.dcp_stabilisers = scientific_logical;
    scientific_problem.fault_location_count = 1;
    scientific_problem.fault_location_partitions = scientific_partitions;
    scientific_problem.s = 1;
    scientific_problem.output_error_count = 1;
    scientific_problem.output_errors = scientific_output;

    assert(comb_cuda_flat_problem_pack(
               &scientific_problem, &scientific_blob,
               &scientific_blob_size) == COMB_CUDA_STATUS_SUCCESS);
    assert(comb_cuda_flat_problem_validate_science(
               scientific_blob, scientific_blob_size) ==
           COMB_CUDA_STATUS_SUCCESS);
    scientific_header = (comb_cuda_flat_problem_t *)scientific_blob;
    scientific_output_f =
        (uint64_t *)((unsigned char *)scientific_blob +
                     scientific_header->output_dcp_f_offset);
    scientific_output_f[0] |= UINT64_C(1) << 1;
    assert(comb_cuda_flat_problem_validate_science(
               scientific_blob, scientific_blob_size) ==
           COMB_CUDA_STATUS_INVALID_PROBLEM);
    scientific_output_f[1] |= UINT64_C(1);
    assert(comb_cuda_flat_problem_validate_science(
               scientific_blob, scientific_blob_size) ==
           COMB_CUDA_STATUS_SUCCESS);

    scientific_x_rows =
        (uint64_t *)((unsigned char *)scientific_blob +
                     scientific_header->x_stabilisers_offset);
    scientific_x_rows[0] = 0;
    assert(comb_cuda_flat_problem_validate_science(
               scientific_blob, scientific_blob_size) ==
           COMB_CUDA_STATUS_INVALID_PROBLEM);

    comb_cuda_flat_problem_free(scientific_blob);
    free(scientific_x_lookup);
    free(scientific_z_lookup);
    vectors_destroy(scientific_x_stabilisers, 1);
    dcps_destroy(scientific_logical, 1);
    dcps_destroy(scientific_output, 1);
  }

  problem.x_lookup_bit_count = 3;
  expect_pack_status(&problem, COMB_CUDA_STATUS_INVALID_PROBLEM);
  problem.x_lookup_bit_count = 4;

  x_stabilisers[0].n = 32;
  expect_pack_status(&problem, COMB_CUDA_STATUS_INVALID_PROBLEM);
  x_stabilisers[0].n = N;

  uint32_t zero_partitions[3] = {2, 0, 1};
  problem.fault_location_partitions = zero_partitions;
  expect_pack_status(&problem, COMB_CUDA_STATUS_INVALID_PROBLEM);
  problem.fault_location_partitions = partitions;

  /* A fully allocated input large enough to overflow W is not practical in a
   * unit test because output_error_count also requires that many validated
   * DCPs.  Corrupted near-overflow DP values are nevertheless rejected above. */

  free(x_lookup);
  free(z_lookup);
  vectors_destroy(x_stabilisers, N_X);
  vectors_destroy(z_stabilisers, N_Z);
  vectors_destroy(x_measurements, 1);
  vectors_destroy(z_measurements, 1);
  dcps_destroy(logical, K);
  dcps_destroy(output, OUTPUT_COUNT);
  return 0;
}
