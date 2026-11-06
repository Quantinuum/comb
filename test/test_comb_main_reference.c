#include "comb_external.h"
#include "comb_rank.h"
#include "comb_reference.h"

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define CHECK(condition)                                                       \
  do {                                                                         \
    if (!(condition)) {                                                        \
      fprintf(stderr, "check failed at %s:%d: %s\n", __FILE__, __LINE__,     \
              #condition);                                                     \
      exit(EXIT_FAILURE);                                                      \
    }                                                                          \
  } while (0)

static uint64_t generator_state = UINT64_C(0xe7037ed1a0b428db);

static uint64_t next_random(void) {
  uint64_t value = generator_state;
  value ^= value >> 12U;
  value ^= value << 25U;
  value ^= value >> 27U;
  generator_state = value;
  return value * UINT64_C(2685821657736338717);
}

static comb_cuda_result_t fresh_result(void) {
  const comb_cuda_result_t value = COMB_CUDA_RESULT_INITIALIZER;
  return value;
}

static comb_cuda_witness_t fresh_witness(uint32_t *faults,
                                         uint32_t capacity) {
  comb_cuda_witness_t value = COMB_CUDA_WITNESS_INITIALIZER;
  value.fault_indices = faults;
  value.fault_capacity = capacity;
  return value;
}

static uint64_t main_cpu_witness_rank(uint32_t location_count,
                                      const uint32_t *partitions, uint32_t s,
                                      const uint32_t *faults) {
  comb_rank_plan_t *plan = NULL;
  uint32_t locations[3] = {0U, 0U, 0U};
  uint32_t offsets[3] = {0U, 0U, 0U};
  uint64_t rank = UINT64_MAX;
  uint32_t selected;
  CHECK(s <= 3U);
  CHECK(comb_rank_plan_create(location_count, partitions, s, &plan) ==
        COMB_RANK_OK);
  for (selected = 0U; selected < s; ++selected) {
    uint32_t location;
    uint32_t prefix = 0U;
    unsigned found = 0U;
    for (location = 0U; location < location_count; ++location) {
      const uint32_t end = prefix + partitions[location];
      if (faults[selected] >= prefix && faults[selected] < end) {
        locations[selected] = location;
        offsets[selected] = faults[selected] - prefix;
        found = 1U;
        break;
      }
      prefix = end;
    }
    CHECK(found != 0U);
    if (selected != 0U) CHECK(locations[selected - 1U] < locations[selected]);
  }
  CHECK(comb_rank_plan_rank(plan, locations, offsets, &rank) == COMB_RANK_OK);
  comb_rank_plan_destroy(plan);
  return rank;
}

static comb_code_t run_differential_case(uint32_t case_index) {
  const uint32_t n = 3U + (case_index & 1U);
  const unsigned measured = case_index >= 12U;
  const uint32_t n_all = n + (measured != 0U ? 2U : 0U);
  const uint32_t n_x = 1U;
  const uint32_t n_z = n - 2U;
  const uint32_t k = 1U;
  const uint32_t location_count = 3U;
  const uint32_t s = (case_index % 3U) == 0U ? 1U : 2U;
  uint32_t partitions[3];
  uint32_t output_count;
  comb_pauli_x_t *x_stabilisers;
  comb_pauli_z_t *z_stabilisers;
  comb_dcp_t *logical_dcps;
  comb_dcp_t *output_errors;
  comb_pauli_x_t *x_measurements = NULL;
  comb_pauli_z_t *z_measurements = NULL;
  comb_bit_array_t x_lookup;
  comb_bit_array_t z_lookup;
  const uint64_t x_lookup_count = UINT64_C(1) << (n - n_x);
  const uint64_t z_lookup_count = UINT64_C(1) << (n - n_z);
  comb_problem_t problem = COMB_PROBLEM_INITIALIZER;
  comb_reference_context_t *reference = NULL;
  comb_cuda_result_t result;
  uint32_t reference_faults[3] = {UINT32_MAX, UINT32_MAX, UINT32_MAX};
  uint32_t main_cpu_faults[3] = {UINT32_MAX, UINT32_MAX, UINT32_MAX};
  comb_cuda_witness_t witness = fresh_witness(reference_faults, 3U);
  comb_code_t main_cpu_code;
  uint32_t index;
  uint32_t bit;

  /*
   * With three locations and s <= 2, this fixture is smaller than main's
   * chunk threshold.  Witness equality is therefore meaningful here.  For a
   * multi-chunk problem main does not specify the globally lowest witness;
   * only its scientific decision and per-candidate behavior are the oracle.
   */

  partitions[0] = 1U + (uint32_t)(next_random() % 3U);
  partitions[1] = 1U + (uint32_t)(next_random() % 3U);
  partitions[2] = 1U + (uint32_t)(next_random() % 3U);
  output_count = partitions[0] + partitions[1] + partitions[2];

  x_stabilisers = comb_calloc_vector((int)n_x, (int)n);
  z_stabilisers = comb_calloc_vector((int)n_z, (int)n);
  logical_dcps = comb_calloc_dcp((int)k, (int)n_all);
  output_errors = comb_calloc_dcp((int)output_count, (int)n_all);
  if (measured != 0U) {
    x_measurements = comb_calloc_vector(1, (int)n_all);
    z_measurements = comb_calloc_vector(1, (int)n_all);
    CHECK(x_measurements != NULL);
    CHECK(z_measurements != NULL);
    comb_set_cvector_aref(x_measurements, 0, (int)(n + 1U), 1);
    comb_set_cvector_aref(z_measurements, 0, (int)n, 1);
  }
  x_lookup = comb_alloc_bit_array((int)x_lookup_count);
  z_lookup = comb_alloc_bit_array((int)z_lookup_count);
  CHECK(x_stabilisers != NULL);
  CHECK(z_stabilisers != NULL);
  CHECK(logical_dcps != NULL);
  CHECK(output_errors != NULL);
  CHECK(x_lookup != NULL);
  CHECK(z_lookup != NULL);

  /* A one-logical-qubit CSS code: Z_0,...,Z_{n-3}, X_{n-2}. */
  comb_set_cvector_aref(x_stabilisers, 0, (int)(n - 2U), 1);
  for (index = 0U; index < n_z; ++index)
    comb_set_cvector_aref(z_stabilisers, (int)index, (int)index, 1);

  /*
   * A code-preserving logical DCP: diagonal action on logical qubit n-1,
   * together with a symmetric coupling to Z-stabilised qubit 0.  The latter
   * is trivial on the +Z_0 code subspace but exercises off-diagonal packing and
   * multiplication in both implementations.
   */
  comb_set_cdcp_F_aref(logical_dcps, 0, (int)(n - 1U),
                       (int)(n - 1U), 1);
  comb_set_cdcp_F_aref(logical_dcps, 0, 0, (int)(n - 1U), 1);
  comb_set_cdcp_F_aref(logical_dcps, 0, (int)(n - 1U), 0, 1);

  for (index = 0U; index < output_count; ++index) {
    uint32_t row;
    for (bit = 0U; bit < n_all; ++bit) {
      comb_set_cdcp_h_x_aref(output_errors, (int)index, (int)bit,
                             (int)((next_random() >> 7U) & 1U));
      comb_set_cdcp_h_z_aref(output_errors, (int)index, (int)bit,
                             (int)((next_random() >> 19U) & 1U));
    }
    for (row = 0U; row < n_all; ++row) {
      uint32_t column;
      for (column = row; column < n_all; ++column) {
        const int value = (int)((next_random() >> 23U) & 1U);
        comb_set_cdcp_F_aref(output_errors, (int)index, (int)row,
                             (int)column, value);
        comb_set_cdcp_F_aref(output_errors, (int)index, (int)column,
                             (int)row, value);
      }
    }
    /* Guarantee both kinds of F term in every generated output DCP. */
    comb_set_cdcp_F_aref(output_errors, (int)index,
                         (int)(index % n), (int)(index % n), 1);
    comb_set_cdcp_F_aref(output_errors, (int)index, 0, (int)(n - 1U), 1);
    comb_set_cdcp_F_aref(output_errors, (int)index, (int)(n - 1U), 0, 1);
  }

  /*
   * Include known coverage shapes without asserting their scientific result:
   * one permissive table, one rejecting X table, then seeded mixed tables.
   */
  for (uint64_t lookup = 0U; lookup < x_lookup_count; ++lookup) {
    unsigned value;
    if (case_index == 0U)
      value = 1U;
    else if (case_index == 1U)
      value = 0U;
    else
      value = (unsigned)((next_random() >> 11U) & 1U);
    comb_set_bit_array_aref(x_lookup, (int)lookup, (int)value);
  }
  for (uint64_t lookup = 0U; lookup < z_lookup_count; ++lookup) {
    const unsigned value = case_index == 0U
                               ? 1U
                               : (unsigned)((next_random() >> 13U) & 1U);
    comb_set_bit_array_aref(z_lookup, (int)lookup, (int)value);
  }

  problem.n_all = n_all;
  problem.n = n;
  problem.n_x = n_x;
  problem.n_z = n_z;
  problem.k = k;
  problem.x_stabilisers = x_stabilisers;
  problem.x_lookup_table = x_lookup;
  problem.x_lookup_bit_count = x_lookup_count;
  problem.z_stabilisers = z_stabilisers;
  problem.z_lookup_table = z_lookup;
  problem.z_lookup_bit_count = z_lookup_count;
  problem.dcp_stabilisers = logical_dcps;
  problem.fault_location_count = location_count;
  problem.fault_location_partitions = partitions;
  problem.s = s;
  problem.output_error_count = output_count;
  problem.output_errors = output_errors;
  problem.x_measurement_count = measured;
  problem.x_measurement_paulis = x_measurements;
  problem.z_measurement_count = measured;
  problem.z_measurement_paulis = z_measurements;

  CHECK(comb_reference_context_create(&problem, &reference) ==
        COMB_CUDA_STATUS_SUCCESS);
  result = fresh_result();
  CHECK(comb_reference_verify(reference, COMB_CUDA_WITNESS_DETERMINISTIC,
                              &result, &witness) ==
        COMB_CUDA_STATUS_SUCCESS);

  main_cpu_code = comb_verify_css_magic_state(
      n_all, n, n_x, n_z, k, x_stabilisers, x_lookup, z_stabilisers, z_lookup,
      logical_dcps, location_count, partitions, s, main_cpu_faults,
      output_count, output_errors, measured, x_measurements,
      measured, z_measurements, 1U);
  CHECK(main_cpu_code == result.verification_code);
  if (main_cpu_code == BAD_FAULT_FOUND) {
    const uint64_t rank = main_cpu_witness_rank(
        location_count, partitions, s, main_cpu_faults);
    CHECK(result.minimum_bad_rank == rank);
    CHECK(witness.rank == rank);
    CHECK(witness.fault_count == s);
    for (index = 0U; index < s; ++index)
      CHECK(reference_faults[index] == main_cpu_faults[index]);
  } else {
    CHECK(main_cpu_code == ALL_FAULTS_GOOD);
    CHECK(result.minimum_bad_rank == COMB_CUDA_RANK_NONE);
    CHECK(witness.fault_count == 0U);
  }

  comb_reference_context_destroy(reference);
  comb_free_bit_array(z_lookup);
  comb_free_bit_array(x_lookup);
  comb_cfree_dcp((int)output_count, output_errors);
  comb_cfree_dcp((int)k, logical_dcps);
  comb_cfree_vector((int)n_z, z_stabilisers);
  comb_cfree_vector((int)n_x, x_stabilisers);
  if (measured != 0U) {
    comb_cfree_vector(1, x_measurements);
    comb_cfree_vector(1, z_measurements);
  }
  return main_cpu_code;
}

int main(void) {
  unsigned saw_good = 0U;
  unsigned saw_bad = 0U;
  uint32_t case_index;
  for (case_index = 0U; case_index < 24U; ++case_index) {
    const comb_code_t decision = run_differential_case(case_index);
    if (decision == ALL_FAULTS_GOOD) saw_good = 1U;
    if (decision == BAD_FAULT_FOUND) saw_bad = 1U;
  }
  CHECK(saw_good != 0U);
  CHECK(saw_bad != 0U);
  puts("comb main CPU/reference: 24 synthetic CSS differential cases passed");
  puts("coverage note: no user-supplied production corpus was available; "
       "this gate does not replace production-corpus differential testing");
  return EXIT_SUCCESS;
}
