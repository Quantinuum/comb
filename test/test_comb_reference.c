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

#define CHECK_STATUS(expression, expected) CHECK((expression) == (expected))

typedef struct {
  comb_problem_t problem;
  uint32_t *partitions;
  comb_bit_word_t *x_lookup;
  comb_bit_word_t *z_lookup;
  comb_dcp_t *outputs;
  comb_dcp_t *logicals;
  comb_vector_t *x_stabilisers;
  comb_vector_t *z_stabilisers;
  comb_vector_t *x_measurements;
  comb_vector_t *z_measurements;
} problem_fixture_t;

static uint64_t random_state = UINT64_C(0xa0761d6478bd642f);

static uint64_t next_random(void) {
  uint64_t value = random_state;
  value ^= value >> 12U;
  value ^= value << 25U;
  value ^= value >> 27U;
  random_state = value;
  return value * UINT64_C(2685821657736338717);
}

static size_t external_word_count(uint64_t bit_count) {
  return (size_t)((bit_count / 8U) + ((bit_count % 8U) != 0U ? 1U : 0U));
}

static comb_bit_word_t *external_bits_allocate(uint64_t bit_count) {
  const size_t words = external_word_count(bit_count);
  if (words == 0U) return NULL;
  return (comb_bit_word_t *)calloc(words, sizeof(comb_bit_word_t));
}

static void external_bit_set(comb_bit_word_t *words, uint64_t bit,
                             unsigned value) {
  const comb_bit_word_t mask = UINT64_C(1) << (bit % 8U);
  CHECK(words != NULL);
  if (value != 0U)
    words[bit / 8U] |= mask;
  else
    words[bit / 8U] &= ~mask;
}

static unsigned external_bit_get(const comb_bit_word_t *words, uint64_t bit) {
  CHECK(words != NULL);
  return (unsigned)((words[bit / 8U] >> (bit % 8U)) & UINT64_C(1));
}

static void vector_array_allocate(comb_vector_t **vectors, uint32_t count,
                                  uint32_t width) {
  uint32_t index;
  *vectors = NULL;
  if (count == 0U) return;
  *vectors = (comb_vector_t *)calloc(count, sizeof(comb_vector_t));
  CHECK(*vectors != NULL);
  for (index = 0U; index < count; ++index) {
    (*vectors)[index].n = width;
    (*vectors)[index].data = external_bits_allocate(width);
    CHECK((*vectors)[index].data != NULL);
  }
}

static void vector_array_free(comb_vector_t *vectors, uint32_t count) {
  uint32_t index;
  if (vectors == NULL) return;
  for (index = 0U; index < count; ++index) free(vectors[index].data);
  free(vectors);
}

static void dcp_array_allocate(comb_dcp_t **dcps, uint32_t count,
                               uint32_t qubits) {
  uint32_t index;
  *dcps = NULL;
  if (count == 0U) return;
  *dcps = (comb_dcp_t *)calloc(count, sizeof(comb_dcp_t));
  CHECK(*dcps != NULL);
  for (index = 0U; index < count; ++index) {
    (*dcps)[index].F.m = qubits;
    (*dcps)[index].F.n = qubits;
    (*dcps)[index].F.data =
        external_bits_allocate((uint64_t)qubits * qubits);
    (*dcps)[index].h_x.n = qubits;
    (*dcps)[index].h_x.data = external_bits_allocate(qubits);
    (*dcps)[index].h_z.n = qubits;
    (*dcps)[index].h_z.data = external_bits_allocate(qubits);
    CHECK((*dcps)[index].F.data != NULL);
    CHECK((*dcps)[index].h_x.data != NULL);
    CHECK((*dcps)[index].h_z.data != NULL);
  }
}

static void dcp_array_free(comb_dcp_t *dcps, uint32_t count) {
  uint32_t index;
  if (dcps == NULL) return;
  for (index = 0U; index < count; ++index) {
    free(dcps[index].F.data);
    free(dcps[index].h_x.data);
    free(dcps[index].h_z.data);
  }
  free(dcps);
}

static void fixture_replace_lookup(problem_fixture_t *fixture, unsigned x,
                                   uint64_t bit_count) {
  comb_bit_word_t **storage = x != 0U ? &fixture->x_lookup : &fixture->z_lookup;
  free(*storage);
  *storage = external_bits_allocate(bit_count);
  CHECK(*storage != NULL);
  if (x != 0U) {
    fixture->problem.x_lookup_table = *storage;
    fixture->problem.x_lookup_bit_count = bit_count;
  } else {
    fixture->problem.z_lookup_table = *storage;
    fixture->problem.z_lookup_bit_count = bit_count;
  }
}

static void fixture_init(problem_fixture_t *fixture, uint32_t qubits,
                         uint32_t output_count) {
  uint64_t lookup_count;
  uint32_t bit;
  memset(fixture, 0, sizeof(*fixture));
  fixture->problem = (comb_problem_t)COMB_PROBLEM_INITIALIZER;
  fixture->partitions = (uint32_t *)malloc(sizeof(uint32_t));
  CHECK(fixture->partitions != NULL);
  fixture->partitions[0] = output_count;
  lookup_count = UINT64_C(1) << qubits;
  fixture->x_lookup = external_bits_allocate(lookup_count);
  fixture->z_lookup = external_bits_allocate(lookup_count);
  CHECK(fixture->x_lookup != NULL);
  CHECK(fixture->z_lookup != NULL);
  for (bit = 0U; bit < lookup_count; ++bit) {
    external_bit_set(fixture->x_lookup, bit, 1U);
    external_bit_set(fixture->z_lookup, bit, 1U);
  }
  dcp_array_allocate(&fixture->outputs, output_count, qubits);
  dcp_array_allocate(&fixture->logicals, qubits, qubits);
  fixture->problem.n_all = qubits;
  fixture->problem.n = qubits;
  fixture->problem.n_x = 0U;
  fixture->problem.n_z = 0U;
  fixture->problem.k = qubits;
  fixture->problem.x_stabilisers = NULL;
  fixture->problem.x_lookup_table = fixture->x_lookup;
  fixture->problem.x_lookup_bit_count = lookup_count;
  fixture->problem.z_stabilisers = NULL;
  fixture->problem.z_lookup_table = fixture->z_lookup;
  fixture->problem.z_lookup_bit_count = lookup_count;
  fixture->problem.dcp_stabilisers = fixture->logicals;
  fixture->problem.fault_location_count = 1U;
  fixture->problem.fault_location_partitions = fixture->partitions;
  fixture->problem.s = 1U;
  fixture->problem.output_error_count = output_count;
  fixture->problem.output_errors = fixture->outputs;
  fixture->problem.x_measurement_count = 0U;
  fixture->problem.x_measurement_paulis = NULL;
  fixture->problem.z_measurement_count = 0U;
  fixture->problem.z_measurement_paulis = NULL;
}

static void fixture_set_logical_count(problem_fixture_t *fixture,
                                      uint32_t count) {
  dcp_array_free(fixture->logicals, fixture->problem.k);
  fixture->logicals = NULL;
  dcp_array_allocate(&fixture->logicals, count, fixture->problem.n_all);
  fixture->problem.k = count;
  fixture->problem.dcp_stabilisers = fixture->logicals;
}

static void fixture_set_x_stabiliser_count(problem_fixture_t *fixture,
                                           uint32_t count) {
  vector_array_free(fixture->x_stabilisers, fixture->problem.n_x);
  fixture->x_stabilisers = NULL;
  vector_array_allocate(&fixture->x_stabilisers, count, fixture->problem.n);
  fixture->problem.n_x = count;
  fixture->problem.x_stabilisers = fixture->x_stabilisers;
  fixture_replace_lookup(fixture, 1U,
                         UINT64_C(1) << (fixture->problem.n - count));
  if (count <= fixture->problem.n - fixture->problem.n_z) {
    fixture_set_logical_count(
        fixture, fixture->problem.n - count - fixture->problem.n_z);
  }
}

static void fixture_set_z_stabiliser_count(problem_fixture_t *fixture,
                                           uint32_t count) {
  vector_array_free(fixture->z_stabilisers, fixture->problem.n_z);
  fixture->z_stabilisers = NULL;
  vector_array_allocate(&fixture->z_stabilisers, count, fixture->problem.n);
  fixture->problem.n_z = count;
  fixture->problem.z_stabilisers = fixture->z_stabilisers;
  fixture_replace_lookup(fixture, 0U,
                         UINT64_C(1) << (fixture->problem.n - count));
  if (count <= fixture->problem.n - fixture->problem.n_x) {
    fixture_set_logical_count(
        fixture, fixture->problem.n - fixture->problem.n_x - count);
  }
}

static void fixture_set_z_measurement_count(problem_fixture_t *fixture,
                                            uint32_t count) {
  vector_array_free(fixture->z_measurements,
                    fixture->problem.z_measurement_count);
  fixture->z_measurements = NULL;
  vector_array_allocate(&fixture->z_measurements, count,
                        fixture->problem.n_all);
  fixture->problem.z_measurement_count = count;
  fixture->problem.z_measurement_paulis = fixture->z_measurements;
}

static void fixture_destroy(problem_fixture_t *fixture) {
  dcp_array_free(fixture->outputs, fixture->problem.output_error_count);
  dcp_array_free(fixture->logicals, fixture->problem.k);
  vector_array_free(fixture->x_stabilisers, fixture->problem.n_x);
  vector_array_free(fixture->z_stabilisers, fixture->problem.n_z);
  vector_array_free(fixture->x_measurements,
                    fixture->problem.x_measurement_count);
  vector_array_free(fixture->z_measurements,
                    fixture->problem.z_measurement_count);
  free(fixture->x_lookup);
  free(fixture->z_lookup);
  free(fixture->partitions);
  memset(fixture, 0, sizeof(*fixture));
}

static void dcp_set_hx(comb_dcp_t *dcp, uint32_t bit, unsigned value) {
  external_bit_set(dcp->h_x.data, bit, value);
}

static void dcp_set_hz(comb_dcp_t *dcp, uint32_t bit, unsigned value) {
  external_bit_set(dcp->h_z.data, bit, value);
}

static void dcp_set_f(comb_dcp_t *dcp, uint32_t row, uint32_t column,
                      unsigned value) {
  external_bit_set(dcp->F.data, (uint64_t)row * dcp->F.n + column, value);
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

static void test_known_good_bad_and_ranges(void) {
  problem_fixture_t fixture;
  comb_reference_context_t *context = NULL;
  comb_cuda_result_t result;
  uint32_t fault = UINT32_MAX;
  comb_cuda_witness_t witness;
  fixture_init(&fixture, 1U, 3U);
  dcp_set_hx(&fixture.outputs[1], 0U, 1U);
  dcp_set_hx(&fixture.outputs[2], 0U, 1U);
  external_bit_set(fixture.x_lookup, 1U, 0U);
  CHECK_STATUS(comb_reference_context_create(&fixture.problem, &context),
               COMB_CUDA_STATUS_SUCCESS);

  result = fresh_result();
  witness = fresh_witness(&fault, 1U);
  CHECK_STATUS(comb_reference_verify(context, COMB_CUDA_WITNESS_DETERMINISTIC,
                                     &result, &witness),
               COMB_CUDA_STATUS_SUCCESS);
  CHECK(result.verification_code == BAD_FAULT_FOUND);
  CHECK(result.total_combination_count == 3U);
  CHECK(result.checked_combination_count == 3U);
  CHECK(result.minimum_bad_rank == 1U);
  CHECK(witness.rank == 1U);
  CHECK(witness.fault_count == 1U);
  CHECK(fault == 1U);

  result = fresh_result();
  witness = fresh_witness(&fault, 1U);
  CHECK_STATUS(comb_reference_verify_range(
                   context, 0U, 1U, COMB_CUDA_WITNESS_FAST_ANY, &result,
                   &witness),
               COMB_CUDA_STATUS_SUCCESS);
  CHECK(result.verification_code == ALL_FAULTS_GOOD);
  CHECK(result.checked_combination_count == 1U);
  CHECK(result.minimum_bad_rank == COMB_CUDA_RANK_NONE);
  CHECK(witness.fault_count == 0U);

  result = fresh_result();
  witness = fresh_witness(&fault, 1U);
  CHECK_STATUS(comb_reference_verify_range(
                   context, 2U, 3U, COMB_CUDA_WITNESS_DETERMINISTIC, &result,
                   &witness),
               COMB_CUDA_STATUS_SUCCESS);
  CHECK(result.verification_code == BAD_FAULT_FOUND);
  CHECK(result.checked_combination_count == 1U);
  CHECK(result.minimum_bad_rank == 2U);
  CHECK(witness.rank == 2U);
  CHECK(fault == 2U);

  result = fresh_result();
  witness = fresh_witness(&fault, 1U);
  CHECK_STATUS(comb_reference_verify_range(
                   context, 1U, 1U, COMB_CUDA_WITNESS_DETERMINISTIC, &result,
                   &witness),
               COMB_CUDA_STATUS_SUCCESS);
  CHECK(result.verification_code == ALL_FAULTS_GOOD);
  CHECK(result.checked_combination_count == 0U);
  CHECK(result.minimum_bad_rank == COMB_CUDA_RANK_NONE);

  result = fresh_result();
  CHECK_STATUS(comb_reference_verify_range(
                   context, 3U, 2U, COMB_CUDA_WITNESS_DETERMINISTIC, &result,
                   NULL),
               COMB_CUDA_STATUS_INVALID_ARGUMENT);
  comb_reference_context_destroy(context);
  fixture_destroy(&fixture);
}

static void test_z_lookup_and_logical_existential(void) {
  problem_fixture_t fixture;
  comb_reference_context_t *context = NULL;
  comb_cuda_result_t result;

  fixture_init(&fixture, 1U, 1U);
  dcp_set_hz(&fixture.outputs[0], 0U, 1U);
  external_bit_set(fixture.z_lookup, 1U, 0U);
  CHECK_STATUS(comb_reference_context_create(&fixture.problem, &context),
               COMB_CUDA_STATUS_SUCCESS);
  result = fresh_result();
  CHECK_STATUS(comb_reference_verify(context, COMB_CUDA_WITNESS_DEFAULT,
                                     &result, NULL),
               COMB_CUDA_STATUS_SUCCESS);
  CHECK(result.verification_code == BAD_FAULT_FOUND);
  comb_reference_context_destroy(context);
  fixture_destroy(&fixture);

  fixture_init(&fixture, 1U, 1U);
  dcp_set_hx(&fixture.outputs[0], 0U, 1U);
  external_bit_set(fixture.x_lookup, 1U, 0U);
  fixture_set_logical_count(&fixture, 1U);
  dcp_set_hx(&fixture.logicals[0], 0U, 1U);
  CHECK_STATUS(comb_reference_context_create(&fixture.problem, &context),
               COMB_CUDA_STATUS_SUCCESS);
  result = fresh_result();
  CHECK_STATUS(comb_reference_verify(context,
                                     COMB_CUDA_WITNESS_DETERMINISTIC, &result,
                                     NULL),
               COMB_CUDA_STATUS_SUCCESS);
  CHECK(result.verification_code == ALL_FAULTS_GOOD);
  CHECK(result.checked_combination_count == 1U);
  comb_reference_context_destroy(context);
  fixture_destroy(&fixture);
}

static void test_context_validation(void) {
  problem_fixture_t fixture;
  comb_reference_context_t *context = NULL;

  fixture_init(&fixture, 2U, 1U);
  dcp_set_f(&fixture.outputs[0], 0U, 1U, 1U);
  CHECK_STATUS(comb_reference_context_create(&fixture.problem, &context),
               COMB_CUDA_STATUS_INVALID_PROBLEM);
  CHECK(context == NULL);
  fixture_destroy(&fixture);

  fixture_init(&fixture, 1U, 1U);
  fixture_set_x_stabiliser_count(&fixture, 1U);
  CHECK_STATUS(comb_reference_context_create(&fixture.problem, &context),
               COMB_CUDA_STATUS_INVALID_PROBLEM);
  CHECK(context == NULL);
  fixture_destroy(&fixture);

  fixture_init(&fixture, 1U, 1U);
  fixture_set_x_stabiliser_count(&fixture, 1U);
  fixture_set_z_stabiliser_count(&fixture, 1U);
  external_bit_set(fixture.x_stabilisers[0].data, 0U, 1U);
  external_bit_set(fixture.z_stabilisers[0].data, 0U, 1U);
  CHECK_STATUS(comb_reference_context_create(&fixture.problem, &context),
               COMB_CUDA_STATUS_INVALID_PROBLEM);
  CHECK(context == NULL);
  fixture_destroy(&fixture);

  fixture_init(&fixture, 2U, 1U);
  fixture_set_x_stabiliser_count(&fixture, 2U);
  external_bit_set(fixture.x_stabilisers[0].data, 1U, 1U);
  external_bit_set(fixture.x_stabilisers[1].data, 0U, 1U);
  CHECK_STATUS(comb_reference_context_create(&fixture.problem, &context),
               COMB_CUDA_STATUS_INVALID_PROBLEM);
  CHECK(context == NULL);
  fixture_destroy(&fixture);

  fixture_init(&fixture, 2U, 1U);
  fixture_set_z_measurement_count(&fixture, 1U);
  external_bit_set(fixture.z_measurements[0].data, 0U, 1U);
  external_bit_set(fixture.z_measurements[0].data, 1U, 1U);
  CHECK_STATUS(comb_reference_context_create(&fixture.problem, &context),
               COMB_CUDA_STATUS_INVALID_PROBLEM);
  CHECK(context == NULL);
  fixture_destroy(&fixture);

  fixture_init(&fixture, 3U, 1U);
  fixture_set_z_measurement_count(&fixture, 2U);
  external_bit_set(fixture.z_measurements[0].data, 0U, 1U);
  external_bit_set(fixture.z_measurements[1].data, 0U, 1U);
  external_bit_set(fixture.z_measurements[1].data, 1U, 1U);
  external_bit_set(fixture.z_measurements[1].data, 2U, 1U);
  CHECK_STATUS(comb_reference_context_create(&fixture.problem, &context),
               COMB_CUDA_STATUS_INVALID_PROBLEM);
  CHECK(context == NULL);
  fixture_destroy(&fixture);

  fixture_init(&fixture, 2U, 1U);
  fixture_set_z_stabiliser_count(&fixture, 1U);
  external_bit_set(fixture.z_stabilisers[0].data, 0U, 1U);
  fixture_set_z_measurement_count(&fixture, 1U);
  external_bit_set(fixture.z_measurements[0].data, 1U, 1U);
  CHECK_STATUS(comb_reference_context_create(&fixture.problem, &context),
               COMB_CUDA_STATUS_SUCCESS);
  CHECK(context != NULL);
  comb_reference_context_destroy(context);
  context = NULL;
  fixture_destroy(&fixture);

  fixture_init(&fixture, 1U, 1U);
  fixture.problem.x_lookup_bit_count = 1U;
  CHECK_STATUS(comb_reference_context_create(&fixture.problem, &context),
               COMB_CUDA_STATUS_INVALID_PROBLEM);
  CHECK(context == NULL);
  fixture_destroy(&fixture);
}

static void test_sequential_measurement_contract(void) {
  problem_fixture_t fixture;
  comb_reference_context_t *context = NULL;
  comb_cuda_result_t result;
  fixture_init(&fixture, 1U, 2U);
  fixture_set_z_measurement_count(&fixture, 1U);
  external_bit_set(fixture.z_measurements[0].data, 0U, 1U);
  dcp_set_hx(&fixture.outputs[1], 0U, 1U);
  CHECK_STATUS(comb_reference_context_create(&fixture.problem, &context),
               COMB_CUDA_STATUS_SUCCESS);
  result = fresh_result();
  CHECK_STATUS(comb_reference_verify(context,
                                     COMB_CUDA_WITNESS_DETERMINISTIC, &result,
                                     NULL),
               COMB_CUDA_STATUS_SUCCESS);
  /* The second error anticommutes with the retained Z measurement. */
  CHECK(result.verification_code == ALL_FAULTS_GOOD);
  CHECK(result.checked_combination_count == 2U);
  comb_reference_context_destroy(context);
  fixture_destroy(&fixture);
}

/*
 * Lock the main CPU implementation's E_1 E_2 convention with a one-qubit
 * pair whose phase vector changes if the factors are reversed.  For
 * E_1=(F=1,h_x=0,h_z=0) and E_2=(F=1,h_x=1,h_z=0), E_1 E_2 has
 * (F=0,h_x=1,h_z=1), while E_2 E_1 has h_z=0.  The lookup rejects only the
 * main-order result.
 */
static void test_main_cpu_fault_composition_order(void) {
  problem_fixture_t fixture;
  comb_reference_context_t *context = NULL;
  comb_cuda_result_t result;
  uint32_t faults[2] = {UINT32_MAX, UINT32_MAX};
  comb_cuda_witness_t witness;

  fixture_init(&fixture, 1U, 2U);
  free(fixture.partitions);
  fixture.partitions = (uint32_t *)malloc(2U * sizeof(uint32_t));
  CHECK(fixture.partitions != NULL);
  fixture.partitions[0] = 1U;
  fixture.partitions[1] = 1U;
  fixture.problem.fault_location_count = 2U;
  fixture.problem.fault_location_partitions = fixture.partitions;
  fixture.problem.s = 2U;
  dcp_set_f(&fixture.outputs[0], 0U, 0U, 1U);
  dcp_set_f(&fixture.outputs[1], 0U, 0U, 1U);
  dcp_set_hx(&fixture.outputs[1], 0U, 1U);
  external_bit_set(fixture.z_lookup, 1U, 0U);

  CHECK_STATUS(comb_reference_context_create(&fixture.problem, &context),
               COMB_CUDA_STATUS_SUCCESS);
  result = fresh_result();
  witness = fresh_witness(faults, 2U);
  CHECK_STATUS(comb_reference_verify(context,
                                     COMB_CUDA_WITNESS_DETERMINISTIC, &result,
                                     &witness),
               COMB_CUDA_STATUS_SUCCESS);
  CHECK(result.verification_code == BAD_FAULT_FOUND);
  CHECK(result.total_combination_count == 1U);
  CHECK(result.checked_combination_count == 1U);
  CHECK(result.minimum_bad_rank == 0U);
  CHECK(witness.rank == 0U);
  CHECK(witness.fault_count == 2U);
  CHECK(faults[0] == 0U);
  CHECK(faults[1] == 1U);

  comb_reference_context_destroy(context);
  fixture_destroy(&fixture);
}

static void test_missing_measurement_is_rejected(void) {
  problem_fixture_t fixture;
  comb_reference_context_t *context = NULL;
  comb_cuda_result_t result;
  fixture_init(&fixture, 2U, 1U);
  fixture.problem.n = 1U;
  fixture_set_logical_count(&fixture, 1U);
  fixture_replace_lookup(&fixture, 1U, 2U);
  fixture_replace_lookup(&fixture, 0U, 2U);
  external_bit_set(fixture.x_lookup, 0U, 1U);
  external_bit_set(fixture.x_lookup, 1U, 1U);
  external_bit_set(fixture.z_lookup, 0U, 1U);
  external_bit_set(fixture.z_lookup, 1U, 1U);
  dcp_set_hz(&fixture.outputs[0], 1U, 1U);
  CHECK_STATUS(comb_reference_context_create(&fixture.problem, &context),
               COMB_CUDA_STATUS_SUCCESS);
  result = fresh_result();
  CHECK_STATUS(comb_reference_verify(context,
                                     COMB_CUDA_WITNESS_DETERMINISTIC, &result,
                                     NULL),
               COMB_CUDA_STATUS_INVALID_PROBLEM);
  comb_reference_context_destroy(context);
  fixture_destroy(&fixture);
}

static void test_seeded_scalar_crosscheck(void) {
  size_t case_index;
  for (case_index = 0U; case_index < 250U; ++case_index) {
    const uint32_t qubits = 1U + (uint32_t)(next_random() % 3U);
    const uint32_t output_count = 1U + (uint32_t)(next_random() % 8U);
    const uint64_t lookup_count = UINT64_C(1) << qubits;
    problem_fixture_t fixture;
    comb_reference_context_t *context = NULL;
    comb_cuda_result_t result;
    uint64_t expected_bad = COMB_CUDA_RANK_NONE;
    uint32_t output;
    uint64_t lookup;
    fixture_init(&fixture, qubits, output_count);
    for (lookup = 0U; lookup < lookup_count; ++lookup) {
      external_bit_set(fixture.x_lookup, lookup,
                       (unsigned)((next_random() >> 3U) & 1U));
      external_bit_set(fixture.z_lookup, lookup,
                       (unsigned)((next_random() >> 7U) & 1U));
    }
    for (output = 0U; output < output_count; ++output) {
      uint64_t hx = 0U;
      uint64_t hz = 0U;
      uint32_t bit;
      for (bit = 0U; bit < qubits; ++bit) {
        const unsigned x = (unsigned)((next_random() >> 11U) & 1U);
        const unsigned z = (unsigned)((next_random() >> 17U) & 1U);
        dcp_set_hx(&fixture.outputs[output], bit, x);
        dcp_set_hz(&fixture.outputs[output], bit, z);
        hx |= (uint64_t)x << bit;
        hz |= (uint64_t)z << bit;
      }
      if (expected_bad == COMB_CUDA_RANK_NONE &&
          (external_bit_get(fixture.x_lookup, hx) == 0U ||
           external_bit_get(fixture.z_lookup, hz) == 0U))
        expected_bad = output;
    }
    CHECK_STATUS(comb_reference_context_create(&fixture.problem, &context),
                 COMB_CUDA_STATUS_SUCCESS);
    result = fresh_result();
    CHECK_STATUS(comb_reference_verify(context,
                                       COMB_CUDA_WITNESS_DETERMINISTIC,
                                       &result, NULL),
                 COMB_CUDA_STATUS_SUCCESS);
    CHECK(result.checked_combination_count == output_count);
    CHECK(result.minimum_bad_rank == expected_bad);
    CHECK(result.verification_code ==
          (expected_bad == COMB_CUDA_RANK_NONE ? ALL_FAULTS_GOOD
                                               : BAD_FAULT_FOUND));
    comb_reference_context_destroy(context);
    fixture_destroy(&fixture);
  }
}

int main(void) {
  test_known_good_bad_and_ranges();
  test_z_lookup_and_logical_existential();
  test_context_validation();
  test_sequential_measurement_contract();
  test_main_cpu_fault_composition_order();
  test_missing_measurement_is_rejected();
  test_seeded_scalar_crosscheck();
  puts("comb_reference: sequential reference checks passed");
  return EXIT_SUCCESS;
}
