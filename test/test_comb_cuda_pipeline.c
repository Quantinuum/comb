#include "comb_cuda_internal.h"
#include "comb_rank.h"
#include "comb_reference.h"
#ifndef COMB_TEST_CUDA_DEVICE
#include "comb_cuda_pipeline.h"
#endif

#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define CHECK(x) do { if (!(x)) { \
  fprintf(stderr, "%s:%d: check failed: %s\n", __FILE__, __LINE__, #x); \
  exit(EXIT_FAILURE); } } while (0)

typedef struct {
  comb_problem_t problem;
  uint32_t partitions[3];
  comb_dcp_t *errors;
  comb_dcp_t *logicals;
  comb_vector_t *xs;
  comb_vector_t *zs;
  comb_vector_t *xm;
  comb_vector_t *zm;
  comb_bit_word_t *xl;
  comb_bit_word_t *zl;
} fixture_t;

static uint64_t random_state = UINT64_C(0x78d97204ee3db625);

static uint64_t random_u64(void) {
  uint64_t x = random_state;
  x ^= x >> 12U;
  x ^= x << 25U;
  x ^= x >> 27U;
  random_state = x;
  return x * UINT64_C(2685821657736338717);
}

static comb_bit_word_t *bits(uint64_t count) {
  comb_bit_word_t *p = calloc((size_t)((count + 7U) / 8U), sizeof(*p));
  CHECK(p != NULL);
  return p;
}

static void bit_set(comb_bit_word_t *p, uint64_t bit, unsigned value) {
  const comb_bit_word_t mask = (comb_bit_word_t)(1U << (bit % 8U));
  if (value != 0U) p[bit / 8U] |= mask;
  else p[bit / 8U] &= (comb_bit_word_t)~mask;
}

static comb_vector_t *vectors(uint32_t count, uint32_t width) {
  comb_vector_t *p = count != 0U ? calloc(count, sizeof(*p)) : NULL;
  CHECK(count == 0U || p != NULL);
  for (uint32_t i = 0; i < count; ++i) {
    p[i].n = width;
    p[i].data = bits(width);
  }
  return p;
}

static comb_dcp_t *dcps(uint32_t count, uint32_t width) {
  comb_dcp_t *p = count != 0U ? calloc(count, sizeof(*p)) : NULL;
  CHECK(count == 0U || p != NULL);
  for (uint32_t i = 0; i < count; ++i) {
    p[i].F.m = p[i].F.n = width;
    p[i].F.data = bits((uint64_t)width * width);
    p[i].h_x.n = p[i].h_z.n = width;
    p[i].h_x.data = bits(width);
    p[i].h_z.data = bits(width);
  }
  return p;
}

static void f_set(comb_dcp_t *p, uint32_t row, uint32_t column,
                  unsigned value) {
  bit_set(p->F.data, (uint64_t)row * p->F.n + column, value);
  bit_set(p->F.data, (uint64_t)column * p->F.n + row, value);
}

static void fixture_init(fixture_t *f, uint32_t all, uint32_t n,
                         uint32_t nx, uint32_t nz, uint32_t xm,
                         uint32_t zm, uint32_t s) {
  memset(f, 0, sizeof(*f));
  f->problem = (comb_problem_t)COMB_PROBLEM_INITIALIZER;
  f->problem.n_all = all;
  f->problem.n = n;
  f->problem.n_x = nx;
  f->problem.n_z = nz;
  f->problem.k = n - nx - nz;
  f->partitions[0] = 2U;
  f->partitions[1] = 1U;
  f->partitions[2] = 3U;
  f->problem.fault_location_count = 3U;
  f->problem.fault_location_partitions = f->partitions;
  f->problem.s = s;
  f->problem.output_error_count = 6U;
  f->problem.output_errors = f->errors = dcps(6U, all);
  f->problem.dcp_stabilisers = f->logicals = dcps(f->problem.k, all);
  f->problem.x_stabilisers = f->xs = vectors(nx, n);
  f->problem.z_stabilisers = f->zs = vectors(nz, n);
  f->problem.x_measurement_count = xm;
  f->problem.x_measurement_paulis = f->xm = vectors(xm, all);
  f->problem.z_measurement_count = zm;
  f->problem.z_measurement_paulis = f->zm = vectors(zm, all);
  f->problem.x_lookup_bit_count = UINT64_C(1) << (n - nx);
  f->problem.z_lookup_bit_count = UINT64_C(1) << (n - nz);
  f->problem.x_lookup_table = f->xl = bits(f->problem.x_lookup_bit_count);
  f->problem.z_lookup_table = f->zl = bits(f->problem.z_lookup_bit_count);
  for (uint64_t i = 0; i < f->problem.x_lookup_bit_count; ++i)
    bit_set(f->xl, i, 1U);
  for (uint64_t i = 0; i < f->problem.z_lookup_bit_count; ++i)
    bit_set(f->zl, i, 1U);
}

static void vectors_free(comb_vector_t *p, uint32_t count) {
  for (uint32_t i = 0; i < count; ++i) free(p[i].data);
  free(p);
}

static void dcps_free(comb_dcp_t *p, uint32_t count) {
  for (uint32_t i = 0; i < count; ++i) {
    free(p[i].F.data);
    free(p[i].h_x.data);
    free(p[i].h_z.data);
  }
  free(p);
}

static void fixture_free(fixture_t *f) {
  dcps_free(f->errors, f->problem.output_error_count);
  dcps_free(f->logicals, f->problem.k);
  vectors_free(f->xs, f->problem.n_x);
  vectors_free(f->zs, f->problem.n_z);
  vectors_free(f->xm, f->problem.x_measurement_count);
  vectors_free(f->zm, f->problem.z_measurement_count);
  free(f->xl);
  free(f->zl);
}

static comb_cuda_witness_t witness_for(uint32_t *indices) {
  comb_cuda_witness_t witness = COMB_CUDA_WITNESS_INITIALIZER;
  witness.fault_capacity = 3U;
  witness.fault_indices = indices;
  return witness;
}

static void check_case(const fixture_t *f, uint32_t case_index) {
  comb_reference_context_t *reference = NULL;
  comb_rank_plan_t *rank_plan = NULL;
  uint64_t total = 0;
  CHECK(comb_reference_context_create(&f->problem, &reference) ==
        COMB_CUDA_STATUS_SUCCESS);
  CHECK(comb_rank_plan_create(3U, f->partitions, f->problem.s, &rank_plan) ==
        COMB_RANK_OK);
  CHECK(comb_rank_plan_total(rank_plan, &total) == COMB_RANK_OK);
#ifdef COMB_TEST_CUDA_DEVICE
  comb_cuda_context_t *backend = NULL;
  comb_cuda_options_t options = COMB_CUDA_OPTIONS_INITIALIZER;
  options.epoch_size = 3U;
  options.scratch_limit_bytes = UINT64_C(16) * 1024U * 1024U;
  CHECK(comb_cuda_context_create(&f->problem, &options, &backend) ==
        COMB_CUDA_STATUS_SUCCESS);
#else
  void *blob = NULL;
  void *prep = NULL;
  uint64_t blob_size = 0;
  uint64_t prep_size = 0;
  uint64_t scratch_bytes = 0;
  CHECK(comb_cuda_flat_problem_pack(&f->problem, &blob, &blob_size) ==
        COMB_CUDA_STATUS_SUCCESS);
  CHECK(comb_cuda_prepare(blob, blob_size, &prep, &prep_size, &scratch_bytes) ==
        COMB_CUDA_STATUS_SUCCESS);
  void *scratch = calloc(1U, (size_t)scratch_bytes);
  CHECK(scratch != NULL);
#endif

  for (uint64_t rank = 0; rank < total; ++rank) {
    comb_cuda_result_t expected = COMB_CUDA_RESULT_INITIALIZER;
    uint32_t expected_faults[3] = {UINT32_MAX, UINT32_MAX, UINT32_MAX};
    comb_cuda_witness_t expected_witness = witness_for(expected_faults);
    CHECK(comb_reference_verify_range(reference, rank, rank + 1U,
              COMB_CUDA_WITNESS_DETERMINISTIC, &expected,
              &expected_witness) == COMB_CUDA_STATUS_SUCCESS);
#ifdef COMB_TEST_CUDA_DEVICE
    comb_cuda_result_t actual = COMB_CUDA_RESULT_INITIALIZER;
    uint32_t actual_faults[3] = {UINT32_MAX, UINT32_MAX, UINT32_MAX};
    comb_cuda_witness_t actual_witness = witness_for(actual_faults);
    CHECK(comb_cuda_verify_range(backend, rank, rank + 1U,
              COMB_CUDA_WITNESS_DETERMINISTIC, &actual,
              &actual_witness) == COMB_CUDA_STATUS_SUCCESS);
    CHECK(actual.verification_code == expected.verification_code);
    CHECK(actual.minimum_bad_rank == expected.minimum_bad_rank);
    CHECK(actual.checked_combination_count == 1U);
    CHECK(actual_witness.fault_count == expected_witness.fault_count);
    for (uint32_t i = 0; i < actual_witness.fault_count; ++i)
      CHECK(actual_faults[i] == expected_faults[i]);
#else
    uint32_t actual_locations[3] = {0U, 0U, 0U};
    uint32_t actual_faults[3] = {0U, 0U, 0U};
    uint32_t expected_locations[3] = {0U, 0U, 0U};
    uint32_t expected_offsets[3] = {0U, 0U, 0U};
    CHECK(comb_gpu_unrank(blob, rank, actual_locations, actual_faults));
    CHECK(comb_rank_plan_unrank(rank_plan, rank, expected_locations,
                                expected_offsets) == COMB_RANK_OK);
    for (uint32_t i = 0; i < f->problem.s; ++i) {
      uint32_t prefix = 0;
      CHECK(actual_locations[i] == expected_locations[i]);
      for (uint32_t j = 0; j < expected_locations[i]; ++j)
        prefix += f->partitions[j];
      CHECK(actual_faults[i] == prefix + expected_offsets[i]);
    }
    int bad = -1;
    const int status = comb_gpu_candidate(blob, prep, rank, scratch, &bad);
    if (status != COMB_CUDA_STATUS_SUCCESS ||
        bad != (expected.verification_code == BAD_FAULT_FOUND)) {
      fprintf(stderr, "case=%u n_all=%u rank=%" PRIu64
              " status=%d bad=%d expected=%d\n", case_index,
              f->problem.n_all, rank, status, bad,
              expected.verification_code == BAD_FAULT_FOUND);
      exit(EXIT_FAILURE);
    }
#endif
  }
#ifdef COMB_TEST_CUDA_DEVICE
  {
    comb_cuda_result_t expected = COMB_CUDA_RESULT_INITIALIZER;
    comb_cuda_result_t actual = COMB_CUDA_RESULT_INITIALIZER;
    CHECK(comb_reference_verify(reference, COMB_CUDA_WITNESS_DETERMINISTIC,
                                 &expected, NULL) == COMB_CUDA_STATUS_SUCCESS);
    CHECK(comb_cuda_verify(backend, COMB_CUDA_WITNESS_DETERMINISTIC,
                            &actual, NULL) == COMB_CUDA_STATUS_SUCCESS);
    CHECK(actual.verification_code == expected.verification_code);
    CHECK(actual.minimum_bad_rank == expected.minimum_bad_rank);
    CHECK(actual.checked_combination_count == total);
    actual = (comb_cuda_result_t)COMB_CUDA_RESULT_INITIALIZER;
    uint32_t fast_faults[3] = {UINT32_MAX, UINT32_MAX, UINT32_MAX};
    comb_cuda_witness_t fast_witness = witness_for(fast_faults);
    CHECK(comb_cuda_verify(backend, COMB_CUDA_WITNESS_FAST_ANY,
                            &actual, &fast_witness) == COMB_CUDA_STATUS_SUCCESS);
    CHECK(actual.verification_code == expected.verification_code);
    CHECK(actual.checked_combination_count > 0U);
    CHECK(actual.checked_combination_count <= total);
    if (actual.verification_code == BAD_FAULT_FOUND) {
      comb_cuda_result_t recheck = COMB_CUDA_RESULT_INITIALIZER;
      uint32_t verified_faults[3] = {UINT32_MAX, UINT32_MAX, UINT32_MAX};
      comb_cuda_witness_t verified_witness = witness_for(verified_faults);
      CHECK(actual.minimum_bad_rank < total);
      CHECK(fast_witness.rank == actual.minimum_bad_rank);
      CHECK(comb_reference_verify_range(reference, actual.minimum_bad_rank,
              actual.minimum_bad_rank + 1U, COMB_CUDA_WITNESS_DETERMINISTIC,
              &recheck, &verified_witness) == COMB_CUDA_STATUS_SUCCESS);
      CHECK(recheck.verification_code == BAD_FAULT_FOUND);
      CHECK(fast_witness.fault_count == verified_witness.fault_count);
      for (uint32_t i = 0; i < fast_witness.fault_count; ++i)
        CHECK(fast_faults[i] == verified_faults[i]);
    } else {
      CHECK(actual.checked_combination_count == total);
      CHECK(fast_witness.fault_count == 0U);
    }
    actual = (comb_cuda_result_t)COMB_CUDA_RESULT_INITIALIZER;
    CHECK(comb_cuda_verify_range(backend, total, total,
            COMB_CUDA_WITNESS_FAST_ANY, &actual, NULL) == COMB_CUDA_STATUS_SUCCESS);
    CHECK(actual.verification_code == ALL_FAULTS_GOOD);
    CHECK(actual.checked_combination_count == 0U);
  }
  comb_cuda_context_destroy(backend);
#else
  free(scratch);
  free(prep);
  comb_cuda_flat_problem_free(blob);
#endif
  comb_rank_plan_destroy(rank_plan);
  comb_reference_context_destroy(reference);
  (void)case_index;
}

/* These DCPs preserve the CSS code's Z-stabilised first qubits and the
 * X-stabilised penultimate qubit.  The code is deliberately small enough
 * for exhaustive rank and logical/coset enumeration. */
static void css_cases(void) {
  for (uint32_t c = 0; c < 80U; ++c) {
    fixture_t f;
    const uint32_t n = 3U + c % 2U;
    fixture_init(&f, n, n, 1U, n - 2U, 0U, 0U, 1U + c % 3U);
    bit_set(f.xs[0].data, n - 2U, 1U);
    for (uint32_t row = 0; row < n - 2U; ++row)
      bit_set(f.zs[row].data, row, 1U);
    f_set(&f.logicals[0], n - 1U, n - 1U, 1U);
    f_set(&f.logicals[0], 0U, n - 1U, 1U);
    for (uint32_t error = 0; error < 6U; ++error) {
      for (uint32_t row = 0; row < n; ++row) {
        bit_set(f.errors[error].h_x.data, row,
                (unsigned)((random_u64() >> 7U) & 1U));
        bit_set(f.errors[error].h_z.data, row,
                (unsigned)((random_u64() >> 9U) & 1U));
        for (uint32_t col = row; col < n; ++col)
          f_set(&f.errors[error], row, col,
                  (unsigned)((random_u64() >> 13U) & 1U));
      }
    }
    if (c != 0U) {
      for (uint64_t i = 0; i < f.problem.x_lookup_bit_count; ++i)
        bit_set(f.xl, i, (unsigned)((random_u64() >> 11U) & 1U));
      for (uint64_t i = 0; i < f.problem.z_lookup_bit_count; ++i)
        bit_set(f.zl, i, (unsigned)((random_u64() >> 17U) & 1U));
    }
    check_case(&f, c);
    fixture_free(&f);
  }
}

static void measurement_cases(void) {
  for (uint32_t c = 0; c < 32U; ++c) {
    fixture_t f;
    fixture_init(&f, 3U, 1U, 0U, 0U, 1U, 1U, 1U + c % 3U);
    bit_set(f.zm[0].data, 1U, 1U);
    bit_set(f.xm[0].data, 2U, 1U);
    f_set(&f.logicals[0], 0U, 0U, 1U);
    for (uint32_t error = 0; error < 6U; ++error) {
      for (uint32_t row = 0; row < 3U; ++row) {
        bit_set(f.errors[error].h_x.data, row,
                (unsigned)((random_u64() >> 3U) & 1U));
        bit_set(f.errors[error].h_z.data, row,
                (unsigned)((random_u64() >> 19U) & 1U));
        for (uint32_t col = row; col < 3U; ++col)
          f_set(&f.errors[error], row, col,
                  (unsigned)((random_u64() >> 29U) & 1U));
      }
    }
    bit_set(f.xl, 1U, c % 2U);
    bit_set(f.zl, 1U, (c / 2U) % 2U);
    check_case(&f, 80U + c);
    fixture_free(&f);
  }
}

static void multiword_cases(void) {
  for (uint32_t c = 0; c < 4U; ++c) {
    fixture_t f;
    fixture_init(&f, 65U, 1U, 0U, 0U, 64U, 0U, 1U + c % 3U);
    for (uint32_t row = 0; row < 64U; ++row)
      bit_set(f.xm[row].data, row + 1U, 1U);
    f_set(&f.logicals[0], 0U, 0U, 1U);
    for (uint32_t error = 0; error < 6U; ++error) {
      bit_set(f.errors[error].h_x.data, 64U, error % 2U);
      bit_set(f.errors[error].h_z.data, 64U, (error / 2U) % 2U);
      bit_set(f.errors[error].h_z.data, 0U, error % 2U);
      f_set(&f.errors[error], 0U, 64U, (error / 2U) % 2U);
      f_set(&f.errors[error], 64U, 64U, (error / 3U) % 2U);
    }
    bit_set(f.zl, 1U, c % 2U);
    check_case(&f, 112U + c);
    fixture_free(&f);
  }
}

#ifdef COMB_TEST_CUDA_DEVICE
/* Exercise multiple kernel blocks, repeated bounded epochs, and early stop.
 * The small exhaustive fixtures above intentionally use tiny rank ranges. */
static void scheduler_cases(void) {
  fixture_t f;
  comb_cuda_context_t *context = NULL;
  comb_cuda_options_t options = COMB_CUDA_OPTIONS_INITIALIZER;
  comb_cuda_result_t result = COMB_CUDA_RESULT_INITIALIZER;
  uint32_t faults[3] = {0U, 0U, 0U};
  comb_cuda_witness_t witness = witness_for(faults);
  fixture_init(&f, 1U, 1U, 0U, 0U, 0U, 0U, 1U);
  dcps_free(f.errors, 6U);
  f.problem.output_error_count = 4096U;
  f.problem.output_errors = f.errors = dcps(4096U, 1U);
  f.partitions[0] = 1500U;
  f.partitions[1] = 1500U;
  f.partitions[2] = 1096U;
  options.epoch_size = 129U;
  options.scratch_limit_bytes = 1024U * 1024U;
  CHECK(comb_cuda_context_create(&f.problem, &options, &context) ==
        COMB_CUDA_STATUS_SUCCESS);
  CHECK(comb_cuda_verify(context, COMB_CUDA_WITNESS_DETERMINISTIC,
                          &result, NULL) == COMB_CUDA_STATUS_SUCCESS);
  CHECK(result.verification_code == ALL_FAULTS_GOOD);
  CHECK(result.checked_combination_count == 4096U);
  comb_cuda_context_destroy(context);
  bit_set(f.xl, 1U, 0U);
  bit_set(f.errors[3000].h_x.data, 0U, 1U);
  CHECK(comb_cuda_context_create(&f.problem, &options, &context) ==
        COMB_CUDA_STATUS_SUCCESS);
  result = (comb_cuda_result_t)COMB_CUDA_RESULT_INITIALIZER;
  CHECK(comb_cuda_verify(context, COMB_CUDA_WITNESS_FAST_ANY,
                          &result, &witness) == COMB_CUDA_STATUS_SUCCESS);
  CHECK(result.verification_code == BAD_FAULT_FOUND);
  CHECK(result.minimum_bad_rank == 3000U);
  CHECK(result.checked_combination_count >= 3001U);
  CHECK(result.checked_combination_count <= 4096U);
  CHECK(witness.rank == 3000U);
  CHECK(witness.fault_count == 1U);
  CHECK(faults[0] == 3000U);
  comb_cuda_context_destroy(context);
  context = NULL;
  options.scratch_limit_bytes = 1U;
  CHECK(comb_cuda_context_create(&f.problem, &options, &context) ==
        COMB_CUDA_STATUS_OUT_OF_MEMORY);
  CHECK(context == NULL);
  fixture_free(&f);
}
#endif

int main(void) {
#ifdef COMB_TEST_CUDA_DEVICE
  /* A CUDA build on a login node can run host checks but must explicitly
   * report that its genuine device test was skipped. */
  fixture_t probe;
  comb_cuda_context_t *context = NULL;
  fixture_init(&probe, 1U, 1U, 0U, 0U, 0U, 0U, 1U);
  comb_cuda_status_t status = comb_cuda_context_create(&probe.problem, NULL,
                                                       &context);
  fixture_free(&probe);
  if (status == COMB_CUDA_STATUS_CUDA_UNAVAILABLE) {
    puts("SKIP: no usable CUDA device");
    return 77;
  }
  CHECK(status == COMB_CUDA_STATUS_SUCCESS);
  comb_cuda_context_destroy(context);
#endif
  css_cases();
  measurement_cases();
  multiword_cases();
#ifdef COMB_TEST_CUDA_DEVICE
  scheduler_cases();
  puts("comb CUDA device: 116 exhaustive differential fixtures passed");
#else
  puts("comb CUDA scientific functions on host: 116 differential fixtures passed");
#endif
  return EXIT_SUCCESS;
}
