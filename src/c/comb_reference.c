#define _POSIX_C_SOURCE 200809L

#include "comb_reference.h"

#include "comb_cuda_internal.h"
#include "comb_dense.h"
#include "comb_rank.h"

#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

struct comb_reference_context {
  void *blob;
  uint64_t blob_size;
  const comb_cuda_flat_problem_t *flat;
  comb_rank_plan_t *rank_plan;
  uint32_t *x_pivots;
  uint32_t *z_pivots;
  uint32_t *x_free_columns;
  uint32_t *z_free_columns;
  uint32_t x_free_count;
  uint32_t z_free_count;
  comb_dense_matrix_t z_sequence;
  comb_dense_matrix_t x_sequence;
  uint64_t logical_variant_count;
  double setup_milliseconds;
};

static double reference_now_milliseconds(void) {
  struct timespec value;
  if (clock_gettime(CLOCK_MONOTONIC, &value) != 0) return 0.0;
  return (double)value.tv_sec * 1000.0 + (double)value.tv_nsec / 1000000.0;
}

static const unsigned char *reference_bytes(
    const comb_reference_context_t *context) {
  return (const unsigned char *)context->blob;
}

static const uint64_t *reference_u64_section(
    const comb_reference_context_t *context, uint64_t offset) {
  return (const uint64_t *)(reference_bytes(context) + offset);
}

static const uint32_t *reference_u32_section(
    const comb_reference_context_t *context, uint64_t offset) {
  return (const uint32_t *)(reference_bytes(context) + offset);
}

static unsigned reference_word_bit(const uint64_t *words, size_t bit) {
  return (unsigned)((words[bit / 64U] >> (bit % 64U)) & UINT64_C(1));
}

static void reference_word_set(uint64_t *words, size_t bit, unsigned value) {
  const uint64_t mask = UINT64_C(1) << (bit % 64U);
  if (value != 0U)
    words[bit / 64U] |= mask;
  else
    words[bit / 64U] &= ~mask;
}

static unsigned reference_parity64(uint64_t value) {
  value ^= value >> 32U;
  value ^= value >> 16U;
  value ^= value >> 8U;
  value ^= value >> 4U;
  value &= UINT64_C(0xf);
  return (unsigned)((UINT64_C(0x6996) >> value) & UINT64_C(1));
}

static unsigned reference_dot(const uint64_t *left, const uint64_t *right,
                              size_t word_count) {
  size_t word;
  unsigned parity = 0U;
  for (word = 0U; word < word_count; ++word)
    parity ^= reference_parity64(left[word] & right[word]);
  return parity & 1U;
}

static unsigned reference_lookup_bit(const uint64_t *lookup,
                                     uint64_t index) {
  return (unsigned)((lookup[index / 64U] >> (index % 64U)) & UINT64_C(1));
}

static comb_cuda_status_t reference_dense_status(comb_dense_status_t status) {
  switch (status) {
    case COMB_DENSE_OK:
      return COMB_CUDA_STATUS_SUCCESS;
    case COMB_DENSE_ALLOCATION_FAILED:
      return COMB_CUDA_STATUS_OUT_OF_MEMORY;
    case COMB_DENSE_OVERFLOW:
      return COMB_CUDA_STATUS_PROBLEM_TOO_LARGE;
    case COMB_DENSE_INVALID_ARGUMENT:
    case COMB_DENSE_DIMENSION_MISMATCH:
    case COMB_DENSE_OUT_OF_RANGE:
    case COMB_DENSE_NONZERO_PADDING:
    case COMB_DENSE_NOT_GENERATED:
    case COMB_DENSE_INCONSISTENT_TABLEAU:
      return COMB_CUDA_STATUS_INVALID_PROBLEM;
    default:
      return COMB_CUDA_STATUS_INTERNAL_ERROR;
  }
}

static comb_cuda_status_t reference_rank_status(comb_rank_status_t status) {
  switch (status) {
    case COMB_RANK_OK:
      return COMB_CUDA_STATUS_SUCCESS;
    case COMB_RANK_ALLOCATION_FAILED:
      return COMB_CUDA_STATUS_OUT_OF_MEMORY;
    case COMB_RANK_OVERFLOW:
      return COMB_CUDA_STATUS_PROBLEM_TOO_LARGE;
    case COMB_RANK_INVALID_ARGUMENT:
    case COMB_RANK_OUT_OF_RANGE:
      return COMB_CUDA_STATUS_INVALID_PROBLEM;
    default:
      return COMB_CUDA_STATUS_INTERNAL_ERROR;
  }
}

static comb_cuda_status_t reference_allocate_u32(uint32_t count,
                                                 uint32_t **values) {
  if (values == NULL) return COMB_CUDA_STATUS_INTERNAL_ERROR;
  *values = NULL;
  if (count == 0U) return COMB_CUDA_STATUS_SUCCESS;
  if ((size_t)count > SIZE_MAX / sizeof(uint32_t))
    return COMB_CUDA_STATUS_PROBLEM_TOO_LARGE;
  *values = (uint32_t *)malloc((size_t)count * sizeof(uint32_t));
  return *values == NULL ? COMB_CUDA_STATUS_OUT_OF_MEMORY
                         : COMB_CUDA_STATUS_SUCCESS;
}

static comb_cuda_status_t reference_extract_rref(
    const uint64_t *rows, uint32_t row_count, uint32_t width,
    uint32_t stride_words, uint32_t **pivots, uint32_t **free_columns,
    uint32_t *free_count) {
  uint32_t row;
  uint32_t column;
  uint32_t previous = 0U;
  unsigned have_previous = 0U;
  comb_cuda_status_t status;
  if (pivots == NULL || free_columns == NULL || free_count == NULL ||
      (row_count != 0U && rows == NULL))
    return COMB_CUDA_STATUS_INTERNAL_ERROR;
  *pivots = NULL;
  *free_columns = NULL;
  *free_count = width - row_count;
  status = reference_allocate_u32(row_count, pivots);
  if (status != COMB_CUDA_STATUS_SUCCESS) return status;
  status = reference_allocate_u32(*free_count, free_columns);
  if (status != COMB_CUDA_STATUS_SUCCESS) {
    free(*pivots);
    *pivots = NULL;
    return status;
  }
  for (row = 0U; row < row_count; ++row) {
    uint32_t leading = width;
    const uint64_t *current = rows + (uint64_t)row * stride_words;
    for (column = 0U; column < width; ++column) {
      if (reference_word_bit(current, column) != 0U) {
        leading = column;
        break;
      }
    }
    if (leading == width ||
        (have_previous != 0U && leading <= previous)) {
      free(*free_columns);
      free(*pivots);
      *free_columns = NULL;
      *pivots = NULL;
      return COMB_CUDA_STATUS_INVALID_PROBLEM;
    }
    for (uint32_t other = 0U; other < row_count; ++other) {
      if (other != row &&
          reference_word_bit(rows + (uint64_t)other * stride_words,
                             leading) != 0U) {
        free(*free_columns);
        free(*pivots);
        *free_columns = NULL;
        *pivots = NULL;
        return COMB_CUDA_STATUS_INVALID_PROBLEM;
      }
    }
    (*pivots)[row] = leading;
    previous = leading;
    have_previous = 1U;
  }
  {
    uint32_t output = 0U;
    uint32_t pivot_index = 0U;
    for (column = 0U; column < width; ++column) {
      if (pivot_index < row_count && (*pivots)[pivot_index] == column)
        ++pivot_index;
      else
        (*free_columns)[output++] = column;
    }
    if (output != *free_count) return COMB_CUDA_STATUS_INTERNAL_ERROR;
  }
  return COMB_CUDA_STATUS_SUCCESS;
}

static comb_cuda_status_t reference_wrap_flat_dcp(
    const comb_reference_context_t *context, unsigned logical,
    uint32_t index, comb_dense_dcp_t *dcp) {
  const comb_cuda_flat_problem_t *flat = context->flat;
  uint64_t f_offset;
  uint64_t hx_offset;
  uint64_t hz_offset;
  uint64_t f_index;
  uint64_t vector_index;
  if (dcp == NULL ||
      (logical != 0U ? index >= flat->k : index >= flat->output_error_count))
    return COMB_CUDA_STATUS_INTERNAL_ERROR;
  f_offset = logical != 0U ? flat->logical_dcp_f_offset
                           : flat->output_dcp_f_offset;
  hx_offset = logical != 0U ? flat->logical_dcp_hx_offset
                            : flat->output_dcp_hx_offset;
  hz_offset = logical != 0U ? flat->logical_dcp_hz_offset
                            : flat->output_dcp_hz_offset;
  f_index = (uint64_t)index * flat->n_all * flat->words_all;
  vector_index = (uint64_t)index * flat->words_all;
  return reference_dense_status(comb_dense_dcp_wrap(
      dcp, flat->n_all, flat->words_all,
      (uint64_t *)(uintptr_t)(reference_u64_section(context, f_offset) +
                              f_index),
      (uint64_t *)(uintptr_t)(reference_u64_section(context, hx_offset) +
                              vector_index),
      (uint64_t *)(uintptr_t)(reference_u64_section(context, hz_offset) +
                              vector_index)));
}

static comb_cuda_status_t reference_validate_dcp_symmetry(
    const comb_reference_context_t *context, unsigned logical,
    uint32_t count) {
  uint32_t index;
  for (index = 0U; index < count; ++index) {
    comb_dense_dcp_t dcp;
    comb_dense_status_t dense_status;
    comb_cuda_status_t status =
        reference_wrap_flat_dcp(context, logical, index, &dcp);
    if (status != COMB_CUDA_STATUS_SUCCESS) return status;
    dense_status = comb_dense_dcp_validate(&dcp, 1U);
    comb_dense_dcp_destroy(&dcp);
    if (dense_status != COMB_DENSE_OK)
      return reference_dense_status(dense_status);
  }
  return COMB_CUDA_STATUS_SUCCESS;
}

static comb_cuda_status_t reference_validate_css_commutation(
    const comb_reference_context_t *context) {
  const comb_cuda_flat_problem_t *flat = context->flat;
  const uint64_t *x_rows =
      reference_u64_section(context, flat->x_stabilisers_offset);
  const uint64_t *z_rows =
      reference_u64_section(context, flat->z_stabilisers_offset);
  uint32_t x;
  uint32_t z;
  for (x = 0U; x < flat->n_x; ++x) {
    for (z = 0U; z < flat->n_z; ++z) {
      if (reference_dot(x_rows + (uint64_t)x * flat->words_code,
                        z_rows + (uint64_t)z * flat->words_code,
                        flat->words_code) != 0U)
        return COMB_CUDA_STATUS_INVALID_PROBLEM;
    }
  }
  return COMB_CUDA_STATUS_SUCCESS;
}

static comb_cuda_status_t reference_build_measurement_sequences(
    comb_reference_context_t *context) {
  const comb_cuda_flat_problem_t *flat = context->flat;
  const size_t z_rows =
      (size_t)flat->z_measurement_count + (size_t)flat->n_z;
  size_t x_rows;
  const uint64_t *z_measurements =
      reference_u64_section(context, flat->z_measurements_offset);
  const uint64_t *x_measurements =
      reference_u64_section(context, flat->x_measurements_offset);
  const uint64_t *z_stabilisers =
      reference_u64_section(context, flat->z_stabilisers_offset);
  comb_dense_status_t dense_status;
  uint32_t row;
  uint32_t word;
  if (z_rows > SIZE_MAX - (size_t)flat->x_measurement_count)
    return COMB_CUDA_STATUS_PROBLEM_TOO_LARGE;
  x_rows = z_rows + (size_t)flat->x_measurement_count;
  dense_status =
      comb_dense_matrix_init(&context->z_sequence, z_rows, flat->n_all);
  if (dense_status != COMB_DENSE_OK) return reference_dense_status(dense_status);
  dense_status =
      comb_dense_matrix_init(&context->x_sequence, x_rows, flat->n_all);
  if (dense_status != COMB_DENSE_OK) return reference_dense_status(dense_status);

  for (row = 0U; row < flat->z_measurement_count; ++row) {
    const uint64_t *source =
        z_measurements + (uint64_t)row * flat->words_all;
    uint64_t *z_destination =
        context->z_sequence.words + (size_t)row *
                                        context->z_sequence.stride_words;
    uint64_t *x_destination =
        context->x_sequence.words + (size_t)row *
                                        context->x_sequence.stride_words;
    if (reference_dot(source, source, flat->words_all) == 0U)
      return COMB_CUDA_STATUS_INVALID_PROBLEM;
    for (word = 0U; word < flat->words_all; ++word) {
      z_destination[word] = source[word];
      x_destination[word] = source[word];
    }
  }
  for (row = 0U; row < flat->n_z; ++row) {
    uint64_t *z_destination =
        context->z_sequence.words +
        ((size_t)flat->z_measurement_count + row) *
            context->z_sequence.stride_words;
    for (word = 0U; word < flat->words_code; ++word)
      z_destination[word] =
          z_stabilisers[(uint64_t)row * flat->words_code + word];
    reference_word_set(
        context->x_sequence.words +
            ((size_t)flat->z_measurement_count + row) *
                context->x_sequence.stride_words,
        context->z_pivots[row], 1U);
  }
  for (row = 0U; row < flat->x_measurement_count; ++row) {
    uint64_t *destination =
        context->x_sequence.words +
        (z_rows + row) * context->x_sequence.stride_words;
    const uint64_t *source =
        x_measurements + (uint64_t)row * flat->words_all;
    for (word = 0U; word < flat->words_all; ++word)
      destination[word] = source[word];
  }
  dense_status = comb_dense_matrix_validate(&context->z_sequence);
  if (dense_status != COMB_DENSE_OK) return reference_dense_status(dense_status);
  dense_status = comb_dense_matrix_validate(&context->x_sequence);
  return reference_dense_status(dense_status);
}

static comb_cuda_status_t reference_context_finish_create(
    comb_reference_context_t *context) {
  const comb_cuda_flat_problem_t *flat = context->flat;
  const uint32_t *partitions =
      reference_u32_section(context, flat->fault_partitions_offset);
  const uint64_t *x_rows =
      reference_u64_section(context, flat->x_stabilisers_offset);
  const uint64_t *z_rows =
      reference_u64_section(context, flat->z_stabilisers_offset);
  uint64_t total = 0U;
  comb_rank_status_t rank_status;
  comb_cuda_status_t status;
  if (flat->k >= 64U) return COMB_CUDA_STATUS_PROBLEM_TOO_LARGE;
  context->logical_variant_count = UINT64_C(1) << flat->k;
  rank_status = comb_rank_plan_create(flat->fault_location_count, partitions,
                                      flat->s, &context->rank_plan);
  if (rank_status != COMB_RANK_OK) return reference_rank_status(rank_status);
  rank_status = comb_rank_plan_total(context->rank_plan, &total);
  if (rank_status != COMB_RANK_OK) return reference_rank_status(rank_status);
  if (total != flat->total_combination_count)
    return COMB_CUDA_STATUS_INVALID_PROBLEM;
  status = reference_extract_rref(
      x_rows, flat->n_x, flat->n, flat->words_code, &context->x_pivots,
      &context->x_free_columns, &context->x_free_count);
  if (status != COMB_CUDA_STATUS_SUCCESS) return status;
  status = reference_extract_rref(
      z_rows, flat->n_z, flat->n, flat->words_code, &context->z_pivots,
      &context->z_free_columns, &context->z_free_count);
  if (status != COMB_CUDA_STATUS_SUCCESS) return status;
  if (context->x_free_count != flat->n - flat->n_x ||
      context->z_free_count != flat->n - flat->n_z ||
      flat->x_lookup_bit_count !=
          (UINT64_C(1) << context->x_free_count) ||
      flat->z_lookup_bit_count !=
          (UINT64_C(1) << context->z_free_count))
    return COMB_CUDA_STATUS_INVALID_PROBLEM;
  status = reference_validate_css_commutation(context);
  if (status != COMB_CUDA_STATUS_SUCCESS) return status;
  status = reference_validate_dcp_symmetry(context, 1U, flat->k);
  if (status != COMB_CUDA_STATUS_SUCCESS) return status;
  status = reference_validate_dcp_symmetry(context, 0U,
                                           flat->output_error_count);
  if (status != COMB_CUDA_STATUS_SUCCESS) return status;
  return reference_build_measurement_sequences(context);
}

comb_cuda_status_t comb_reference_context_create(
    const comb_problem_t *problem, comb_reference_context_t **context) {
  comb_reference_context_t *created;
  comb_cuda_status_t status;
  double started;
  if (context == NULL) return COMB_CUDA_STATUS_INVALID_ARGUMENT;
  *context = NULL;
  started = reference_now_milliseconds();
  created = (comb_reference_context_t *)calloc(1, sizeof(*created));
  if (created == NULL) return COMB_CUDA_STATUS_OUT_OF_MEMORY;
  status = comb_cuda_flat_problem_pack(problem, &created->blob,
                                       &created->blob_size);
  if (status != COMB_CUDA_STATUS_SUCCESS) {
    comb_reference_context_destroy(created);
    return status;
  }
  status = comb_cuda_flat_problem_validate_science(created->blob,
                                                   created->blob_size);
  if (status != COMB_CUDA_STATUS_SUCCESS) {
    comb_reference_context_destroy(created);
    return status;
  }
  created->flat = (const comb_cuda_flat_problem_t *)created->blob;
  status = reference_context_finish_create(created);
  if (status != COMB_CUDA_STATUS_SUCCESS) {
    comb_reference_context_destroy(created);
    return status;
  }
  created->setup_milliseconds = reference_now_milliseconds() - started;
  *context = created;
  return COMB_CUDA_STATUS_SUCCESS;
}

void comb_reference_context_destroy(comb_reference_context_t *context) {
  if (context == NULL) return;
  comb_dense_matrix_destroy(&context->z_sequence);
  comb_dense_matrix_destroy(&context->x_sequence);
  comb_rank_plan_destroy(context->rank_plan);
  free(context->x_pivots);
  free(context->z_pivots);
  free(context->x_free_columns);
  free(context->z_free_columns);
  comb_cuda_flat_problem_free(context->blob);
  free(context);
}

static comb_cuda_status_t reference_dcp_copy(comb_dense_dcp_t *destination,
                                             const comb_dense_dcp_t *source) {
  size_t row;
  comb_dense_status_t dense_status;
  dense_status = comb_dense_dcp_validate(destination, 0U);
  if (dense_status != COMB_DENSE_OK) return reference_dense_status(dense_status);
  dense_status = comb_dense_dcp_validate(source, 0U);
  if (dense_status != COMB_DENSE_OK) return reference_dense_status(dense_status);
  if (destination->qubits != source->qubits)
    return COMB_CUDA_STATUS_INTERNAL_ERROR;
  for (row = 0U; row < destination->qubits; ++row) {
    memcpy(destination->f.words + row * destination->f.stride_words,
           source->f.words + row * source->f.stride_words,
           destination->word_count * sizeof(uint64_t));
  }
  if (destination->word_count != 0U) {
    memcpy(destination->h_x, source->h_x,
           destination->word_count * sizeof(uint64_t));
    memcpy(destination->h_z, source->h_z,
           destination->word_count * sizeof(uint64_t));
  }
  return COMB_CUDA_STATUS_SUCCESS;
}

static comb_cuda_status_t reference_tableau_from_dcp(
    const comb_dense_dcp_t *dcp, comb_dense_tableau_t *tableau) {
  size_t row;
  comb_dense_status_t dense_status = comb_dense_tableau_init(
      tableau, dcp->qubits, dcp->qubits + 1U, dcp->qubits);
  if (dense_status != COMB_DENSE_OK) return reference_dense_status(dense_status);
  for (row = 0U; row < dcp->qubits; ++row) {
    memcpy(tableau->z.words + row * tableau->z.stride_words,
           dcp->f.words + row * dcp->f.stride_words,
           dcp->word_count * sizeof(uint64_t));
    reference_word_set(tableau->x.words + row * tableau->x.stride_words, row,
                       1U);
    tableau->delta[row] = (uint8_t)reference_word_bit(
        dcp->f.words + row * dcp->f.stride_words, row);
    tableau->epsilon[row] =
        (uint8_t)reference_word_bit(dcp->h_z, row);
  }
  dense_status = comb_dense_tableau_validate(tableau);
  return reference_dense_status(dense_status);
}

static uint64_t reference_compress_columns(const uint64_t *row,
                                           const uint32_t *columns,
                                           uint32_t count) {
  uint64_t result = 0U;
  uint32_t bit;
  for (bit = 0U; bit < count; ++bit) {
    if (reference_word_bit(row, columns[bit]) != 0U)
      result |= UINT64_C(1) << bit;
  }
  return result;
}

static int reference_row_uses_only_columns(const uint64_t *row,
                                           uint32_t width,
                                           const uint32_t *columns,
                                           uint32_t count) {
  uint32_t bit;
  for (bit = 0U; bit < width; ++bit) {
    uint32_t column;
    unsigned allowed = 0U;
    if (reference_word_bit(row, bit) == 0U) continue;
    for (column = 0U; column < count; ++column) {
      if (columns[column] == bit) {
        allowed = 1U;
        break;
      }
    }
    if (allowed == 0U) return 0;
  }
  return 1;
}

static comb_cuda_status_t reference_variant_good(
    const comb_reference_context_t *context, const comb_dense_dcp_t *error,
    unsigned *good) {
  const comb_cuda_flat_problem_t *flat = context->flat;
  const uint64_t *x_stabilisers =
      reference_u64_section(context, flat->x_stabilisers_offset);
  const uint64_t *x_lookup =
      reference_u64_section(context, flat->x_lookup_offset);
  const uint64_t *z_lookup =
      reference_u64_section(context, flat->z_lookup_offset);
  uint64_t *x_error = NULL;
  comb_dense_tableau_t tableau;
  uint64_t *generators = NULL;
  size_t z_rank = 0U;
  comb_cuda_status_t status = COMB_CUDA_STATUS_SUCCESS;
  comb_dense_status_t dense_status;
  uint32_t row;
  memset(&tableau, 0, sizeof(tableau));
  if (good == NULL) return COMB_CUDA_STATUS_INTERNAL_ERROR;
  *good = 0U;
  x_error = (uint64_t *)malloc((size_t)flat->words_all * sizeof(uint64_t));
  if (x_error == NULL) return COMB_CUDA_STATUS_OUT_OF_MEMORY;
  memcpy(x_error, error->h_x,
         (size_t)flat->words_all * sizeof(uint64_t));

  for (row = 0U; row < flat->z_measurement_count; ++row) {
    const uint64_t *measurement =
        context->z_sequence.words +
        (size_t)row * context->z_sequence.stride_words;
    if (reference_dot(x_error, measurement, flat->words_all) != 0U) {
      *good = 1U;
      goto cleanup;
    }
  }
  for (row = 0U; row < flat->n_x; ++row) {
    if (reference_word_bit(x_error, context->x_pivots[row]) != 0U) {
      uint32_t word;
      for (word = 0U; word < flat->words_code; ++word)
        x_error[word] ^=
            x_stabilisers[(uint64_t)row * flat->words_code + word];
    }
  }
  {
    const uint64_t reduced = reference_compress_columns(
        x_error, context->x_free_columns, context->x_free_count);
    if (reduced >= flat->x_lookup_bit_count) {
      status = COMB_CUDA_STATUS_INTERNAL_ERROR;
      goto cleanup;
    }
    if (reference_lookup_bit(x_lookup, reduced) == 0U) goto cleanup;
  }

  status = reference_tableau_from_dcp(error, &tableau);
  if (status != COMB_CUDA_STATUS_SUCCESS) goto cleanup;
  {
    unsigned rejected = 0U;
    dense_status = comb_dense_tableau_postselect(
        &tableau, &context->z_sequence, COMB_DENSE_MEASURE_Z, &rejected);
    if (dense_status != COMB_DENSE_OK) {
      status = reference_dense_status(dense_status);
      goto cleanup;
    }
    if (rejected != 0U) {
      *good = 1U;
      goto cleanup;
    }
    dense_status = comb_dense_tableau_postselect(
        &tableau, &context->x_sequence, COMB_DENSE_MEASURE_X, &rejected);
    if (dense_status != COMB_DENSE_OK) {
      status = reference_dense_status(dense_status);
      goto cleanup;
    }
    if (rejected != 0U) {
      *good = 1U;
      goto cleanup;
    }
  }
  dense_status = comb_dense_tableau_reduce_coset(&tableau, &z_rank);
  if (dense_status != COMB_DENSE_OK) {
    status = reference_dense_status(dense_status);
    goto cleanup;
  }
  if (z_rank >= 64U) {
    status = COMB_CUDA_STATUS_PROBLEM_TOO_LARGE;
    goto cleanup;
  }
  if (z_rank != 0U) {
    if (z_rank > SIZE_MAX / sizeof(uint64_t)) {
      status = COMB_CUDA_STATUS_PROBLEM_TOO_LARGE;
      goto cleanup;
    }
    generators = (uint64_t *)malloc(z_rank * sizeof(uint64_t));
    if (generators == NULL) {
      status = COMB_CUDA_STATUS_OUT_OF_MEMORY;
      goto cleanup;
    }
    for (size_t generator = 0U; generator < z_rank; ++generator) {
      const uint64_t *generator_row =
          tableau.z.words + generator * tableau.z.stride_words;
      if (!reference_row_uses_only_columns(
              generator_row, flat->n_all, context->z_free_columns,
              context->z_free_count)) {
        status = COMB_CUDA_STATUS_INVALID_PROBLEM;
        goto cleanup;
      }
      generators[generator] = reference_compress_columns(
          generator_row, context->z_free_columns, context->z_free_count);
    }
  }
  {
    uint64_t factor = 0U;
    size_t state_row;
    uint64_t coset_element;
    const uint64_t coset_count = UINT64_C(1) << z_rank;
    for (state_row = z_rank; state_row < flat->n_all; ++state_row) {
      if (tableau.epsilon[state_row] != 0U) {
        size_t leading = flat->n_all;
        size_t column;
        const uint64_t *x_row =
            tableau.x.words + state_row * tableau.x.stride_words;
        for (column = 0U; column < flat->n_all; ++column) {
          if (reference_word_bit(x_row, column) != 0U) {
            leading = column;
            break;
          }
        }
        if (leading == flat->n_all) {
          status = COMB_CUDA_STATUS_INVALID_PROBLEM;
          goto cleanup;
        }
        {
          unsigned found = 0U;
          for (uint32_t compressed = 0U;
               compressed < context->z_free_count; ++compressed) {
            if (context->z_free_columns[compressed] == leading) {
              factor |= UINT64_C(1) << compressed;
              found = 1U;
              break;
            }
          }
          if (found == 0U) {
            status = COMB_CUDA_STATUS_INVALID_PROBLEM;
            goto cleanup;
          }
        }
      }
    }
    for (coset_element = 0U; coset_element < coset_count; ++coset_element) {
      uint64_t reduced = factor;
      size_t generator;
      for (generator = 0U; generator < z_rank; ++generator) {
        if (((coset_element >> generator) & UINT64_C(1)) != 0U)
          reduced ^= generators[generator];
      }
      if (reduced >= flat->z_lookup_bit_count) {
        status = COMB_CUDA_STATUS_INTERNAL_ERROR;
        goto cleanup;
      }
      if (reference_lookup_bit(z_lookup, reduced) == 0U) goto cleanup;
    }
  }
  *good = 1U;

cleanup:
  free(generators);
  comb_dense_tableau_destroy(&tableau);
  free(x_error);
  return status;
}

static comb_cuda_status_t reference_rank_fault_indices(
    const comb_reference_context_t *context, uint64_t rank,
    uint32_t *locations, uint32_t *offsets, uint32_t *fault_indices) {
  const uint64_t *prefixes =
      reference_u64_section(context, context->flat->fault_prefixes_offset);
  comb_rank_status_t rank_status = comb_rank_plan_unrank(
      context->rank_plan, rank, locations, offsets);
  uint32_t index;
  if (rank_status != COMB_RANK_OK) return reference_rank_status(rank_status);
  for (index = 0U; index < context->flat->s; ++index) {
    const uint64_t fault = prefixes[locations[index]] + offsets[index];
    if (fault >= context->flat->output_error_count)
      return COMB_CUDA_STATUS_INTERNAL_ERROR;
    fault_indices[index] = (uint32_t)fault;
  }
  return COMB_CUDA_STATUS_SUCCESS;
}

static comb_cuda_status_t reference_tuple_bad(
    const comb_reference_context_t *context, uint64_t rank,
    uint32_t *locations, uint32_t *offsets, uint32_t *fault_indices,
    unsigned *bad) {
  comb_dense_dcp_t combined;
  comb_dense_dcp_t error;
  comb_dense_dcp_t view;
  comb_cuda_status_t status;
  comb_dense_status_t dense_status;
  uint32_t selected;
  uint64_t logical_mask;
  memset(&combined, 0, sizeof(combined));
  memset(&error, 0, sizeof(error));
  memset(&view, 0, sizeof(view));
  if (bad == NULL) return COMB_CUDA_STATUS_INTERNAL_ERROR;
  *bad = 0U;
  status = reference_rank_fault_indices(context, rank, locations, offsets,
                                        fault_indices);
  if (status != COMB_CUDA_STATUS_SUCCESS) return status;
  dense_status = comb_dense_dcp_init(&combined, context->flat->n_all);
  if (dense_status != COMB_DENSE_OK) return reference_dense_status(dense_status);
  dense_status = comb_dense_dcp_init(&error, context->flat->n_all);
  if (dense_status != COMB_DENSE_OK) {
    status = reference_dense_status(dense_status);
    goto cleanup;
  }
  status = reference_wrap_flat_dcp(context, 0U, fault_indices[0], &view);
  if (status != COMB_CUDA_STATUS_SUCCESS) goto cleanup;
  status = reference_dcp_copy(&combined, &view);
  comb_dense_dcp_destroy(&view);
  if (status != COMB_CUDA_STATUS_SUCCESS) goto cleanup;
  for (selected = 1U; selected < context->flat->s; ++selected) {
    status =
        reference_wrap_flat_dcp(context, 0U, fault_indices[selected], &view);
    if (status != COMB_CUDA_STATUS_SUCCESS) goto cleanup;
    dense_status = comb_dense_dcp_multiply(&combined, &view, &combined);
    comb_dense_dcp_destroy(&view);
    if (dense_status != COMB_DENSE_OK) {
      status = reference_dense_status(dense_status);
      goto cleanup;
    }
  }
  for (logical_mask = 0U; logical_mask < context->logical_variant_count;
       ++logical_mask) {
    unsigned good = 0U;
    uint32_t logical;
    status = reference_dcp_copy(&error, &combined);
    if (status != COMB_CUDA_STATUS_SUCCESS) goto cleanup;
    for (logical = 0U; logical < context->flat->k; ++logical) {
      if (((logical_mask >> logical) & UINT64_C(1)) != 0U) {
        status = reference_wrap_flat_dcp(context, 1U, logical, &view);
        if (status != COMB_CUDA_STATUS_SUCCESS) goto cleanup;
        dense_status = comb_dense_dcp_multiply(&error, &view, &error);
        comb_dense_dcp_destroy(&view);
        if (dense_status != COMB_DENSE_OK) {
          status = reference_dense_status(dense_status);
          goto cleanup;
        }
      }
    }
    status = reference_variant_good(context, &error, &good);
    if (status != COMB_CUDA_STATUS_SUCCESS) goto cleanup;
    if (good != 0U) {
      *bad = 0U;
      goto cleanup;
    }
  }
  *bad = 1U;

cleanup:
  comb_dense_dcp_destroy(&view);
  comb_dense_dcp_destroy(&error);
  comb_dense_dcp_destroy(&combined);
  return status;
}

static comb_cuda_status_t reference_validate_outputs(
    const comb_reference_context_t *context, comb_cuda_witness_mode_t mode,
    comb_cuda_result_t *result, comb_cuda_witness_t *witness) {
  if (context == NULL || result == NULL ||
      result->struct_size < sizeof(*result) ||
      result->abi_version != COMB_CUDA_ABI_VERSION ||
      (mode != COMB_CUDA_WITNESS_DEFAULT &&
       mode != COMB_CUDA_WITNESS_DETERMINISTIC &&
       mode != COMB_CUDA_WITNESS_FAST_ANY))
    return COMB_CUDA_STATUS_INVALID_ARGUMENT;
  if (witness != NULL) {
    if (witness->struct_size < sizeof(*witness) ||
        witness->abi_version != COMB_CUDA_ABI_VERSION ||
        witness->fault_capacity < context->flat->s ||
        witness->fault_indices == NULL)
      return COMB_CUDA_STATUS_INVALID_ARGUMENT;
  }
  return COMB_CUDA_STATUS_SUCCESS;
}

comb_cuda_status_t comb_reference_verify_range(
    comb_reference_context_t *context, uint64_t rank_begin,
    uint64_t rank_end, comb_cuda_witness_mode_t mode,
    comb_cuda_result_t *result, comb_cuda_witness_t *witness) {
  uint32_t *locations = NULL;
  uint32_t *offsets = NULL;
  uint32_t *fault_indices = NULL;
  uint64_t rank;
  uint64_t minimum_bad = COMB_CUDA_RANK_NONE;
  comb_cuda_status_t status;
  double started;
  status = reference_validate_outputs(context, mode, result, witness);
  if (status != COMB_CUDA_STATUS_SUCCESS) return status;
  if (rank_begin > rank_end ||
      rank_end > context->flat->total_combination_count)
    return COMB_CUDA_STATUS_INVALID_ARGUMENT;
  result->verification_code = OTHER_ERROR;
  result->reserved0 = 0U;
  result->total_combination_count = context->flat->total_combination_count;
  result->checked_combination_count = 0U;
  result->minimum_bad_rank = COMB_CUDA_RANK_NONE;
  result->setup_milliseconds = context->setup_milliseconds;
  result->verification_milliseconds = 0.0;
  if (witness != NULL) {
    witness->rank = COMB_CUDA_RANK_NONE;
    witness->fault_count = 0U;
  }
  status = reference_allocate_u32(context->flat->s, &locations);
  if (status != COMB_CUDA_STATUS_SUCCESS) return status;
  status = reference_allocate_u32(context->flat->s, &offsets);
  if (status != COMB_CUDA_STATUS_SUCCESS) goto cleanup;
  status = reference_allocate_u32(context->flat->s, &fault_indices);
  if (status != COMB_CUDA_STATUS_SUCCESS) goto cleanup;
  started = reference_now_milliseconds();
  for (rank = rank_begin; rank < rank_end; ++rank) {
    unsigned bad = 0U;
    status = reference_tuple_bad(context, rank, locations, offsets,
                                 fault_indices, &bad);
    if (status != COMB_CUDA_STATUS_SUCCESS) {
      result->verification_code = INTERNAL_ERROR;
      result->verification_milliseconds =
          reference_now_milliseconds() - started;
      goto cleanup;
    }
    ++result->checked_combination_count;
    if (bad != 0U && minimum_bad == COMB_CUDA_RANK_NONE) minimum_bad = rank;
  }
  result->minimum_bad_rank = minimum_bad;
  result->verification_code = minimum_bad == COMB_CUDA_RANK_NONE
                                  ? ALL_FAULTS_GOOD
                                  : BAD_FAULT_FOUND;
  result->verification_milliseconds =
      reference_now_milliseconds() - started;
  if (witness != NULL && minimum_bad != COMB_CUDA_RANK_NONE) {
    status = reference_rank_fault_indices(context, minimum_bad, locations,
                                          offsets, witness->fault_indices);
    if (status != COMB_CUDA_STATUS_SUCCESS) {
      result->verification_code = INTERNAL_ERROR;
      goto cleanup;
    }
    witness->rank = minimum_bad;
    witness->fault_count = context->flat->s;
  }

cleanup:
  free(fault_indices);
  free(offsets);
  free(locations);
  return status;
}

comb_cuda_status_t comb_reference_verify(
    comb_reference_context_t *context, comb_cuda_witness_mode_t mode,
    comb_cuda_result_t *result, comb_cuda_witness_t *witness) {
  if (context == NULL) return COMB_CUDA_STATUS_INVALID_ARGUMENT;
  return comb_reference_verify_range(context, 0U,
                                     context->flat->total_combination_count,
                                     mode, result, witness);
}
