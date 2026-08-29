#include "comb_dense.h"

#include <limits.h>
#include <stdlib.h>
#include <string.h>

static int comb_dense_checked_multiply(size_t left, size_t right,
                                       size_t *result) {
  if (result == NULL) return 0;
  if (left != 0U && right > SIZE_MAX / left) return 0;
  *result = left * right;
  return 1;
}

static size_t comb_dense_words_for_bits(size_t bits) {
  return (bits / 64U) + ((bits % 64U) != 0U ? 1U : 0U);
}

static uint64_t comb_dense_tail_mask(size_t columns) {
  const size_t remainder = columns % 64U;
  if (remainder == 0U) return UINT64_MAX;
  return (UINT64_C(1) << remainder) - UINT64_C(1);
}

static unsigned comb_dense_parity64(uint64_t value) {
  value ^= value >> 32U;
  value ^= value >> 16U;
  value ^= value >> 8U;
  value ^= value >> 4U;
  value &= UINT64_C(0xf);
  return (unsigned)((UINT64_C(0x6996) >> value) & UINT64_C(1));
}

static unsigned comb_dense_word_bit(const uint64_t *words, size_t bit) {
  return (unsigned)((words[bit / 64U] >> (bit % 64U)) & UINT64_C(1));
}

static void comb_dense_word_set(uint64_t *words, size_t bit, unsigned value) {
  const uint64_t mask = UINT64_C(1) << (bit % 64U);
  if (value != 0U)
    words[bit / 64U] |= mask;
  else
    words[bit / 64U] &= ~mask;
}

comb_dense_status_t comb_dense_matrix_layout(size_t rows, size_t columns,
                                             size_t *stride_words,
                                             size_t *word_count) {
  size_t stride;
  size_t count;
  if (stride_words == NULL || word_count == NULL)
    return COMB_DENSE_INVALID_ARGUMENT;
  stride = comb_dense_words_for_bits(columns);
  if (!comb_dense_checked_multiply(rows, stride, &count))
    return COMB_DENSE_OVERFLOW;
  *stride_words = stride;
  *word_count = count;
  return COMB_DENSE_OK;
}

comb_dense_status_t comb_dense_matrix_init(comb_dense_matrix_t *matrix,
                                           size_t rows, size_t columns) {
  size_t stride;
  size_t count;
  comb_dense_status_t status;
  if (matrix == NULL) return COMB_DENSE_INVALID_ARGUMENT;
  memset(matrix, 0, sizeof(*matrix));
  status = comb_dense_matrix_layout(rows, columns, &stride, &count);
  if (status != COMB_DENSE_OK) return status;
  if (count > SIZE_MAX / sizeof(uint64_t)) return COMB_DENSE_OVERFLOW;
  if (count != 0U) {
    matrix->words = (uint64_t *)calloc(count, sizeof(uint64_t));
    if (matrix->words == NULL) return COMB_DENSE_ALLOCATION_FAILED;
  }
  matrix->rows = rows;
  matrix->columns = columns;
  matrix->stride_words = stride;
  matrix->owns_words = 1U;
  return COMB_DENSE_OK;
}

comb_dense_status_t comb_dense_matrix_wrap(comb_dense_matrix_t *matrix,
                                           size_t rows, size_t columns,
                                           size_t stride_words,
                                           uint64_t *words) {
  size_t canonical_stride;
  size_t ignored_count;
  size_t storage_count;
  comb_dense_status_t status;
  if (matrix == NULL) return COMB_DENSE_INVALID_ARGUMENT;
  memset(matrix, 0, sizeof(*matrix));
  status = comb_dense_matrix_layout(rows, columns, &canonical_stride,
                                    &ignored_count);
  if (status != COMB_DENSE_OK) return status;
  if (stride_words < canonical_stride)
    return COMB_DENSE_DIMENSION_MISMATCH;
  if (!comb_dense_checked_multiply(rows, stride_words, &storage_count))
    return COMB_DENSE_OVERFLOW;
  if (storage_count != 0U && words == NULL)
    return COMB_DENSE_INVALID_ARGUMENT;
  matrix->rows = rows;
  matrix->columns = columns;
  matrix->stride_words = stride_words;
  matrix->words = words;
  matrix->owns_words = 0U;
  status = comb_dense_matrix_validate(matrix);
  if (status != COMB_DENSE_OK) memset(matrix, 0, sizeof(*matrix));
  return status;
}

comb_dense_status_t comb_dense_matrix_validate(
    const comb_dense_matrix_t *matrix) {
  size_t canonical_stride;
  size_t ignored_count;
  size_t storage_count;
  size_t row;
  size_t word;
  comb_dense_status_t status;
  if (matrix == NULL) return COMB_DENSE_INVALID_ARGUMENT;
  status = comb_dense_matrix_layout(matrix->rows, matrix->columns,
                                    &canonical_stride, &ignored_count);
  if (status != COMB_DENSE_OK) return status;
  if (matrix->stride_words < canonical_stride)
    return COMB_DENSE_DIMENSION_MISMATCH;
  if (!comb_dense_checked_multiply(matrix->rows, matrix->stride_words,
                                   &storage_count))
    return COMB_DENSE_OVERFLOW;
  if (storage_count != 0U && matrix->words == NULL)
    return COMB_DENSE_INVALID_ARGUMENT;
  for (row = 0U; row < matrix->rows; ++row) {
    const uint64_t *row_words =
        matrix->words + (row * matrix->stride_words);
    if (canonical_stride != 0U && (matrix->columns % 64U) != 0U &&
        (row_words[canonical_stride - 1U] &
         ~comb_dense_tail_mask(matrix->columns)) != 0U)
      return COMB_DENSE_NONZERO_PADDING;
    for (word = canonical_stride; word < matrix->stride_words; ++word) {
      if (row_words[word] != 0U) return COMB_DENSE_NONZERO_PADDING;
    }
  }
  return COMB_DENSE_OK;
}

void comb_dense_matrix_mask_tail(comb_dense_matrix_t *matrix) {
  size_t canonical_stride;
  size_t row;
  size_t word;
  if (matrix == NULL || matrix->words == NULL) return;
  canonical_stride = comb_dense_words_for_bits(matrix->columns);
  for (row = 0U; row < matrix->rows; ++row) {
    uint64_t *row_words = matrix->words + (row * matrix->stride_words);
    if (canonical_stride != 0U && (matrix->columns % 64U) != 0U)
      row_words[canonical_stride - 1U] &=
          comb_dense_tail_mask(matrix->columns);
    for (word = canonical_stride; word < matrix->stride_words; ++word)
      row_words[word] = 0U;
  }
}

void comb_dense_matrix_destroy(comb_dense_matrix_t *matrix) {
  if (matrix == NULL) return;
  if (matrix->owns_words != 0U) free(matrix->words);
  memset(matrix, 0, sizeof(*matrix));
}

comb_dense_status_t comb_dense_matrix_get(const comb_dense_matrix_t *matrix,
                                          size_t row, size_t column,
                                          unsigned *value) {
  comb_dense_status_t status;
  if (value == NULL) return COMB_DENSE_INVALID_ARGUMENT;
  status = comb_dense_matrix_validate(matrix);
  if (status != COMB_DENSE_OK) return status;
  if (row >= matrix->rows || column >= matrix->columns)
    return COMB_DENSE_OUT_OF_RANGE;
  *value = comb_dense_word_bit(
      matrix->words + (row * matrix->stride_words), column);
  return COMB_DENSE_OK;
}

comb_dense_status_t comb_dense_matrix_set(comb_dense_matrix_t *matrix,
                                          size_t row, size_t column,
                                          unsigned value) {
  comb_dense_status_t status;
  if (value > 1U) return COMB_DENSE_INVALID_ARGUMENT;
  status = comb_dense_matrix_validate(matrix);
  if (status != COMB_DENSE_OK) return status;
  if (row >= matrix->rows || column >= matrix->columns)
    return COMB_DENSE_OUT_OF_RANGE;
  comb_dense_word_set(matrix->words + (row * matrix->stride_words), column,
                      value);
  return COMB_DENSE_OK;
}

comb_dense_status_t comb_dense_matrix_row_parity(
    const comb_dense_matrix_t *left, size_t left_row,
    const comb_dense_matrix_t *right, size_t right_row, unsigned *parity) {
  size_t word_count;
  size_t word;
  unsigned result = 0U;
  comb_dense_status_t status;
  if (parity == NULL) return COMB_DENSE_INVALID_ARGUMENT;
  status = comb_dense_matrix_validate(left);
  if (status != COMB_DENSE_OK) return status;
  status = comb_dense_matrix_validate(right);
  if (status != COMB_DENSE_OK) return status;
  if (left->columns != right->columns)
    return COMB_DENSE_DIMENSION_MISMATCH;
  if (left_row >= left->rows || right_row >= right->rows)
    return COMB_DENSE_OUT_OF_RANGE;
  word_count = comb_dense_words_for_bits(left->columns);
  for (word = 0U; word < word_count; ++word) {
    const uint64_t left_word =
        left->words[(left_row * left->stride_words) + word];
    const uint64_t right_word =
        right->words[(right_row * right->stride_words) + word];
    result ^= comb_dense_parity64(left_word & right_word);
  }
  *parity = result & 1U;
  return COMB_DENSE_OK;
}

comb_dense_status_t comb_dense_matrix_row_xor(
    comb_dense_matrix_t *destination, size_t destination_row,
    const comb_dense_matrix_t *source, size_t source_row) {
  size_t word_count;
  size_t word;
  comb_dense_status_t status;
  status = comb_dense_matrix_validate(destination);
  if (status != COMB_DENSE_OK) return status;
  status = comb_dense_matrix_validate(source);
  if (status != COMB_DENSE_OK) return status;
  if (destination->columns != source->columns)
    return COMB_DENSE_DIMENSION_MISMATCH;
  if (destination_row >= destination->rows || source_row >= source->rows)
    return COMB_DENSE_OUT_OF_RANGE;
  word_count = comb_dense_words_for_bits(destination->columns);
  for (word = 0U; word < word_count; ++word)
    destination->words[(destination_row * destination->stride_words) + word] ^=
        source->words[(source_row * source->stride_words) + word];
  if (word_count != 0U && (destination->columns % 64U) != 0U)
    destination->words[(destination_row * destination->stride_words) +
                       word_count - 1U] &=
        comb_dense_tail_mask(destination->columns);
  return COMB_DENSE_OK;
}

comb_dense_status_t comb_dense_matrix_row_swap(comb_dense_matrix_t *matrix,
                                               size_t first_row,
                                               size_t second_row) {
  size_t word;
  comb_dense_status_t status = comb_dense_matrix_validate(matrix);
  if (status != COMB_DENSE_OK) return status;
  if (first_row >= matrix->rows || second_row >= matrix->rows)
    return COMB_DENSE_OUT_OF_RANGE;
  for (word = 0U; word < matrix->stride_words; ++word) {
    uint64_t *first =
        &matrix->words[(first_row * matrix->stride_words) + word];
    uint64_t *second =
        &matrix->words[(second_row * matrix->stride_words) + word];
    const uint64_t temporary = *first;
    *first = *second;
    *second = temporary;
  }
  return COMB_DENSE_OK;
}

comb_dense_status_t comb_dense_tableau_init(comb_dense_tableau_t *tableau,
                                            size_t qubits,
                                            size_t row_capacity,
                                            size_t row_count) {
  comb_dense_status_t status;
  if (tableau == NULL || row_count > row_capacity)
    return COMB_DENSE_INVALID_ARGUMENT;
  memset(tableau, 0, sizeof(*tableau));
  status = comb_dense_matrix_init(&tableau->x, row_capacity, qubits);
  if (status != COMB_DENSE_OK) return status;
  status = comb_dense_matrix_init(&tableau->z, row_capacity, qubits);
  if (status != COMB_DENSE_OK) {
    comb_dense_matrix_destroy(&tableau->x);
    return status;
  }
  if (row_capacity > SIZE_MAX / sizeof(uint8_t)) {
    comb_dense_tableau_destroy(tableau);
    return COMB_DENSE_OVERFLOW;
  }
  tableau->qubits = qubits;
  tableau->row_capacity = row_capacity;
  tableau->row_count = row_count;
  tableau->owns_phase = 1U;
  if (row_capacity != 0U) {
    tableau->delta = (uint8_t *)calloc(row_capacity, sizeof(uint8_t));
    tableau->epsilon = (uint8_t *)calloc(row_capacity, sizeof(uint8_t));
    if (tableau->delta == NULL || tableau->epsilon == NULL) {
      comb_dense_tableau_destroy(tableau);
      return COMB_DENSE_ALLOCATION_FAILED;
    }
  }
  return COMB_DENSE_OK;
}

comb_dense_status_t comb_dense_tableau_wrap(
    comb_dense_tableau_t *tableau, size_t qubits, size_t row_capacity,
    size_t row_count, size_t stride_words, uint64_t *x_words,
    uint64_t *z_words, uint8_t *delta, uint8_t *epsilon) {
  comb_dense_status_t status;
  if (tableau == NULL || row_count > row_capacity ||
      (row_capacity != 0U && (delta == NULL || epsilon == NULL)))
    return COMB_DENSE_INVALID_ARGUMENT;
  memset(tableau, 0, sizeof(*tableau));
  status = comb_dense_matrix_wrap(&tableau->x, row_capacity, qubits,
                                  stride_words, x_words);
  if (status != COMB_DENSE_OK) return status;
  status = comb_dense_matrix_wrap(&tableau->z, row_capacity, qubits,
                                  stride_words, z_words);
  if (status != COMB_DENSE_OK) {
    memset(tableau, 0, sizeof(*tableau));
    return status;
  }
  tableau->qubits = qubits;
  tableau->row_capacity = row_capacity;
  tableau->row_count = row_count;
  tableau->delta = delta;
  tableau->epsilon = epsilon;
  tableau->owns_phase = 0U;
  status = comb_dense_tableau_validate(tableau);
  if (status != COMB_DENSE_OK) memset(tableau, 0, sizeof(*tableau));
  return status;
}

comb_dense_status_t comb_dense_tableau_validate(
    const comb_dense_tableau_t *tableau) {
  size_t row;
  comb_dense_status_t status;
  if (tableau == NULL) return COMB_DENSE_INVALID_ARGUMENT;
  if (tableau->row_count > tableau->row_capacity)
    return COMB_DENSE_DIMENSION_MISMATCH;
  if (tableau->x.rows != tableau->row_capacity ||
      tableau->z.rows != tableau->row_capacity ||
      tableau->x.columns != tableau->qubits ||
      tableau->z.columns != tableau->qubits)
    return COMB_DENSE_DIMENSION_MISMATCH;
  status = comb_dense_matrix_validate(&tableau->x);
  if (status != COMB_DENSE_OK) return status;
  status = comb_dense_matrix_validate(&tableau->z);
  if (status != COMB_DENSE_OK) return status;
  if (tableau->row_capacity != 0U &&
      (tableau->delta == NULL || tableau->epsilon == NULL))
    return COMB_DENSE_INVALID_ARGUMENT;
  for (row = 0U; row < tableau->row_capacity; ++row) {
    if (tableau->delta[row] > 1U || tableau->epsilon[row] > 1U)
      return COMB_DENSE_INCONSISTENT_TABLEAU;
  }
  return COMB_DENSE_OK;
}

void comb_dense_tableau_destroy(comb_dense_tableau_t *tableau) {
  if (tableau == NULL) return;
  comb_dense_matrix_destroy(&tableau->x);
  comb_dense_matrix_destroy(&tableau->z);
  if (tableau->owns_phase != 0U) {
    free(tableau->delta);
    free(tableau->epsilon);
  }
  memset(tableau, 0, sizeof(*tableau));
}

comb_dense_status_t comb_dense_tableau_row_swap(
    comb_dense_tableau_t *tableau, size_t first_row, size_t second_row) {
  uint8_t temporary;
  comb_dense_status_t status = comb_dense_tableau_validate(tableau);
  if (status != COMB_DENSE_OK) return status;
  if (first_row >= tableau->row_count || second_row >= tableau->row_count)
    return COMB_DENSE_OUT_OF_RANGE;
  status = comb_dense_matrix_row_swap(&tableau->x, first_row, second_row);
  if (status != COMB_DENSE_OK) return status;
  status = comb_dense_matrix_row_swap(&tableau->z, first_row, second_row);
  if (status != COMB_DENSE_OK) return status;
  temporary = tableau->delta[first_row];
  tableau->delta[first_row] = tableau->delta[second_row];
  tableau->delta[second_row] = temporary;
  temporary = tableau->epsilon[first_row];
  tableau->epsilon[first_row] = tableau->epsilon[second_row];
  tableau->epsilon[second_row] = temporary;
  return COMB_DENSE_OK;
}

comb_dense_status_t comb_dense_tableau_row_add(
    const comb_dense_tableau_t *source, size_t source_row,
    comb_dense_tableau_t *destination, size_t destination_row) {
  const uint8_t delta_source =
      (source != NULL && source->delta != NULL && source_row < source->row_count)
          ? source->delta[source_row]
          : 0U;
  uint8_t delta_destination;
  unsigned a2_t_u_a1;
  comb_dense_status_t status;
  status = comb_dense_tableau_validate(source);
  if (status != COMB_DENSE_OK) return status;
  status = comb_dense_tableau_validate(destination);
  if (status != COMB_DENSE_OK) return status;
  if (source->qubits != destination->qubits)
    return COMB_DENSE_DIMENSION_MISMATCH;
  if (source_row >= source->row_count ||
      destination_row >= destination->row_count)
    return COMB_DENSE_OUT_OF_RANGE;
  delta_destination = destination->delta[destination_row];
  status = comb_dense_matrix_row_parity(&source->x, source_row,
                                        &destination->z, destination_row,
                                        &a2_t_u_a1);
  if (status != COMB_DENSE_OK) return status;
  status = comb_dense_matrix_row_xor(&destination->x, destination_row,
                                     &source->x, source_row);
  if (status != COMB_DENSE_OK) return status;
  status = comb_dense_matrix_row_xor(&destination->z, destination_row,
                                     &source->z, source_row);
  if (status != COMB_DENSE_OK) return status;
  destination->epsilon[destination_row] =
      (uint8_t)(destination->epsilon[destination_row] ^
                source->epsilon[source_row]);
  destination->delta[destination_row] =
      (uint8_t)(delta_destination ^ delta_source);
  destination->epsilon[destination_row] =
      (uint8_t)(destination->epsilon[destination_row] ^
                (uint8_t)(delta_source & delta_destination) ^
                (uint8_t)a2_t_u_a1);
  return COMB_DENSE_OK;
}

static unsigned comb_dense_matrix_bit_unchecked(
    const comb_dense_matrix_t *matrix, size_t row, size_t column) {
  return comb_dense_word_bit(matrix->words + (row * matrix->stride_words),
                             column);
}

comb_dense_status_t comb_dense_tableau_form_identity(
    comb_dense_tableau_t *tableau, unsigned *negative_identity) {
  size_t column;
  size_t row_start = 0U;
  size_t target;
  comb_dense_status_t status;
  if (negative_identity == NULL) return COMB_DENSE_INVALID_ARGUMENT;
  status = comb_dense_tableau_validate(tableau);
  if (status != COMB_DENSE_OK) return status;
  if (tableau->qubits == SIZE_MAX ||
      tableau->row_count != tableau->qubits + 1U)
    return COMB_DENSE_DIMENSION_MISMATCH;
  for (column = 0U; column < tableau->qubits; ++column) {
    size_t row;
    size_t found = SIZE_MAX;
    for (row = row_start; row < tableau->row_count; ++row) {
      if (comb_dense_matrix_bit_unchecked(&tableau->x, row, column) != 0U) {
        if (found == SIZE_MAX)
          found = row;
        else {
          status = comb_dense_tableau_row_add(tableau, found, tableau, row);
          if (status != COMB_DENSE_OK) return status;
        }
      }
    }
    if (found != SIZE_MAX) {
      status = comb_dense_tableau_row_swap(tableau, found, row_start);
      if (status != COMB_DENSE_OK) return status;
      ++row_start;
    }
    found = SIZE_MAX;
    for (row = row_start; row < tableau->row_count; ++row) {
      if (comb_dense_matrix_bit_unchecked(&tableau->z, row, column) != 0U) {
        if (found == SIZE_MAX)
          found = row;
        else {
          status = comb_dense_tableau_row_add(tableau, found, tableau, row);
          if (status != COMB_DENSE_OK) return status;
        }
      }
    }
    if (found != SIZE_MAX) {
      status = comb_dense_tableau_row_swap(tableau, found, row_start);
      if (status != COMB_DENSE_OK) return status;
      ++row_start;
    }
  }
  target = tableau->qubits;
  for (column = 0U; column < tableau->qubits; ++column) {
    if (comb_dense_matrix_bit_unchecked(&tableau->x, target, column) != 0U ||
        comb_dense_matrix_bit_unchecked(&tableau->z, target, column) != 0U)
      return COMB_DENSE_NOT_GENERATED;
  }
  if (tableau->delta[target] != 0U)
    return COMB_DENSE_INCONSISTENT_TABLEAU;
  *negative_identity = tableau->epsilon[target];
  return COMB_DENSE_OK;
}

comb_dense_status_t comb_dense_tableau_reduce_coset(
    comb_dense_tableau_t *tableau, size_t *z_rank) {
  size_t row_start = 0U;
  size_t column;
  size_t row;
  comb_dense_status_t status;
  if (z_rank == NULL) return COMB_DENSE_INVALID_ARGUMENT;
  status = comb_dense_tableau_validate(tableau);
  if (status != COMB_DENSE_OK) return status;
  if (tableau->row_count != tableau->qubits)
    return COMB_DENSE_DIMENSION_MISMATCH;
  for (column = 0U; column < tableau->qubits; ++column) {
    size_t found = SIZE_MAX;
    for (row = row_start; row < tableau->row_count; ++row) {
      if (comb_dense_matrix_bit_unchecked(&tableau->z, row, column) != 0U) {
        if (found == SIZE_MAX)
          found = row;
        else {
          status = comb_dense_tableau_row_add(tableau, found, tableau, row);
          if (status != COMB_DENSE_OK) return status;
        }
      }
    }
    if (found != SIZE_MAX) {
      status = comb_dense_tableau_row_swap(tableau, found, row_start);
      if (status != COMB_DENSE_OK) return status;
      ++row_start;
    }
  }
  *z_rank = row_start;
  for (column = 0U; column < tableau->qubits; ++column) {
    size_t found = SIZE_MAX;
    for (row = row_start; row < tableau->row_count; ++row) {
      if (comb_dense_matrix_bit_unchecked(&tableau->x, row, column) != 0U) {
        if (found == SIZE_MAX)
          found = row;
        else {
          status = comb_dense_tableau_row_add(tableau, found, tableau, row);
          if (status != COMB_DENSE_OK) return status;
        }
      }
    }
    if (found != SIZE_MAX) {
      status = comb_dense_tableau_row_swap(tableau, found, row_start);
      if (status != COMB_DENSE_OK) return status;
      ++row_start;
    }
  }
  for (row = 0U; row < tableau->row_count; ++row) {
    size_t leading = SIZE_MAX;
    for (column = 0U; column < tableau->qubits; ++column) {
      if (comb_dense_matrix_bit_unchecked(&tableau->x, row, column) != 0U) {
        leading = column;
        break;
      }
    }
    if (leading != SIZE_MAX) {
      size_t row_above = row;
      while (row_above != 0U) {
        --row_above;
        if (comb_dense_matrix_bit_unchecked(&tableau->x, row_above,
                                            leading) != 0U) {
          status =
              comb_dense_tableau_row_add(tableau, row, tableau, row_above);
          if (status != COMB_DENSE_OK) return status;
        }
      }
    }
  }
  return COMB_DENSE_OK;
}

static void comb_dense_tableau_clear_row(comb_dense_tableau_t *tableau,
                                         size_t row) {
  memset(tableau->x.words + (row * tableau->x.stride_words), 0,
         tableau->x.stride_words * sizeof(uint64_t));
  memset(tableau->z.words + (row * tableau->z.stride_words), 0,
         tableau->z.stride_words * sizeof(uint64_t));
  tableau->delta[row] = 0U;
  tableau->epsilon[row] = 0U;
}

static void comb_dense_tableau_set_measurement(
    comb_dense_tableau_t *tableau, size_t row,
    const comb_dense_matrix_t *measurements, size_t measurement_row,
    comb_dense_measurement_basis_t basis) {
  comb_dense_matrix_t *part =
      basis == COMB_DENSE_MEASURE_X ? &tableau->x : &tableau->z;
  const size_t word_count = comb_dense_words_for_bits(tableau->qubits);
  size_t word;
  comb_dense_tableau_clear_row(tableau, row);
  for (word = 0U; word < word_count; ++word)
    part->words[(row * part->stride_words) + word] =
        measurements->words[(measurement_row * measurements->stride_words) +
                            word];
}

comb_dense_status_t comb_dense_tableau_postselect(
    comb_dense_tableau_t *tableau,
    const comb_dense_matrix_t *measurements,
    comb_dense_measurement_basis_t basis, unsigned *rejected) {
  size_t measurement;
  comb_dense_status_t status;
  if (rejected == NULL ||
      (basis != COMB_DENSE_MEASURE_X && basis != COMB_DENSE_MEASURE_Z))
    return COMB_DENSE_INVALID_ARGUMENT;
  *rejected = 0U;
  status = comb_dense_tableau_validate(tableau);
  if (status != COMB_DENSE_OK) return status;
  status = comb_dense_matrix_validate(measurements);
  if (status != COMB_DENSE_OK) return status;
  if (tableau->row_count != tableau->qubits ||
      tableau->row_capacity <= tableau->qubits ||
      measurements->columns != tableau->qubits)
    return COMB_DENSE_DIMENSION_MISMATCH;
  for (measurement = 0U; measurement < measurements->rows; ++measurement) {
    const comb_dense_matrix_t *opposite =
        basis == COMB_DENSE_MEASURE_X ? &tableau->z : &tableau->x;
    size_t pivot = SIZE_MAX;
    size_t row;
    for (row = 0U; row < tableau->qubits; ++row) {
      unsigned anticommutes;
      status = comb_dense_matrix_row_parity(measurements, measurement,
                                            opposite, row, &anticommutes);
      if (status != COMB_DENSE_OK) return status;
      if (anticommutes != 0U) {
        pivot = row;
        break;
      }
    }
    if (pivot != SIZE_MAX) {
      for (row = pivot + 1U; row < tableau->qubits; ++row) {
        unsigned anticommutes;
        status = comb_dense_matrix_row_parity(measurements, measurement,
                                              opposite, row, &anticommutes);
        if (status != COMB_DENSE_OK) return status;
        if (anticommutes != 0U) {
          status =
              comb_dense_tableau_row_add(tableau, pivot, tableau, row);
          if (status != COMB_DENSE_OK) return status;
        }
      }
      comb_dense_tableau_set_measurement(tableau, pivot, measurements,
                                         measurement, basis);
    } else {
      unsigned negative_identity;
      const size_t appended_row = tableau->qubits;
      comb_dense_tableau_set_measurement(tableau, appended_row, measurements,
                                         measurement, basis);
      tableau->row_count = tableau->qubits + 1U;
      status =
          comb_dense_tableau_form_identity(tableau, &negative_identity);
      tableau->row_count = tableau->qubits;
      if (status != COMB_DENSE_OK) return status;
      if (negative_identity != 0U) {
        *rejected = 1U;
        return COMB_DENSE_OK;
      }
    }
  }
  return COMB_DENSE_OK;
}

static int comb_dense_vector_tail_is_zero(const uint64_t *vector,
                                          size_t qubits) {
  const size_t word_count = comb_dense_words_for_bits(qubits);
  if (word_count == 0U) return 1;
  if (vector == NULL) return 0;
  if ((qubits % 64U) == 0U) return 1;
  return (vector[word_count - 1U] & ~comb_dense_tail_mask(qubits)) == 0U;
}

comb_dense_status_t comb_dense_dcp_init(comb_dense_dcp_t *dcp,
                                        size_t qubits) {
  size_t bytes;
  comb_dense_status_t status;
  if (dcp == NULL) return COMB_DENSE_INVALID_ARGUMENT;
  memset(dcp, 0, sizeof(*dcp));
  status = comb_dense_matrix_init(&dcp->f, qubits, qubits);
  if (status != COMB_DENSE_OK) return status;
  dcp->qubits = qubits;
  dcp->word_count = comb_dense_words_for_bits(qubits);
  dcp->owns_vectors = 1U;
  if (!comb_dense_checked_multiply(dcp->word_count, sizeof(uint64_t), &bytes)) {
    comb_dense_dcp_destroy(dcp);
    return COMB_DENSE_OVERFLOW;
  }
  if (bytes != 0U) {
    dcp->h_x = (uint64_t *)calloc(dcp->word_count, sizeof(uint64_t));
    dcp->h_z = (uint64_t *)calloc(dcp->word_count, sizeof(uint64_t));
    if (dcp->h_x == NULL || dcp->h_z == NULL) {
      comb_dense_dcp_destroy(dcp);
      return COMB_DENSE_ALLOCATION_FAILED;
    }
  }
  return COMB_DENSE_OK;
}

comb_dense_status_t comb_dense_dcp_wrap(comb_dense_dcp_t *dcp,
                                        size_t qubits, size_t f_stride_words,
                                        uint64_t *f_words, uint64_t *h_x,
                                        uint64_t *h_z) {
  comb_dense_status_t status;
  if (dcp == NULL) return COMB_DENSE_INVALID_ARGUMENT;
  memset(dcp, 0, sizeof(*dcp));
  status = comb_dense_matrix_wrap(&dcp->f, qubits, qubits, f_stride_words,
                                  f_words);
  if (status != COMB_DENSE_OK) return status;
  dcp->qubits = qubits;
  dcp->word_count = comb_dense_words_for_bits(qubits);
  dcp->h_x = h_x;
  dcp->h_z = h_z;
  dcp->owns_vectors = 0U;
  status = comb_dense_dcp_validate(dcp, 0U);
  if (status != COMB_DENSE_OK) memset(dcp, 0, sizeof(*dcp));
  return status;
}

comb_dense_status_t comb_dense_dcp_validate(const comb_dense_dcp_t *dcp,
                                            unsigned require_symmetric_f) {
  size_t row;
  size_t column;
  comb_dense_status_t status;
  if (dcp == NULL) return COMB_DENSE_INVALID_ARGUMENT;
  if (dcp->f.rows != dcp->qubits || dcp->f.columns != dcp->qubits ||
      dcp->word_count != comb_dense_words_for_bits(dcp->qubits))
    return COMB_DENSE_DIMENSION_MISMATCH;
  status = comb_dense_matrix_validate(&dcp->f);
  if (status != COMB_DENSE_OK) return status;
  if (dcp->word_count != 0U && (dcp->h_x == NULL || dcp->h_z == NULL))
    return COMB_DENSE_INVALID_ARGUMENT;
  if (!comb_dense_vector_tail_is_zero(dcp->h_x, dcp->qubits) ||
      !comb_dense_vector_tail_is_zero(dcp->h_z, dcp->qubits))
    return COMB_DENSE_NONZERO_PADDING;
  if (require_symmetric_f != 0U) {
    for (row = 0U; row < dcp->qubits; ++row) {
      for (column = row + 1U; column < dcp->qubits; ++column) {
        if (comb_dense_matrix_bit_unchecked(&dcp->f, row, column) !=
            comb_dense_matrix_bit_unchecked(&dcp->f, column, row))
          return COMB_DENSE_INVALID_ARGUMENT;
      }
    }
  }
  return COMB_DENSE_OK;
}

void comb_dense_dcp_mask_tail(comb_dense_dcp_t *dcp) {
  if (dcp == NULL) return;
  comb_dense_matrix_mask_tail(&dcp->f);
  if (dcp->word_count != 0U && (dcp->qubits % 64U) != 0U) {
    const uint64_t mask = comb_dense_tail_mask(dcp->qubits);
    if (dcp->h_x != NULL) dcp->h_x[dcp->word_count - 1U] &= mask;
    if (dcp->h_z != NULL) dcp->h_z[dcp->word_count - 1U] &= mask;
  }
}

void comb_dense_dcp_destroy(comb_dense_dcp_t *dcp) {
  if (dcp == NULL) return;
  comb_dense_matrix_destroy(&dcp->f);
  if (dcp->owns_vectors != 0U) {
    free(dcp->h_x);
    free(dcp->h_z);
  }
  memset(dcp, 0, sizeof(*dcp));
}

comb_dense_status_t comb_dense_dcp_multiply(comb_dense_dcp_t *result,
                                            const comb_dense_dcp_t *c2,
                                            const comb_dense_dcp_t *c1) {
  comb_dense_dcp_t temporary;
  size_t row;
  size_t column;
  size_t word;
  comb_dense_status_t status;
  if (result == NULL) return COMB_DENSE_INVALID_ARGUMENT;
  status = comb_dense_dcp_validate(result, 0U);
  if (status != COMB_DENSE_OK) return status;
  status = comb_dense_dcp_validate(c2, 0U);
  if (status != COMB_DENSE_OK) return status;
  status = comb_dense_dcp_validate(c1, 0U);
  if (status != COMB_DENSE_OK) return status;
  if (result->qubits != c2->qubits || result->qubits != c1->qubits)
    return COMB_DENSE_DIMENSION_MISMATCH;
  status = comb_dense_dcp_init(&temporary, result->qubits);
  if (status != COMB_DENSE_OK) return status;
  for (row = 0U; row < result->qubits; ++row) {
    for (word = 0U; word < result->word_count; ++word) {
      temporary.f.words[(row * temporary.f.stride_words) + word] =
          c2->f.words[(row * c2->f.stride_words) + word] ^
          c1->f.words[(row * c1->f.stride_words) + word];
    }
  }
  for (word = 0U; word < result->word_count; ++word) {
    temporary.h_x[word] = c2->h_x[word] ^ c1->h_x[word];
    temporary.h_z[word] = c2->h_z[word] ^ c1->h_z[word];
  }
  for (row = 0U; row < result->qubits; ++row) {
    unsigned phase = 0U;
    for (column = 0U; column < result->qubits; ++column) {
      phase ^=
          comb_dense_matrix_bit_unchecked(&c1->f, row, column) &
          comb_dense_word_bit(c2->h_x, column);
    }
    phase ^=
        comb_dense_matrix_bit_unchecked(&c1->f, row, row) &
        comb_dense_matrix_bit_unchecked(&c2->f, row, row);
    if (phase != 0U)
      comb_dense_word_set(temporary.h_z, row,
                          comb_dense_word_bit(temporary.h_z, row) ^ 1U);
  }
  comb_dense_dcp_mask_tail(&temporary);
  for (row = 0U; row < result->qubits; ++row) {
    memset(result->f.words + (row * result->f.stride_words), 0,
           result->f.stride_words * sizeof(uint64_t));
    memcpy(result->f.words + (row * result->f.stride_words),
           temporary.f.words + (row * temporary.f.stride_words),
           result->word_count * sizeof(uint64_t));
  }
  if (result->word_count != 0U) {
    memcpy(result->h_x, temporary.h_x,
           result->word_count * sizeof(uint64_t));
    memcpy(result->h_z, temporary.h_z,
           result->word_count * sizeof(uint64_t));
  }
  comb_dense_dcp_destroy(&temporary);
  return COMB_DENSE_OK;
}
