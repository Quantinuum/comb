#include "comb_dense.h"

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

#define CHECK_OK(expression) CHECK((expression) == COMB_DENSE_OK)

enum {
  RANDOM_ROW_CASES = 10000,
  RANDOM_PHASE_CASES = 10000,
  RANDOM_DCP_CASES = 10000,
  RANDOM_COSET_CASES = 2000,
  MAX_REFERENCE_QUBITS = 140,
  MAX_REFERENCE_ROWS = 16,
  MAX_DCP_QUBITS = 17
};

static uint64_t random_state = UINT64_C(0xd1b54a32d192ed03);

static uint64_t next_random(void) {
  uint64_t value = random_state;
  value ^= value >> 12U;
  value ^= value << 25U;
  value ^= value >> 27U;
  random_state = value;
  return value * UINT64_C(2685821657736338717);
}

static unsigned random_bit(void) { return (unsigned)(next_random() & 1U); }

static unsigned matrix_bit(const comb_dense_matrix_t *matrix, size_t row,
                           size_t column) {
  unsigned value = 2U;
  CHECK_OK(comb_dense_matrix_get(matrix, row, column, &value));
  CHECK(value <= 1U);
  return value;
}

static void test_layout_and_validation(void) {
  size_t stride = 99U;
  size_t count = 99U;
  uint64_t storage[6] = {0U, 0U, 0U, 0U, 0U, 0U};
  comb_dense_matrix_t matrix;
  CHECK_OK(comb_dense_matrix_layout(3U, 65U, &stride, &count));
  CHECK(stride == 2U);
  CHECK(count == 6U);
  CHECK(comb_dense_matrix_layout(SIZE_MAX, 65U, &stride, &count) ==
        COMB_DENSE_OVERFLOW);
  CHECK(comb_dense_matrix_wrap(&matrix, 3U, 65U, 1U, storage) ==
        COMB_DENSE_DIMENSION_MISMATCH);
  storage[1] = UINT64_C(2);
  CHECK(comb_dense_matrix_wrap(&matrix, 3U, 65U, 2U, storage) ==
        COMB_DENSE_NONZERO_PADDING);
  storage[1] = UINT64_C(1);
  CHECK_OK(comb_dense_matrix_wrap(&matrix, 3U, 65U, 2U, storage));
  CHECK(comb_dense_matrix_set(&matrix, 3U, 0U, 1U) ==
        COMB_DENSE_OUT_OF_RANGE);
  CHECK(comb_dense_matrix_set(&matrix, 0U, 65U, 1U) ==
        COMB_DENSE_OUT_OF_RANGE);
  CHECK(comb_dense_matrix_set(&matrix, 0U, 0U, 2U) ==
        COMB_DENSE_INVALID_ARGUMENT);
  comb_dense_matrix_destroy(&matrix);
}

static void test_width_boundaries(void) {
  const size_t widths[] = {31U, 32U, 33U, 63U, 64U, 65U, 127U};
  size_t width_index;
  for (width_index = 0U;
       width_index < sizeof(widths) / sizeof(widths[0]); ++width_index) {
    const size_t width = widths[width_index];
    const size_t words = (width / 64U) + ((width % 64U) != 0U ? 1U : 0U);
    comb_dense_matrix_t matrix;
    size_t bit;
    CHECK_OK(comb_dense_matrix_init(&matrix, 2U, width));
    for (bit = 0U; bit < width; ++bit) {
      CHECK_OK(comb_dense_matrix_set(&matrix, 0U, bit,
                                     (unsigned)((bit * 7U + 3U) & 1U)));
      CHECK_OK(comb_dense_matrix_set(&matrix, 1U, bit,
                                     (unsigned)((bit * 11U + 1U) & 1U)));
    }
    CHECK_OK(comb_dense_matrix_validate(&matrix));
    if ((width % 64U) != 0U) {
      const uint64_t invalid_bit = UINT64_C(1) << (width % 64U);
      matrix.words[words - 1U] |= invalid_bit;
      CHECK(comb_dense_matrix_validate(&matrix) ==
            COMB_DENSE_NONZERO_PADDING);
      comb_dense_matrix_mask_tail(&matrix);
      CHECK_OK(comb_dense_matrix_validate(&matrix));
    }
    comb_dense_matrix_destroy(&matrix);
  }
}

static void test_random_rows(void) {
  size_t iteration;
  for (iteration = 0U; iteration < RANDOM_ROW_CASES; ++iteration) {
    const size_t width = 1U + (size_t)(next_random() % 140U);
    uint8_t left[MAX_REFERENCE_QUBITS] = {0U};
    uint8_t right[MAX_REFERENCE_QUBITS] = {0U};
    comb_dense_matrix_t matrix;
    unsigned expected_parity = 0U;
    unsigned actual_parity = 2U;
    size_t bit;
    CHECK_OK(comb_dense_matrix_init(&matrix, 2U, width));
    for (bit = 0U; bit < width; ++bit) {
      left[bit] = (uint8_t)random_bit();
      right[bit] = (uint8_t)random_bit();
      CHECK_OK(comb_dense_matrix_set(&matrix, 0U, bit, left[bit]));
      CHECK_OK(comb_dense_matrix_set(&matrix, 1U, bit, right[bit]));
      expected_parity ^=
          (unsigned)((unsigned)left[bit] & (unsigned)right[bit]);
    }
    CHECK_OK(comb_dense_matrix_row_parity(&matrix, 0U, &matrix, 1U,
                                          &actual_parity));
    CHECK(actual_parity == expected_parity);
    CHECK_OK(comb_dense_matrix_row_xor(&matrix, 1U, &matrix, 0U));
    for (bit = 0U; bit < width; ++bit)
      CHECK(matrix_bit(&matrix, 1U, bit) ==
            ((unsigned)left[bit] ^ (unsigned)right[bit]));
    CHECK_OK(comb_dense_matrix_row_swap(&matrix, 0U, 1U));
    for (bit = 0U; bit < width; ++bit)
      CHECK(matrix_bit(&matrix, 0U, bit) ==
            ((unsigned)left[bit] ^ (unsigned)right[bit]));
    comb_dense_matrix_destroy(&matrix);
  }
}

static void test_random_phase_row_add(void) {
  size_t iteration;
  for (iteration = 0U; iteration < RANDOM_PHASE_CASES; ++iteration) {
    const size_t qubits = 1U + (size_t)(next_random() % 140U);
    uint8_t source_x[MAX_REFERENCE_QUBITS] = {0U};
    uint8_t source_z[MAX_REFERENCE_QUBITS] = {0U};
    uint8_t destination_x[MAX_REFERENCE_QUBITS] = {0U};
    uint8_t destination_z[MAX_REFERENCE_QUBITS] = {0U};
    const unsigned delta_source = random_bit();
    const unsigned delta_destination = random_bit();
    const unsigned epsilon_source = random_bit();
    const unsigned epsilon_destination = random_bit();
    unsigned dot = 0U;
    unsigned expected_epsilon;
    comb_dense_tableau_t tableau;
    size_t bit;
    CHECK_OK(comb_dense_tableau_init(&tableau, qubits, 2U, 2U));
    for (bit = 0U; bit < qubits; ++bit) {
      source_x[bit] = (uint8_t)random_bit();
      source_z[bit] = (uint8_t)random_bit();
      destination_x[bit] = (uint8_t)random_bit();
      destination_z[bit] = (uint8_t)random_bit();
      CHECK_OK(comb_dense_matrix_set(&tableau.x, 0U, bit, source_x[bit]));
      CHECK_OK(comb_dense_matrix_set(&tableau.z, 0U, bit, source_z[bit]));
      CHECK_OK(
          comb_dense_matrix_set(&tableau.x, 1U, bit, destination_x[bit]));
      CHECK_OK(
          comb_dense_matrix_set(&tableau.z, 1U, bit, destination_z[bit]));
      dot ^=
          (unsigned)((unsigned)source_x[bit] & (unsigned)destination_z[bit]);
    }
    tableau.delta[0] = (uint8_t)delta_source;
    tableau.delta[1] = (uint8_t)delta_destination;
    tableau.epsilon[0] = (uint8_t)epsilon_source;
    tableau.epsilon[1] = (uint8_t)epsilon_destination;
    expected_epsilon = epsilon_destination ^ epsilon_source ^
                       (delta_source & delta_destination) ^ dot;
    CHECK_OK(comb_dense_tableau_row_add(&tableau, 0U, &tableau, 1U));
    CHECK(tableau.delta[1] ==
          (uint8_t)(delta_source ^ delta_destination));
    CHECK(tableau.epsilon[1] == (uint8_t)expected_epsilon);
    for (bit = 0U; bit < qubits; ++bit) {
      CHECK(matrix_bit(&tableau.x, 1U, bit) ==
            ((unsigned)source_x[bit] ^ (unsigned)destination_x[bit]));
      CHECK(matrix_bit(&tableau.z, 1U, bit) ==
            ((unsigned)source_z[bit] ^ (unsigned)destination_z[bit]));
      CHECK(matrix_bit(&tableau.x, 0U, bit) == (unsigned)source_x[bit]);
      CHECK(matrix_bit(&tableau.z, 0U, bit) == (unsigned)source_z[bit]);
    }
    comb_dense_tableau_destroy(&tableau);
  }
}

typedef struct {
  size_t qubits;
  size_t rows;
  uint8_t x[MAX_REFERENCE_ROWS][MAX_REFERENCE_QUBITS];
  uint8_t z[MAX_REFERENCE_ROWS][MAX_REFERENCE_QUBITS];
  uint8_t delta[MAX_REFERENCE_ROWS];
  uint8_t epsilon[MAX_REFERENCE_ROWS];
} reference_tableau_t;

static void reference_row_add(reference_tableau_t *tableau, size_t source,
                              size_t destination) {
  const unsigned delta_source = tableau->delta[source];
  const unsigned delta_destination = tableau->delta[destination];
  unsigned dot = 0U;
  size_t bit;
  for (bit = 0U; bit < tableau->qubits; ++bit)
    dot ^= (unsigned)tableau->x[source][bit] &
           (unsigned)tableau->z[destination][bit];
  for (bit = 0U; bit < tableau->qubits; ++bit) {
    tableau->x[destination][bit] ^= tableau->x[source][bit];
    tableau->z[destination][bit] ^= tableau->z[source][bit];
  }
  tableau->epsilon[destination] =
      (uint8_t)((unsigned)tableau->epsilon[destination] ^
                (unsigned)tableau->epsilon[source] ^
                (delta_source & delta_destination) ^ dot);
  tableau->delta[destination] =
      (uint8_t)(delta_source ^ delta_destination);
}

static void reference_row_swap(reference_tableau_t *tableau, size_t first,
                               size_t second) {
  size_t bit;
  uint8_t temporary;
  for (bit = 0U; bit < tableau->qubits; ++bit) {
    temporary = tableau->x[first][bit];
    tableau->x[first][bit] = tableau->x[second][bit];
    tableau->x[second][bit] = temporary;
    temporary = tableau->z[first][bit];
    tableau->z[first][bit] = tableau->z[second][bit];
    tableau->z[second][bit] = temporary;
  }
  temporary = tableau->delta[first];
  tableau->delta[first] = tableau->delta[second];
  tableau->delta[second] = temporary;
  temporary = tableau->epsilon[first];
  tableau->epsilon[first] = tableau->epsilon[second];
  tableau->epsilon[second] = temporary;
}

static size_t reference_reduce_coset(reference_tableau_t *tableau) {
  size_t row_start = 0U;
  size_t column;
  size_t row;
  for (column = 0U; column < tableau->qubits; ++column) {
    size_t found = SIZE_MAX;
    for (row = row_start; row < tableau->rows; ++row) {
      if (tableau->z[row][column] != 0U) {
        if (found == SIZE_MAX)
          found = row;
        else
          reference_row_add(tableau, found, row);
      }
    }
    if (found != SIZE_MAX) {
      reference_row_swap(tableau, found, row_start);
      ++row_start;
    }
  }
  {
    const size_t z_rank = row_start;
    for (column = 0U; column < tableau->qubits; ++column) {
      size_t found = SIZE_MAX;
      for (row = row_start; row < tableau->rows; ++row) {
        if (tableau->x[row][column] != 0U) {
          if (found == SIZE_MAX)
            found = row;
          else
            reference_row_add(tableau, found, row);
        }
      }
      if (found != SIZE_MAX) {
        reference_row_swap(tableau, found, row_start);
        ++row_start;
      }
    }
    for (row = 0U; row < tableau->rows; ++row) {
      size_t leading = SIZE_MAX;
      for (column = 0U; column < tableau->qubits; ++column) {
        if (tableau->x[row][column] != 0U) {
          leading = column;
          break;
        }
      }
      if (leading != SIZE_MAX) {
        size_t above = row;
        while (above != 0U) {
          --above;
          if (tableau->x[above][leading] != 0U)
            reference_row_add(tableau, row, above);
        }
      }
    }
    return z_rank;
  }
}

static void test_random_coset_reduction(void) {
  size_t iteration;
  for (iteration = 0U; iteration < RANDOM_COSET_CASES; ++iteration) {
    const size_t qubits = 1U + (size_t)(next_random() % 8U);
    comb_dense_tableau_t tableau;
    reference_tableau_t reference;
    size_t expected_z_rank;
    size_t actual_z_rank = SIZE_MAX;
    size_t row;
    size_t bit;
    memset(&reference, 0, sizeof(reference));
    reference.qubits = qubits;
    reference.rows = qubits;
    CHECK_OK(comb_dense_tableau_init(&tableau, qubits, qubits, qubits));
    for (row = 0U; row < qubits; ++row) {
      for (bit = 0U; bit < qubits; ++bit) {
        reference.x[row][bit] = (uint8_t)random_bit();
        reference.z[row][bit] = (uint8_t)random_bit();
        CHECK_OK(comb_dense_matrix_set(&tableau.x, row, bit,
                                       reference.x[row][bit]));
        CHECK_OK(comb_dense_matrix_set(&tableau.z, row, bit,
                                       reference.z[row][bit]));
      }
      reference.delta[row] = (uint8_t)random_bit();
      reference.epsilon[row] = (uint8_t)random_bit();
      tableau.delta[row] = reference.delta[row];
      tableau.epsilon[row] = reference.epsilon[row];
    }
    expected_z_rank = reference_reduce_coset(&reference);
    CHECK_OK(comb_dense_tableau_reduce_coset(&tableau, &actual_z_rank));
    CHECK(actual_z_rank == expected_z_rank);
    for (row = 0U; row < qubits; ++row) {
      CHECK(tableau.delta[row] == reference.delta[row]);
      CHECK(tableau.epsilon[row] == reference.epsilon[row]);
      for (bit = 0U; bit < qubits; ++bit) {
        CHECK(matrix_bit(&tableau.x, row, bit) == reference.x[row][bit]);
        CHECK(matrix_bit(&tableau.z, row, bit) == reference.z[row][bit]);
      }
    }
    comb_dense_tableau_destroy(&tableau);
  }
}

static void set_dcp_f(comb_dense_dcp_t *dcp, size_t row, size_t column,
                      unsigned value) {
  CHECK_OK(comb_dense_matrix_set(&dcp->f, row, column, value));
}

static unsigned dcp_vector_bit(const uint64_t *vector, size_t bit) {
  return (unsigned)((vector[bit / 64U] >> (bit % 64U)) & UINT64_C(1));
}

static void set_dcp_vector_bit(uint64_t *vector, size_t bit, unsigned value) {
  const uint64_t mask = UINT64_C(1) << (bit % 64U);
  if (value != 0U)
    vector[bit / 64U] |= mask;
  else
    vector[bit / 64U] &= ~mask;
}

static void test_random_dcp_multiplication(void) {
  size_t iteration;
  for (iteration = 0U; iteration < RANDOM_DCP_CASES; ++iteration) {
    const size_t qubits = 1U + (size_t)(next_random() % MAX_DCP_QUBITS);
    uint8_t f1[MAX_DCP_QUBITS][MAX_DCP_QUBITS] = {{0U}};
    uint8_t f2[MAX_DCP_QUBITS][MAX_DCP_QUBITS] = {{0U}};
    uint8_t hx1[MAX_DCP_QUBITS] = {0U};
    uint8_t hx2[MAX_DCP_QUBITS] = {0U};
    uint8_t hz1[MAX_DCP_QUBITS] = {0U};
    uint8_t hz2[MAX_DCP_QUBITS] = {0U};
    comb_dense_dcp_t c1;
    comb_dense_dcp_t c2;
    comb_dense_dcp_t result;
    size_t row;
    size_t column;
    CHECK_OK(comb_dense_dcp_init(&c1, qubits));
    CHECK_OK(comb_dense_dcp_init(&c2, qubits));
    CHECK_OK(comb_dense_dcp_init(&result, qubits));
    for (row = 0U; row < qubits; ++row) {
      hx1[row] = (uint8_t)random_bit();
      hx2[row] = (uint8_t)random_bit();
      hz1[row] = (uint8_t)random_bit();
      hz2[row] = (uint8_t)random_bit();
      set_dcp_vector_bit(c1.h_x, row, hx1[row]);
      set_dcp_vector_bit(c2.h_x, row, hx2[row]);
      set_dcp_vector_bit(c1.h_z, row, hz1[row]);
      set_dcp_vector_bit(c2.h_z, row, hz2[row]);
      for (column = row; column < qubits; ++column) {
        const uint8_t first = (uint8_t)random_bit();
        const uint8_t second = (uint8_t)random_bit();
        f1[row][column] = first;
        f1[column][row] = first;
        f2[row][column] = second;
        f2[column][row] = second;
        set_dcp_f(&c1, row, column, first);
        set_dcp_f(&c1, column, row, first);
        set_dcp_f(&c2, row, column, second);
        set_dcp_f(&c2, column, row, second);
      }
    }
    CHECK_OK(comb_dense_dcp_validate(&c1, 1U));
    CHECK_OK(comb_dense_dcp_validate(&c2, 1U));
    CHECK_OK(comb_dense_dcp_multiply(&result, &c2, &c1));
    for (row = 0U; row < qubits; ++row) {
      unsigned expected_hz =
          (unsigned)hz2[row] ^ (unsigned)hz1[row] ^
          ((unsigned)f1[row][row] & (unsigned)f2[row][row]);
      for (column = 0U; column < qubits; ++column) {
        const unsigned expected_f =
            (unsigned)f2[row][column] ^ (unsigned)f1[row][column];
        CHECK(matrix_bit(&result.f, row, column) == expected_f);
        expected_hz ^=
            (unsigned)f1[row][column] & (unsigned)hx2[column];
      }
      CHECK(dcp_vector_bit(result.h_x, row) ==
            ((unsigned)hx2[row] ^ (unsigned)hx1[row]));
      CHECK(dcp_vector_bit(result.h_z, row) == expected_hz);
    }
    if ((iteration % 97U) == 0U) {
      CHECK_OK(comb_dense_dcp_multiply(&c2, &c2, &c1));
      for (row = 0U; row < qubits; ++row) {
        unsigned expected_hz =
            (unsigned)hz2[row] ^ (unsigned)hz1[row] ^
            ((unsigned)f1[row][row] & (unsigned)f2[row][row]);
        for (column = 0U; column < qubits; ++column)
          expected_hz ^=
              (unsigned)f1[row][column] & (unsigned)hx2[column];
        CHECK(dcp_vector_bit(c2.h_z, row) == expected_hz);
      }
    }
    comb_dense_dcp_destroy(&result);
    comb_dense_dcp_destroy(&c2);
    comb_dense_dcp_destroy(&c1);
  }
}

static void test_dcp_order_fixture(void) {
  comb_dense_dcp_t c1;
  comb_dense_dcp_t c2;
  comb_dense_dcp_t forward;
  comb_dense_dcp_t reverse;
  CHECK_OK(comb_dense_dcp_init(&c1, 2U));
  CHECK_OK(comb_dense_dcp_init(&c2, 2U));
  CHECK_OK(comb_dense_dcp_init(&forward, 2U));
  CHECK_OK(comb_dense_dcp_init(&reverse, 2U));
  set_dcp_f(&c1, 0U, 1U, 1U);
  set_dcp_f(&c1, 1U, 0U, 1U);
  set_dcp_vector_bit(c2.h_x, 1U, 1U);
  CHECK_OK(comb_dense_dcp_multiply(&forward, &c2, &c1));
  CHECK_OK(comb_dense_dcp_multiply(&reverse, &c1, &c2));
  CHECK(dcp_vector_bit(forward.h_z, 0U) == 1U);
  CHECK(dcp_vector_bit(reverse.h_z, 0U) == 0U);
  comb_dense_dcp_destroy(&reverse);
  comb_dense_dcp_destroy(&forward);
  comb_dense_dcp_destroy(&c2);
  comb_dense_dcp_destroy(&c1);
}

static void test_postselection_fixtures(void) {
  comb_dense_tableau_t state;
  comb_dense_matrix_t measurement;
  unsigned rejected = 2U;

  /* |0>: +X is random, so it replaces the anticommuting +Z stabiliser. */
  CHECK_OK(comb_dense_tableau_init(&state, 1U, 2U, 1U));
  CHECK_OK(comb_dense_matrix_init(&measurement, 1U, 1U));
  CHECK_OK(comb_dense_matrix_set(&state.z, 0U, 0U, 1U));
  CHECK_OK(comb_dense_matrix_set(&measurement, 0U, 0U, 1U));
  CHECK_OK(comb_dense_tableau_postselect(&state, &measurement,
                                         COMB_DENSE_MEASURE_X, &rejected));
  CHECK(rejected == 0U);
  CHECK(matrix_bit(&state.x, 0U, 0U) == 1U);
  CHECK(matrix_bit(&state.z, 0U, 0U) == 0U);
  comb_dense_matrix_destroy(&measurement);
  comb_dense_tableau_destroy(&state);

  /* -Z forces rejection when postselecting +Z. */
  CHECK_OK(comb_dense_tableau_init(&state, 1U, 2U, 1U));
  CHECK_OK(comb_dense_matrix_init(&measurement, 1U, 1U));
  CHECK_OK(comb_dense_matrix_set(&state.z, 0U, 0U, 1U));
  CHECK_OK(comb_dense_matrix_set(&measurement, 0U, 0U, 1U));
  state.epsilon[0] = 1U;
  CHECK_OK(comb_dense_tableau_postselect(&state, &measurement,
                                         COMB_DENSE_MEASURE_Z, &rejected));
  CHECK(rejected == 1U);
  comb_dense_matrix_destroy(&measurement);
  comb_dense_tableau_destroy(&state);

  /* X and Z are independent, so the last row cannot be reduced to I. */
  CHECK_OK(comb_dense_tableau_init(&state, 1U, 2U, 2U));
  CHECK_OK(comb_dense_matrix_set(&state.x, 0U, 0U, 1U));
  CHECK_OK(comb_dense_matrix_set(&state.z, 1U, 0U, 1U));
  CHECK(comb_dense_tableau_form_identity(&state, &rejected) ==
        COMB_DENSE_NOT_GENERATED);
  comb_dense_tableau_destroy(&state);
}

int main(void) {
  test_layout_and_validation();
  test_width_boundaries();
  test_random_rows();
  test_random_phase_row_add();
  test_random_coset_reduction();
  test_random_dcp_multiplication();
  test_dcp_order_fixture();
  test_postselection_fixtures();
  puts("comb_dense: all deterministic and randomized checks passed");
  return EXIT_SUCCESS;
}
