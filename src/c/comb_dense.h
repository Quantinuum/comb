#ifndef COMB_DENSE_H
#define COMB_DENSE_H

/*
 * Dependency-free dense GF(2) correctness primitives for Comb.
 *
 * This module is deliberately independent of M4RI and CUDA.  It is intended
 * to be a small, inspectable CPU oracle for differential tests of the packed
 * GPU implementation.  It covers packed matrices, phase-aware stabiliser row
 * operations, deterministic tableau elimination, sequential pure-X/pure-Z
 * postselection, and DCP multiplication.  Fault enumeration, lookup-table
 * policy, and the complete CSS verification criterion are outside its scope.
 *
 * Packed rows use bit i of words[i / 64], with bit i % 64 within that word.
 * Every row is required to have zero padding bits above its logical width.
 */

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
  COMB_DENSE_OK = 0,
  COMB_DENSE_INVALID_ARGUMENT = 1,
  COMB_DENSE_OVERFLOW = 2,
  COMB_DENSE_ALLOCATION_FAILED = 3,
  COMB_DENSE_DIMENSION_MISMATCH = 4,
  COMB_DENSE_OUT_OF_RANGE = 5,
  COMB_DENSE_NONZERO_PADDING = 6,
  COMB_DENSE_NOT_GENERATED = 7,
  COMB_DENSE_INCONSISTENT_TABLEAU = 8
} comb_dense_status_t;

typedef struct {
  size_t rows;
  size_t columns;
  size_t stride_words;
  uint64_t *words;
  unsigned owns_words;
} comb_dense_matrix_t;

/* Compute the canonical packed layout, with checked size arithmetic. */
comb_dense_status_t comb_dense_matrix_layout(size_t rows, size_t columns,
                                             size_t *stride_words,
                                             size_t *word_count);

/* Allocate a zero-filled, canonically strided matrix. */
comb_dense_status_t comb_dense_matrix_init(comb_dense_matrix_t *matrix,
                                           size_t rows, size_t columns);

/*
 * Wrap non-owning storage.  stride_words may exceed the canonical stride, but
 * logical rows must not overlap and logical tail bits must already be zero.
 */
comb_dense_status_t comb_dense_matrix_wrap(comb_dense_matrix_t *matrix,
                                           size_t rows, size_t columns,
                                           size_t stride_words,
                                           uint64_t *words);

comb_dense_status_t comb_dense_matrix_validate(
    const comb_dense_matrix_t *matrix);
void comb_dense_matrix_mask_tail(comb_dense_matrix_t *matrix);
void comb_dense_matrix_destroy(comb_dense_matrix_t *matrix);

comb_dense_status_t comb_dense_matrix_get(const comb_dense_matrix_t *matrix,
                                          size_t row, size_t column,
                                          unsigned *value);
comb_dense_status_t comb_dense_matrix_set(comb_dense_matrix_t *matrix,
                                          size_t row, size_t column,
                                          unsigned value);
comb_dense_status_t comb_dense_matrix_row_parity(
    const comb_dense_matrix_t *left, size_t left_row,
    const comb_dense_matrix_t *right, size_t right_row, unsigned *parity);
comb_dense_status_t comb_dense_matrix_row_xor(
    comb_dense_matrix_t *destination, size_t destination_row,
    const comb_dense_matrix_t *source, size_t source_row);
comb_dense_status_t comb_dense_matrix_row_swap(comb_dense_matrix_t *matrix,
                                               size_t first_row,
                                               size_t second_row);

typedef struct {
  size_t qubits;
  size_t row_count;
  size_t row_capacity;
  comb_dense_matrix_t x;
  comb_dense_matrix_t z;
  uint8_t *delta;
  uint8_t *epsilon;
  unsigned owns_phase;
} comb_dense_tableau_t;

/* Allocate row_capacity rows and expose the first row_count rows. */
comb_dense_status_t comb_dense_tableau_init(comb_dense_tableau_t *tableau,
                                            size_t qubits,
                                            size_t row_capacity,
                                            size_t row_count);

/* Wrap non-owning storage.  x_words and z_words use the same stride. */
comb_dense_status_t comb_dense_tableau_wrap(
    comb_dense_tableau_t *tableau, size_t qubits, size_t row_capacity,
    size_t row_count, size_t stride_words, uint64_t *x_words,
    uint64_t *z_words, uint8_t *delta, uint8_t *epsilon);

comb_dense_status_t comb_dense_tableau_validate(
    const comb_dense_tableau_t *tableau);
void comb_dense_tableau_destroy(comb_dense_tableau_t *tableau);

comb_dense_status_t comb_dense_tableau_row_swap(
    comb_dense_tableau_t *tableau, size_t first_row, size_t second_row);

/*
 * Multiply destination Pauli row by source Pauli row using Comb's exact
 * Dehaene-De Moor phase convention.  In particular, source-x dot
 * destination-z, delta(source), and delta(destination) are read before either
 * packed row is changed.
 */
comb_dense_status_t comb_dense_tableau_row_add(
    const comb_dense_tableau_t *source, size_t source_row,
    comb_dense_tableau_t *destination, size_t destination_row);

/*
 * Reproduce form_identity_ge using lowest-index pivots.  The active tableau
 * must contain qubits + 1 rows; the last row is reduced to +/-I.  If it cannot
 * be reduced to I, COMB_DENSE_NOT_GENERATED is returned.
 */
comb_dense_status_t comb_dense_tableau_form_identity(
    comb_dense_tableau_t *tableau, unsigned *negative_identity);

/*
 * Reproduce set_coset_to_stabiliser_state's deterministic Z-then-X
 * elimination and X back-elimination.  The tableau must have qubits active
 * rows.  z_rank receives the number of independent Z parts.
 */
comb_dense_status_t comb_dense_tableau_reduce_coset(
    comb_dense_tableau_t *tableau, size_t *z_rank);

typedef enum {
  COMB_DENSE_MEASURE_X = 0,
  COMB_DENSE_MEASURE_Z = 1
} comb_dense_measurement_basis_t;

/*
 * Sequentially postselect pure-X or pure-Z measurements onto +1.  This is the
 * direct, measurement-by-measurement oracle corresponding to dcp.c, not the
 * optimized batched-Z transform.  tableau must have qubits active rows and at
 * least qubits + 1 rows of capacity.  rejected is set when -1 is forced.
 */
comb_dense_status_t comb_dense_tableau_postselect(
    comb_dense_tableau_t *tableau,
    const comb_dense_matrix_t *measurements,
    comb_dense_measurement_basis_t basis, unsigned *rejected);

typedef struct {
  size_t qubits;
  size_t word_count;
  comb_dense_matrix_t f;
  uint64_t *h_x;
  uint64_t *h_z;
  unsigned owns_vectors;
} comb_dense_dcp_t;

comb_dense_status_t comb_dense_dcp_init(comb_dense_dcp_t *dcp,
                                        size_t qubits);
comb_dense_status_t comb_dense_dcp_wrap(comb_dense_dcp_t *dcp,
                                        size_t qubits, size_t f_stride_words,
                                        uint64_t *f_words, uint64_t *h_x,
                                        uint64_t *h_z);
comb_dense_status_t comb_dense_dcp_validate(const comb_dense_dcp_t *dcp,
                                            unsigned require_symmetric_f);
void comb_dense_dcp_mask_tail(comb_dense_dcp_t *dcp);
void comb_dense_dcp_destroy(comb_dense_dcp_t *dcp);

/*
 * Compute c2 * c1.  The operation is alias-safe.  The phase vector is
 * h_z2 + h_z1 + F1 h_x2 + diag(F1) diag(F2), over GF(2), matching dcp.c.
 */
comb_dense_status_t comb_dense_dcp_multiply(comb_dense_dcp_t *result,
                                            const comb_dense_dcp_t *c2,
                                            const comb_dense_dcp_t *c1);

#ifdef __cplusplus
}
#endif

#endif /* COMB_DENSE_H */
