#ifndef COMB_H
#define COMB_H

#include <m4ri/m4ri.h>

#include "comb_external.h"
#include "lwlog.h"
#define MAX_SIZE_COUNT 10240
// if TIMING_MODE is set to an integer, then repeat until TIMING_MODE faults
// have been iterated over
#ifdef TIMING_MODE
#define BAD_FAULT_FOUND ALL_FAULTS_GOOD
#endif

// Conventions:
// vectors are stored as column matrices by default (n x 1) matrices.
// However, in certain cases, it is much more preferable to store them as row
// vectors. In such cases, they are /always/ labelled with '_transpose' to avoid
// confusion, even in cases where they are the only form of storage. Matrices
// will considered to have their rows consisting of Paulis (or parts of Paulis).

typedef struct {  // Stores a DCP operator as a tableau.
  mzd_t *F;    // the upper right block element of c (not c-bar). n x n matrix
  mzd_t *h_x;  // n x 1 matrix
  mzd_t *h_x_transpose;  // 1 x n matrix (row vector).
  mzd_t *h_z;            // n x 1 matrix
  mzd_t *F_diagonal;     // diagonal part of F. n x 1 matrix.
} dcp_t;

typedef struct {
  mzd_t *x_parts;       // n x n - rows are the x parts of paulis
  mzd_t *z_parts;       // n x n - rows are the z parts of paulis
  mzd_t *delta_vector;  // n x 1 - entries are the imaginary bits of the vectors
  // idea: delta_vector is determined by x_parts and z_parts as it must always
  // be hermitian...
  mzd_t *epsilon_vector;  // n x 1 - entries are the sign bits of the vectors.
} stabiliser_state_t;     // A stabiliser state

typedef struct {
  mzd_t *factor_transpose;  // factor: 1 x n matrix (called transpose to make
                            // clear its a row vector)
  mzd_t *set_generators;    // k x n matrix (each row a generator)
} pauli_z_coset_t;
// Useful macros

#define DCP_N(t) ((t)->F->nrows)
#define STATE_N(t) ((((t)->x_parts->nrows)))
#define COSET_K(t) ((((t)->set_generators->nrows)))
#define COSET_N(t) ((((t)->set_generators->ncols)))

// convert.c
void convert_permutation(mzd_t *permutation, uint32_t n_z, uint32_t n,
                         uint32_t *rrefs);
void convert_dcp(dcp_t *dcp, comb_dcp_t *value);

// dcp.c
dcp_t *dcp_mul(dcp_t *result, dcp_t *c2, dcp_t *c1);
void dcp_mulf_right(dcp_t* c2, dcp_t* c1); // Computes c2 * c1 and stores the result in c2, effectively multiplying c2 on the right by c1.
void dcp_mulf_left(dcp_t* c2, dcp_t* c1); // Computes c2 * c1 and stores the result in c1, effectively multiplying c1 on the left by c2.
void set_stabiliser_state_to_dcp(stabiliser_state_t *state, dcp_t *dcp);
int postselect_on_measuring_plus_one_x(stabiliser_state_t *state,
                                       mzd_t *x_measurement_matrix);
int postselect_on_measuring_plus_one_z(stabiliser_state_t *state,
                                       mzd_t *x_measurement_matrix);
int _postselect_on_x(stabiliser_state_t *state, mzd_t *x_parts, mzd_t *z_parts,
                     mzd_t *x_measurement_matrix, mzd_t *delta_bits,
                     mzd_t *epsilon_bits, char *x_measurement_string);
void set_coset_to_stabiliser_state(pauli_z_coset_t *coset,
                                   stabiliser_state_t *state);

// gaussian_elimination.c
void swap_rows(stabiliser_state_t *state, int i, int j);
void add_row(stabiliser_state_t *state, int sourcerow, int destrow);
void add_row_to(stabiliser_state_t *source_state, int sourcerow,
                stabiliser_state_t *dest_state, int destrow);
void set_row_to(stabiliser_state_t *source_state, int sourcerow,
                stabiliser_state_t *dest_state, int destrow);
int form_identity_ge(stabiliser_state_t *state);
int verify_coset_ge(int r, stabiliser_state_t *state);

// init_free.c
dcp_t *init_dcp(uint32_t n);
stabiliser_state_t *init_state(int n);
void init_dcp_to_size(dcp_t *tableau, uint32_t n);
pauli_z_coset_t *init_pauli_z_coset(uint32_t n, uint32_t k);
void free_dcp_array(uint32_t count, dcp_t *array);
void free_dcp(dcp_t *dcp);
void free_pauli_z_coset(pauli_z_coset_t *coset);
void set_dcp_to(dcp_t *a, dcp_t *b);

// iteration.c
int incf_bit_vector(comb_bit_array_t v, int n);
int ft_iteration_combination_incf(int k, uint32_t *sizes, int n,
                                  int *first_balls, int *last_balls,
                                  int *current_combination);

// print.c
int do_print_pauli(mzd_t *x_parts, mzd_t *z_parts, mzd_t *delta_vector,
                   mzd_t *epsilon_vector, int row_number, char print_newline_p,
                   int counter);
char *print_row_vector_as_pauli(mzd_t *row_vector, int row,
                                const char *character);
char *print_dcp(dcp_t *dcp, char print_newline_p);
char *print_stabiliser_state(stabiliser_state_t *state, char print_newline_p);
char *print_pauli_z_coset(pauli_z_coset_t *coset, char print_newline_p);

// utils.c
void fix_dcp(dcp_t *dcp);
BIT row_product(mzd_t *x, int i, mzd_t *z, int j);
void expand_state(stabiliser_state_t *state);
void contract_state(stabiliser_state_t *state);
void clear_state(stabiliser_state_t *state);
void set_row(mzd_t *x_parts, mzd_t *z_parts, int j, mzd_t *x_measurement_matrix,
             int i, mzd_t *epsilon_bits, mzd_t *delta_bits);
void permute_columns(mzd_t *data, unsigned mask);
void convert_row(mzd_t *M, int row, comb_vector_t *vector);
void add_and_column_vector(mzd_t *result, mzd_t *F1, mzd_t *F2);

// inline functions

inline word mzd_data(const mzd_t *mzd) {  // assumes mzd is a row vector
  return mzd->data[0];
}

inline uint32_t gray_bit_to_switch(uint32_t n) {
  uint32_t n_add_one = n + 1;
  uint32_t diff = n ^ n_add_one ^ (n >> 1) ^ (n_add_one >> 1);
  return n ? 31 - __builtin_clz(diff) : 0;
}

#endif  // COMB_H
