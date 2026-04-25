#include "comb_external.h"
#include <m4ri/m4ri.h>
#include "lwlog.h"
#define MAX_SIZE_COUNT 1024
// if TIMING_MODE is set, then repeat until all faults good
// #define TIMING_MODE 1000000
#ifdef TIMING_MODE
#define BAD_FAULT_FOUND ALL_FAULTS_GOOD
#endif

// Conventions:
// vectors are stored as column matrices by default (n x 1) matrices.
// However, in certain cases, it is much more preferable to store them as row vectors. In such cases, they are /always/ labelled with '_transpose' to avoid confusion, even in cases where they are the only form of storage.
// Matrices will considered to have their rows consisting of Paulis (or parts of Paulis).

typedef struct { // Stores a DCP operator as a tableau.
  mzd_t* F; // the upper right block element of c (not c-bar). n x n matrix
  mzd_t* h_x; // n x 1 matrix
  mzd_t* h_x_transpose; // 1 x n matrix (row vector).
  mzd_t* h_z; // n x 1 matrix
  mzd_t* F_diagonal; // diagonal part of F. n x 1 matrix.
} dcp_t;

typedef struct {
  mzd_t* x_parts; // n x n - rows are the x parts of paulis
  mzd_t* z_parts; // n x n - rows are the z parts of paulis
  mzd_t* delta_vector; // n x 1 - entries are the imaginary bits of the vectors
  // idea: delta_vector is determined by x_parts and z_parts as it must always be hermitian...
  mzd_t* epsilon_vector; // n x 1 - entries are the sign bits of the vectors.
} stabiliser_state_t; // A stabiliser state

typedef struct {
  mzd_t* factor_transpose; // factor: 1 x n matrix (called transpose to make clear its a row vector)
  mzd_t* set_generators; // k x n matrix (each row a generator)
} pauli_z_coset_t;
// Useful macros

#define DCP_N(t) ((t)->F->nrows)
#define STATE_N(t) ((((t)->x_parts->nrows)))
#define COSET_K(t) ((((t)->set_generators->nrows)))
#define COSET_N(t) ((((t)->set_generators->ncols)))

inline uint32_t gray_bit_to_switch(uint32_t n)
{
  uint32_t n_add_one = n + 1;
  uint32_t diff = n ^ n_add_one ^ (n >> 1) ^ (n_add_one >> 1);
  return n ? 31 - __builtin_clz(diff) : 0;
}

// Creation
dcp_t* init_dcp(uint32_t n);
pauli_z_coset_t* init_pauli_z_coset(uint32_t n, uint32_t k);
// Freeing
void free_dcp(dcp_t* dcp);
void free_pauli_z_coset(pauli_z_coset_t* coset);
// Setting
// ? what to put here. they are structs so not that much point?
void set_dcp_to(dcp_t *a, dcp_t* b); // assigns b to a: a = b.
// Printing
char* print_dcp(dcp_t *dcp, char print_newline_p);
char* print_pauli_z_coset(pauli_z_coset_t *coset, char print_newline_p);
char* print_stabiliser_state(stabiliser_state_t *coset, char print_newline_p);
// more convenient creation
void add_and_column_vector(mzd_t* result, mzd_t* F1,mzd_t* F2);
// Converting
void convert_dcp(dcp_t* dcp, comb_dcp_t* comb);

char* print_row_vector_as_pauli(mzd_t* row_vector, int row, const char* character);
// for algorithm:
// 1. Multiplies DCP operators
dcp_t *dcp_mul(dcp_t *result, dcp_t *c1, dcp_t *c2);
// 2. Sets a stabiliser state to the diagonal part of a DCP
void set_stabiliser_state_to_dcp(stabiliser_state_t* state, dcp_t* dcp);
// 3. Manipulate stabiliser states more easily
void swap_rows(stabiliser_state_t* state, int i, int j);
void add_row(stabiliser_state_t* state, int sourcerow, int destrow);
void left_multiply(stabiliser_state_t* state, mzd_t* M, stabiliser_state_t* new_state);

// 3. Postselects on a stabiliser state (destructive). Can be used to canonicalise too...
int postselect_on_measuring_plus_one_x(stabiliser_state_t* state, mzd_t* x_measurement_matrix);
int postselect_on_measuring_plus_one_z(stabiliser_state_t* state, mzd_t* z_measurement_matrix);

// 4. Applies the diagonal clifford structure theorem to express a clifford as a linear combination of a set of paulis.
void set_coset_to_stabiliser_state(pauli_z_coset_t *coset, stabiliser_state_t *state);
// iteration

int ft_iteration_combination_incf(int k, uint32_t *sizes, int n, int *first_balls,
                                  int *last_balls, int *current_combination);

void init_dcp_to_size(dcp_t* tableau, uint32_t n);

inline word mzd_data(const mzd_t* mzd){ // assumes mzd is a row vector
  return mzd->data[0];
}
void  permute_qubits(pauli_z_coset_t* coset,int n, mzd_t* c, mzd_t* temp, mzd_t* big_c);
// conversion
void convert_dcp(dcp_t* dcp, comb_dcp_t* value);
void convert_permutation(mzd_t* permutation, uint32_t n, uint32_t n_z, uint32_t *rrefs);
void convert_row(mzd_t* M, int row, comb_vector_t* vector);

void free_dcp_array(uint32_t count, dcp_t* array);

BIT row_product(mzd_t* x, int i, mzd_t* z, int j);

void fix_dcp(dcp_t* dcp);

stabiliser_state_t* init_state(int n);

void expand_state(stabiliser_state_t* state);

void set_row(mzd_t* x_parts, mzd_t* z_parts, int j, mzd_t* x_measurement_matrix, int i,mzd_t* epsilon_bits, mzd_t* delta_bits);
void contract_state(stabiliser_state_t* state);
unsigned compress4(unsigned x, unsigned m);
void permute_columns(mzd_t* data, unsigned mask);

int form_identity_ge(stabiliser_state_t* state);
int verify_coset_ge(int r, stabiliser_state_t* state);
