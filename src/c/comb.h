#include <m4ri/m4ri.h>

#include "comb_external.h"
#include "lwlog.h"

#include <math.h>
#include <stdio.h>
#include <time.h>
#include <omp.h>
#ifdef COMB_MPI
#include <mpi.h>
#endif

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
#define min(a,b)                                \
  ({ __typeof__ (a) _a = (a);                   \
    __typeof__ (b) _b = (b);                    \
    _a < _b ? _a : _b; })

typedef struct {
  unsigned long long faults_already_checked; // faults already checked.
  int* current_combination; // current fault combination.
  int* final_combination_to_check; // final fault combination to check. shouldn't go beyond this
  uint32_t s;// number of faults in the combination.
} chunk_t; // A chunk of a fault tolerance computation
// chunk.c
void init_chunk(chunk_t* chunk, unsigned long long faults_to_check_count, int* starting_combination, int s);
int chunk_finished_p(chunk_t* chunk);
chunk_t* partition_into_chunks(uint32_t fault_location_count, uint32_t *fault_location_partitions, uint32_t s, size_t minimum_chunk_size, size_t* chunk_count);
void describe_chunks(chunk_t* chunks, size_t chunk_count);

void init_first_and_last_balls(int* first_balls, int* last_balls, uint32_t fault_location_count, uint32_t *fault_location_partitions);
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

// fault_checking.c
int dcp_fault_bad_p_css(dcp_t *error, mzd_t* C_z, mzd_t *C_z_transpose, mzd_t *C_x, mzd_t *C_x_transpose, mzd_t *P_x,
    mzd_t *z_measurements_for_x,
    unsigned x_permutation,  // permutation mask for x canonicalisation
    unsigned z_permutation,  // permutation mask for z canonicalisation
    comb_bit_array_t x_lookup_table, comb_bit_array_t z_lookup_table,
    mzd_t *x_stabilisers, uint32_t *rref_x, int n_all, int n, int n_x, int n_z,
    stabiliser_state_t *state, pauli_z_coset_t *coset, mzd_t *D, mzp_t *Q,
                        int row_ops_count, int *row_ops_a, int *row_ops_b);

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

unsigned long long combi(int n,int k);
void unrank(int s, int* combination, unsigned long long n);
unsigned long long rank(int s, int* combination);
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

// serialise.c

// Serialise the run to output_file so it can be subsequently run later by e.g. a MPI implementation.
int serialise_run_to_file(int fd,    uint32_t n_all,  ///< number of qubits in code + all ancillas
    uint32_t n,      // number of qubits in code
    uint32_t n_x,    // number of independent x stabiliser generators in code
    uint32_t n_z,    // number of independent z stabiliser generators in code
    uint32_t k,  // number of logical qubits in code. equivalently, number of
                 // logical dcp operators that stabilise the code
    comb_pauli_x_t *_x_stabilisers,   // n-x of these
    comb_bit_array_t x_lookup_table,  // x lookup table. of size 2^(n-n_x)
    comb_pauli_z_t *_z_stabilisers, comb_bit_array_t z_lookup_table,
    comb_dcp_t *_dcp_stabilisers,   // k of these
    uint32_t fault_location_count,  // number of fault locations
    uint32_t
        *fault_location_partitions,  // at least fault_location_count of these.
                                     // number of faults at each location
    uint32_t s,  // number of faults to which circuit should be resilient to
    uint32_t output_error_count,   // number of output errors that result from
                                   // propagating all the faults. should be at
                                   // least sum(fault_location_partitions)
    comb_dcp_t *_output_errors,    // at least output_error_count of these
    uint32_t x_measurement_count,  // number of x measurement paulis
    comb_pauli_x_t
        *_x_measurement_paulis,    // postselections made in x basis. at
                                   // least x_measurement_count of these
    uint32_t z_measurement_count,  // number of z measurement paulis
    comb_pauli_z_t *_z_measurement_paulis);

// Opens a previously serialised run from a file and checks it using comb_verify_css_magic_state
int deserialise_run_from_file(char *input_file);

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
int stabilisers_rref_p(comb_vector_t *operators, uint32_t n, uint32_t m,
                       uint32_t *coefficients, unsigned *permutation_mask);
// inline functions

inline word mzd_data(const mzd_t *mzd) {  // assumes mzd is a row vector
  return mzd->data[0];
}

inline uint32_t gray_bit_to_switch(uint32_t n) {
  uint32_t n_add_one = n + 1;
  uint32_t diff = n ^ n_add_one ^ (n >> 1) ^ (n_add_one >> 1);
  return n ? 31 - __builtin_clz(diff) : 0;
}
#ifdef COMB_MPI

#define REQUEST_CHUNK_TAG 1
#define RECEIVE_CHUNK_TAG 2
#define NOTIFY_CHUNK_GOOD 3

int mpi_get_chunk_to_process(int world_rank);
void mpi_notify_primary_chunk_good(int chunk_num);
size_t mpi_get_free_chunk(comb_bit_array_t free_chunks_index_array, int chunk_count);
void mpi_primary(int chunk_count);
void mpi_secondary(int* chunks_checked, unsigned long long* fault_counter, uint32_t n_all,  ///< number of qubits in code + all ancillas
    uint32_t n,      // number of qubits in code
    uint32_t n_x,    // number of independent x stabiliser generators in code
    uint32_t n_z,    // number of independent z stabiliser generators in code
    uint32_t k,  // number of logical qubits in code. equivalently, number of
                 // logical dcp operators that stabilise the code
    mzd_t *x_stabilisers,             // n-x of these
    comb_bit_array_t x_lookup_table,  // x lookup table. of size 2^(n-n_x)
    mzd_t *z_stabilisers, comb_bit_array_t z_lookup_table,
    dcp_t *dcp_stabilisers,         // k of these
    uint32_t fault_location_count,  // number of fault locations
    uint32_t
        *fault_location_partitions,  // at least fault_location_count of these.
                                     // number of faults at each location
    uint32_t s,  // number of faults to which circuit should be resilient to
    uint32_t *bad_fault_array,     // where to store indices of bad_faults
    uint32_t output_error_count,   // number of output errors that result from
                                   // propagating all the faults. should be at
                                   // least sum(fault_location_partitions)
    dcp_t *output_errors,          // at least output_error_count of these
    uint32_t x_measurement_count,  // number of x measurement paulis
    mzd_t *x_measurement_paulis,   // postselections made in x basis. at least
                                   // x_measurement_count of these
    uint32_t z_measurement_count,  // number of z measurement paulis
    mzd_t *z_measurement_paulis, uint32_t *rref_x, uint32_t *rref_z,
                   unsigned x_permutation, unsigned z_permutation, chunk_t* chunks);
#endif

comb_code_t __do_verify_css_magic_state(
    uint32_t n_all,  ///< number of qubits in code + all ancillas
    uint32_t n,      // number of qubits in code
    uint32_t n_x,    // number of independent x stabiliser generators in code
    uint32_t n_z,    // number of independent z stabiliser generators in code
    uint32_t k,  // number of logical qubits in code. equivalently, number of
                 // logical dcp operators that stabilise the code
    mzd_t *x_stabilisers,             // n-x of these
    comb_bit_array_t x_lookup_table,  // x lookup table. of size 2^(n-n_x)
    mzd_t *z_stabilisers, comb_bit_array_t z_lookup_table,
    dcp_t *dcp_stabilisers,         // k of these
    uint32_t fault_location_count,  // number of fault locations
    uint32_t
        *fault_location_partitions,  // at least fault_location_count of these.
                                     // number of faults at each location
    uint32_t s,  // number of faults to which circuit should be resilient to
    uint32_t *bad_fault_array,     // where to store indices of bad_faults
    uint32_t output_error_count,   // number of output errors that result from
                                   // propagating all the faults. should be at
                                   // least sum(fault_location_partitions)
    dcp_t *output_errors,          // at least output_error_count of these
    uint32_t x_measurement_count,  // number of x measurement paulis
    mzd_t *x_measurement_paulis,   // postselections made in x basis. at least
                                   // x_measurement_count of these
    uint32_t z_measurement_count,  // number of z measurement paulis
    mzd_t *z_measurement_paulis, uint32_t *rref_x, uint32_t *rref_z,
    unsigned x_permutation, unsigned z_permutation,
    chunk_t* chunk);
