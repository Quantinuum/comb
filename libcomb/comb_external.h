// External interface
#include <stdint.h>
#include <stdlib.h>
// Error codes

typedef enum {ALL_FAULTS_GOOD = 0,
              BAD_FAULT_FOUND = 1,
              INVALID_INPUT = 2,
              PROBLEM_TOO_LARGE = 3,
              INTERNAL_ERROR = 4,
              OTHER_ERROR = 5} comb_code_t;

#define SetBit(A, k) (A[(k) / sizeof(comb_bit_word_t)] |= (1U << ((k) % sizeof(comb_bit_word_t))))
#define ClearBit(A, k) (A[(k) / sizeof(comb_bit_word_t)] &= ~(1U << ((k) % sizeof(comb_bit_word_t))))
#define TestBit(A, k) (A[(k) / sizeof(comb_bit_word_t)] & (1U << ((k) % sizeof(comb_bit_word_t))))
// A bit vector is an array of words

typedef uint64_t comb_bit_word_t;
typedef comb_bit_word_t* comb_bit_array_t; // maybe not needed?

//#define MAX_BIT_VECTOR_LENGTH 64
//#define MAX_BIT_VECTOR_WORD_LENGTH (MAX_BIT_VECTOR_LENGTH / sizeof(comb_bit_word_t))
//#define MAX_BIT_MATRIX_LENGTH (MAX_BIT_VECTOR_LENGTH * MAX_BIT_VECTOR_LENGTH )
//#define MAX_BIT_MATRIX_WORD_LENGTH (MAX_BIT_VECTOR_LENGTH / sizeof(comb_bit_word_t))

typedef struct { // n x 1 vector
  uint32_t n;
  comb_bit_array_t data;
} comb_vector_t;

typedef struct { // m x n matrix
  uint32_t m;
  uint32_t n;
  comb_bit_array_t data;
} comb_matrix_t;

typedef struct {
  comb_matrix_t F;
  comb_vector_t h_x;
  comb_vector_t h_z;
} comb_dcp_t;

typedef struct {
  uint32_t n;
  uint32_t* indices;
} comb_permutation_t;

typedef comb_vector_t comb_pauli_x_t;
typedef comb_vector_t comb_pauli_z_t;

comb_bit_array_t comb_alloc_bit_array(int n);
comb_vector_t *comb_alloc_vector(int n);
comb_matrix_t *comb_alloc_matrix(int m, int n);
comb_dcp_t *comb_alloc_dcp(int n);
comb_permutation_t *comb_alloc_permutation(int n);

void comb_set_vector_aref(comb_vector_t* vector, int i, int value);
void comb_set_matrix_aref(comb_matrix_t* matrix, int i, int j, int value);

int comb_vector_aref(comb_vector_t* vector, int i);
int comb_matrix_aref(comb_matrix_t* matrix, int i, int j);

void comb_free_bit_array(comb_bit_array_t pointer);
void comb_free_vector(comb_vector_t* pointer);
void comb_free_matrix(comb_matrix_t* pointer);
void comb_free_dcp(comb_dcp_t* pointer);
void comb_free_permutation(comb_permutation_t pointer);

comb_code_t comb_verify_css_magic_state( // describing the magic state
                                           uint32_t n_all, // number of qubits in code + all ancillas
                                           uint32_t n, // number of qubits in code
                                           uint32_t n_x, // number of independent x stabiliser generators in code
                                           uint32_t n_z, // number of independent z stabiliser generators in code
                                           uint32_t k, // number of logical qubits in code. equivalently, number of logical dcp operators that stabilise the code
                                           comb_pauli_x_t* x_stabilisers, // n-x of these
                                           comb_permutation_t* x_permutation, //permutation that puts x_stabilisers in rref form
                                           comb_bit_array_t x_lookup_table, // x lookup table. of size 2^(n-n_x)
                                           comb_pauli_z_t* z_stabilisers,
                                           comb_permutation_t* z_permutation, // permutation that puts z_stabilisers in rref form
                                           comb_bit_array_t z_lookup_table,
                                           comb_dcp_t* dcp_stabilisers, // k of these

                                           // other parameters to the problem
                                           int fault_location_count, // number of fault locations
                                           int *fault_location_partitions, // at least fault_location_count of these. number of faults at each location
                                           int s, // number of faults to which circuit should be resilient to
                                           uint32_t *bad_fault_array, // where to store indices of bad faults. must have enough space for s elements.
                                           uint32_t output_error_count, // number of output errors that result from propagating all the faults. should be at least sum(fault_location_partitions)
                                           comb_dcp_t* output_errors, // at least output_error_count of these
                                           uint32_t x_measurement_count, // number of x measurement paulis
                                           comb_pauli_x_t* x_measurement_paulis, // postselections made in x basis. at least x_measurement_count of these
                                           uint32_t z_measurement_count, // number of z measurement paulis
                                           comb_pauli_z_t* z_measurement_paulis); // postselections made in z basis
