// Putative new interface for libcomb

#ifdef __cplusplus
extern "C" {
#endif

// need fixed width integers
#include <stdint.h>
#include <stdlib.h>
#include <m4ri/m4ri.h>
// Error codes

  typedef enum {ALL_FAULTS_GOOD = 0,
                BAD_FAULT_FOUND = 1,
                INVALID_INPUT = 2,
                PROBLEM_TOO_LARGE = 3,
                INTERNAL_ERROR = 4,
                OTHER_ERROR = 5} comb_code_t;

  // A bit vector is an array of words

  typedef uint64_t comb_bit_word_t;
  typedef comb_bit_word_t* comb_bit_array_t; // maybe not needed?


#define MAX_SIZE_COUNT 1024

#define SetBit(A, k) (A[(k) / sizeof(comb_bit_word_t)] |= (1U << ((k) % sizeof(comb_bit_word_t))))
#define ClearBit(A, k) (A[(k) / sizeof(comb_bit_word_t)] &= ~(1U << ((k) % sizeof(comb_bit_word_t))))
#define TestBit(A, k) (A[(k) / sizeof(comb_bit_word_t)] & (1U << ((k) % sizeof(comb_bit_word_t))))

  // External interface

  typedef struct {
    size_t pauli_length; // number of qubits in the pauli
    // convention is: bits 0 - n-1 z parts
    comb_bit_array_t data; // bit array. must have enough space to store pauli_length bits
  } comb_pauli_z_t;

  typedef struct {
    size_t pauli_length; // number of qubits in the pauli
    // convention is: bits 0 - n-1 x parts
    comb_bit_array_t data; // bit array. must have enough space to store pauli_length bits
  } comb_pauli_x_t;

  typedef struct {
    size_t pauli_length; // number of qubits in the pauli
    // convention is: bits 0 - n-1 z parts, n - 2n-1 x parts, 2n is the `phase'` bit, 2n+1 sign bit
    comb_bit_array_t x_data; // bit array. must have enough space to store pauli_length bits
    comb_bit_array_t z_data; // bit array. must have enough space to store pauli_length bits
    char sign_bit; // the sign and phase bits of the pauli... potentially unused?
  } comb_pauli_t;

  typedef struct {
    uint32_t term_count;
    comb_pauli_z_t* terms;
    int32_t* coefficients;
  } comb_pauli_z_sum_t; // A sum of Z paulis.

  typedef struct {
    comb_pauli_z_sum_t exponent_sum;
    comb_pauli_x_t pauli_x_part;
  } comb_dcp_t; // A DCP operator (represented as an exponential of a Pauli sum)

  typedef struct {
    uint32_t n;
    uint32_t* indices;
  } comb_permutation_t;
  // If you want, you can supply a function to be called on each pauli after canonicalisation

  typedef int (*comb_ft_iteration_xz_evaluator)(comb_pauli_t* output_error_combination);
  typedef int (*comb_ft_iteration_x_evaluator)(comb_pauli_x_t* output_error_combination_x);
  typedef int (*comb_ft_iteration_z_evaluator)(comb_pauli_z_t* output_error_combination_z);

  // Specifies magic states of CSS codes prepared by the CNOT + T gate set.
  // It includes lookup tables to check the weight of operators
  typedef struct {
    size_t n; // number of qubits
    size_t n_x; // number of independent x stabiliser generators
    size_t n_z; // number of independent z stabiliser generators
    size_t k; // number of dcp stabilisers
    comb_pauli_x_t* x_stabilisers; // n-x of these. they are also in rref; after a permutation, they are in the appropriate diagonal form in the final n-x entries
    size_t* x_permutation; // permutation to apply to qubits to get into canonical form needed for x weight lookup
    comb_bit_array_t x_lookup_table; // z weight lookup table. of size 2^(n - n_x)
    comb_pauli_z_t* z_stabilisers; // n-z of these. they are also in rref; after a permutation, they are in the appropriate diagonal form in the final n-z entries
    size_t* z_permutation;  // permutation to apply to qubits to get into canonical form needed for z weight lookup
    comb_bit_array_t z_lookup_table; // z weight lookup table. of size 2^(n - n_z)
    comb_dcp_t* dcp_stabilisers; // k of these
    comb_ft_iteration_x_evaluator x_evaluator; // optional. if not NULL, it is called instead of dealing with lookup tables.
    comb_ft_iteration_z_evaluator z_evaluator; // optional. if not NULL, it is called instead of dealing with lookup tables.
  } comb_css_magic_state_t;

  // for pauli eigenstates of css codes. Assumed to be of X type
  typedef struct {
    uint32_t n; // number of qubits in the code
    uint32_t n_x; // number of x stabilisers in the code
    comb_pauli_x_t* stabilisers; // n_x of these, in rref
    comb_bit_array_t weight_lookup_table; // of size 2^(n-n_x)
  } comb_css_pauli_eigenstate_t; // one half of a css code (assumed to be x)

  typedef struct {
    uint32_t n; // number of qubits in the code
    comb_pauli_t* stabilisers; // n of these, in a rref canonical form
    comb_bit_array_t weight_lookup_table;
    comb_ft_iteration_xz_evaluator xz_evaluator;
  } comb_pauli_eigenstate_t; // a general stabiliser code.

  // Functions

  /* Verifies the fault tolerance of a circuit that prepares a magic state of a CSS code with a CNOT + T circuit by checking that all possible combinations of faults are detected or sufficiently low weight
     `n' is the number of fault locations.
     `m' is a pointer to an array of size at least n, and m[i] is the number of faults in the i-th location.
     `k' is the number of faults to which the circuit should be resilient to.
     `bad_fault_array` is a pointer to an array of at least k elements, to which the indices of any bad faults will be written if any are found.
     `state' is the magic state being prepared.
     `output_error_count` is the number of propagated output errors. It must be at least m[0]+...m[s-1] else an error is returned.
     `output_errors` is a pointer to at least `output_error_count` number of errors.
     `x_measurement_count' is the number of independent postselections made on ancilla measurements of purely X type
     `x_measurement_paulis' is a pointer to at least `x_measurement_count` measurement paulis.
     `z_measurement_count' is the number of independent postselections made on ancilla measurements of purely Z type
     `z_measurement_paulis' is a pointer to at least `z_measurement_count` measurement paulis.
     The x_measurement_paulis and z_measurement_paulis must have entirely disjoint support, else an error is returned.
  */
  comb_code_t comb_verify_css_magic_state(int n, int *m, int k,
                                          size_t *bad_fault_array,
                                          comb_css_magic_state_t* state,
                                          size_t output_error_count, comb_dcp_t* output_errors,
                                          size_t x_measurement_count, comb_pauli_x_t* x_measurement_paulis,
                                          size_t z_measurement_count, comb_pauli_z_t* z_measurement_paulis);

  comb_code_t comb_verify_css_pauli_eigenstate(int s, int *m, int k,
                                               size_t *bad_fault_array,
                                               comb_css_pauli_eigenstate_t* state,
                                               size_t output_error_count, comb_pauli_z_t* output_errors,
                                               size_t x_measurement_count, comb_pauli_x_t* x_measurement_paulis);

  comb_code_t comb_verify_pauli_eigenstate(int s, int *m, int k,
                                           size_t *bad_fault_array,
                                           comb_pauli_eigenstate_t* state,
                                           size_t output_error_count, comb_pauli_t* output_errors,
                                           size_t x_measurement_count, comb_pauli_x_t* x_measurement_paulis,
                                           size_t z_measurement_count, comb_pauli_z_t* z_measurement_paulis);

  // Some legacy methods here.. they are faster for some things



  // idea. do rref in c? maybe a todo later. people should do it themselves...
  // it goes against the philosophy of the library

  // Internal definitions - not necessary to expose these in CFFI

  // use m4rie for matrix multiplication
  // no need to implement addition / multiplication etc.
  typedef struct { // A Clifford tableau following Dehaene and Demoor.
    mzd_t* c_bar; // 2n + 1 square matrix
    mzd_t* h_bar; // 2n + 1 column vector
  } clifford_t;

  typedef struct { // Stores a DCP operator as a tableau.
    mzd_t* F; // the upper right block element of c (not c-bar). n x n matrix
    mzd_t* h_x; // signs of images of z_i
    mzd_t* h_z;
  } dcp_t;

  typedef struct {
    char sign_bit; // sign bit
    mzd_t* data; // 2n + 1 column vector
  } pauli_t; // A Pauli following Dehaene and Demoor.

  typedef struct {
    int n;
    pauli_t** generators;
  } stabiliser_state_t; // A stabiliser state

  typedef struct {
    pauli_t* factor; // factor
    size_t k;
    pauli_t** set_generators; // array of k generators
  } pauli_z_coset_t;
  // Useful macros

  uint32_t gray_bit_to_switch(uint32_t n);
#define DCP_N(t) ((t)->F->nrows)
#define PAULI_N(t) ((((t)->data->nrows)-1)/2)
#define STABILISER_N(t) ((((t)->data->nrows)-1)/2)

  // Creation
  clifford_t* init_clifford(size_t n);
  dcp_t* init_dcp(size_t n);
  pauli_t* init_pauli(size_t n);
  pauli_z_coset_t* init_pauli_z_coset(size_t n, size_t k);
  dcp_t* dcp_pauli_x(size_t i, size_t n);
  dcp_t* dcp_pauli_z(size_t i, size_t n);
  dcp_t* dcp_cz(size_t i, size_t j,size_t n);
  dcp_t* dcp_s(size_t i, size_t n);
  // Freeing
  void free_clifford(clifford_t* clifford);
  void free_pauli(pauli_t* pauli);
  void free_dcp(dcp_t* dcp);
  void free_pauli_z_coset(pauli_z_coset_t* coset);
  // Setting
  // ? what to put here. they are structs so not that much point?
  void set_dcp_to(dcp_t *a, dcp_t* b); // assigns b to a: a = b.
  // Printing
  void print_dcp(dcp_t *dcp, char print_newline_p);
  void print_pauli(pauli_t *pauli, char print_newline_p);
  void print_pauli_z_coset(pauli_z_coset_t *coset, char print_newline_p);
  void print_stabiliser_state(stabiliser_state_t *coset, char print_newline_p);
  // more convenient creation
  pauli_t *pauli_x(size_t i, size_t n);
  pauli_t *pauli_z(size_t i, size_t n);

  mzd_t* diag(mzd_t* matrix);
  void add_and(mzd_t* result, mzd_t* F1,mzd_t* F2);

  // Manipulating
  int pauli_commute_p(pauli_t *p1, pauli_t *p2);
  // Converting
  void convert_dcp(dcp_t* dcp, comb_dcp_t* comb);
  // Multiplies Paulis
  pauli_t *pauli_mul(pauli_t *result, pauli_t *p1, pauli_t *p2);
  // Conjugates Clifford via Pauli
  pauli_t *dcp_conjugate_pauli(pauli_t *result, dcp_t *clifford, pauli_t *pauli);
  // for algorithm:
  // 1. Multiplies DCP operators
  dcp_t *dcp_mul(dcp_t *result, dcp_t *c1, dcp_t *c2);

  pauli_t* dcp_image_x(dcp_t *dcp, size_t i);
  // 2. converts a diagonal clifford to a stabiliser state (operator state duality)
  stabiliser_state_t* dcp_to_stabiliser_state(dcp_t* dcp);

  // 3. Postselects on a stabiliser state (destructive). Can be used to canonicalise too...
  int postselect_on_measuring_plus_one(stabiliser_state_t* state, pauli_t* measurement_operator);

  // 4. Applies the diagonal clifford structure theorem to express a clifford as a linear combination of a set of paulis.
  pauli_z_coset_t *diagonal_clifford_to_coset(stabiliser_state_t *state);

  // Postselects the diagonal clifford on measuring p to be +1. returns non-zero if the result vanishes.
  int postselect_diagonal_clifford(dcp_t* result, dcp_t* clifford, pauli_t *p);
  /* // Given a diagonal Clifford C, returns C|+> */
  /* stabiliser_state_t *clifford_to_stabiliser_state(clifford_t *clifford); */
  /* // Measures a Pauli on a stabiliser state. Returns non-zero if state vanishes. destructive */
  /* int measure_state(stabiliser_state_t* state, pauli_t *measurement_pauli); // destructive! */
  /* // Turns a stabiliser state corresponding to a diagonal Clifford to a coset. */
  /* pauli_z_coset_t *stabiliser_state_to_coset(stabiliser_state_t *state); */

  int postselect_pauli_z_coset(pauli_z_coset_t* result, pauli_z_coset_t* coset, int k, pauli_t **measurement_operators);

  // iteration

  int ft_iteration_combination_incf(int k, int *sizes, int n, int *first_balls,
                                    int *last_balls, int *current_combination);

  void init_dcp_to_size(dcp_t* tableau, size_t n);

  int dcp_fault_bad_p_css(dcp_t *error, comb_css_magic_state_t* state, pauli_t* p, int count);

  uint32_t mzd_data(mzd_t* mzd);
  void  permute_qubits(pauli_z_coset_t* coset,int n, size_t* c);
  void set_to_pauli(pauli_t* a, pauli_t* b);
  void init_pauli_to_size(pauli_t* a, int size);
  void set_dcp_to_transvection(dcp_t* dcp, pauli_t* z_pauli);
  void convert_pauli(pauli_t* temp_pauli,comb_pauli_z_t* pauli);
#ifdef __cplusplus
}
#endif
