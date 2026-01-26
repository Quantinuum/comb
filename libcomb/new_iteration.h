// Putative new interface for libcomb

#ifdef __cplusplus
extern "C" {
#endif

// need fixed width integers
#include <stdint.h>

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
     `s' is the number of fault locations.
     `m' is a pointer to an array of size at least s, and m[i] is the number of faults in the i-th location.
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
  comb_code_t comb_verify_css_magic_state(int s, int *m, int k,
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

#ifdef __cplusplus
}
#endif
