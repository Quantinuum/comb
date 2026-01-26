#ifdef __cplusplus
extern "C" {
#endif


// need fixed width integers
#include <stdint.h>

// error codes

#define ALL_FAULTS_GOOD 0
#define BAD_FAULT_FOUND 1
#define PROBLEM_TOO_LARGE 2
#define INTERNAL_ERROR 3
#define OTHER_ERROR 4

  // maximum sizes

#define MAX_PAULI_SIZE 4
#define MAX_TERM_COUNT 128
#define MAX_SIZE_COUNT 1024


// A bit vector is an array of uint32_t. Here are some convenience methods for accessing them

#define SetBit(A, k) (A[(k) / 32U] |= (1U << ((k) % 32U)))
#define ClearBit(A, k) (A[(k) / 32U] &= ~(1U << ((k) % 32U)))
#define TestBit(A, k) (A[(k) / 32U] & (1U << ((k) % 32U)))

  // common data types

  // A pauli of z type
  typedef struct {
    uint32_t pauli_length; // number of qubits in the pauli
    // convention is: bits 0 - n-1 z parts
    uint32_t data[MAX_PAULI_SIZE]; // bit array
  } ctc_pauli_z_t;

  // A pauli of x type
  typedef struct {
    uint32_t pauli_length; // number of qubits in the pauli
    // convention is: bits 0 - n-1 x parts
    uint32_t data[MAX_PAULI_SIZE]; // bit array of size n (uint32_t array of size n/32)
  } ctc_pauli_x_t;

  typedef struct {
    int32_t real;
  int32_t imaginary;
  } ctc_gaussian_integer_t;

  // A sum of z paulis
  typedef struct {
    uint32_t term_count;
    ctc_pauli_z_t* terms;
    ctc_gaussian_integer_t* coefficients;
  } ctc_pauli_z_sum_t; // problem: need complex coefficients (gaussian integers should suffice)

  typedef struct {
    uint32_t diagonal_operators_count; // number of diagonal operators in the exponent
    ctc_pauli_z_t diagonal_operators[MAX_TERM_COUNT]; // a pointer to an array of size diagonal_operators_count
    int32_t coefficients[MAX_TERM_COUNT]; // a pointer to an array of size diagonal_operators_count
    ctc_pauli_x_t pauli_x_part; // the x part of the DCP operator
  } ctc_dcp_t;



  // typedef struct {
    // a struct that contains enough information to canonicalise/compactify an operator in a css code.
    // weight finder just a list of paulis and which index to multiply them by.
    // or: could demand that all output errors are already canonicalised.. this could work
    //  } ctc_weight_finder;

  // The signature of a function that evaluates whether a given Pauli is good or not. the arguments are pointers to the bit arrays that make up the pauli.
  // Ideally this function should not take very long to execute (lookup tables etc.)
  typedef int (*ctc_ft_iteration_xz_evaluator)(ctc_pauli_x_t* output_error_combination_x_part_data,
                                               ctc_pauli_z_t* output_error_combination_z_part_data);
  typedef int (*ctc_ft_iteration_x_evaluator)(ctc_pauli_x_t* output_error_combination_x_part_data);
  typedef int (*ctc_ft_iteration_z_evaluator)(ctc_pauli_z_t* output_error_combination_z_part_data);
  /**
   Iterates over all combinations of exactly `k' faults, such that no two of
   them belong to the same fault location, computes the product of their output errors and verifies if
   the Pauli terms in the expansion of the combination are valid according to some user-supplied function `evaluator'. `s' names the number
   of fault locations. `m' is an array such that m[i] denotes the number of
   faults at the i-th location. `m' must have enough room for s elements. `k'
   denotes the number of faults that are combined to make a given fault
   combination. `propagated_output_errors' is an array such that propagated_output_errors[i] is
   the DCP error that results at the end of the circuit from propagating the i-th
   fault, where faults are labelled lexicographically according to their
   location in `m`. `propagated_output_errors' must have enough room for m[0] + ... +
   m[s-1] elements. `bad_fault_array' stores the indices of the faults that make
   up a bad combination if one is found. `bad_fault_array' must have enough room
   for k elements. `evaluator' is the function that is called to evaluate the goodness of a Pauli error in each
   output fault combination.
*/
  // TODO Come up with a way to specify logical Cliffords that stabilise the operator
  int ctc_ft_iteration_iterate(uint32_t s, uint32_t *m, uint32_t k, ctc_dcp_t* propagated_output_errors,
                               uint32_t *bad_fault_array,
                               uint32_t logical_clifford_count, ctc_dcp_t* logical_cliffords,
                               ctc_ft_iteration_xz_evaluator evaluator);


  int ctc_ft_iteration_iterate_css(uint32_t s, int *m, uint32_t k, uint32_t output_error_count,
                                   ctc_dcp_t* propagated_output_errors,
                                   uint32_t *bad_fault_array,
                                   uint32_t logical_clifford_count, ctc_dcp_t* logical_cliffords,
                                   ctc_ft_iteration_x_evaluator x_evaluator,
                                   ctc_ft_iteration_z_evaluator z_evaluator,
                                   uint32_t x_measurement_count, uint32_t* x_measurement_indices,
                                   uint32_t z_measurement_count, uint32_t* z_measurement_indices,
                                   int canonicaliser_count, ctc_pauli_x_t* canonicalisers);

  void add_term_to_sum(ctc_pauli_z_sum_t* z_sum,ctc_pauli_z_t* current_term, ctc_gaussian_integer_t coefficient,
                       uint32_t x_measurement_count, uint32_t* x_measurement_indices,
                       uint32_t z_measurement_count, uint32_t* z_measurement_indices);
  // extra definitions

  // A number of the form n + im, also called a gaussian integer.
  // Included here for completeness, but possibly not needed
  //typedef struct {
    //uint32_t imaginary_part;
    // uint32_t real_part;
    //} ctc_gaussian_integer_t;

  // General pauli (with phase and sign). Probably will not be used but put here for completeness.
  //typedef struct {
  //  uint32_t pauli_length; // number of qubits in the pauli. note that the bit vector actually has length 2n + 2
  //  // convention is (possibly subject to change!): bit 0 sign bit, bit 1 phase bit, bits 2 - n+1 x parts, bits n+2 - 2n+1 z parts
  //  uint32_t* data; // pointer to a bit array of size 2n+2 (uint32_t array of size 2n+2/32)
  //} ctc_pauli_t;

  // cffi them?

  // gray codes
  uint32_t integer_to_gray(uint32_t integer);
  uint32_t gray_to_integer(uint32_t integer);
  uint32_t gray_bit_to_switch(uint32_t n);



void ctc_print_pauli_x(ctc_pauli_x_t* operator);
void ctc_print_pauli_z(ctc_pauli_z_t* operator);
void ctc_print_dcp(ctc_dcp_t* error);
void ctc_print_pauli_sum(ctc_pauli_z_sum_t* sum);

  int ctc_pauli_anticommute_p(ctc_pauli_x_t* p1, ctc_pauli_z_t* p2);

void ctc_canonicalise_operator(ctc_dcp_t* operator);

void ctc_store_dcp_mul(ctc_dcp_t* p1, ctc_dcp_t* p2, ctc_dcp_t* result);
void ctc_store_pauli_x_mul(ctc_pauli_x_t* p1, ctc_pauli_x_t* p2, ctc_pauli_x_t* result);
void ctc_store_pauli_z_mul(ctc_pauli_z_t* p1, ctc_pauli_z_t* p2, ctc_pauli_z_t* result);

void ctc_postselect_x(ctc_dcp_t* operator);
void ctc_postselect_z(ctc_dcp_t* operator);

  void set_mul_complex(ctc_gaussian_integer_t* a,
                       ctc_gaussian_integer_t* b,
                       ctc_gaussian_integer_t* result);

  void set_coefficient_to(ctc_gaussian_integer_t* identity_part, ctc_gaussian_integer_t* matrix_part, int32_t coeff);
ctc_pauli_z_sum_t* allocate_pauli_z_sum(uint32_t term_count);
  void free_pauli_z_sum(ctc_pauli_z_sum_t* sum);

int ft_iteration_combination_incf(int k, int *sizes, int n, int *first_balls,
                                  int *last_balls, int *current_combination);

  int pauli_z_equal_p(ctc_pauli_z_t* a, ctc_pauli_z_t* b);

  void canonicalise_pauli_x(ctc_pauli_x_t* operator, int k, ctc_pauli_x_t* canonicalisers);
  void canonicalise_pauli_z(ctc_pauli_z_t* operator, int k, ctc_pauli_z_t* canonicalisers);

#ifdef __cplusplus
}
#endif
