// Putative new interface for libcomb

#ifdef __cplusplus
extern "C" {
#endif

// need fixed width integers
#include "comb_external.h"
#include <m4ri/m4ri.h>
#include "lwlog.h"
#define MAX_SIZE_COUNT 1024
  // Functions

  // idea. do rref in c? maybe a todo later. people should do it themselves...
  // it goes against the philosophy of the library
  // also cant be done as it precludes reuse of lookup tables

  // These are all internal definitions

  // use m4ri for matrix multiplication
  // no need to implement addition / multiplication etc.

  typedef struct { // Stores a DCP operator as a tableau.
    mzd_t* F; // the upper right block element of c (not c-bar). n x n matrix
    mzd_t* h_x; // n x 1 matrix
    mzd_t* h_x_transpose; // 1 x n matrix
    mzd_t* h_z; // n x 1 matrix
    mzd_t* F_diagonal; // diagonal part of F. n x 1 matrix.
  } dcp_t;

  typedef struct {
    char sign_bit; // sign bit
    mzd_t* data; // 2n + 1 x 1 matrix (column vector)
  } pauli_t; // A Pauli following Dehaene and Demoor.

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

uint32_t gray_bit_to_switch(uint32_t n);

#define DCP_N(t) ((t)->F->nrows)
#define PAULI_N(t) ((((t)->data->nrows)-1)/2)
#define STABILISER_N(t) ((((t)->data->nrows)-1)/2)
#define STABILISER_K(t) ((((t)->data->ncols)))

  // Creation
  dcp_t* init_dcp(uint32_t n);
  pauli_t* init_pauli(uint32_t n);
  pauli_z_coset_t* init_pauli_z_coset(uint32_t n, uint32_t k);
  // Freeing
  void free_pauli(pauli_t* pauli);
  void free_dcp(dcp_t* dcp);
  void free_pauli_z_coset(pauli_z_coset_t* coset);
  // Setting
  // ? what to put here. they are structs so not that much point?
  void set_dcp_to(dcp_t *a, dcp_t* b); // assigns b to a: a = b.
  // Printing
  char* print_dcp(dcp_t *dcp, char print_newline_p);
  char* print_pauli(pauli_t *pauli, char print_newline_p);
  char* print_pauli_z_coset(pauli_z_coset_t *coset, char print_newline_p);
  char* print_stabiliser_state(stabiliser_state_t *coset, char print_newline_p);
  // more convenient creation
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

  pauli_t* dcp_image_x(dcp_t *dcp, uint32_t i);
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

  int ft_iteration_combination_incf(int k, uint32_t *sizes, int n, int *first_balls,
                                    int *last_balls, int *current_combination);

  void init_dcp_to_size(dcp_t* tableau, uint32_t n);


  uint32_t mzd_data(mzd_t* mzd);
  void  permute_qubits(pauli_z_coset_t* coset,int n, mzd_t* c, mzd_t* temp, mzd_t* big_c);
  void set_to_pauli(pauli_t* a, pauli_t* b);
  void init_pauli_to_size(pauli_t* a, int size);
  void set_dcp_to_transvection(dcp_t* dcp, pauli_t* z_pauli);
  void convert_pauli(pauli_t* temp_pauli,comb_pauli_z_t* pauli);
  // conversion
  void convert_dcp(dcp_t* dcp, comb_dcp_t* value);
  void convert_pauli_x(pauli_t* pauli, comb_pauli_x_t* value);
  void convert_pauli_z(pauli_t* pauli, comb_pauli_z_t* value);
  void convert_permutation(mzd_t* permutation, uint32_t n, uint32_t n_z, uint32_t *rrefs);

  void free_pauli_array(uint32_t count, pauli_t* array);

  void free_dcp_array(uint32_t count, dcp_t* array);
  void *incf_coset(pauli_z_coset_t *coset, int index);
#ifdef __cplusplus
}
#endif
