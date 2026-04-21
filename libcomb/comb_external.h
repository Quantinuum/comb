/**
 * \file comb_external.h
 *
 * \brief External interface to Comb.
 *
 * It is intended that exposing this interface will be sufficient for invoking Comb from another language via FFI.
 *
 * \author Selwyn Simsek
 *
 */

#include <stdint.h>
#include <stdlib.h>

/**
 * The possible outcomes from checking a circuit for fault tolerance.
 */
typedef enum {ALL_FAULTS_GOOD = 0,
              BAD_FAULT_FOUND = 1,
              INVALID_INPUT = 2,
              PROBLEM_TOO_LARGE = 3,
              INTERNAL_ERROR = 4,
              OTHER_ERROR = 5} comb_code_t;

/// A word; assumed to always be 64 bits wide.
typedef uint64_t comb_bit_word_t;
/// A bit array is an array of words.
typedef comb_bit_word_t* comb_bit_array_t;

/// Sets the k-th bit of a bit array A to 1.
#define SetBit(A, k) (A[(k) / sizeof(comb_bit_word_t)] |= (1U << ((k) % sizeof(comb_bit_word_t))))
/// Clears the k-th bit of a bit array A (sets it to 0).
#define ClearBit(A, k) (A[(k) / sizeof(comb_bit_word_t)] &= ~(1U << ((k) % sizeof(comb_bit_word_t))))
/// Evaluates to 0 if the k-th bit of the bit array A is 0, and non-zero otherwise.
#define TestBit(A, k) (A[(k) / sizeof(comb_bit_word_t)] & (1U << ((k) % sizeof(comb_bit_word_t))))

/**
 * A bit vector.
 */
typedef struct {
  uint32_t n; /// The length of the vector
  /**
   * Bit array that stores the entries of the vector.
   * Must have enough space for n bits.
   */
  comb_bit_array_t data;
} comb_vector_t;

/**
 * A bit matrix of fixed size.
 */
typedef struct {
  uint32_t m; ///< The number of rows of the matrix.
  uint32_t n; ///< The number of columns of the matrix.
  /**
   * A bit array that stores the entries of the array in row-major order.
   * Must have enough room to store \f$ mn \f$ bits.
   */
  comb_bit_array_t data;
} comb_matrix_t;

/**
 * The product of a diagonal Clifford and a Pauli (DCP).
 */
typedef struct {
  /**
   * Stores the Z parts of the images of the \f$X_i\f$.
   * Must be a \f$ n \times n\f$ dimensional symmetric matrix.
   */
  comb_matrix_t F;
  /**
   * Stores the signs of the images of the Paulis \f$Z_i\f$ by conjugation.
   * i.e. the i-th bit is set to 1 if the Pauli \f$Z_i\f$ changes sign and 0 otherwise.
   * Must be a \f$ n \f$-dimensional vector
   */
  comb_vector_t h_x;
  /**
   * Stores the sign bits of the images of the Paulis \f$X_i\f$ by conjugation.
   * Must be a \f$ n \f$-dimensional vector
   */
  comb_vector_t h_z;
} comb_dcp_t;


typedef comb_vector_t comb_pauli_x_t; ///< A Pauli of X type
typedef comb_vector_t comb_pauli_z_t; ///< A Pauli of Z type

// Allocate bit arrays, arrays of vectors, dcps and permutations
/// Allocate a bit array.
comb_bit_array_t comb_alloc_bit_array(int n /*!< The number of bits in the array.*/);
/// Allocate an array of comb_vector_t.
comb_vector_t *comb_calloc_vector(int count /*!< The number of vectors to allocate.*/, int n /*!< The length of each vector.*/);
/// Allocate an array of comb_dcp_t.
comb_dcp_t *comb_calloc_dcp(int count, int n);

// Convenience methods for manipulating arrays of dcps, vectors and bit arrays

/// Sets a bit of a vector in an array of vectors.
void comb_set_cvector_aref(comb_vector_t* vector, int n /*!<The index of the vector in the array.*/, int i /*!< The index of the bit to set.*/, int value /*!< The new value of the bit.*/);

/// Sets an entry of the matrix F of a DCP operator in an array of DCP operators.
void comb_set_cdcp_F_aref(comb_dcp_t* dcp, int n  /*!<The index of the vector in the array.*/, int i /*!< The row of the entry to set.*/, int j /*!< The column of the entry to set.*/, int value /*!< The new value of the bit.*/);

/// Sets an entry of the vector \f$h_x\f$ of a DCP operator in an array of DCP operators.
void comb_set_cdcp_h_x_aref(comb_dcp_t* dcp, int n  /*!<The index of the vector in the array.*/, int i /*!< The index of the bit to set.*/, int value /*!< The new value of the bit.*/);

/// Sets an entry of the vector \f$h_z\f$ of a DCP operator in an array of DCP operators.
void comb_set_cdcp_h_z_aref(comb_dcp_t* dcp, int n  /*!<The index of the vector in the array.*/, int i /*!< The index of the bit to set.*/, int value /*!< The new value of the bit.*/);

/// Sets a bit in a bit array.
void comb_set_bit_array_aref(comb_bit_array_t array, int i /*!< The index of the bit to set.*/, int value /*!< The new value of the bit.*/);

/// Frees a bit array previously created with comb_alloc_bit_array.
void comb_free_bit_array(comb_bit_array_t pointer /*!< The pointer to previously allocated storage.*/);
/// Frees an array of comb_vector_t previously created with comb_calloc_vector.
void comb_cfree_vector(int count /*!< The length of the array.*/, comb_vector_t* pointer /*!< The pointer to previously allocated storage.*/);
/// Frees an array of comb_dcp_t previously created with comb_calloc_dcp.
void comb_cfree_dcp(int count/*!< The length of the array.*/, comb_dcp_t* pointer /*!< The pointer to previously allocated storage.*/);

/**
 * Determines whether a circuit that prepares the magic state of a CSS code is fault tolerant or not.
 * One might invoke this function for all \f$1 \leq s \leq \frac{d}{2}\f$, choosing the appropriate lookup tables for each invocation.
 * @return A code indicating the result of the verification. If BAD_FAULT_FOUND, then bad_fault_array is set to the indices of the bad faults in output_error_array.
 */
comb_code_t comb_verify_css_magic_state( // 1. Describing the magic state
                                           uint32_t n_all, ///< The number of qubits that comprise the code and all the ancillas.
                                           uint32_t n, ///< The number of qubits in the code.
                                           uint32_t n_x, ///< The number of independent stabiliser generators of X type in the code.
                                           uint32_t n_z, ///< The number of independent stabiliser generators of Z type in the code.
                                           uint32_t k, ///< The number of logical qubits in the code.
                                           comb_pauli_x_t* x_stabilisers, ///< An pointer to the stabiliser generators of X type. They must be in reduced row echelon form. There must be n_x stabiliser generators.
                                           comb_bit_array_t x_lookup_table, ///< A bit array that is the lookup table that determines whether the X parts of the postselected errors are good or not. The lookup table must be of size \f$2^{n-n_x} \f$ bits.
                                           comb_pauli_z_t* z_stabilisers, ///< A pointer to the stabiliser generators of Z type. They must be in reduced row echelon form. There must be n_z stabiliser generators.
                                           comb_bit_array_t z_lookup_table,///< A bit array that is the lookup table that determines whether the Z parts of the postselected errors are good or not. The lookup table must be of size \f$2^{n-n_z} \f$ bits.
                                           comb_dcp_t* dcp_stabilisers, ///< A pointer to the logical DCP operators that stabilise the magic state and thus characterise it uniquely. There must be k such operators.
                                           // 2. Other parameters
                                           uint32_t fault_location_count, ///< The number of fault locations in the circuit.
                                           uint32_t *fault_location_partitions, ///< A pointer to the numbers of faults at each fault location. There must be at least fault_location_count of them.
                                           uint32_t s, ///< The number of faults that make up each fault combination.
                                           uint32_t *bad_fault_array, ///< A pointer to where indices of bad faults can be stored. It must have enough space for s elements.
                                           uint32_t output_error_count, ///< The number of output errors that result from propagating all the faults. This number should be at least the sum of all the fault_location_partitions.
                                           comb_dcp_t* output_errors, ///< A pointer to the DCP output errors that result from propagating each fault to the end of the circuit. There must be at least output_error_count errors.
                                           uint32_t x_measurement_count, ///< The number of ancilla measurements made in the X basis.
                                           comb_pauli_x_t* x_measurement_paulis, ///< A pointer to the measurement operators corresponding to postselected measurements made in the X basis. There must be at least x_measurement_count of these
                                           uint32_t z_measurement_count, ///< The number of ancilla measurements made in the Z basis.
                                           comb_pauli_z_t* z_measurement_paulis ///< A pointer to the measurement operators corresponding to postselected measurements made in the Z basis. There must be at least z_measurement_count of these
                                         );
