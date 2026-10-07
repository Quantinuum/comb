#include "comb.h"

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
    uint32_t *bad_fault_array,     // where to store indices of bad_faults
    uint32_t output_error_count,   // number of output errors that result from
                                   // propagating all the faults. should be at
                                   // least sum(fault_location_partitions)
    comb_dcp_t *_output_errors,    // at least output_error_count of these
    uint32_t x_measurement_count,  // number of x measurement paulis
    comb_pauli_x_t
        *_x_measurement_paulis,    // postselections made in x basis. at
                                   // least x_measurement_count of these
    uint32_t z_measurement_count,  // number of z measurement paulis
                          comb_pauli_z_t *_z_measurement_paulis){
  // TODO
  return 0; // success
}

int deserialise_run_from_file(char *input_file){
  // TODO
  return 0; // success
}
