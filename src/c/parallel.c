#include "comb.h"

#include <math.h>
#include <stdio.h>
#include <time.h>
#include <omp.h>

// Routines to help with parallelism with OpenMP

comb_code_t serial_check_chunks(size_t chunk_count, chunk_t* chunks,
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
                                unsigned x_permutation, unsigned z_permutation){
  for(int i = 0; i < chunk_count; i++){
    comb_code_t chunk_return_code = __do_verify_css_magic_state(
                                                                n_all, n, n_x, n_z, k, x_stabilisers, x_lookup_table, z_stabilisers,
                                                                z_lookup_table, dcp_stabilisers, fault_location_count,
                                                                fault_location_partitions, s, bad_fault_array, output_error_count,
                                                                output_errors, x_measurement_count, x_measurement_paulis,
                                                                z_measurement_count, z_measurement_paulis, rref_x, rref_z, x_permutation,
                                                                z_permutation, chunks+i, 0);
    if(chunk_return_code != ALL_FAULTS_GOOD)return chunk_return_code;

  }
  return ALL_FAULTS_GOOD;
}


comb_code_t openmp_check_chunks(size_t chunk_count, chunk_t* chunks,  uint32_t n_all,  ///< number of qubits in code + all ancillas
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
                                unsigned x_permutation, unsigned z_permutation){
  comb_bit_array_t array = comb_alloc_bit_array(chunk_count);
  size_t chunks_checked_so_far = 0;
  unsigned long long faults_checked_so_far = 0;
  int bad_fault_foundp = 0;
  comb_code_t return_code;


  #pragma omp parallel
  {
    int continue_with_loop = 1;
    while(continue_with_loop){
      int id_to_scan = -1;
#pragma omp critical
      {
        for(int i = 0; i < chunk_count && (id_to_scan == -1) && !bad_fault_foundp; i++){
          if(!TestBit(array,i)){
            id_to_scan=i;
            SetBit(array,i);
          }
        }
      }
      if(id_to_scan == -1){
        // no free block found
        continue_with_loop = 0;
        printf("thread %d; no free chunk found, quitting\n",omp_get_thread_num());
      }
      else{
        time_t currentTime;
        time(&currentTime);
        printf("thread %d: started scanning %d [%ld / %ld (%.2f%%) chunks checked so far] %s",omp_get_thread_num(),id_to_scan, chunks_checked_so_far , chunk_count,
               ((100.0 *(chunks_checked_so_far ))/ chunk_count),ctime(&currentTime));

        comb_code_t chunk_return_code = __do_verify_css_magic_state(
                                                                    n_all, n, n_x, n_z, k, x_stabilisers, x_lookup_table, z_stabilisers,
                                                                    z_lookup_table, dcp_stabilisers, fault_location_count,
                                                                    fault_location_partitions, s, bad_fault_array, output_error_count,
                                                                    output_errors, x_measurement_count, x_measurement_paulis,
                                                                    z_measurement_count, z_measurement_paulis, rref_x, rref_z, x_permutation,
                                                                    z_permutation, chunks+ id_to_scan, 0);
        #pragma omp atomic
        chunks_checked_so_far++;
        #pragma omp atomic
        faults_checked_so_far+=chunks[id_to_scan].faults_already_checked;
        printf("thread: %d: checked chunk %d, %llu faults checked so far across all threads\n",omp_get_thread_num(),id_to_scan,faults_checked_so_far);
        if(chunk_return_code != ALL_FAULTS_GOOD){
          bad_fault_foundp=1;
          return_code = chunk_return_code;
        }
      }
    }
  }

  free(array);
  return return_code;
}
