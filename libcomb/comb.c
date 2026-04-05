#ifdef __cplusplus
extern "C" {
#endif
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h> 
#include "comb.h"

  comb_code_t comb_verify_css_magic_state( uint32_t n_all, // number of qubits in code + all ancillas
                                           uint32_t n, // number of qubits in code
                                           uint32_t n_x, // number of independent x stabiliser generators in code
                                           uint32_t n_z, // number of independent z stabiliser generators in code
                                           uint32_t k, // number of logical qubits in code. equivalently, number of logical dcp operators that stabilise the code
                                           comb_pauli_x_t* _x_stabilisers, // n-x of these
                                           comb_permutation_t* _x_permutation, //permutation that puts x_stabilisers in rref form
                                           comb_bit_array_t x_lookup_table, // x lookup table. of size 2^(n-n_x)
                                           comb_pauli_z_t* _z_stabilisers,
                                           comb_permutation_t* _z_permutation, // permutation that puts z_stabilisers in rref form
                                           comb_bit_array_t z_lookup_table,
                                           comb_dcp_t* _dcp_stabilisers, // k of these

                                           uint32_t fault_location_count, // number of fault locations
                                           uint32_t *fault_location_partitions, // at least fault_location_count of these. number of faults at each location
                                           uint32_t s, // number of faults to which circuit should be resilient to
                                           uint32_t *bad_fault_array, // where to store indices of bad_faults
                                           uint32_t output_error_count, // number of output errors that result from propagating all the faults. should be at least sum(fault_location_partitions)
                                           comb_dcp_t* _output_errors, // at least output_error_count of these
                                           uint32_t x_measurement_count, // number of x measurement paulis
                                           comb_pauli_x_t* _x_measurement_paulis, // postselections made in x basis. at least x_measurement_count of these
                                           uint32_t z_measurement_count, // number of z measurement paulis
                                           comb_pauli_z_t* _z_measurement_paulis){
    lwlog_info("start timing \n");
    // Start timing
    int fault_counter=0;
    struct timespec start, finish;
    clock_gettime( CLOCK_MONOTONIC, &start );
    lwlog_debug("n_all=%d, n=%d, n_x = %d, n_z = %d, k = %d,fault_location_count = %d, s = %d, output_error_count = %d, x_measurement_count = %d, z_measurement_count = %d",n_all,n,n_x,n_z,k,fault_location_count,s,output_error_count,x_measurement_count,z_measurement_count);
    lwlog_info("allocating objects");
    dcp_t* dcp_stabilisers = malloc(k * sizeof(dcp_t));
    dcp_t* output_errors = malloc(output_error_count * sizeof(dcp_t));
    pauli_t* x_measurement_paulis = malloc(x_measurement_count * sizeof(pauli_t));
    pauli_t* z_measurement_paulis = malloc(z_measurement_count * sizeof(pauli_t));
    pauli_t* x_stabilisers = malloc(n_x * sizeof(pauli_t));
    pauli_t* z_stabilisers = malloc(n_z * sizeof(pauli_t));
    mzp_t* x_permutation = mzp_init(n);
    mzp_t* z_permutation = mzp_init(n);

    lwlog_info("converting objects");
    printf("stabilisers: \n");
    for(int i = 0; i < n_x; i++){
      convert_pauli_x(&(x_stabilisers[i]),&(_x_stabilisers[i]));
      print_pauli(&(x_stabilisers[i]),1);
    }
    for(int i = 0; i < n_z; i++){
      convert_pauli_z(&(z_stabilisers[i]),&(_z_stabilisers[i]));
      print_pauli(&(z_stabilisers[i]),1);
    }
    for(int i = 0; i < k; i++){
      convert_dcp(&(dcp_stabilisers[i]),&(_dcp_stabilisers[i]));
      printf("logical stabiliser %d: ",i);
      print_dcp(&(dcp_stabilisers[i]),1);
    }
    for(int i = 0; i < output_error_count; i++){
      convert_dcp(&(output_errors[i]),&(_output_errors[i]));
      printf("output error %d: ",i);
      print_dcp(&(output_errors[i]),1);
    }
    printf("measurements: \n");
    for(int i = 0; i < x_measurement_count; i++){
      convert_pauli_x(&(x_measurement_paulis[i]),&(_x_measurement_paulis[i]));
      print_pauli(&(x_measurement_paulis[i]),1);
    }
    for(int i = 0; i < z_measurement_count; i++){
      convert_pauli_z(&(z_measurement_paulis[i]),&(_z_measurement_paulis[i]));
      print_pauli(&(z_measurement_paulis[i]),1);
    }
    printf("x lookup table: ");
    for(int i = 0; i < (1 << (n-n_x)); i++){
      printf("%d",TestBit(x_lookup_table,i));
        }
    printf("\n");
    printf("z lookup table: ");
    for(int i = 0; i < (1 << (n-n_z)); i++){
      printf("%d",TestBit(z_lookup_table,i));
    }
    printf("\n");
    convert_permutation(x_permutation,_x_permutation);
    convert_permutation(z_permutation,_z_permutation);
    mzp_print(x_permutation);
    mzp_print(z_permutation);

  // Stop timing
  clock_gettime( CLOCK_MONOTONIC, &finish );
  double seconds = (finish.tv_sec - start.tv_sec) + 1.0e-9 * (finish.tv_nsec - start.tv_nsec);
  printf("Took %f seconds.\n", seconds);
  // TODO free everything

  return ALL_FAULTS_GOOD;
  }

  void convert_dcp(dcp_t* dcp, comb_dcp_t* value){}
  void convert_pauli_x(pauli_t* pauli, comb_pauli_x_t* value)
  {
    int n = value->n;
    pauli->data=mzd_init(1+(2*n),1);
    pauli->sign_bit=0;
    for(int i = 0; i < n;i++){
      //printf("bit %d=%d\n",i,TestBit(value->data,i));
      if(TestBit(value->data,i)){
        //printf("writing bit %d\n",i);
      mzd_write_bit(pauli->data,i+n,0,1);
      }
    }
  }
  void convert_pauli_z(pauli_t* pauli, comb_pauli_z_t* value){
    int n = value->n;
    pauli->data=mzd_init(1+(2*n),1);
    pauli->sign_bit=0;
    for(int i = 0; i < n;i++){
      //printf("bit %d=%d\n",i,TestBit(value->data,i));
      if(TestBit(value->data,i)){
        //printf("writing bit %d\n",i);
        mzd_write_bit(pauli->data,i,0,1);
      }
    }
  }
  void convert_permutation(mzp_t* permutation, comb_permutation_t* value){}


#ifdef __cplusplus
}
#endif
