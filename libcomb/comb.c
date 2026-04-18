#ifdef __cplusplus
extern "C" {
#endif
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "comb.h"


  int dcp_fault_bad_p_css(dcp_t *error,
                          mzd_t* x_measurement_matrix,
                          mzd_t* z_measurement_matrix,
                          mzd_t *z_measurements_for_x,
                          unsigned x_permutation, // permutation mask for x canonicalisation
                          unsigned z_permutation, // permutation mask for z canonicalisation
                          comb_bit_array_t x_lookup_table,
                          comb_bit_array_t z_lookup_table,
                          mzd_t* x_stabilisers,
                          uint32_t* rref_x,
                          int n_all, int n, int n_x,
                          stabiliser_state_t* state,
                          pauli_z_coset_t* coset){
    lwlog_debug("Testing %s",print_dcp(error,0));
    mzd_t* x_error = error->h_x_transpose; // dont need to transpose to get x error
    lwlog_debug("X error is %s",print_row_vector_as_pauli(x_error,0,"X"));
    for(int i = 0; i < z_measurements_for_x->nrows; i++){
      if(row_product(x_error,0,z_measurements_for_x,i)){
        lwlog_debug("X error anticommutes with Z measurement %d %s, so whole error is good",i,
                    print_row_vector_as_pauli(z_measurements_for_x,i,"Z"));
        return 0;
      }
      else{
        lwlog_debug("X error commutes with Z measurement %d %s, continuing...",i,
                    print_row_vector_as_pauli(z_measurements_for_x,i,"Z"));
      }
    }
    // need to canonicalise the x error.
    for(int i = 0; i < x_stabilisers->nrows; i++){
      if(mzd_read_bit(x_error,0,rref_x[i])){ // eliminate by multiplying by a stabiliser
        lwlog_debug("Multiplying by %s to canonicalise...",print_row_vector_as_pauli(x_stabilisers,i,"X"));
        mzd_combine(x_error,0,0,x_error,0,0,x_stabilisers,i,0);
        lwlog_debug("New X error is %s",print_row_vector_as_pauli(x_error,0,"X"));
      }
    }
    // apply the permutation to x.
    lwlog_debug("Error before permutation and eliminating ancillas is %s (mask %x)",print_row_vector_as_pauli(x_error,0,"X"),x_permutation);
    x_error->data[0]&=((1<<(n)) - 1);
    lwlog_debug("Error after eliminating ancilla terms is %s",print_row_vector_as_pauli(x_error,0,"X"));
    permute_columns(x_error,x_permutation);
    lwlog_debug("Error after permutation is %s",print_row_vector_as_pauli(x_error,0,"X"));
    int reduced_x_error = mzd_data(x_error);
    lwlog_debug("Reduced x error is %x",reduced_x_error);
    if(!TestBit(x_lookup_table,reduced_x_error)){
      lwlog_debug("x error undetected and high weight, so whole error is bad");
      return 1;
    }
    // Now, test the z part. Start by turning it into a state
    set_stabiliser_state_to_dcp(state, error);
    // Do first the Z and then the X postselections. If state is postselected can return 0 (good).
    if(postselect_on_measuring_plus_one_z(state,z_measurement_matrix)){
      lwlog_debug("error detected in z postselections, so error is good");
      return 0;
    }
    if(postselect_on_measuring_plus_one_x(state,x_measurement_matrix)){
      lwlog_debug("error detected in x postselections, so error is good");
      return 0;
    }
    // error is now canonicalised. turn it into a coset and apply the z permutation to it:
    set_coset_to_stabiliser_state(coset,state);
    permute_columns(coset->factor_transpose, z_permutation);
    permute_columns(coset->set_generators, z_permutation);
    // iterate over all elements in the coset and verify that each one is good. if not, return bad
    for(int i = 0; i < (1 << COSET_K(coset)); i++){
      int reduced_z_error = mzd_data(coset->factor_transpose);
      if(!TestBit(z_lookup_table,reduced_z_error)){
        lwlog_debug("z error in coset is high weight, so error is bad");
        return 1;
      }
      if((i+1)!=(1 << COSET_K(coset)))incf_coset(coset,i);
    }
    lwlog_debug("all errors in coset were low weight, so error is good");
    return 0;
  }

  int stabilisers_rref_p(comb_vector_t* operators,uint32_t n, uint32_t m, uint32_t* coefficients, unsigned* permutation_mask){
  // ensures the operators are in rref. n operators on m qubits.
  int current_leading_index=-1;
  for(int row = 0; row < n; row++){
    comb_vector_t* current_operator = &(operators[row]);
    int new_leading_index = m; // new_leading_index = n means that the row contains no entries at all.
    for(int column = 0; column < m && (new_leading_index == m); column ++){
      if(TestBit(current_operator->data,column))new_leading_index=column;
    }
    if(new_leading_index < current_leading_index && !((new_leading_index=m) && (current_leading_index==m)))
      return 0; // new leading index must be to the right of the old one
    // verify if new leading index has no entries in other rows
    for(int row_2 = 0; row_2 < n; row_2++){
      if(row_2!=row){
        if(TestBit(operators[row_2].data,new_leading_index))return 0; // leading index shares entries in other rows
      }
    }
    current_leading_index=new_leading_index;
    coefficients[row]=current_leading_index;
    *permutation_mask |= (1 << current_leading_index);
    lwlog_info("rref: coefficients[%d]=%d",row,current_leading_index);
  }
  *permutation_mask = ~*permutation_mask;
  return 1;
}

  comb_code_t __do_verify_css_magic_state( uint32_t n_all, ///< number of qubits in code + all ancillas
                                           uint32_t n, // number of qubits in code
                                           uint32_t n_x, // number of independent x stabiliser generators in code
                                           uint32_t n_z, // number of independent z stabiliser generators in code
                                           uint32_t k, // number of logical qubits in code. equivalently, number of logical dcp operators that stabilise the code
                                           mzd_t* x_stabilisers, // n-x of these
                                           comb_bit_array_t x_lookup_table, // x lookup table. of size 2^(n-n_x)
                                           mzd_t* z_stabilisers,
                                           comb_bit_array_t z_lookup_table,
                                           dcp_t* dcp_stabilisers, // k of these
                                           uint32_t fault_location_count, // number of fault locations
                                           uint32_t *fault_location_partitions, // at least fault_location_count of these. number of faults at each location
                                           uint32_t s, // number of faults to which circuit should be resilient to
                                           uint32_t *bad_fault_array, // where to store indices of bad_faults
                                           uint32_t output_error_count, // number of output errors that result from propagating all the faults. should be at least sum(fault_location_partitions)
                                           dcp_t* output_errors, // at least output_error_count of these
                                           uint32_t x_measurement_count, // number of x measurement paulis
                                           mzd_t* x_measurement_paulis, // postselections made in x basis. at least x_measurement_count of these
                                           uint32_t z_measurement_count, // number of z measurement paulis
                                           mzd_t* z_measurement_paulis,
                                           uint32_t* rref_x,
                                           uint32_t* rref_z,
                                           unsigned x_permutation,
                                           unsigned z_permutation,
                                           unsigned long long* fault_counter
                                           ){


    lwlog_debug("setting up the iteration");
    int current_combination[MAX_SIZE_COUNT] = {0};
    int first_balls[MAX_SIZE_COUNT];
    int last_balls[MAX_SIZE_COUNT];
    dcp_t* current_dcp_output_combination = init_dcp(n_all);
    dcp_t* temp = init_dcp(n_all);
    dcp_t* error_with_stabilisers = init_dcp(n_all);
    stabiliser_state_t* state = init_state(n_all);
    pauli_z_coset_t* coset = init_pauli_z_coset(n_all, n_all);
    mzd_t* z_measurement_matrix = mzd_init(z_measurement_count + n_z, n_all);  //int measurement_sequence_count = (x_measurement_count) + (2*z_measurement_count) + (2*n_z);
    mzd_t* x_measurement_matrix = mzd_init(z_measurement_count + n_z + x_measurement_count, n_all);
    mzd_t* z_measurements_for_x = mzd_init(z_measurement_count, n_all);
    for(int i = 0; i < z_measurement_count; i++){
      // WARNING! wont work for parities
      mzd_copy_row(z_measurement_matrix, i ,
                   z_measurement_paulis, i);
      mzd_copy_row(x_measurement_matrix, i ,
                   z_measurement_paulis, i);
    }
    for(int i = 0; i < n_z; i++){
      mzd_write_bit(x_measurement_matrix, i  + z_measurement_count, rref_z[i], 1);
      mzd_copy_row(z_measurement_matrix, i  + z_measurement_count,
                   z_stabilisers, i);
    }
    for(int i = 0; i < x_measurement_count; i++){
      mzd_copy_row(x_measurement_matrix, i + z_measurement_count + n_z,
                   x_measurement_paulis, i);
    }
    // now copy z measurements for x
    for(int i = 0; i < z_measurement_count; i++){
      mzd_copy_row(z_measurements_for_x,i,z_measurement_paulis,i);
    }
    for (int i = 0; i < k; i++) {
      current_combination[i] = 0;
      for (int j = 0; j < i; j++) {
        current_combination[i] += fault_location_partitions[j];
      }
    }
    for (int i = 0; i < fault_location_count; i++) {
      first_balls[i] = 0;
      for (int j = 0; j < i; j++) {
        first_balls[i] += fault_location_partitions[j];
      }
    }
    for (int i = 0; i < fault_location_count; i++) {
      last_balls[i] = 0;
      for (int j = 0; j < i + 1; j++) {
        last_balls[i] += fault_location_partitions[j];
      }
      last_balls[i]--;
    }
      do {
        ++*fault_counter;
        lwlog_debug("checking fault number %lld\n",*fault_counter);
        // Multiply DCP errors together. TODO verify that order is correct
        set_dcp_to(current_dcp_output_combination,&(output_errors[current_combination[0]]));
        for(int i = 1 ; i < s ; i++){
          dcp_mul(temp, current_dcp_output_combination, &(output_errors[current_combination[i]]));
          set_dcp_to(current_dcp_output_combination, temp);
      }
      // output errors multiplied. now, multiply by all combinations of logical operators and test them:
      // TODO they all commute (i think) so a gray code should be fine here really
      // loop over all combinations of dcp stabilisers
      int good_combination_found=0;
      for(uint32_t dcp_index=0;dcp_index < (1 << k)  && !good_combination_found;
          dcp_index++){
        set_dcp_to(error_with_stabilisers,current_dcp_output_combination);
        //
        for(int i = 0 ; i < k && !good_combination_found; i++){
          if(dcp_index & (1 << i)){
            // multiply the output error by dcp stabilisers
            lwlog_info("multiplying error by logical stabiliser");
            dcp_mul(temp,error_with_stabilisers,&(dcp_stabilisers[i]));
            set_dcp_to(error_with_stabilisers,temp);
          }
        }
        // now check to see if it is good or not.
        if(!dcp_fault_bad_p_css(error_with_stabilisers, x_measurement_matrix, z_measurement_matrix, z_measurements_for_x, x_permutation, z_permutation, x_lookup_table, z_lookup_table, x_stabilisers, rref_x, n_all, n,  n_x, state, coset)){
          // a good one was found after multiplying by a choice of logicals, set the flag
          good_combination_found=1;
          //break outer;
        }
        //error_with_stabilisers=0; //error_with_stabilisers is now bad after the destructive operation
      }
      if (!good_combination_found) {
        lwlog_debug("no good combination was found: fault combination was bad");
        // set the indices and indicate a bad fault was found
        for (int i = 0; i < s; i++) {
          bad_fault_array[i] = current_combination[i];
        }
        return BAD_FAULT_FOUND;
      }
  } while (ft_iteration_combination_incf(s, fault_location_partitions, fault_location_count, first_balls, last_balls,
                                         current_combination));
      return ALL_FAULTS_GOOD; // got to end of loop, so all faults are good
  }
  comb_code_t comb_verify_css_magic_state( uint32_t n_all, ///< number of qubits in code + all ancillas
                                           uint32_t n, // number of qubits in code
                                           uint32_t n_x, // number of independent x stabiliser generators in code
                                           uint32_t n_z, // number of independent z stabiliser generators in code
                                           uint32_t k, // number of logical qubits in code. equivalently, number of logical dcp operators that stabilise the code
                                           comb_pauli_x_t* _x_stabilisers, // n-x of these
                                           comb_bit_array_t x_lookup_table, // x lookup table. of size 2^(n-n_x)
                                           comb_pauli_z_t* _z_stabilisers,
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

    init_temps(n_all);
    lwlog_info("start timing \n");
    // Start timing
    struct timespec start, finish;
    clock_gettime( CLOCK_MONOTONIC, &start );
    lwlog_info("n_all=%d, n=%d, n_x = %d, n_z = %d, k = %d,fault_location_count = %d, s = %d, output_error_count = %d, x_measurement_count = %d, z_measurement_count = %d",n_all,n,n_x,n_z,k,fault_location_count,s,output_error_count,x_measurement_count,z_measurement_count);
    lwlog_info("allocating and converting objects");
    uint32_t rref_x[1024];
    uint32_t rref_z[1024];
    unsigned x_permutation=0;
    unsigned z_permutation=0;
    if(!stabilisers_rref_p(_x_stabilisers,n_x,n,rref_x,&x_permutation)){
      lwlog_err("x stabilisers not in rref");
      return INVALID_INPUT;
    }
    if(!stabilisers_rref_p(_z_stabilisers,n_z,n,rref_z,&z_permutation)){
      lwlog_err("z stabilisers not in rref");
      return INVALID_INPUT;
    }
    dcp_t* dcp_stabilisers = malloc(k * sizeof(dcp_t));
    dcp_t* output_errors = malloc(output_error_count * sizeof(dcp_t));
    mzd_t* x_measurement_paulis = mzd_init(x_measurement_count,n_all);
    mzd_t* z_measurement_paulis = mzd_init(z_measurement_count,n_all);
    mzd_t* x_stabilisers = mzd_init(n_x,n);
    mzd_t* z_stabilisers = mzd_init(n_z,n);
    unsigned long long fault_counter=0;
    lwlog_info("converting objects");
    for(int i = 0; i < n_x; i++){
      lwlog_debug("converting x stabiliser %d",i);
      convert_row(x_stabilisers,i,&(_x_stabilisers[i]));
      lwlog_debug("x_stabilisers[%d]=%s",i,print_row_vector_as_pauli(x_stabilisers,i,"X"));
    }
    for(int i = 0; i < n_z; i++){
      lwlog_debug("converting z stabiliser %d",i);
      convert_row(z_stabilisers,i,&(_z_stabilisers[i]));
      lwlog_debug("z_stabilisers[%d]=%s",i,print_row_vector_as_pauli(z_stabilisers,i,"X"));
    }
    for(int i = 0; i < k; i++){
      lwlog_debug("converting logical stabiliser %d",i);
      convert_dcp(&(dcp_stabilisers[i]),&(_dcp_stabilisers[i]));
      lwlog_debug("logical stabiliser %d: %s",i,print_dcp(&(dcp_stabilisers[i]),0));
    }
    lwlog_debug("converting %d output errors",output_error_count);
    for(int i = 0; i < output_error_count; i++){
      convert_dcp(&(output_errors[i]),&(_output_errors[i]));
    }
    for(int i = 0; i < x_measurement_count; i++){
      lwlog_debug("converting x measurement %d",i);
      convert_row(x_measurement_paulis,i,&(_x_measurement_paulis[i]));
      lwlog_debug("x_measurement_paulis[%d]=%s",i,print_row_vector_as_pauli(x_measurement_paulis,i,"X"));
    }
    for(int i = 0; i < z_measurement_count; i++){
      lwlog_debug("converting z measurement %d",i);
      convert_row(z_measurement_paulis,i,&(_z_measurement_paulis[i]));
      lwlog_debug("z_measurement_paulis[%d]=%s",i,print_row_vector_as_pauli(z_measurement_paulis,i,"Z"));
    }

    for(int i = 0; i < fault_location_count ; i++){
      lwlog_debug("fault location partition %d: %d\n",i,fault_location_partitions[i]);
    }
    lwlog_debug("all objects converted");
    comb_code_t return_code = __do_verify_css_magic_state(n_all, n, n_x, n_z, k, x_stabilisers, x_lookup_table, z_stabilisers, z_lookup_table,
                                                   dcp_stabilisers, fault_location_count, fault_location_partitions, s, bad_fault_array,
                                                   output_error_count, output_errors, x_measurement_count, x_measurement_paulis,
                                                   z_measurement_count, z_measurement_paulis, rref_x, rref_z, x_permutation, z_permutation, &fault_counter);
  // Stop timing
  clock_gettime( CLOCK_MONOTONIC, &finish );
  double seconds = (finish.tv_sec - start.tv_sec) + 1.0e-9 * (finish.tv_nsec - start.tv_nsec);
  printf("Took %f seconds, %lld faults checked.\n", seconds, fault_counter);
  lwlog_debug("freeing memory");
  free_dcp_array(k,dcp_stabilisers);
  free_dcp_array(output_error_count,output_errors);
  mzd_free(x_measurement_paulis);
  mzd_free(z_measurement_paulis);
  mzd_free(x_stabilisers);
  mzd_free(z_stabilisers);
  return return_code;
  }
#ifdef __cplusplus
}
#endif
