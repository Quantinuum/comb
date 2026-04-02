#ifdef __cplusplus
extern "C" {
#endif
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h> 
#include "comb.h"

  comb_code_t comb_verify_css_magic_state(int n, int *sizes, int k,
                                          size_t *bad_fault_array,
                                          comb_css_magic_state_t* state,
                                          size_t output_error_count, comb_dcp_t* output_errors,
                                          size_t x_measurement_count, comb_pauli_x_t* x_measurement_paulis,
                                          size_t z_measurement_count, comb_pauli_z_t* z_measurement_paulis){
    // Start timing
    int fault_counter=0;
    struct timespec start, finish;
    clock_gettime( CLOCK_MONOTONIC, &start );

    // Convert the output errors to dcp.
    dcp_t* dcp_output_errors = malloc(output_error_count * sizeof(dcp_t));
    for(int i = 0 ; i < output_error_count ; i++){
      convert_dcp(&(dcp_output_errors[i]),&(output_errors[i])); // TODO write convert_dcp
    }
    dcp_t* logical_stabilisers = malloc(state->k * sizeof(dcp_t));
    for(int i = 0 ; i < state->k ; i++){
      convert_dcp(&(logical_stabilisers[i]),&(state->dcp_stabilisers[i]));
    }

    // Check that the number of faults is not too large
    uint32_t fault_count = 0;
    for(int i = 0;i<n; i++){
      fault_count+=sizes[i];
    }
    if (n > MAX_SIZE_COUNT) {
      printf("problem too large: %i > %i\n", n, MAX_SIZE_COUNT);
      return PROBLEM_TOO_LARGE;
    }



    // Allocate space for multiplying dcp errors together in the iteration
    int qubit_count = output_errors[0].pauli_x_part.pauli_length;
    dcp_t current_dcp_output_combination;
    dcp_t temp;
    dcp_t error_with_stabilisers;
    init_dcp_to_size(&temp,qubit_count);
    init_dcp_to_size(&current_dcp_output_combination,qubit_count);
    init_dcp_to_size(&error_with_stabilisers,qubit_count);
    // Construct the sequence of measurements that postselect and canonicalise all dcp errors
    int measurement_sequence_count = 0;
    pauli_t* measurement_sequence = malloc(measurement_sequence_count*sizeof(pauli_t));
    // TODO construct the sequence
    // Set up the iteration
    int current_combination[MAX_SIZE_COUNT];
    int first_balls[MAX_SIZE_COUNT];
    int last_balls[MAX_SIZE_COUNT];
    for (int i = 0; i < k; i++) {
      current_combination[i] = 0;
      for (int j = 0; j < i; j++) {
        current_combination[i] += sizes[j];
      }
    }
    for (int i = 0; i < n; i++) {
      first_balls[i] = 0;
      for (int j = 0; j < i; j++) {
        first_balls[i] += sizes[j];
      }
    }
    for (int i = 0; i < n; i++) {
      last_balls[i] = 0;
      for (int j = 0; j < i + 1; j++) {
        last_balls[i] += sizes[j];
      }
      last_balls[i]--;
    }
    // do the iteration
    unsigned long long c = 0;
    do {
      if(fault_counter%10000 == 0)printf("checking fault number %d\n",fault_counter++);
      // Multiply DCP errors together. TODO verify that order is correct
      set_dcp_to(&current_dcp_output_combination,&(dcp_output_errors[current_combination[0]]));
      for(int i = 1 ; i < k ; i++){
        dcp_mul(&temp,&current_dcp_output_combination,&(dcp_output_errors[current_combination[i]]));
        current_dcp_output_combination=temp;
      }
      // output errors multiplied. now, multiply by all combinations of logical operators and test them:
      // TODO they all commute (i think) so a gray code should be fine here really
      // loop over all combinations of dcp stabilisers
      int good_combination_found=0;
      for(uint32_t dcp_index=0;dcp_index < (1 << state->k)  && !good_combination_found;
          dcp_index++){
        set_dcp_to(&error_with_stabilisers,&current_dcp_output_combination);
        //
        for(int i = 0 ; i < state->k && !good_combination_found; i++){
          if(dcp_index & (1 << i)){
            // multiply the output error by dcp stabilisers
            dcp_mul(&temp,&error_with_stabilisers,&(logical_stabilisers[i]));
            set_dcp_to(&error_with_stabilisers,&temp);
          }

        }
        // now check to see if it is good or not.
        if(!dcp_fault_bad_p_css(&error_with_stabilisers,state,measurement_sequence,measurement_sequence_count)){
          // a good one was found after multiplying by a choice of logicals, set the flag
          good_combination_found=1;
          //break outer;
        }
        //error_with_stabilisers=0; //error_with_stabilisers is now bad after the destructive operation
      }
      if (!good_combination_found) {
        // need to return a fault
        for (int i = 0; i < k; i++) {
          bad_fault_array[i] = current_combination[i];
        }
        fflush(stdout);
        clock_gettime( CLOCK_MONOTONIC, &finish );
        double seconds = (finish.tv_sec - start.tv_sec) + 1.0e-9 * (finish.tv_nsec - start.tv_nsec);
        printf("Took %f seconds.\n", seconds);
        return BAD_FAULT_FOUND;
      }

    c++;
  } while (ft_iteration_combination_incf(k, sizes, n, first_balls, last_balls,
                                         current_combination));
  printf("%llu faults processed.\n", c);
  fflush(stdout);

  // Stop timing
  clock_gettime( CLOCK_MONOTONIC, &finish );
  double seconds = (finish.tv_sec - start.tv_sec) + 1.0e-9 * (finish.tv_nsec - start.tv_nsec);
  printf("Took %f seconds.\n", seconds);

  return ALL_FAULTS_GOOD;
  }
  int dcp_fault_bad_p_css(dcp_t *error, comb_css_magic_state_t* magic_state, pauli_t* p, int count){
    // First, check the X part of the error, as this is easier to deal with.
    // Does it anticommute with any Z basis postselections? If so, the error is always detected, so the entire error is good
    // Is it low weight? If so, it is good, but we must check the rest of the error. If it is not low weight, then the whole error is bad
    // Now we do the more expensive Z check.
    // Convert the dcp into a stabiliser state:
    stabiliser_state_t* state = dcp_to_stabiliser_state(error);
    // Do all postselections and canonicalisations by measuring all the Paulis on the state.
    // If the state ever vanishes, the error is good, so return.
    for(int i =0 ; i< count; i++){
      if(postselect_on_measuring_plus_one(state,&(p[i])))return 0; // vanishes so return good
    }
    // State now postselected and canonicalised. Expand the resulting state as a coset
    pauli_z_coset_t *coset = diagonal_clifford_to_coset(state);
    // Apply the qubit permutation:
    permute_qubits(coset,magic_state->n,magic_state->z_permutation);

    pauli_t current_element;
    pauli_t temp;
    set_to_pauli(&current_element,coset->factor);
    init_pauli_to_size(&temp,PAULI_N(coset->factor));
    // Iterate over every element in the coset using a gray code
    for(int i = 0 ; i < (1 << coset->k); i++){
      // If any term is bad, the error is bad, otherwise it is good.
      pauli_mul(&temp,&current_element,coset->set_generators[gray_bit_to_switch(i)]);
      if(!TestBit(magic_state->z_lookup_table,mzd_data(temp.data))) {
        // it's bad so return, otherwise check the other ones
        return 1;
      }
      current_element=temp;
    }
  // all elements in the coset were good so return
    return 0;
  }
  uint32_t gray_bit_to_switch(uint32_t n)
{
  uint32_t n_add_one = n + 1;
  uint32_t diff = n ^ n_add_one ^ (n >> 1) ^ (n_add_one >> 1);
  return n ? 31 - __builtin_clz(diff) : 0;
}

uint32_t mzd_data(mzd_t* mzd){
  int n = mzd->nrows;
  int m = mzd->ncols;
  uint32_t elt=0;
  int count = 0;
  for(int i = 0; i < n; i++){
    for(int j = 0; j < m; j++){
      if(mzd_read_bit(mzd,i,j))elt|=(1<< count);
      count++;
    }
  }
  return elt;
}


  void  permute_qubits(pauli_z_coset_t* coset,int n, size_t* c){
  mzp_t* perm = mzp_init(2*n+1);
  for(int i =0; i < n; i++)
    perm->values[i]=c[i]; // only need to permute zs
  mzd_apply_p_left(coset->factor->data,perm);
  for(int i =0 ; i< coset->k; i++)
    mzd_apply_p_left(coset->set_generators[i]->data,perm);
}

  void set_to_pauli(pauli_t* a, pauli_t* b){
    a->sign_bit=b->sign_bit;
    mzd_copy(a->data,b->data);
  }
  void init_pauli_to_size(pauli_t* a, int n){
    a->sign_bit=0;
    a->data=mzd_init(1+(2*n),1);
  }

  void convert_dcp(dcp_t* dcp, comb_dcp_t* comb){
    // TODO do something
    int n = comb->pauli_x_part.pauli_length;
    dcp_t* temp1 = init_dcp(n);
    dcp_t* temp2 = init_dcp(n);
    pauli_t* temp_pauli = init_pauli(n);
    for(int i =0; i < comb->exponent_sum.term_count; i ++){
      int coefficient = comb->exponent_sum.coefficients[i];
      comb_pauli_z_t* pauli= &(comb->exponent_sum.terms[i]);
      while(coefficient < 0)coefficient +=8;
      while(coefficient > 8)coefficient -=8;
      for(int j = 0; j < coefficient; j++){
        convert_pauli(temp_pauli,pauli);
        set_dcp_to_transvection(temp1,temp_pauli);
        dcp_mul(temp2, temp1, dcp);
        set_dcp_to(dcp,temp2);
      }
    }

    // finally write the x part
    for(int i =0 ; i< n; i++){
      mzd_write_bit(dcp->h_z,i,0,TestBit(comb->pauli_x_part.data,i));
    }
    free_dcp(temp1);
    free_dcp(temp2);
    free_pauli(temp_pauli);

  }


  void convert_pauli(pauli_t* temp_pauli,comb_pauli_z_t* pauli){
    int n = PAULI_N(temp_pauli);
    for(int i =0 ; i< n; i++)mzd_write_bit(temp_pauli->data, i,0,TestBit(pauli->data,i));
  }
#ifdef __cplusplus
}
#endif
