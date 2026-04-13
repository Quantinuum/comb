#ifdef __cplusplus
extern "C" {
#endif
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h> 
#include "comb.h"
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
    lwlog_debug("in mzd_data, returning %d",elt);
    return elt;
  }

  uint32_t mzd_x_data(mzd_t* mzd){
    int n = mzd->nrows;
    int m = mzd->ncols;
    uint32_t elt=0;
    int count = 0;
    for(int i = ((n-1)/2); i < n; i++){
      for(int j = 0; j < m; j++){
        if(mzd_read_bit(mzd,i,j))elt|=(1<< count);
        count++;
      }
    }
    lwlog_debug("in mzd_x_data, returning %d",elt);
    return elt;
  }

  void  permute_qubits(pauli_z_coset_t* coset,int n, mzd_t* c){
    lwlog_debug("in permute qubits");
    mzd_t* temp = mzd_init(2*n+1,1);
    mzd_t* big_c = mzd_init((2*n)+1,(2*n)+1);
    mzd_set_ui(big_c,1);
    mzd_copy(big_c,c);
    mzd_mul_naive(temp,big_c,coset->factor->data);
    mzd_copy(coset->factor->data,temp);
    for(int i =0 ; i< coset->k; i++){
      mzd_mul_naive(temp,big_c,coset->set_generators[i]->data);
      mzd_copy(coset->set_generators[i]->data,temp);
    }
    mzd_free(temp);
    mzd_free(big_c);
  }

  void  permute_x_pauli(pauli_t* pauli, int n, int n_all, mzd_t* permutation){
    lwlog_debug("in permute qubits");
    mzd_t* temp = mzd_init(2*n_all+1,1);
    mzd_t* big_c = mzd_init((2*n_all)+1,(2*n_all)+1);
    mzd_set_ui(big_c,1);
    for(int i = 0; i < n; i ++){
      for(int j = 0; j < n; j++){
        mzd_write_bit(big_c,i+n_all,j+n_all,mzd_read_bit(permutation,i,j));
      }
    }
    mzd_mul_naive(temp,big_c,pauli->data);
    mzd_copy(pauli->data,temp);
    mzd_free(temp);
    mzd_free(big_c);
  }
  void set_to_pauli(pauli_t* a, pauli_t* b){
    //lwlog_debug("in set_to_pauli");
    //print_pauli(a,1);
    //print_pauli(b,1);
    a->sign_bit=b->sign_bit;
    mzd_copy(a->data,b->data);
  }
  void init_pauli_to_size(pauli_t* a, int n){
    a->sign_bit=0;
    a->data=mzd_init(1+(2*n),1);
  }

  void* init_to_pauli_x (pauli_t* p, int i, int n_all){
    p->data=mzd_init((n_all*2) + 1,1);
    mzd_write_bit(p->data,i+ n_all,0,1);
  }
  void* init_to_pauli (pauli_t* p, int n_all){
    p->data=mzd_init((n_all*2) + 1,1);
  }

  int dcp_fault_bad_p_css(dcp_t *error, pauli_t* p, int count, int n, mzd_t* z_permutation, comb_bit_array_t z_lookup_table, mzd_t* x_permutation, comb_bit_array_t x_lookup_table,pauli_t* x_stabilisers, uint32_t n_x, uint32_t* rref_x, uint32_t z_measurement_count, pauli_t* z_measurement_paulis, int n_all){
    lwlog_debug("Testing %s",print_dcp(error,0));
    //print_dcp(error,1);
    // TODO First, check the X part of the error, as this is easier to deal with.
    // Does it anticommute with any Z basis postselections? If so, the error is always detected, so the entire error is good
    pauli_t* x_part = init_pauli(n_all);//mzd_init(1+(2*n),1);
    for(int i = 0; i < n_all; i++){
      mzd_write_bit(x_part->data,i+n_all,0,mzd_read_bit(error->h_x,i,0));
    }
    lwlog_debug("x_part pauli %s before measurement",print_pauli(x_part,0));
    for(int j = 0; j < z_measurement_count; j++){
      if(!pauli_commute_p(&(z_measurement_paulis[j]),x_part))return 0; // detected
    }
    pauli_t* x_part_residual = init_pauli(n);//mzd_init(1+(2*n),1);
    for(int i = 0; i < n; i++){
      mzd_write_bit(x_part_residual->data,i+n,0,mzd_read_bit(x_part->data,i+n_all,0));
    }
    pauli_t* x_temp = init_pauli(n);//mzd_init(1+(2*n),1);
      // now put in canonical form
    for(int i = 0; i < n_x; i++){
      if(mzd_read_bit(x_part_residual->data,n+rref_x[i],0)){
        pauli_mul(x_temp, x_part_residual,&(x_stabilisers[i]));
        set_to_pauli(x_part_residual,x_temp);
      }
    }
    lwlog_debug("x_part pauli %s after canonicalisation",print_pauli(x_part_residual,0));
    // now permute it
    permute_x_pauli(x_part_residual,n,n,x_permutation);
    lwlog_debug("x_part pauli %s after zeroing and permutation",print_pauli(x_part_residual,0));
    // Is it low weight? If so, it is good, but we must check the rest of the error. If it is not low weight, then the whole error is bad
    if(!TestBit(x_lookup_table,mzd_x_data(x_part_residual->data)))
      {
        lwlog_debug("x_part is bad, return 1",print_pauli(x_part_residual,0));
        return 1;
      }
    lwlog_debug("x_part is good, continue to z part",print_pauli(x_part_residual,0));
    // Now we do the more expensive Z check.
    // Convert the dcp into a stabiliser state:
    stabiliser_state_t* state = dcp_to_stabiliser_state(error);
    // Do all postselections and canonicalisations by measuring all the Paulis on the state.
    // If the state ever vanishes, the error is good, so return.
    for(int i =0 ; i< count; i++){
      lwlog_debug("Postselecting measurement %d: %s",i,print_pauli(&(p[i]),0));
      if(postselect_on_measuring_plus_one(state,&(p[i]))){
        lwlog_debug("Postselecting measurement %d vanished, returning 0",i);
      return 0; // vanishes so return good
      }
      lwlog_debug("Finished postselecting %d",i);
    }
    lwlog_debug("state after postselecting is %s",print_stabiliser_state(state,0));
    // State now postselected and canonicalised. Expand the resulting state as a coset
    pauli_z_coset_t *coset = diagonal_clifford_to_coset(state);
    lwlog_debug("Diagonal clifford to coset %s",print_pauli_z_coset(coset,0));
    // Apply the qubit permutation:
    permute_qubits(coset,n,z_permutation);
    lwlog_debug("Permuted coset: %s",print_pauli_z_coset(coset,0));
    pauli_t current_element;
    pauli_t temp;
    init_pauli_to_size(&current_element,PAULI_N(coset->factor));
    set_to_pauli(&current_element,coset->factor);
    init_pauli_to_size(&temp,PAULI_N(coset->factor));
    // Iterate over every element in the coset using a gray code
    for(int i = 0 ; i < (1 << coset->k); i++){
      lwlog_debug("iterating %d: %s in gray code",i,print_pauli(&current_element,0));
      // If any term is bad, the error is bad, otherwise it is good.
      if(!TestBit(z_lookup_table,mzd_data(current_element.data))) {
        //it's bad so return, otherwise check the other ones
        lwlog_debug("%s is bad, returning 1",print_pauli(&current_element,0));
        return 1;
       }
      else{
        lwlog_debug("%s was good",print_pauli(&current_element,0));
      }

      if((i+1)!=(1<<coset->k)){
        pauli_mul(&temp,&current_element,coset->set_generators[gray_bit_to_switch(i)]);
        set_to_pauli(&current_element,&temp);
      }
    }
    lwlog_debug("Finished gray code loop, returning 0");
    mzd_free(temp.data);
    mzd_free(current_element.data);
    lwlog_debug("freed data");
    return 0;
  // all elements in the coset were good so return
    return 0;
  }
  int stabilisers_rref_p(comb_vector_t* operators,uint32_t n, uint32_t m, uint32_t* coefficients){
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
    lwlog_info("rref: coefficients[%d]=%d",row,current_leading_index);
  }
  return 1;
}

  comb_code_t comb_verify_css_magic_state( uint32_t n_all, // number of qubits in code + all ancillas
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
    lwlog_info("start timing \n");
    // Start timing
    int fault_counter=0;
    struct timespec start, finish;
    clock_gettime( CLOCK_MONOTONIC, &start );
    lwlog_info("n_all=%d, n=%d, n_x = %d, n_z = %d, k = %d,fault_location_count = %d, s = %d, output_error_count = %d, x_measurement_count = %d, z_measurement_count = %d",n_all,n,n_x,n_z,k,fault_location_count,s,output_error_count,x_measurement_count,z_measurement_count);
    lwlog_info("allocating and converting objects");
    uint32_t rref_x[1024];
    uint32_t rref_z[1024];
    if(!stabilisers_rref_p(_x_stabilisers,n_x,n,rref_x)){
      lwlog_err("x stabilisers not in rref");
      return INVALID_INPUT;
    }
    if(!stabilisers_rref_p(_z_stabilisers,n_z,n,rref_z)){
      lwlog_err("z stabilisers not in rref");
      return INVALID_INPUT;
    }
    dcp_t* dcp_stabilisers = malloc(k * sizeof(dcp_t));
    dcp_t* output_errors = malloc(output_error_count * sizeof(dcp_t));
    pauli_t* x_measurement_paulis = malloc(x_measurement_count * sizeof(pauli_t));
    pauli_t* z_measurement_paulis = malloc(z_measurement_count * sizeof(pauli_t));
    pauli_t* x_stabilisers = malloc(n_x * sizeof(pauli_t));
    pauli_t* z_stabilisers = malloc(n_z * sizeof(pauli_t));
    mzd_t* x_permutation = mzd_init(n,n);
    mzd_t* z_permutation = mzd_init(n,n);

    lwlog_info("converting objects");
    for(int i = 0; i < n_x; i++){
      lwlog_debug("converting x stabiliser %d",i);
      convert_pauli_x(&(x_stabilisers[i]),&(_x_stabilisers[i]));
      //print_pauli(&(x_stabilisers[i]),1);
    }
    for(int i = 0; i < n_z; i++){
      lwlog_debug("converting z stabiliser %d",i);
      convert_pauli_z(&(z_stabilisers[i]),&(_z_stabilisers[i]));
      //print_pauli(&(z_stabilisers[i]),1);
    }
    for(int i = 0; i < k; i++){
      lwlog_debug("converting logical stabiliser %d",i);
      convert_dcp(&(dcp_stabilisers[i]),&(_dcp_stabilisers[i]));
      lwlog_debug("logical stabiliser %d: %s",i,print_dcp(&(dcp_stabilisers[i]),0));
      //print_dcp(&(dcp_stabilisers[i]),1);
    }
    lwlog_debug("converting %d output errors",output_error_count);
    for(int i = 0; i < output_error_count; i++){
      convert_dcp(&(output_errors[i]),&(_output_errors[i]));
      //printf("output error %d: ",i);
      //print_dcp(&(output_errors[i]),1);
    }
    for(int i = 0; i < x_measurement_count; i++){
      lwlog_debug("converting x measurement %d",i);
      convert_pauli_x(&(x_measurement_paulis[i]),&(_x_measurement_paulis[i]));
      //print_pauli(&(x_measurement_paulis[i]),1);
    }
    for(int i = 0; i < z_measurement_count; i++){
      lwlog_debug("converting z measurement %d",i);
      convert_pauli_z(&(z_measurement_paulis[i]),&(_z_measurement_paulis[i]));
      //print_pauli(&(z_measurement_paulis[i]),1);
    }
    lwlog_debug("x lookup table: ");
    for(int i = 0; i < (1 << (n-n_x)); i++){
      //lwlog_debug("%d",TestBit(x_lookup_table,i)?1:0);
        }
    lwlog_debug("z lookup table: ");
    for(int i = 0; i < (1 << (n-n_z)); i++){
      //lwlog_debug("%d",TestBit(z_lookup_table,i)?1:0);
    }
    convert_permutation(x_permutation,n_x,n,rref_x);
    convert_permutation(z_permutation,n_z,n,rref_z);
    //mzp_print(x_permutation);
    //mzp_print(z_permutation);
    for(int i = 0; i < fault_location_count ; i++){
      //lwlog_debug("fault location partition %d: %d\n",i,fault_location_partitions[i]);
    }
    lwlog_debug("all objects converted");

    // Allocate space for multiplying dcp errors together in the iteration
    dcp_t current_dcp_output_combination;
    dcp_t temp;
    dcp_t error_with_stabilisers;
    init_dcp_to_size(&temp,n_all);
    init_dcp_to_size(&current_dcp_output_combination,n_all);
    init_dcp_to_size(&error_with_stabilisers,n_all);

    // Construct the sequence of measurements that postselect and canonicalise all dcp errors
    int measurement_sequence_count = (x_measurement_count) + (2*z_measurement_count) + (2*n_z);
    pauli_t* measurement_sequence = malloc(measurement_sequence_count*sizeof(pauli_t));
    for(int i = 0; i < x_measurement_count; i++){
      // postselect on x measurements
      measurement_sequence[i]=x_measurement_paulis[i];
    }
    for(int i = 0; i < z_measurement_count; i++){
      // postselect on z measurements
      int index=0;
      for(int j = 0; j < n_all; j++){
        if(mzd_read_bit(z_measurement_paulis[i].data,j,0))index=j;
      }
      measurement_sequence[x_measurement_count + (2*i)]     = z_measurement_paulis[i];
      //measurement_sequence[x_measurement_count + (2*i) + 1] = pauli_x(index,n_all);
      init_to_pauli_x(&(measurement_sequence[x_measurement_count + (2*i) + 1]),index,n_all);
    }
    for(int i = 0; i < n_z; i++){
      int index = rref_z[i];

      init_to_pauli(&(measurement_sequence[x_measurement_count + 2*z_measurement_count + (2*i)]),n_all);
      for(int j = 0; j < PAULI_N(&(z_stabilisers[i])); j++){
        mzd_write_bit(measurement_sequence[x_measurement_count + 2*z_measurement_count + (2*i)].data,j,0,mzd_read_bit(z_stabilisers[i].data,j,0));
        }
      init_to_pauli_x(&(measurement_sequence[x_measurement_count + 2*z_measurement_count + (2*i) + 1]),index,n_all);
    }
    // TODO construct the sequence
    // Check that the number of faults is not too large
    uint32_t fault_count = 0;
    for(int i = 0;i<n; i++){
      fault_count+=fault_location_partitions[i];
    }
    if (fault_location_count > MAX_SIZE_COUNT) {
      printf("problem too large: %i > %i\n", fault_location_count, MAX_SIZE_COUNT);
      return PROBLEM_TOO_LARGE;
    }

    lwlog_debug("setting up the iteration");
    int current_combination[MAX_SIZE_COUNT];
    int first_balls[MAX_SIZE_COUNT];
    int last_balls[MAX_SIZE_COUNT];
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
    // do the iteration
    unsigned long long c = 0;
    do {
      lwlog_debug("checking fault number %d\n",fault_counter++);
      // Multiply DCP errors together. TODO verify that order is correct
      set_dcp_to(&current_dcp_output_combination,&(output_errors[current_combination[0]]));
      for(int i = 1 ; i < s ; i++){
        dcp_mul(&temp,&current_dcp_output_combination,&(output_errors[current_combination[i]]));
        set_dcp_to(&current_dcp_output_combination,&temp);
      }
      // output errors multiplied. now, multiply by all combinations of logical operators and test them:
      // TODO they all commute (i think) so a gray code should be fine here really
      // loop over all combinations of dcp stabilisers
      int good_combination_found=0;
      for(uint32_t dcp_index=0;dcp_index < (1 << k)  && !good_combination_found;
          dcp_index++){
        set_dcp_to(&error_with_stabilisers,&current_dcp_output_combination);
        //
        for(int i = 0 ; i < k && !good_combination_found; i++){
          if(dcp_index & (1 << i)){
            // multiply the output error by dcp stabilisers
            lwlog_info("multiplying error by logical stabiliser");
            dcp_mul(&temp,&error_with_stabilisers,&(dcp_stabilisers[i]));
            set_dcp_to(&error_with_stabilisers,&temp);
          }
          // (dcp_t *error, pauli_t* p, int cozunt, int n, uint32_t* z_permutation, comb_bit_array_t z_lookup_table)
        }
        // now check to see if it is good or not.
        if(!dcp_fault_bad_p_css(&error_with_stabilisers,measurement_sequence,measurement_sequence_count,n,z_permutation,z_lookup_table,x_permutation,x_lookup_table,x_stabilisers,n_x,rref_x,z_measurement_count,z_measurement_paulis,n_all)){
          // a good one was found after multiplying by a choice of logicals, set the flag
          good_combination_found=1;
          //break outer;
        }
        //error_with_stabilisers=0; //error_with_stabilisers is now bad after the destructive operation
      }
      if (!good_combination_found) {
        lwlog_debug("no good combination was found: fault combination was bad");
        // need to return a fault
        for (int i = 0; i < s; i++) {
          bad_fault_array[i] = current_combination[i];
        }
        fflush(stdout);
        clock_gettime( CLOCK_MONOTONIC, &finish );
        double seconds = (finish.tv_sec - start.tv_sec) + 1.0e-9 * (finish.tv_nsec - start.tv_nsec);
        printf("Took %f seconds.\n", seconds);
        return BAD_FAULT_FOUND;
      }

    c++;
  } while (ft_iteration_combination_incf(s, fault_location_partitions, fault_location_count, first_balls, last_balls,
                                         current_combination));
  lwlog_info("%llu faults processed.\n", c);
  // Stop timing
  clock_gettime( CLOCK_MONOTONIC, &finish );
  double seconds = (finish.tv_sec - start.tv_sec) + 1.0e-9 * (finish.tv_nsec - start.tv_nsec);
  printf("Took %f seconds.\n", seconds);
  //lwlog_debug("TODO freeing everything");


  return ALL_FAULTS_GOOD;
  }

#ifdef __cplusplus
}
#endif
