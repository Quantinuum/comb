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


  void  permute_qubits(pauli_z_coset_t* coset,int n, uint32_t* c){
    lwlog_debug("in permute qubits");
    mzp_t* perm = mzp_init(2*n+1);
    for(int i =0; i < n; i++)
      perm->values[i]=c[i]; // only need to permute zs
    mzd_apply_p_left(coset->factor->data,perm);
    for(int i =0 ; i< coset->k; i++)
      mzd_apply_p_left(coset->set_generators[i]->data,perm);
  }

  void set_to_pauli(pauli_t* a, pauli_t* b){
    lwlog_debug("in set_to_pauli");
    print_pauli(a,1);
    print_pauli(b,1);
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

  int dcp_fault_bad_p_css(dcp_t *error, pauli_t* p, int count, int n, uint32_t* z_permutation, comb_bit_array_t z_lookup_table, uint32_t* x_permutation, comb_bit_array_t x_lookup_table,pauli_t* x_stabilisers, uint32_t n_x){
    lwlog_debug("in dcp_fault_bad_p_css:");
    print_dcp(error,1);
    // TODO First, check the X part of the error, as this is easier to deal with.
    // Does it anticommute with any Z basis postselections? If so, the error is always detected, so the entire error is good
    // Is it low weight? If so, it is good, but we must check the rest of the error. If it is not low weight, then the whole error is bad
    // Now we do the more expensive Z check.
    // Convert the dcp into a stabiliser state:
    stabiliser_state_t* state = dcp_to_stabiliser_state(error);
    // Do all postselections and canonicalisations by measuring all the Paulis on the state.
    // If the state ever vanishes, the error is good, so return.
    for(int i =0 ; i< count; i++){
      lwlog_debug("postselecting %d",i);
      if(postselect_on_measuring_plus_one(state,&(p[i]))){
        lwlog_debug("postselecting %d vanished, so return good",i);
      return 0; // vanishes so return good
      }
      lwlog_debug("finished postselecting %d",i);
    }
    // State now postselected and canonicalised. Expand the resulting state as a coset
    lwlog_debug("before diagonal clifford to coset");
    pauli_z_coset_t *coset = diagonal_clifford_to_coset(state);
    lwlog_debug("after diagonal clifford to coset");
    print_pauli_z_coset(coset,1);
    // Apply the qubit permutation:
    permute_qubits(coset,n,z_permutation);

    lwlog_debug("qubits permuted");
    pauli_t current_element;
    pauli_t temp;
    init_pauli_to_size(&current_element,PAULI_N(coset->factor));
    set_to_pauli(&current_element,coset->factor);
    init_pauli_to_size(&temp,PAULI_N(coset->factor));
    lwlog_debug("n is %d",PAULI_N(coset->factor));
    // Iterate over every element in the coset using a gray code
    for(int i = 0 ; i < (1 << coset->k); i++){
      lwlog_debug("iterating %d in gray code",i);
      // If any term is bad, the error is bad, otherwise it is good.
      print_pauli(&current_element,1);
      if(!TestBit(z_lookup_table,mzd_data(temp.data))) {
        //it's bad so return, otherwise check the other ones
        lwlog_debug("returning 1");
        return 1;
       }
    

      if((i+1)!=(1<<coset->k)){
        pauli_mul(&temp,&current_element,coset->set_generators[gray_bit_to_switch(i)]);
        set_to_pauli(&current_element,&temp);
      }
    }
    lwlog_debug("finished loop");
    return 0;
    mzd_free(temp.data);
    mzd_free(current_element.data);
    lwlog_debug("freed data");
  // all elements in the coset were good so return
    return 0;
  }


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
    lwlog_info("allocating and converting objects");
    dcp_t* dcp_stabilisers = malloc(k * sizeof(dcp_t));
    dcp_t* output_errors = malloc(output_error_count * sizeof(dcp_t));
    pauli_t* x_measurement_paulis = malloc(x_measurement_count * sizeof(pauli_t));
    pauli_t* z_measurement_paulis = malloc(z_measurement_count * sizeof(pauli_t));
    pauli_t* x_stabilisers = malloc(n_x * sizeof(pauli_t));
    pauli_t* z_stabilisers = malloc(n_z * sizeof(pauli_t));
    mzp_t* x_permutation = mzp_init(n);
    mzp_t* z_permutation = mzp_init(n);

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
      //printf("logical stabiliser %d: ",i);
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
      print_pauli(&(z_measurement_paulis[i]),1);
    }
    lwlog_debug("x lookup table: ");
    for(int i = 0; i < (1 << (n-n_x)); i++){
      //lwlog_debug("%d",TestBit(x_lookup_table,i)?1:0);
        }
    lwlog_debug("z lookup table: ");
    for(int i = 0; i < (1 << (n-n_z)); i++){
      //lwlog_debug("%d",TestBit(z_lookup_table,i)?1:0);
    }
    convert_permutation(x_permutation,_x_permutation);
    convert_permutation(z_permutation,_z_permutation);
    mzp_print(x_permutation);
    mzp_print(z_permutation);
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
      int index = -1;
      for(int j = n-n_z; j < n && (index== -1); j++){
        int candidate_index=_z_permutation->indices[j];
        if(mzd_read_bit(z_stabilisers[i].data,candidate_index,0)){
          index=candidate_index;
        }
      }

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
      for(int i = 1 ; i < k ; i++){
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
            dcp_mul(&temp,&error_with_stabilisers,&(dcp_stabilisers[i]));
            set_dcp_to(&error_with_stabilisers,&temp);
          }
          // (dcp_t *error, pauli_t* p, int count, int n, uint32_t* z_permutation, comb_bit_array_t z_lookup_table)
        }
        // now check to see if it is good or not.
        if(!dcp_fault_bad_p_css(&error_with_stabilisers,measurement_sequence,measurement_sequence_count,n,_z_permutation->indices,z_lookup_table,_x_permutation->indices,x_lookup_table,x_stabilisers,n_x)){
          // a good one was found after multiplying by a choice of logicals, set the flag
          good_combination_found=1;
          //break outer;
        }
        //error_with_stabilisers=0; //error_with_stabilisers is now bad after the destructive operation
      }
      if (!good_combination_found) {
        lwlog_debug("no good combination was found: fault combination was bad");
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
  } while (ft_iteration_combination_incf(k, fault_location_partitions, fault_location_count, first_balls, last_balls,
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
