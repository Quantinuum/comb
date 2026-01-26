#include "cnot_plus_t_iteration.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h> 
#ifdef __cplusplus
extern "C" {
#endif

  void ctc_print_pauli_x(ctc_pauli_x_t* operator){
    //printf("#G(");
    for(int i =0;i<(operator->pauli_length<128?operator->pauli_length:128);i++){
      printf("%c",TestBit(operator->data,i)?'X':'I');
    }
    //printf(")");
  }



  void ctc_print_pauli_z(ctc_pauli_z_t* operator){
    //printf("#G(");
    for(int i =0;i<(operator->pauli_length<128?operator->pauli_length:128);i++){
      printf("%c",TestBit(operator->data,i)?'Z':'I');
    }
    //printf(")");
  }
  void ctc_print_dcp(ctc_dcp_t* operator){
    printf("exp[ i pi/4(");
    for(int i = 0; i < operator->diagonal_operators_count; i++){
      printf("%d.",operator->coefficients[i]);
      ctc_print_pauli_z(&(operator->diagonal_operators[i]));
      if(i < (operator->diagonal_operators_count - 1))printf("+");
    }
    printf(")]");
      ctc_print_pauli_x(&(operator->pauli_x_part));
  }


  int ctc_pauli_anticommute_p(ctc_pauli_x_t* p1, ctc_pauli_z_t* p2){
    uint result = 0;
    for(int i=0;i<p1->pauli_length;i++){
      if(TestBit(p1->data,i) && TestBit(p2->data,i))result++;
    }
    return (result & 1);
  }

  void ctc_store_pauli_x_product(ctc_pauli_x_t* p1, ctc_pauli_x_t* p2, ctc_pauli_x_t* result){
    result->pauli_length = p1->pauli_length;// + p2->pauli_length;
    for(int i=0;i<MAX_PAULI_SIZE;i++){
      result->data[i]=p1->data[i]^p2->data[i];
    }
  }

  void ctc_store_dcp_mul(ctc_dcp_t* p1, ctc_dcp_t* p2, ctc_dcp_t* result){
    // Multiplies p1 and p2 together and stores the result.
    // Does not do cancellations on the sum in result - this must be done separately.
    result->diagonal_operators_count = p1->diagonal_operators_count + p2->diagonal_operators_count;
    for(int i = 0; i< p1->diagonal_operators_count; i++){
      result->diagonal_operators[i]=p1->diagonal_operators[i];
      result->coefficients[i]=p1->coefficients[i];
    }
    for(int i = 0; i< p2->diagonal_operators_count; i++){
      result->diagonal_operators[i+p1->diagonal_operators_count]=p2->diagonal_operators[i];
      if(ctc_pauli_anticommute_p(&(p1->pauli_x_part),&(p2->diagonal_operators[i]))){
        //printf("anticommutes\n");
      result->coefficients[i + p1->diagonal_operators_count]=-p2->coefficients[i];
      }
      else
       result->coefficients[i + p1->diagonal_operators_count]=p2->coefficients[i];
    }
    ctc_store_pauli_x_product(&(p1->pauli_x_part),&(p2->pauli_x_part),&(result->pauli_x_part));
  }
  int x_pauli_detected_by_z_measurements(ctc_pauli_x_t* pauli,
                                     uint32_t z_measurement_count,
                                     uint32_t* z_measurement_indices){
    int x_detected_p=0;
    for(int i=0;i<z_measurement_count;i++){
      //printf("seeing if z[%d] is measured\n",z_measurement_indices[i]);
      if(TestBit(pauli->data, z_measurement_indices[i]))x_detected_p=1;
    }
    return x_detected_p;
  }
  int z_pauli_detected_by_x_measurements(ctc_pauli_z_t* pauli, uint32_t x_measurement_count,
                                         uint32_t* x_measurement_indices){
    int z_detected_p=0;
    for(int i=0;i<x_measurement_count;i++){
      //printf("seeing if x[%d] is measured\n",x_measurement_indices[i]);
      if(TestBit(pauli->data, x_measurement_indices[i]))z_detected_p=1;
    }
    return z_detected_p;
  }
  ctc_pauli_z_sum_t* allocate_pauli_z_sum(uint32_t term_count){
    ctc_pauli_z_sum_t* sum = malloc(sizeof(ctc_pauli_z_sum_t));
    sum->term_count=0;
    sum->coefficients=malloc(sizeof(ctc_gaussian_integer_t)*term_count);
    sum->terms=malloc(sizeof(ctc_pauli_z_t)*term_count);
    return sum;
  }
  void free_pauli_z_sum(ctc_pauli_z_sum_t* sum){
    free(sum->terms);
    free(sum->coefficients);
    free(sum);
  }

  int pauli_z_equal_p(ctc_pauli_z_t* a, ctc_pauli_z_t* b){
    if(a->pauli_length != b->pauli_length)return 0;
    for(int i = 0;i<(a->pauli_length);i++){
      if(TestBit(a->data, i) != TestBit(b->data, i))return 0;
    }
    return 1;
  }

  void add_term_to_sum(ctc_pauli_z_sum_t* z_sum,ctc_pauli_z_t* current_term, ctc_gaussian_integer_t coefficient,
                       uint32_t x_measurement_count, uint32_t* x_measurement_indices,
                       uint32_t z_measurement_count, uint32_t* z_measurement_indices){
    ctc_pauli_z_t erased_term=*current_term;
    // first, erase the measurement bits
    for(int i = 0;i<x_measurement_count;i++){
      ClearBit(erased_term.data,x_measurement_indices[i]);
      erased_term.pauli_length=fmin(erased_term.pauli_length,x_measurement_indices[i]);
    }
    for(int i = 0;i<z_measurement_count;i++){
      ClearBit(erased_term.data,z_measurement_indices[i]);
      erased_term.pauli_length=fmin(erased_term.pauli_length,z_measurement_indices[i]);
    }
    // should do canonicalisation here
    for(int i = 0;i<z_sum->term_count; i++){
      if(pauli_z_equal_p(&erased_term, &(z_sum->terms[i]))){
        z_sum->coefficients[i].real+=coefficient.real;
        z_sum->coefficients[i].imaginary+=coefficient.imaginary;
          z_sum->terms[i]=erased_term;
          return;
      }
    }
    int i = z_sum->term_count++;
    z_sum->coefficients[i]=coefficient;
    z_sum->terms[i]=erased_term;
  }

  void set_mul_complex(ctc_gaussian_integer_t* a,
                                         ctc_gaussian_integer_t* b,
                                         ctc_gaussian_integer_t* result){
    result->real = (a->real* b->real) - (a->imaginary * b->imaginary);
    result->imaginary= (a->real * b->imaginary) + (a->imaginary * b->real);
  }

  void ctc_print_pauli_sum(ctc_pauli_z_sum_t* sum){
    for(int i = 0;i < sum->term_count; i++){
      if(sum->coefficients[i].real != 0 || sum->coefficients[i].imaginary !=0){
      printf("(%d + %d i)",sum->coefficients[i].real,sum->coefficients[i].imaginary);
      ctc_print_pauli_z(&(sum->terms[i]));
      if((i+1)!=sum->term_count)printf("+");
      }
    }
  }
#define CANONICALISE_DEBUG 0
  void canonicalise_pauli_x(ctc_pauli_x_t* operator, int k, ctc_pauli_x_t* canonicalisers){
    canonicalise_pauli_z((ctc_pauli_z_t*)operator, k, (ctc_pauli_z_t*)canonicalisers);
  }
  void canonicalise_pauli_z(ctc_pauli_z_t* operator, int k, ctc_pauli_z_t* canonicalisers){
    uint32_t n = canonicalisers[0].pauli_length;
    //printf("canonicalising ");
    //ctc_print_pauli_z(operator);
    //printf("\n");
    for(int i = 0; i < k ; i++){
      //printf("canonicaliser %d ",i);
      //ctc_print_pauli_z(&canonicalisers[i]);
      //printf("\n");
      int m = (n - k + i);
      //ctc_print_pauli_z(&canonicalisers[i]);
      //printf("\n");
#if CANONICALISE_DEBUG
      for(int j = n-k; j< n; j++){
        if(TestBit(canonicalisers[i].data,j)){
          if(j!=m)printf("canonicalisers wrong \n");
        }
        else{
          if(j == m )printf("canonicalisers wrong \n");
        }
      }
#endif
      if(TestBit(operator->data,m)){
        // operator has support on qubit m; eliminate by multiplying by the canonicaliser
        ctc_pauli_z_t result = {0};
        //printf("multiplying by ");
        //ctc_print_pauli_z(&(canonicalisers[i]));
        //printf("\n");
        ctc_store_pauli_z_mul(operator,&(canonicalisers[i]),&result);
        
        for(int s = 0 ;s < MAX_PAULI_SIZE; s++){
          operator->data[s]=result.data[s];
        }

        //printf("result is ");
        //ctc_print_pauli_z(operator);
        //printf("\n");
      }
    }
    #if CANONICALISE_DEBUG
    for(int i = n - k; i < n ; i++){
      if(TestBit(operator->data,i)){
        printf("operator not canonicalised correctly ");
        ctc_print_pauli_z(operator);
        printf("\n");
      }
    }
    #endif
    //for(int i = n; i<32*MAX_PAULI_SIZE;i++){ ClearBit(operator->data,i);} // think not necessary
  }

  void set_coefficient_to(ctc_gaussian_integer_t* identity_part, ctc_gaussian_integer_t*  matrix_part, int32_t c){
    // sets a to exp(i c pi/4)
    switch( c % 8){
    case 0:
    case 8:
    case -8:
      identity_part->real=1;identity_part->imaginary=0;
      matrix_part->real=0;matrix_part->imaginary=0;
      break;
    case 1:
    case -7:
      identity_part->real=1;identity_part->imaginary=0;
      matrix_part->real=0;matrix_part->imaginary=1;
      break;
    case 2:
    case -6:
      identity_part->real=0;identity_part->imaginary=0;
      matrix_part->real=0;matrix_part->imaginary=2;
      break;
    case 3:
    case -5:
      identity_part->real=-2;identity_part->imaginary=0;
      matrix_part->real=0;matrix_part->imaginary=2;
      break;
    case 4:
    case -4:
      identity_part->real=0;identity_part->imaginary=0;
      matrix_part->real=-4;matrix_part->imaginary=0;
      break;
    case 5:
    case -3:
      identity_part->real=-4;identity_part->imaginary=0;
      matrix_part->real=0;matrix_part->imaginary=-4;
      break;
    case 6:
    case -2:
      identity_part->real=0;identity_part->imaginary=0;
      matrix_part->real=0;matrix_part->imaginary=-8;
      break;
    case 7:
    case -1:
      identity_part->real=8;identity_part->imaginary=0;
      matrix_part->real=0;matrix_part->imaginary=-8;
      break;
    default:
      printf("invalid coefficient %d\n",c);
    }
  }
  void ctc_store_pauli_z_mul(ctc_pauli_z_t* p1, ctc_pauli_z_t* p2, ctc_pauli_z_t* result){
    result->pauli_length=p1->pauli_length;
    //if((p1->pauli_length) != (p2->pauli_length))printf("pauli lengths not equal!\n");
    for(int i =0;i<MAX_PAULI_SIZE;i++){
      result->data[i]=p1->data[i]^p2->data[i];
    }
  }

  int dcp_fault_bad_p_css(ctc_dcp_t* output_error_combination,
                          ctc_ft_iteration_x_evaluator x_evaluator,
                          ctc_ft_iteration_z_evaluator z_evaluator,
                          uint32_t x_measurement_count, uint32_t* x_measurement_indices,
                          uint32_t z_measurement_count, uint32_t* z_measurement_indices,
                          int s, ctc_pauli_x_t* canonicalisers){
    //
    // first see if the x part is undetected by a z measurement and high weight
    //printf("testing output error combination\n");
    //ctc_print_dcp(output_error_combination);
    //printf("\n");
    uint32_t data_length = 100000000;
    for(int i = 0;i<x_measurement_count;i++){
      data_length=fmin(data_length,x_measurement_indices[i]);
    }
    for(int i = 0;i<z_measurement_count;i++){
      data_length=fmin(data_length,z_measurement_indices[i]);
    }
    //printf("checking x part ");
    // ctc_print_pauli_x(&(output_error_combination->pauli_x_part));
    // printf("\n");
    if(!x_pauli_detected_by_z_measurements(&(output_error_combination->pauli_x_part),
                                           z_measurement_count, z_measurement_indices)){
      ctc_pauli_x_t trimmed_error=output_error_combination->pauli_x_part;
      trimmed_error.pauli_length=data_length;
      //printf("trimmed_error before canonicalisation: ");
      //ctc_print_pauli_x(&trimmed_error);
      //printf("\n");
      canonicalise_pauli_x(&trimmed_error,s,canonicalisers); // need to canonicalise the error before evaluating
      //printf("checking trimmed_error: ");
      //ctc_print_pauli_x(&trimmed_error);
      //printf("\n");
      if(x_evaluator(&trimmed_error)){
        //printf("x part of ");
        //ctc_print_dcp(output_error_combination);
        //printf(" was undetected and high weight \n");
        return 1;
      }
      // x is undetected and low weight, so check z terms
    uint32_t pauli_length = output_error_combination->diagonal_operators[0].pauli_length;
    uint32_t sum_term_count=(1 << output_error_combination->diagonal_operators_count);
    // now need only worry about the z parts.
    ctc_pauli_z_sum_t* pauli_z_sum = allocate_pauli_z_sum(sum_term_count);
    // iterates over each term in the sum, checks if its detected, if not adds it to a new sum.
    for(int i = 0;i<sum_term_count;i++){
      ctc_pauli_z_t current_term_in_sum={0};
      ctc_pauli_z_t temp={0};
      current_term_in_sum.pauli_length=pauli_length;
      ctc_gaussian_integer_t current_coefficient={1,0};
      // set the current term
      for(int j = 0;j<output_error_combination->diagonal_operators_count;j++){
        ctc_gaussian_integer_t identity_coefficient;
        ctc_gaussian_integer_t matrix_coefficient;
        ctc_gaussian_integer_t temp_complex_result;
        set_coefficient_to(&identity_coefficient,&matrix_coefficient,output_error_combination->coefficients[j]);
        if(i & (1 << j)){
          ctc_store_pauli_z_mul(&current_term_in_sum,&(output_error_combination->diagonal_operators[j]),&temp);
          current_term_in_sum=temp;
          set_mul_complex(&current_coefficient,
                          &matrix_coefficient,
                          &temp_complex_result);
          current_coefficient=temp_complex_result;
        }
        else{
          set_mul_complex(&current_coefficient,
                          &identity_coefficient,
                          &temp_complex_result);
          current_coefficient=temp_complex_result;
        }
      } // canonicalise the z pauli before adding it to the sum in order to ensure cancellations
      canonicalise_pauli_z(&current_term_in_sum,s,(ctc_pauli_z_t*)canonicalisers);
      if(!z_pauli_detected_by_x_measurements(&current_term_in_sum,x_measurement_count,x_measurement_indices)){
        add_term_to_sum(pauli_z_sum,&current_term_in_sum,current_coefficient,
                        x_measurement_count,x_measurement_indices,z_measurement_count,z_measurement_indices);
      }

    }
    //printf("checking pauli sum: ");
    //ctc_print_pauli_sum(pauli_z_sum);
    //printf("\n");
    for(int i = 0;i<pauli_z_sum->term_count; i++){
      if(pauli_z_sum->coefficients[i].real != 0 || pauli_z_sum->coefficients[i].imaginary != 0 ){
        if(z_evaluator(&(pauli_z_sum->terms[i]))){
          //printf("z term ");
          //ctc_print_pauli_z(&(pauli_z_sum->terms[i]));
          //printf(" in the sum was found to be bad\n");
          free_pauli_z_sum(pauli_z_sum);
          return 1;
        }
      }
    }
    //printf("returning \n");
    free_pauli_z_sum(pauli_z_sum);
    //printf("sum freed\n");
    // every term in the sum was good so return 0
    return 0;
    // after all terms are added to the sum, the new sum is canonicalised i.e. sorted in place and each term in turn is checked to see if it passes the evaluator or not.
    // returns 1 if bad, 0 if not
    //printf("checking output error combination: ");
    //ctc_print_dcp(output_error_combination);
    //printf("\n");
    //if(x_evaluator(output_error_combination->pauli_x_part.data))return 1; // the x part was bad
    // TODO
    // iterate over sum with gray code? or an approach thats easier to implement correctly
    // gray code definitely does the terms, but not the coefficients
    // coefficients from a popcnt of the gray code?
    // move to unary encoding for dcp fault. (what does this mean?)
    //---
    // alternative:
    // TODO analyse these phases more. identities?


    return 0;
    }
    else{
      // x part was detected, so fault must be good
      return 0;
    }
  }
  int ctc_ft_iteration_iterate_css(uint32_t n, int *sizes, uint32_t k,
                                   uint32_t output_error_count,
                                   ctc_dcp_t* propagated_output_errors,
                                   uint32_t *bad_fault_array,
                                   uint32_t logical_clifford_count, ctc_dcp_t* logical_cliffords,
                                   ctc_ft_iteration_x_evaluator x_evaluator,
                                   ctc_ft_iteration_z_evaluator z_evaluator,
                                   uint32_t x_measurement_count, uint32_t* x_measurement_indices,
                                   uint32_t z_measurement_count, uint32_t* z_measurement_indices,
                                   int canonicaliser_count, ctc_pauli_x_t* canonicalisers) {
    int fault_counter=0;
    uint32_t fault_count = 0;
    struct timespec start, finish;
    clock_gettime( CLOCK_MONOTONIC, &start );
    for(int i = 0;i<n; i++){
      fault_count+=sizes[i];
    }
  if (n > MAX_SIZE_COUNT) {
    printf("problem too large: %i > %i\n", n, MAX_SIZE_COUNT);
    return PROBLEM_TOO_LARGE;
  }
  for(int i = 0; i < output_error_count; i++){
    //printf("%d/%d: ",i,output_error_count);
    //ctc_print_dcp(&(propagated_output_errors[i]));
    //printf("\n");
  }
  int current_combination[MAX_SIZE_COUNT]={0};
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
  unsigned long long c = 0;
  ctc_dcp_t current_output_error;
  ctc_dcp_t temp;
  do {
    if(fault_counter%10000 == 0)printf("checking fault number %d\n",fault_counter++);
    // print_array(k,current_combination);
    current_output_error=propagated_output_errors[current_combination[0]];
    for (int b = 1; b < k; b++) {
      //current_fault ^= propagated_faults[current_combination[b]];
      ctc_store_dcp_mul(&(propagated_output_errors[current_combination[b]]),&current_output_error,&temp);
      current_output_error=temp;
    }
    //printf("current output error: ");
    //ctc_print_dcp(&current_output_error);
    //printf("\n");
    int good_combination_found=0;

    // loop over multiplying by dcp errors
    for(uint32_t dcp_index=0;dcp_index < (1 << logical_clifford_count)  && !good_combination_found;
             dcp_index++){
      ctc_dcp_t multiplied_output_error=current_output_error;
      for(int i = 0 ; i < logical_clifford_count && !good_combination_found; i++){
        if(dcp_index & (1 << i)){
          //printf("before multiplying: ");
          //ctc_print_dcp( &multiplied_output_error);
          //printf("\n");
          //printf("multiplying by: ");
          //ctc_print_dcp( &(logical_cliffords[i]));
          //printf("\n");
          ctc_store_dcp_mul(&multiplied_output_error, &(logical_cliffords[i]), &temp);
          multiplied_output_error=temp;
          //printf("multiplied output error: ");
          //ctc_print_dcp(&multiplied_output_error);
          //printf("\n");
        }
      }
      //printf("multiplied output error: ");
      //ctc_print_dcp(&multiplied_output_error);
      //printf("\n");
      // check the multiplied error
      if(!dcp_fault_bad_p_css(&multiplied_output_error,x_evaluator,z_evaluator,x_measurement_count, x_measurement_indices, z_measurement_count, z_measurement_indices,canonicaliser_count,canonicalisers)){
        // a good one was found after multiplying by a choice of logicals
        good_combination_found=1;
        //break outer;
      }
    }


    if (!good_combination_found) {


      //printf("bad fault found: ");
      //print_array(k, current_combination);
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
  clock_gettime( CLOCK_MONOTONIC, &finish );
  double seconds = (finish.tv_sec - start.tv_sec) + 1.0e-9 * (finish.tv_nsec - start.tv_nsec);
  printf("Took %f seconds.\n", seconds);
  return ALL_FAULTS_GOOD;
}
  // gray codes (from wikipedia)
  uint32_t integer_to_gray(uint32_t n){
    return n ^ (n >> 1);
  }
  uint32_t gray_to_integer(uint32_t num){
    uint32_t mask = num;
    while(mask){
      mask = mask >> 1;
      num = num ^ mask;
    }
    return num;
  }

  uint32_t gray_bit_to_switch(uint32_t n){
    uint32_t n_add_one = n + 1;
    uint32_t diff = n ^ (n + 1) ^ (n >> 1) ^ (n_add_one >> 1);
    return n ? 31 - __builtin_clz(diff) : 0;
  }

#ifdef __cplusplus
}
#endif
