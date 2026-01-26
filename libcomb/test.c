#include "cnot_plus_t_iteration.h"
#include <stdio.h>

ctc_pauli_x_t test_x(){
  ctc_pauli_x_t a;
  a.pauli_length=8;
  a.data[0]=0b00100101;
  return a;
}

ctc_pauli_z_t test_z(){
  ctc_pauli_z_t a;
  a.pauli_length=64;
  a.data[0]=0b00100101;
  a.data[1]=0b11100111;
  return a;
}

ctc_pauli_z_t test_z2(){
  ctc_pauli_z_t a;
  a.pauli_length=6;
  a.data[0]=0b111111;
  return a;
}
ctc_pauli_z_t test_z3(){
  ctc_pauli_z_t a;
  a.pauli_length=6;
  a.data[0]=0b101011;
  return a;
}

ctc_dcp_t test_dcp(){
  ctc_dcp_t a;
  a.pauli_x_part=test_x();
  a.diagonal_operators_count=2;
  a.diagonal_operators[0]=test_z2();
  a.diagonal_operators[1]=test_z3();
  a.coefficients[0]=-1;
  a.coefficients[1]=1;
  return a;
}
void test_iterate(){
  uint32_t n = 3;
  int sizes[]={1,1,1};
  uint32_t k = 1;
  ctc_dcp_t propagated_output_errors[3];
  uint32_t bad_fault_array[3];
  uint32_t logical_clifford_count=0;
  ctc_dcp_t logical_cliffords[0];
  ctc_ft_iteration_x_evaluator x_evaluator;
  ctc_ft_iteration_z_evaluator z_evaluator;
  // ctc_ft_iteration_iterate_css(n,sizes,k,propagated_output_errors,bad_fault_array,logical_clifford_count
  //                           ,logical_cliffords,x_evaluator,z_evaluator);
}

int main(){
  ctc_pauli_x_t a = test_x();
  ctc_pauli_z_t b = test_z();
  ctc_dcp_t c = test_dcp();
  ctc_print_pauli_x(&a);
  printf("\n");
  ctc_print_pauli_z(&b);
  printf("\n");
  ctc_print_dcp(&c);
  printf("\n");
  //test_iterate();
  return 0;
}
