#include "lwlog.h"
#include "comb.h"

int main(int argv, char** args){
  init_temps(3);
  dcp_t* dcp = init_dcp(3);
  dcp->F = mzd_from_str(3,3,args[1]);
  dcp->h_x = mzd_from_str(3,1,args[2]);
  dcp->h_z = mzd_from_str(3,1,args[3]);
  fix_dcp(dcp);

  stabiliser_state_t* state = init_state(3);
  stabiliser_state_t* new_state = init_state(3);
  pauli_z_coset_t* coset= init_pauli_z_coset(3,3);
  set_stabiliser_state_to_dcp(state,dcp);
  mzp_t *p = mzp_init(3);
  mzp_t *q = mzp_init(3);
  set_coset_to_stabiliser_state(coset,state,temp_n_n_1,temp_n_n_2,p,q,temp_n_n_plus_one,temp_n_plus_one_n,new_state);
  for(int i = 0; i < 1 << COSET_K(coset); i++){
    fprintf(stderr,"%s",print_pauli_z_coset(coset,1));
    incf_coset(coset,i);
  }
  mzd_print(coset->factor_transpose);
  printf("\n");
  mzd_print(coset->set_generators);
  printf("\n");

  /* mzd_t *A = mzd_from_str(1,5,args[1]); */
  /* mzd_print(A); */
  /* printf("row stride: %ld\nhigh bitmask: %ld\n",A->rowstride,A->high_bitmask); */
  /* for(int i = 0; i < A->rowstride * A->nrows; i+= A->rowstride){ */
  /*   printf("%d: %ld\n",i,A->data[i]); */
  /* } */
  return 0;
}
