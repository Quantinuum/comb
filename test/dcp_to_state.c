#include "lwlog.h"
#include "comb.h"

int main(int argv, char** args){
//  init_temps();
  dcp_t* dcp = init_dcp(3);
  dcp->F = mzd_from_str(3,3,args[1]);
  dcp->h_x = mzd_from_str(3,1,args[2]);
  dcp->h_z = mzd_from_str(3,1,args[3]);
  fix_dcp(dcp);

  stabiliser_state_t* state = init_state(3);
  set_stabiliser_state_to_dcp(state,dcp);
  fprintf(stderr,"%s",print_stabiliser_state(state,1));
  mzd_print(state->x_parts);
  printf("\n");
  mzd_print(state->z_parts);
  printf("\n");
  mzd_print(state->x_parts);
  printf("\n");
  mzd_print(state->delta_vector);
  printf("\n");
  mzd_print(state->epsilon_vector);
  printf("\n");

  /* mzd_t *A = mzd_from_str(1,5,args[1]); */
  /* mzd_print(A); */
  /* printf("row stride: %ld\nhigh bitmask: %ld\n",A->rowstride,A->high_bitmask); */
  /* for(int i = 0; i < A->rowstride * A->nrows; i+= A->rowstride){ */
  /*   printf("%d: %ld\n",i,A->data[i]); */
  /* } */
  return 0;
}
