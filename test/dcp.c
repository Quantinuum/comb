#include "lwlog.h"
#include "comb.h"

int main(int argv, char** args){
//  init_temps();
  dcp_t* dcp = init_dcp(3);
  dcp->F = mzd_from_str(3,3,args[1]);
  dcp->h_x = mzd_from_str(3,1,args[2]);
  dcp->h_z = mzd_from_str(3,1,args[3]);
  fix_dcp(dcp);
  fprintf(stderr,"%s",print_dcp(dcp,1));

  mzd_print(dcp->F);
  printf("\n");
  mzd_print(dcp->h_x);
  printf("\n");
  mzd_print(dcp->h_x_transpose);
  printf("\n");
  mzd_print(dcp->h_z);
  printf("\n");
  mzd_print(dcp->F_diagonal);
  printf("\n");

  /* mzd_t *A = mzd_from_str(1,5,args[1]); */
  /* mzd_print(A); */
  /* printf("row stride: %ld\nhigh bitmask: %ld\n",A->rowstride,A->high_bitmask); */
  /* for(int i = 0; i < A->rowstride * A->nrows; i+= A->rowstride){ */
  /*   printf("%d: %ld\n",i,A->data[i]); */
  /* } */
  return 0;
}
