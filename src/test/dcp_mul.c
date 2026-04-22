#include "lwlog.h"
#include "comb.h"

int main(int argv, char** args){
  init_temps();
  dcp_t* dcp1 = init_dcp(3);
  dcp1->F = mzd_from_str(3,3,args[1]);
  dcp1->h_x = mzd_from_str(3,1,args[2]);
  dcp1->h_z = mzd_from_str(3,1,args[3]);


  dcp_t* dcp2 = init_dcp(3);
  dcp2->F = mzd_from_str(3,3,args[4]);
  dcp2->h_x = mzd_from_str(3,1,args[5]);
  dcp2->h_z = mzd_from_str(3,1,args[6]);

  dcp_t* dcp3 = init_dcp(3);
  fix_dcp(dcp1);
  fix_dcp(dcp2);
  dcp_mul(dcp3,dcp1,dcp2);
  fprintf(stderr,"%s",print_dcp(dcp1,1));
  fprintf(stderr,"%s",print_dcp(dcp2,1));
  fprintf(stderr,"%s",print_dcp(dcp3,1));

  mzd_print(dcp3->F);
  printf("\n");
  mzd_print(dcp3->h_x);
  printf("\n");
  mzd_print(dcp3->h_x_transpose);
  printf("\n");
  mzd_print(dcp3->h_z);
  printf("\n");
  mzd_print(dcp3->F_diagonal);
  printf("\n");

  /* mzd_t *A = mzd_from_str(1,5,args[1]); */
  /* mzd_print(A); */
  /* printf("row stride: %ld\nhigh bitmask: %ld\n",A->rowstride,A->high_bitmask); */
  /* for(int i = 0; i < A->rowstride * A->nrows; i+= A->rowstride){ */
  /*   printf("%d: %ld\n",i,A->data[i]); */
  /* } */
  return 0;
}
