#include "lwlog.h"
#include "comb.h"

int main(int argv, char** args){
  mzd_t *A = mzd_from_str(1,5, args[1]);
  mzd_print(A);
  printf("row stride: %ld\nhigh bitmask: %ld\n",A->rowstride,A->high_bitmask);
  for(int i = 0; i < A->rowstride * A->nrows; i+= A->rowstride){
    printf("%d: %ld\n",i,A->data[i]);
  }
  return 0;
}
