#include "comb.h"

int main(int argc, char *argv[]){
  #ifdef COMB_MPI
  int provided;
  MPI_Init_thread(&argc, &argv,MPI_THREAD_MULTIPLE,&provided);
  if(provided!=MPI_THREAD_MULTIPLE){
    printf("MPI_THREAD_MULTIPLE not supported");
    MPI_Abort(MPI_COMM_WORLD,INTERNAL_ERROR);
  }
  #endif
  if(argc <= 1){
    printf("Usage: comb input-file\n");
     #ifdef COMB_MPI
    MPI_Finalize();
    #endif
    return 1;
  }
  
  int result = deserialise_run_from_file(argv[1]);
  #ifdef COMB_MPI
    MPI_Finalize();
    #endif
  return result;
}
