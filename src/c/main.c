#include "comb.h"


int main(int argc, char *argv[]){
  if(argc <= 1){
    printf("Usage: comb input_file\n");
    return 1;
  }
  return deserialise_run_from_file(argv[1]);
}
