#include "comb.h"


int main(int argc, char *argv[]){
  if(argc <= 2){
    printf("Usage: comb input-file thread-count\n");
    return 1;
  }
  uint32_t thread_count = (uint32_t)strtol(argv[2], NULL, 10);
  return deserialise_run_from_file(argv[1],thread_count);
}
