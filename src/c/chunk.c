#include "comb.h"

void init_chunk(chunk_t* chunk, unsigned long long faults_to_check_count, int* starting_combination, int s){
  chunk->faults_already_checked=0;
  for(int i = 0; i < s; i++){
    chunk->current_combination[i]=starting_combination[i];
  }
  chunk->s=s;
}

void describe_chunks(chunk_t* chunks, size_t chunk_count){
  printf("%ld chunks\n",chunk_count);
  fflush(stdout);
  for(size_t i = 0; i < chunk_count; i++){
    printf("chunks[%ld].faults_already_checked=%lld\n",i,chunks[i].faults_already_checked);
    for(int j = 0 ; j < chunks[i].s ; j++)
    printf("chunks[%ld].combination[%d]=%d\n",i,j,chunks[i].current_combination[j]);

    for(int j = 0 ; j < chunks[i].s ; j++)
      printf("chunks[%ld].final_combination_to_check[%d]=%d\n",i,j,chunks[i].final_combination_to_check[j]);
  }
  
}

// Partitions a fault tolerance problem into a set of chunks that can be iterated over independently.
// Returns an array of size that is stored in chunk_count.
// Each chunk has size at least minimum_chunk_size, except when it is the last one.
// Chunks may have more than minimum_chunk_size faults in them.
chunk_t* partition_into_chunks(uint32_t fault_location_count, uint32_t *fault_location_partitions, uint32_t s, size_t minimum_chunk_size, size_t* chunk_count){
  
    int* current_fault_location_combination = (int*)calloc(s,sizeof(int));
    int* first_balls = (int*)calloc(fault_location_count,sizeof(int));
    int* last_balls = (int*)calloc(fault_location_count,sizeof(int));
    init_first_and_last_balls(first_balls, last_balls, fault_location_count, fault_location_partitions);
    unsigned long long N = combi(fault_location_count,s);
    if(!minimum_chunk_size)minimum_chunk_size=1;
    *chunk_count=1+(size_t)(N/minimum_chunk_size);
    chunk_t* return_chunks=(chunk_t*)calloc(*chunk_count,sizeof(chunk_t));
    for(size_t i = 0; i < *chunk_count; i++){
      unsigned long long starting_fault_locations_rank = i * minimum_chunk_size;
      unsigned long long ending_fault_locations_rank = min(((i + 1) * minimum_chunk_size) - 1,N-1);
      // store the fault location combination
      unrank(s,current_fault_location_combination,starting_fault_locations_rank);
      return_chunks[i].current_combination=(int*)calloc(s,sizeof(int));
      for(int j = 0; j < s; j++){
        return_chunks[i].current_combination[j]=first_balls[current_fault_location_combination[j]];
      }
      // now do the ending combination
      return_chunks[i].final_combination_to_check=(int*)calloc(s,sizeof(int));
      unrank(s,current_fault_location_combination,ending_fault_locations_rank);
      for(int j = 0; j < s; j++){
        return_chunks[i].final_combination_to_check[j]=last_balls[current_fault_location_combination[j]];
      }
      return_chunks[i].s=s;
    }
    free(first_balls);
    free(last_balls);
    free(current_fault_location_combination);
    return return_chunks;
  }
  void init_first_and_last_balls(int* first_balls, int* last_balls, uint32_t fault_location_count, uint32_t *fault_location_partitions){
    for (int i = 0; i < fault_location_count; i++) {
      first_balls[i] = 0;
      for (int j = 0; j < i; j++) {
        first_balls[i] += fault_location_partitions[j];
      }
    }
    for (int i = 0; i < fault_location_count; i++) {
      last_balls[i] = 0;
      for (int j = 0; j < i + 1; j++) {
        last_balls[i] += fault_location_partitions[j];
      }
      last_balls[i]--;
    }
  }
int chunk_finished_p(chunk_t* chunk){
  for(int i = 0; i < chunk->s; i++){
    if(chunk->current_combination[i]!=chunk->final_combination_to_check[i])
      return 0; // can't have finished
  }
  // must have finished
  return 1;
}
