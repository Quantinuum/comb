#include "comb.h"

void init_chunk(chunk_t* chunk, unsigned long long faults_to_check_count, int* starting_combination, int s){
  chunk->faults_to_check_count = faults_to_check_count;
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
    printf("chunks[%ld].faults_to_check_count=%lld\n",i,chunks[i].faults_to_check_count);
    printf("chunks[%ld].faults_already_checked=%lld\n",i,chunks[i].faults_already_checked);
    for(int j = 0 ; j < chunks[i].s ; j++)
    printf("chunks[%ld].combination[%d]=%d\n",i,j,chunks[i].current_combination[j]);
  }
}

// Partitions a fault tolerance problem into a set of chunks that can be iterated over independently.
// Returns an array of size that is stored in chunk_count.
// Each chunk has size at least minimum_chunk_size, except when it is the last one.
// Chunks may have more than minimum_chunk_size faults in them.
chunk_t* partition_into_chunks(uint32_t fault_location_count, uint32_t *fault_location_partitions, uint32_t s, size_t minimum_chunk_size, size_t* chunk_count){
  // naive implementation to start off with
  int current_combination[MAX_SIZE_COUNT] = {0};
  int first_balls[MAX_SIZE_COUNT];
  int last_balls[MAX_SIZE_COUNT];

  for (int i = 0; i < s; i++) {
    current_combination[i] = 0;
    for (int j = 0; j < i; j++) {
      current_combination[i] += fault_location_partitions[j];
    }
  }
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
  lwlog_debug("in partition chunks");
  size_t chunk_delta = 1;
  *chunk_count=chunk_delta;
  chunk_t* return_chunks = (chunk_t*)calloc(*chunk_count,sizeof(chunk_t));
  size_t current_chunk_index=0;
  lwlog_debug("allocated return_chunks");
  // initialise the first chunk to the beginning.
  return_chunks[current_chunk_index].current_combination=calloc(s,sizeof(int));
  for(int i = 0; i < s; i++){
    return_chunks[current_chunk_index].current_combination[i]=current_combination[i];
  }
  return_chunks[current_chunk_index].s=s;
  do{
    // check to see if current chunk has enough space for another fault combination
    if(return_chunks[current_chunk_index].faults_to_check_count < minimum_chunk_size){
      // it does; so simply incf and continue
      return_chunks[current_chunk_index].faults_to_check_count++;
    }
    else{
      // it doesn't, so make the current combination the first one of a new chunk.
      // first ensure that there are enough chunks in the array
      current_chunk_index++;
      if(current_chunk_index >= *chunk_count){
        chunk_t* new_chunks = realloc(return_chunks,(*chunk_count + chunk_delta)
                                      * (sizeof(chunk_t)));
        if(!new_chunks)return 0; // no space
        return_chunks=new_chunks;
        *chunk_count += chunk_delta;
      }
      // initialise the new chunk
      return_chunks[current_chunk_index].faults_to_check_count = 0;
      return_chunks[current_chunk_index].faults_already_checked = 0;
      return_chunks[current_chunk_index].current_combination=calloc(s,sizeof(int));
      for(int i = 0; i < s; i++){
        return_chunks[current_chunk_index].current_combination[i]=current_combination[i];
      }
      return_chunks[current_chunk_index].s=s;
      // now continue
    }
  } while(ft_iteration_combination_incf(s, fault_location_partitions,
                                       fault_location_count, first_balls,
                                        last_balls, current_combination));
  return return_chunks;
}
