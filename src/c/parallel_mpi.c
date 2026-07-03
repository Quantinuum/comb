#ifdef USE_MPI
#include "comb.h"

#include <math.h>
#include <stdio.h>
#include <time.h>
#include <omp.h>

// Routines to help with parallelism and MPI
comb_code_t mpi_check_chunks(size_t chunk_count, chunk_t* chunks,
                               uint32_t n_all,  ///< number of qubits in code + all ancillas
    uint32_t n,      // number of qubits in code
    uint32_t n_x,    // number of independent x stabiliser generators in code
    uint32_t n_z,    // number of independent z stabiliser generators in code
    uint32_t k,  // number of logical qubits in code. equivalently, number of
                 // logical dcp operators that stabilise the code
    mzd_t *x_stabilisers,             // n-x of these
    comb_bit_array_t x_lookup_table,  // x lookup table. of size 2^(n-n_x)
    mzd_t *z_stabilisers, comb_bit_array_t z_lookup_table,
    dcp_t *dcp_stabilisers,         // k of these
    uint32_t fault_location_count,  // number of fault locations
    uint32_t
        *fault_location_partitions,  // at least fault_location_count of these.
                                     // number of faults at each location
    uint32_t s,  // number of faults to which circuit should be resilient to
    uint32_t *bad_fault_array,     // where to store indices of bad_faults
    uint32_t output_error_count,   // number of output errors that result from
                                   // propagating all the faults. should be at
                                   // least sum(fault_location_partitions)
    dcp_t *output_errors,          // at least output_error_count of these
    uint32_t x_measurement_count,  // number of x measurement paulis
    mzd_t *x_measurement_paulis,   // postselections made in x basis. at least
                                   // x_measurement_count of these
    uint32_t z_measurement_count,  // number of z measurement paulis
    mzd_t *z_measurement_paulis, uint32_t *rref_x, uint32_t *rref_z,
                             unsigned x_permutation, unsigned z_permutation){

  // A single threaded MPI implementation.
  // Initialize the MPI environment
  MPI_Init(NULL, NULL);

  // Get the number of processes
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  // Get the rank of the process
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  // Get the name of the processor
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  int name_len;
  MPI_Get_processor_name(processor_name, &name_len);

  comb_bit_array_t array;
  int number;
  if (world_rank == 0) {
    array = comb_alloc_bit_array(chunk_count);
    //MPI_Send(&number, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
  } else if (world_rank == 1) {
    //MPI_Recv(&number, 1, MPI_INT, 0, 0, MPI_COMM_WORLD,
    //         MPI_STATUS_IGNORE);
    //printf("Process 1 received number %d from process 0\n",
    //       number);
  }
  mpi_info_t info;
  info.array = array;
  info.chunk_count = chunk_count;
  int chunk_id = -1;
  while((chunk_id = mpi_blocking_get_free_chunk(array, chunk_count)) != -1){
    
    comb_code_t chunk_return_code = __do_verify_css_magic_state(n_all, n, n_x, n_z, k, x_stabilisers, x_lookup_table, z_stabilisers, z_lookup_table, dcp_stabilisers, fault_location_count,fault_location_partitions, s, bad_fault_array, output_error_count, output_errors, x_measurement_count, x_measurement_paulis,  z_measurement_count, z_measurement_paulis, rref_x, rref_z, x_permutation, z_permutation, chunks+ chunk_id, &info);
    if(chunk_return_code != ALL_FAULTS_GOOD){
      if(chunk_return_code == BAD_FAULT_FOUND){
        mpi_broadcast_bad_fault_found(bad_fault_array,s);
      }
      else mpi_broadcast_other_error_code(chunk_return_code);
    }

    
  }



  // clean up mpi
  MPI_Finalize();

}

#define REQUEST_CHUNK_TAG 42

int mpi_blocking_get_free_chunk(comb_bit_array_t array, size_t chunk_count){
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  if(world_rank == 0){
    // we are on the master process so just get one
    for(int i = 0 ; i < chunk_count; i++){
      if(!TestBit(array,i)){
        SetBit(array,i);
        return i;
      }
    }
    // none free, so return -1
    return -1;
  }
  // not on the master process so request a new chunk
  int elt = 42;
  printf("sending a chunk request to master from process %d\n", world_rank);
  MPI_Send(NULL, 0, MPI_INT, 0, REQUEST_CHUNK_TAG, MPI_COMM_WORLD);
  // now block until receiving the ID from the master process.
  MPI_Recv(&elt, 1, MPI_INT, 0, REQUEST_CHUNK_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  printf("receiving a chunk number %d\n", elt);
  return elt; // elt now the ID of the new block
}
void mpi_do_master_work(comb_bit_array_t array, size_t chunk_count, MPI_Request **requests){
  // Assumes we are on the master thread. Loops through messages doing the work.
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  if(*requests == NULL){ // need to initialise the requests
    *requests = (MPI_Request*)calloc(world_size, sizeof(MPI_Request));
    for(int i = 1; i < world_size ; i++){
      MPI_Irecv(NULL, 0, MPI_INT, i, REQUEST_CHUNK_TAG, MPI_COMM_WORLD, *requests + i);
    }
  }

  for(int i = 1; i < world_size; i++){
    // check for message
    int flag;
    MPI_Test(*requests + i, &flag, MPI_STATUS_IGNORE);
    if(flag){
      printf("received a chunk request from process %d\n",i);
      // one of the non blocking requests completed; set up another one
      MPI_Irecv(NULL, 0, MPI_INT, i, REQUEST_CHUNK_TAG, MPI_COMM_WORLD, *requests + i);
      // and now find a free chunk
      int free_chunk = -1;

      for(int j = 0; j < chunk_count; j++){
        if(!TestBit(array,j)){
          SetBit(array,j);
          free_chunk = j;
        }
      }
      // and send it
      printf("sending chunk %d to process %d\n",free_chunk, i);
      MPI_Send(&free_chunk, 1, MPI_INT, i, REQUEST_CHUNK_TAG, MPI_COMM_WORLD);
    }
  }
}

void mpi_broadcast_bad_fault_found(int *bad_fault_array, int s){
  // TODO replace with something better
  for(int i = 0; i < s; i++){
    printf("bad fault found: %d\n",bad_fault_array[i]);
  }
  exit(1);
}
void mpi_broadcast_other_error_code(comb_code_t code){
  printf("other error code: %d\n",code);
  exit(1);
}
#endif
