#ifdef COMB_MPI
#include "comb.h"


int mpi_get_chunk_to_process(int world_rank){
  // To be called from a secondary process. Sends a message requesting a new chunk, and blocks until one is sent from the primary process.
  // This means that secondary processes can never die, the primary process must handle termination.
  int chunk_to_process;
  MPI_Status status;
  fprintf(stderr,"Requesting chunk from primary process; from (%d/%d)\n",world_rank,omp_get_thread_num());
  MPI_Sendrecv(NULL, 0, MPI_INT, 0, REQUEST_CHUNK_TAG,
               &chunk_to_process, 1, MPI_INT, 0, RECEIVE_CHUNK_TAG,
               MPI_COMM_WORLD, &status);
  fprintf(stderr,"Received chunk %d from primary process; from (%d/%d)\n",chunk_to_process,world_rank,omp_get_thread_num());
  return chunk_to_process;
}

void mpi_notify_primary_chunk_good(int chunk_num){
  // Notifies the primary process that chunk_num was found to be good.
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    fprintf(stderr,"About to notify primary that chunk %d was found good; from (%d/%d)\n",chunk_num,world_rank,omp_get_thread_num());
  MPI_Send(&chunk_num, 1, MPI_INT, 0, NOTIFY_CHUNK_GOOD, MPI_COMM_WORLD);
}

size_t mpi_get_free_chunk(comb_bit_array_t free_chunks_index_array, int chunk_count){
  for(int i = 0; i < chunk_count; i++){
    if(!TestBit(free_chunks_index_array,i)){
      return i;
    }
  }
  return -1;
}

void mpi_primary(int chunk_count){
  
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  // we are on the master thread on process 0
      // communicate with the other processes, giving them work as appropriate
      // an array that has 1s where a chunk has not yet started to be checked.
      comb_bit_array_t free_chunks_index_array = comb_alloc_bit_array(chunk_count);
      comb_bit_array_t chunks_found_good_array = comb_alloc_bit_array(chunk_count);
      // 1. See if there is a message requesting a new chunk. Respond with a free chunk if one exists, marking it as not free anymore.
        MPI_Request chunk_to_process_request;
        int chunk_to_process_code = MPI_Irecv(0, 0, MPI_INT, MPI_ANY_SOURCE,
                             REQUEST_CHUNK_TAG, MPI_COMM_WORLD, &chunk_to_process_request);
        // 2. See if there are any messages indicating that a chunk was confirmed to be good. If there are, mark the corresponding entries in chunks_found_good_array.
        MPI_Request chunk_found_good_request;
        int chunk_found_good_index;
        int chunk_found_good_code =  MPI_Irecv(&chunk_found_good_index, 1, MPI_INT, MPI_ANY_SOURCE,
                              NOTIFY_CHUNK_GOOD, MPI_COMM_WORLD, &chunk_found_good_request);
      while(mpi_get_free_chunk(chunks_found_good_array,chunk_count) != -1) { // we are not done until all chunks are confirmed to be good
        int flag=0;
        MPI_Status status;
        MPI_Test(&chunk_to_process_request,&flag,&status);
        if(flag){
          // a chunk was requested; send it back to the receiver
          int source=status.MPI_SOURCE;
            fprintf(stderr,"From primary: Process %d requested a chunk\n", source);
          // now send a free chunk
            int free_chunk_index=mpi_get_free_chunk(free_chunks_index_array,chunk_count);
            // only send one back if its available
            if(free_chunk_index != -1)SetBit(free_chunks_index_array,free_chunk_index);
            fprintf(stderr,"From primary: Sending free chunk %d to %d, marking it as not free\n", free_chunk_index, source);
            MPI_Send(&free_chunk_index, 1, MPI_INT, source, RECEIVE_CHUNK_TAG, MPI_COMM_WORLD);
          // update the request
        int chunk_to_process_code = MPI_Irecv(0, 0, MPI_INT, MPI_ANY_SOURCE,
                             REQUEST_CHUNK_TAG, MPI_COMM_WORLD, &chunk_to_process_request);
        }
        MPI_Test(&chunk_found_good_request,&flag,&status);
        if(flag){
          // a chunk was found good; mark it
          fprintf(stderr,"From primary: chunk %d was found good; marking it\n",chunk_found_good_index);
          SetBit(chunks_found_good_array,chunk_found_good_index);
          // make the request again
          int chunk_found_good_code =  MPI_Irecv(&chunk_found_good_index, 1, MPI_INT, MPI_ANY_SOURCE,
                                                 NOTIFY_CHUNK_GOOD, MPI_COMM_WORLD, &chunk_found_good_request);
        }
      }
      // send -1 to every waiting child just in case there is some secondary process waiting
      int chunk_index=-1;
      for(int i = 0; i < world_size; i++)
        for(int j = 0; j < omp_get_num_threads(); j++)MPI_Send(&chunk_index,1,MPI_INT,i,RECEIVE_CHUNK_TAG, MPI_COMM_WORLD);
      fprintf(stderr,"From primary: All faults good!\n");
}

void mpi_secondary(int* chunks_checked, unsigned long long* fault_counter, uint32_t n_all,  ///< number of qubits in code + all ancillas
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
    unsigned x_permutation, unsigned z_permutation, chunk_t* chunks){
  
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  int chunk_num = -1;
      while((chunk_num=mpi_get_chunk_to_process(world_rank))!=-1){
        comb_code_t chunk_return_code = __do_verify_css_magic_state(
      n_all, n, n_x, n_z, k, x_stabilisers, x_lookup_table, z_stabilisers,
      z_lookup_table, dcp_stabilisers, fault_location_count,
      fault_location_partitions, s, bad_fault_array, output_error_count,
      output_errors, x_measurement_count, x_measurement_paulis,
      z_measurement_count, z_measurement_paulis, rref_x, rref_z, x_permutation,
      z_permutation, chunks+chunk_num); // check the chunk
        #pragma omp criticial
        { // need atomic increments of the fault counter
          *fault_counter+=chunks[chunk_num].faults_already_checked;
          (*chunks_checked)++;
        }
        if(chunk_return_code!=ALL_FAULTS_GOOD){
          if(chunk_return_code == BAD_FAULT_FOUND){
            for(int i = 0; i < s; i++){
              fprintf(stderr,"bad_fault_array[%d]=%d\n",i,bad_fault_array[i]);
            }
          }
          MPI_Abort(MPI_COMM_WORLD,chunk_return_code);
        }
        else mpi_notify_primary_chunk_good(chunk_num);
      }
}

#endif
