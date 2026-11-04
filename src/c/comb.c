#include "comb.h"

comb_code_t __do_verify_css_magic_state(
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
    unsigned x_permutation, unsigned z_permutation,
    chunk_t* chunk) {
  lwlog_debug("setting up the iteration");

  int first_balls[MAX_SIZE_COUNT];
  int last_balls[MAX_SIZE_COUNT];
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
  dcp_t *current_dcp_output_combination = init_dcp(n_all);
  dcp_t *error_with_stabilisers = init_dcp(n_all);
  stabiliser_state_t *state = init_state(n_all);
  pauli_z_coset_t *coset = init_pauli_z_coset(n_all, n_all);
  mzd_t *C_z = mzd_init(z_measurement_count + n_z, n_all);
  mzd_t *C_x = mzd_init(z_measurement_count + n_z, n_all);
  mzd_t *P_x = mzd_init(x_measurement_count, n_all);
  mzd_t *D =
      mzd_init(n_all - (z_measurement_count + n_z), z_measurement_count + n_z);
  mzd_t *z_measurements_for_x = mzd_init(z_measurement_count, n_all);
  for (int i = 0; i < z_measurement_count; i++) {
    // WARNING! wont work for parities. need a matrix solve routine to figure
    // out which x operators anticommute with the z measurement paulis
    mzd_copy_row(C_z, i, z_measurement_paulis, i);
    mzd_copy_row(C_x, i, z_measurement_paulis, i);
  }
  for (int i = 0; i < n_z; i++) {
    mzd_write_bit(C_x, i + z_measurement_count, rref_z[i], 1);
    mzd_copy_row(C_z, i + z_measurement_count, z_stabilisers, i);
  }
  for (int i = 0; i < x_measurement_count; i++) {
    mzd_copy_row(P_x, i, x_measurement_paulis, i);
  }
  // now copy z measurements for x
  for (int i = 0; i < z_measurement_count; i++) {
    mzd_copy_row(z_measurements_for_x, i, z_measurement_paulis, i);
  }
  // Calculate A
  mzd_t *tmp = mzd_copy(0, C_z);
  mzd_t *A_transpose =
      mzd_kernel_left_pluq(tmp, 0);  // C_z A^T = 0 by definition
  mzd_t *A = mzd_transpose(0, A_transpose);
  // free to do row ops on A as row ops preserve A C_Z^T = 0.
  mzd_echelonize_naive(A, 1);  // 1 for full rref (Gauss-Jordan elimination)
  mzp_t *Q = mzp_init(n_all);
  mzp_t *P = mzp_init(n_all - (n_z + z_measurement_count));
  mzd_pluq(A, P, Q, 0);  // now the columns of A have been reordered in order to
                         // later on do the row ops in place.
  mzd_free(A_transpose);
  mzd_free(tmp);
  int row_ops_count = 0;
  for (int i = 0; i < (n_all - (n_z + z_measurement_count)); i++) {
    for (int j = (n_all - (n_z + z_measurement_count)); j < n_all; j++) {
      if (mzd_read_bit(A, i, j)) {
        row_ops_count++;
      }
    }
  }
  int *row_ops_a = malloc(row_ops_count * sizeof(int));
  int *row_ops_b = malloc(row_ops_count * sizeof(int));
  int count = 0;
  for (int i = 0; i < (n_all - (n_z + z_measurement_count)); i++) {
    for (int j = (n_all - (n_z + z_measurement_count)); j < n_all; j++) {
      if (mzd_read_bit(A, i, j)) {
        row_ops_a[count] = i;
        row_ops_b[count] = j;
        count++;
      }
    }
  }
  mzd_t *C_z_transpose = mzd_transpose(0, C_z);
  mzd_t *C_x_transpose = mzd_transpose(0, C_x);

    do {
      chunk->faults_already_checked++;
            // Multiply DCP errors together. TODO verify that order is correct
      set_dcp_to(current_dcp_output_combination,
                 &(output_errors[chunk->current_combination[0]]));
      for (int i = 1; i < s; i++) { // later faults occur earlier (to the left) in the product.
        #if 0
        if(current_combination[i] <= chunk->current_combination[i-1]){
          printf("invalid combination: %d <= %d, fault counter %lld\n", chunk->current_combination[i],
                 chunk->current_combination[i-1], *fault_counter);
          for(int j = 0 ; j <s; j++){
            printf("comb[%d]=%d, ",j,chunk->current_combination[j]);
          }
          exit(1);
        }
        #endif
	dcp_mulf_right(current_dcp_output_combination, &(output_errors[chunk->current_combination[i]]));
      }
      // output errors multiplied. now, multiply by all combinations of logical
      // operators and test them:
      // TODO they all commute (i think) so a gray code should be fine here
      // really loop over all combinations of dcp stabilisers
      int good_combination_found = 0;
      for (uint32_t dcp_index = 0;
           dcp_index < (1 << k) && !good_combination_found; dcp_index++) {
        set_dcp_to(error_with_stabilisers, current_dcp_output_combination);
        //
        for (int i = 0; i < k && !good_combination_found; i++) {
          if (dcp_index & (1 << i)) {
            // multiply the output error on the right by dcp stabilisers
            lwlog_info("multiplying error by logical stabiliser");
            //dcp_mulf_right(error_with_stabilisers, &(dcp_stabilisers[i]));
            dcp_mulf_left(&(dcp_stabilisers[i]), error_with_stabilisers);
	  }
        }
        // now check to see if it is good or not.
        if (!dcp_fault_bad_p_css(error_with_stabilisers, C_z, C_z_transpose, C_x, C_x_transpose,
                                 P_x, z_measurements_for_x, x_permutation,
                                 z_permutation, x_lookup_table, z_lookup_table,
                                 x_stabilisers, rref_x, n_all, n, n_x, n_z,
                                 state, coset, D, Q, row_ops_count, row_ops_a,
                                 row_ops_b)) {
          // a good one was found after multiplying by a choice of logicals, set
          // the flag
          good_combination_found = 1;
        }
      }
      if (!good_combination_found) {
        lwlog_debug("no good combination was found: fault combination was bad");
        // set the indices and indicate a bad fault was found
        for (int i = 0; i < s; i++) {
          bad_fault_array[i] = chunk->current_combination[i];
        }
        return BAD_FAULT_FOUND;
      }
      //lwlog_debug("about to incf on fault number %lld",*fault_counter);
      //#endif
    } while (!chunk_finished_p(chunk) && ft_iteration_combination_incf(s, fault_location_partitions,
                                           fault_location_count, first_balls,
                                           last_balls, chunk->current_combination)

    );
  return ALL_FAULTS_GOOD;  // got to end of loop, so all faults are good
}
comb_code_t comb_verify_css_magic_state(
    uint32_t n_all,  ///< number of qubits in code + all ancillas
    uint32_t n,      // number of qubits in code
    uint32_t n_x,    // number of independent x stabiliser generators in code
    uint32_t n_z,    // number of independent z stabiliser generators in code
    uint32_t k,  // number of logical qubits in code. equivalently, number of
                 // logical dcp operators that stabilise the code
    comb_pauli_x_t *_x_stabilisers,   // n-x of these
    comb_bit_array_t x_lookup_table,  // x lookup table. of size 2^(n-n_x)
    comb_pauli_z_t *_z_stabilisers, comb_bit_array_t z_lookup_table,
    comb_dcp_t *_dcp_stabilisers,   // k of these
    uint32_t fault_location_count,  // number of fault locations
    uint32_t
        *fault_location_partitions,  // at least fault_location_count of these.
                                     // number of faults at each location
    uint32_t s,  // number of faults to which circuit should be resilient to
    uint32_t *bad_fault_array,     // where to store indices of bad_faults
    uint32_t output_error_count,   // number of output errors that result from
                                   // propagating all the faults. should be at
                                   // least sum(fault_location_partitions)
    comb_dcp_t *_output_errors,    // at least output_error_count of these
    uint32_t x_measurement_count,  // number of x measurement paulis
    comb_pauli_x_t
        *_x_measurement_paulis,    // postselections made in x basis. at
                                   // least x_measurement_count of these
    uint32_t z_measurement_count,  // number of z measurement paulis
    comb_pauli_z_t *_z_measurement_paulis, uint32_t should_serialise_p) {
  // serialise run
  if(should_serialise_p){
  char name[] = "comb-runXXXXXX.run";
  int fd = mkstemps(name,4);
  int serialise_success = serialise_run_to_file(fd,n_all,n,n_x,n_z,k,_x_stabilisers,x_lookup_table,
                                                _z_stabilisers,z_lookup_table,_dcp_stabilisers,fault_location_count,
                                                fault_location_partitions,s,output_error_count,
                                                _output_errors,x_measurement_count,_x_measurement_paulis,
                                                z_measurement_count,_z_measurement_paulis);
  if(serialise_success!=0) return SERIALIZATION_ERROR;
  else printf("saved run to: %s\n",name);
  }
  lwlog_info("start timing \n");

  struct timespec start, finish;
  clock_gettime(CLOCK_MONOTONIC, &start);
  lwlog_info(
      "n_all=%d, n=%d, n_x = %d, n_z = %d, k = %d,fault_location_count "
      "= %d, s = %d, output_error_count = %d, x_measurement_count = %d, "
      "z_measurement_count = %d",
      n_all, n, n_x, n_z, k, fault_location_count, s, output_error_count,
      x_measurement_count, z_measurement_count);
  lwlog_info("allocating and converting objects");
  uint32_t rref_x[1024];
  uint32_t rref_z[1024];
  unsigned x_permutation = 0;
  unsigned z_permutation = 0;
  if (!stabilisers_rref_p(_x_stabilisers, n_x, n, rref_x, &x_permutation)) {
    lwlog_err("x stabilisers not in rref");
    return INVALID_INPUT;
  }
  if (!stabilisers_rref_p(_z_stabilisers, n_z, n, rref_z, &z_permutation)) {
    lwlog_err("z stabilisers not in rref");
    return INVALID_INPUT;
  }
  dcp_t *dcp_stabilisers = (dcp_t *)malloc(k * sizeof(dcp_t));
  dcp_t *output_errors = (dcp_t *)malloc(output_error_count * sizeof(dcp_t));
  mzd_t *x_measurement_paulis = mzd_init(x_measurement_count, n_all);
  mzd_t *z_measurement_paulis = mzd_init(z_measurement_count, n_all);
  mzd_t *x_stabilisers = mzd_init(n_x, n);
  mzd_t *z_stabilisers = mzd_init(n_z, n);
  int chunks_checked=0;
  unsigned long long fault_counter = 0;
  lwlog_info("converting objects");
  for (int i = 0; i < n_x; i++) {
    lwlog_debug("converting x stabiliser %d", i);
    convert_row(x_stabilisers, i, &(_x_stabilisers[i]));
    lwlog_debug("x_stabilisers[%d]=%s", i,
                print_row_vector_as_pauli(x_stabilisers, i, "X"));
  }
  for (int i = 0; i < n_z; i++) {
    lwlog_debug("converting z stabiliser %d", i);
    convert_row(z_stabilisers, i, &(_z_stabilisers[i]));
    lwlog_debug("z_stabilisers[%d]=%s", i,
                print_row_vector_as_pauli(z_stabilisers, i, "Z"));
  }
  for (int i = 0; i < k; i++) {
    lwlog_debug("converting logical stabiliser %d", i);
    convert_dcp(&(dcp_stabilisers[i]), &(_dcp_stabilisers[i]));
    lwlog_debug("logical stabiliser %d: %s", i,
                print_dcp(&(dcp_stabilisers[i]), 0));
  }
  lwlog_debug("converting %d output errors", output_error_count);
  for (int i = 0; i < output_error_count; i++) {
    convert_dcp(&(output_errors[i]), &(_output_errors[i]));
  }
  for (int i = 0; i < x_measurement_count; i++) {
    lwlog_debug("converting x measurement %d", i);
    convert_row(x_measurement_paulis, i, &(_x_measurement_paulis[i]));
    lwlog_debug("x_measurement_paulis[%d]=%s", i,
                print_row_vector_as_pauli(x_measurement_paulis, i, "X"));
  }
  for (int i = 0; i < z_measurement_count; i++) {
    lwlog_debug("converting z measurement %d", i);
    convert_row(z_measurement_paulis, i, &(_z_measurement_paulis[i]));
    lwlog_debug("z_measurement_paulis[%d]=%s", i,
                print_row_vector_as_pauli(z_measurement_paulis, i, "Z"));
  }

  for (int i = 0; i < fault_location_count; i++) {
    lwlog_debug("fault location partition %d: %d\n", i,
                fault_location_partitions[i]);
  }
  lwlog_debug("all objects converted");
  size_t chunk_count=0;
  int delta = 10000 / (15^s);
  chunk_t* chunks = partition_into_chunks(fault_location_count, fault_location_partitions,s,delta,&chunk_count);
  describe_chunks(chunks,chunk_count);
  comb_code_t return_code = ALL_FAULTS_GOOD;

  #ifndef COMB_MPI
  #pragma omp parallel for
  for(int i = 0; i < chunk_count; i++){
    printf("Verifying chunk [%d/%ld]\n",i,chunk_count);
  comb_code_t chunk_return_code = __do_verify_css_magic_state(
      n_all, n, n_x, n_z, k, x_stabilisers, x_lookup_table, z_stabilisers,
      z_lookup_table, dcp_stabilisers, fault_location_count,
      fault_location_partitions, s, bad_fault_array, output_error_count,
      output_errors, x_measurement_count, x_measurement_paulis,
      z_measurement_count, z_measurement_paulis, rref_x, rref_z, x_permutation,
      z_permutation, chunks+i);
  printf("Verified chunk [%d/%ld], %lld faults checked.\n",i,chunk_count,chunks[i].faults_already_checked);
  if(chunk_return_code!=ALL_FAULTS_GOOD){
    return_code= chunk_return_code; // TODO cleanup and break here
  }
  fault_counter+=chunks[i].faults_already_checked;

  }
  #else
  // MPI
    // Get the number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

#pragma omp parallel shared(chunks_checked,fault_counter)
  {
    if(world_size==1 && omp_get_num_threads() == 1){ // exit if serial; need at least two parallel workers
      printf("only one process and thread: use serial build instead\n");
      fflush(stdout);
      MPI_Abort(MPI_COMM_WORLD, INTERNAL_ERROR);
    }
    if(omp_get_thread_num() == 0 && world_rank == 0)mpi_primary(chunk_count);
    else mpi_secondary(&chunks_checked, &fault_counter, n_all, n, n_x, n_z, k, x_stabilisers,
                       x_lookup_table, z_stabilisers,
      z_lookup_table, dcp_stabilisers, fault_location_count,
      fault_location_partitions, s, bad_fault_array, output_error_count,
      output_errors, x_measurement_count, x_measurement_paulis,
      z_measurement_count, z_measurement_paulis, rref_x, rref_z, x_permutation,
      z_permutation,chunks);
  }

  #endif
  // Stop timing
 cleanup:
  clock_gettime(CLOCK_MONOTONIC, &finish);
  double seconds = (finish.tv_sec - start.tv_sec) +
                   1.0e-9 * (finish.tv_nsec - start.tv_nsec);
  #ifdef COMB_MPI
    printf("Process %d: %d chunks checked, Took %f seconds, %lld fault combinations checked.\n", world_rank, chunks_checked, seconds,
         fault_counter);
  #else
  printf("Took %f seconds, %lld fault combinations checked.\n", seconds,
         fault_counter);
  #endif
  lwlog_debug("freeing memory");
  free_dcp_array(k, dcp_stabilisers);
  free_dcp_array(output_error_count, output_errors);
  mzd_free(x_measurement_paulis);
  mzd_free(z_measurement_paulis);
  mzd_free(x_stabilisers);
  mzd_free(z_stabilisers);
  return return_code;
}
