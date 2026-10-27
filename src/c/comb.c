#include "comb.h"

#include <math.h>
#include <stdio.h>
#include <time.h>
#include <omp.h>
int dcp_fault_bad_p_css(
                        dcp_t *error, mzd_t* C_z, mzd_t *C_z_transpose, mzd_t *C_x, mzd_t *C_x_transpose, mzd_t *P_x,
    mzd_t *z_measurements_for_x,
    unsigned x_permutation,  // permutation mask for x canonicalisation
    unsigned z_permutation,  // permutation mask for z canonicalisation
    comb_bit_array_t x_lookup_table, comb_bit_array_t z_lookup_table,
    mzd_t *x_stabilisers, uint32_t *rref_x, int n_all, int n, int n_x, int n_z,
    stabiliser_state_t *state, pauli_z_coset_t *coset, mzd_t *D, mzp_t *Q,
    int row_ops_count, int *row_ops_a, int *row_ops_b) {
  lwlog_debug("Testing %s", print_dcp(error, 0));
  mzd_t *x_error =
      error->h_x_transpose;  // dont need to transpose to get x error
  lwlog_debug("X error is %s", print_row_vector_as_pauli(x_error, 0, "X"));
  for (int i = 0; i < z_measurements_for_x->nrows; i++) {
    if (row_product(x_error, 0, z_measurements_for_x, i)) {
      lwlog_debug(
          "X error anticommutes with Z measurement %d %s, so whole "
          "error is good",
          i, print_row_vector_as_pauli(z_measurements_for_x, i, "Z"));
      return 0;
    } else {
      lwlog_debug("X error commutes with Z measurement %d %s, continuing...", i,
                  print_row_vector_as_pauli(z_measurements_for_x, i, "Z"));
    }
  }
  // need to canonicalise the x error.
  for (int i = 0; i < x_stabilisers->nrows; i++) {
    if (mzd_read_bit(x_error, 0,
                     rref_x[i])) {  // eliminate by multiplying by a stabiliser
      lwlog_debug("Multiplying by %s to canonicalise...",
                  print_row_vector_as_pauli(x_stabilisers, i, "X"));
      mzd_combine(x_error, 0, 0, x_error, 0, 0, x_stabilisers, i, 0);
      lwlog_debug("New X error is %s",
                  print_row_vector_as_pauli(x_error, 0, "X"));
    }
  }
  // apply the permutation to x.
  lwlog_debug(
      "Error before permutation and eliminating ancillas is %s (mask %x)",
      print_row_vector_as_pauli(x_error, 0, "X"), x_permutation);
  x_error->data[0] &= ((1UL << (n)) - 1);
  lwlog_debug("Error after eliminating ancilla terms is %s",
              print_row_vector_as_pauli(x_error, 0, "X"));
  permute_columns(x_error, x_permutation);
  lwlog_debug("Error after permutation is %s",
              print_row_vector_as_pauli(x_error, 0, "X"));
  int reduced_x_error = mzd_data(x_error);
  lwlog_debug("Reduced x error is %x", reduced_x_error);
  if (!TestBit(x_lookup_table, reduced_x_error)) {
    lwlog_debug("x error undetected and high weight, so whole error is bad");
    return 1;
  }
  // Now, test the z part. Start by turning it into a state
  set_stabiliser_state_to_dcp(state, error);
  // Instead of doing all the Z and then X measurements separately, we combine
  // them (batched gaussian elimination):
  // 1. Row ops on S_i using A to get S'_i
  int k = C_z_transpose->ncols;
  lwlog_debug("before row ops on S_i, state is %s",
              print_stabiliser_state(state, 0));
  // In order to do the row ops in place, we first apply a permutation
  // calculated beforehand via the PLUQ decomposition of A:
  mzd_apply_p_left(state->x_parts, Q);
  mzd_apply_p_left(state->z_parts, Q);
  mzd_apply_p_left(state->delta_vector, Q);
  mzd_apply_p_left(state->epsilon_vector, Q);
  // We can now do the row ops in place, specified by the parameters row_ops*
  for (int i = 0; i < row_ops_count; i++) {
    add_row_to(state, row_ops_b[i], state, row_ops_a[i]);
  }
  lwlog_debug("after row ops on S_i, state is %s",
              print_stabiliser_state(state, 0));
  // 2. Calculate D using matrix multiplication
  // D is whether the j-th x measurement anticommutes with the i-th S_i
  // shrink state->z_parts for the multiplication; want to consider only the
  // first n_all - k stabs.
  state->z_parts->nrows = (n_all - k);
  // D <- Z_parts . C_x_transpose
  _mzd_mul_m4rm(D,state->z_parts, C_x_transpose, 0, 1);
  //_mzd_mul_naive(D, state->z_parts, C_x,  1);  // C_x k x n, state->Z_parts is n - k x n
  // mzd_t* D_transpose = mzd_transpose(0,D); // calculate D^T as well, maybe it
  // will be needed instead.
  //  3. Calculate S''_i via matrix multiplication D C_Z
  //  i.e. multiply the S''_i = S'_i D_ij C_z

  // Z_parts <- Z_parts + D C_z
  _mzd_mul_m4rm(state->z_parts,D,C_z,0,0);

  //_mzd_mul_naive(state->z_parts, D, C_z_transpose, 0);  // z parts is n-k x n, D is n-k x k, C_z is k x n
  state->z_parts->nrows = (n_all);
  // 4. Add the X measurements as stabilisers by filling the lower k elts of
  // dest_state with C_x
  for (int i = (n_all - k); i < n_all; i++) {
    mzd_row_clear_offset(state->x_parts, i, 0);
    mzd_combine(state->x_parts, i, 0, state->x_parts, i, 0, C_x,
                i - (n_all - k), 0);
    mzd_row_clear_offset(state->z_parts, i, 0);
    mzd_row_clear_offset(state->epsilon_vector, i, 0);
    mzd_row_clear_offset(state->delta_vector, i, 0);
  }
  lwlog_debug("after batch ge, state is %s", print_stabiliser_state(state, 0));

  if (postselect_on_measuring_plus_one_x(state, P_x)) {
    lwlog_debug("error detected in x postselections, so error is good");
    return 0;
  }
  // error is now canonicalised. turn it into a coset and apply the z
  // permutation to it:
  set_coset_to_stabiliser_state(coset, state);
  lwlog_debug("coset before permuting is %s", print_pauli_z_coset(coset, 0));
  permute_columns(coset->factor_transpose, z_permutation);
  permute_columns(coset->set_generators, z_permutation);
  lwlog_debug("coset after permuting is %s", print_pauli_z_coset(coset, 0));
  // iterate over all elements in the coset and verify that each one is good. if
  // not, return bad
  for (int i = 0; i < (1 << COSET_K(coset)); i++) {
    lwlog_debug("Testing Z error %d in coset %s", i,
                print_row_vector_as_pauli(coset->factor_transpose, 0, "Z"));
    int reduced_z_error = mzd_data(coset->factor_transpose);
    lwlog_debug("number is %d", reduced_z_error);
#if 0
    // This was intended to be a test that the reduced z error is not too large.
    // However, this can result from a missing measurement, so it does not indicate a bug with the code per se.
    // In the future, this case should be dealt with by sanitising the inputs better. TODO important
    // Until this is fixed, the behaviour of Comb is undefined if measurements are missing.
    if(reduced_z_error >= (1 << (n - n_z))){
      lwlog_err("reduced error too large for lookup table, should not have happened, exiting");
      exit(1);
    }
#endif
    if (!TestBit(z_lookup_table, reduced_z_error)) {
      lwlog_debug("z error in coset is high weight, so error is bad");
      return 1;
    }
    if ((i + 1) !=
        (1 << COSET_K(coset)))  // multiply the coset factor by a generator of
                                // the vector space according to a gray code
      coset->factor_transpose->data[0] ^=
          mzd_row(coset->set_generators, gray_bit_to_switch(i))[0];
  }
  lwlog_debug("all errors in coset were low weight, so error is good");
  return 0;
}

int stabilisers_rref_p(comb_vector_t *operators, uint32_t n, uint32_t m,
                       uint32_t *coefficients, unsigned *permutation_mask) {
  // ensures the operators are in rref. n operators on m qubits.
  int current_leading_index = -1;
  for (int row = 0; row < n; row++) {
    comb_vector_t *current_operator = &(operators[row]);
    int new_leading_index = m;  // new_leading_index = n means that the row
                                // contains no entries at all.
    for (int column = 0; column < m && (new_leading_index == m); column++) {
      if (TestBit(current_operator->data, column)) new_leading_index = column;
    }
    if (new_leading_index < current_leading_index &&
        !((new_leading_index = m) && (current_leading_index == m)))
      return 0;  // new leading index must be to the right of the old one
    // verify if new leading index has no entries in other rows
    for (int row_2 = 0; row_2 < n; row_2++) {
      if (row_2 != row) {
        if (TestBit(operators[row_2].data, new_leading_index))
          return 0;  // leading index shares entries in other rows
      }
    }
    current_leading_index = new_leading_index;
    coefficients[row] = current_leading_index;
    *permutation_mask |= (1 << current_leading_index);
    lwlog_info("rref: coefficients[%d]=%d", row, current_leading_index);
  }
  *permutation_mask = ~*permutation_mask;
  return 1;
}

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
      // chunk->faults_to_check_count--;
      chunk->faults_already_checked++;
#ifndef DRY_RUN
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
        //dcp_mulf_left(&(output_errors[current_combination[i]]), current_dcp_output_combination);
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
        free(row_ops_a);
        free(row_ops_b);
        return BAD_FAULT_FOUND;
      }
      //lwlog_debug("about to incf on fault number %lld",*fault_counter);
#endif
    } while (!chunk_finished_p(chunk) && ft_iteration_combination_incf(s, fault_location_partitions,
                                           fault_location_count, first_balls,
                                           last_balls, chunk->current_combination)

    );
        free(row_ops_a);
        free(row_ops_b);
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
    comb_pauli_z_t *_z_measurement_paulis, uint32_t thread_count) {
  // serialise run
  char name[] = "comb-runXXXXXX.run";
  int fd = mkstemps(name,4);
  int serialise_success = serialise_run_to_file(fd,n_all,n,n_x,n_z,k,_x_stabilisers,x_lookup_table,
                                                _z_stabilisers,z_lookup_table,_dcp_stabilisers,fault_location_count,
                                                fault_location_partitions,s,output_error_count,
                                                _output_errors,x_measurement_count,_x_measurement_paulis,
                                                z_measurement_count,_z_measurement_paulis);
  if(serialise_success!=0) return SERIALIZATION_ERROR;
  else printf("saved run to: %s\n",name);
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
  // Stop timing
 cleanup:
  clock_gettime(CLOCK_MONOTONIC, &finish);
  double seconds = (finish.tv_sec - start.tv_sec) +
                   1.0e-9 * (finish.tv_nsec - start.tv_nsec);
  printf("Took %f seconds, %lld fault combinations checked.\n", seconds,
         fault_counter);
  lwlog_debug("freeing memory");
  free_dcp_array(k, dcp_stabilisers);
  free_dcp_array(output_error_count, output_errors);
  mzd_free(x_measurement_paulis);
  mzd_free(z_measurement_paulis);
  mzd_free(x_stabilisers);
  mzd_free(z_stabilisers);
  return return_code;
}

