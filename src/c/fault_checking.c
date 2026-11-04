#include "comb.h"

int dcp_fault_bad_p_css(dcp_t *error, mzd_t* C_z, mzd_t *C_z_transpose, mzd_t *C_x, mzd_t *C_x_transpose, mzd_t *P_x,
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
