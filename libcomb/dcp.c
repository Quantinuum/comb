#ifdef __cplusplus
extern "C" {
#endif
#include "comb.h"

  // 1. Multiply DCPs together
  dcp_t *dcp_mul(dcp_t *result, dcp_t *c2, dcp_t *c1){ // allocation free. does not check sizes
    if(!result)result = init_dcp(DCP_N(c1));
    mzd_add(result->F_diagonal,c2->F_diagonal,c1->F_diagonal); // F parts add
    mzd_add(result->F,c2->F,c1->F); // F parts add
    mzd_add(result->h_x,c2->h_x,c1->h_x); // h_x
    mzd_add(result->h_x_transpose,c2->h_x_transpose,c1->h_x_transpose); // h_x_transpose
    mzd_add(result->h_z,c2->h_z,c1->h_z); // h_z
    mzd_addmul_naive(result->h_z,c1->F,c2->h_x); // TODO replace with transpose
    add_and_column_vector(result->h_z,c1->F_diagonal,c2->F_diagonal);
    return result;
  } // doesnt transpose. need add_and to work on column matrices though. keep track of h_x and h_x transpose? should be able to avoid all transposing when multiplying dcp operators.

  // 2. Turn a DCP into a stabiliser state
  void set_stabiliser_state_to_dcp(stabiliser_state_t* state, dcp_t* dcp){
    lwlog_debug("setting stabiliser state to dcp");
    mzd_copy(state->epsilon_vector,dcp->h_z);
    mzd_set_ui(state->x_parts,1);
    mzd_copy(state->z_parts,dcp->F);
    mzd_copy(state->delta_vector,dcp->F_diagonal); //? think it works
  } // no transposing here either... note that as F is symmetric, and x the identity, state can choose which transposition convention it wants. It probably shouldn't, but either way this shows that this method will never ever transpose.

  // 3. Carry out measurements on a stabiliser state.
  // First, we'll need to be able to reorder the generators according to some specified row operations.
  void left_multiply(stabiliser_state_t* state, mzd_t* M /*n x n*/, stabiliser_state_t* new_state){
    // Performs elementary row operations on the stabiliser generators of a state by multiplying on the left by a matrix M and stores the results in a different state.
    mzd_mul_naive(temp_n_n_1,state->delta_vector,mzd_transpose(0,state->delta_vector)); // delta delta^T
    mzd_mul_naive(temp_n_n_2,state->x_parts,mzd_transpose(0,state->z_parts)); //SX S^Z^T
    // do the multiplications
    mzd_mul_naive(new_state->x_parts,M,state->x_parts);
    mzd_mul_naive(new_state->z_parts,M,state->z_parts);
    mzd_mul_naive(new_state->delta_vector,M,state->delta_vector);
    mzd_mul_naive(new_state->epsilon_vector,M,state->epsilon_vector);
    // TODO fix the epsilons using the delta and the a storage (important)
  } // lots of transposing.
  void swap_rows(stabiliser_state_t* state, int i, int j){
    mzd_row_swap(state->x_parts,i,j);
    mzd_row_swap(state->z_parts,i,j);
    mzd_row_swap(state->delta_vector,i,j);
    mzd_row_swap(state->epsilon_vector,i,j);
  }
  void add_row(stabiliser_state_t* state, int sourcerow, int destrow){
    // Adds rows and stores them in destrow.
    BIT delta1 = mzd_read_bit(state->delta_vector,sourcerow,0);
    BIT delta2 = mzd_read_bit(state->delta_vector,destrow,0);

    //char epsilon1 = mzd_read_bit(state->epsilon_vector,i,0);
    //char epsilon2 = mzd_read_bit(state->epsilon_vector,j,0);

    mzd_row_add(state->x_parts, sourcerow, destrow);
    mzd_row_add(state->z_parts, sourcerow, destrow);
    mzd_row_add(state->delta_vector, sourcerow, destrow);
    mzd_row_add(state->epsilon_vector, sourcerow, destrow);
    // correct epsilon
    BIT a2TUa1 = row_product(state->x_parts,sourcerow,state->z_parts,destrow); // this is symmetric under x <-> z for stabiliser states actually
    BIT epsilon_bit = mzd_read_bit(state->epsilon_vector,destrow,0);
    mzd_write_bit(state->epsilon_vector,destrow,0,epsilon_bit ^ (delta1 & delta2) ^ a2TUa1);
  }


  // Postselect a stabiliser state on a set of X measurements.
  int postselect_on_measuring_plus_one_x(stabiliser_state_t* state, mzd_t* x_measurement_matrix /* n_x x n; rows are paulis*/
                                         ){
    //int n_x = x_measurement_matrix->nrows;
    // int n = z_measurement_matrix->ncols;
    //  for(int i = 0; i < n_x; i++){
      // mzd_t* x_measurement = mzd_init(1,n); // row vector
      
      // }//
    //    return 0;
    // TODO
    return 0;
  }
  // Postselect a stabiliser state on a set of Z measurements.
  int postselect_on_measuring_plus_one_z(stabiliser_state_t* state, mzd_t* z_measurement_matrix /* n_z x n*/
                                         ){
    // TODO
    return 0;
  }
  // 4. Turn a stabiliser state into a coset. Idea: reduce size of stabiliser state by means of the permutation first? you know it will have zero support on ancilla qubits. what's the point in doing too much pluq. you can do it simply on the reduced data indices.

    void set_coset_to_stabiliser_state(pauli_z_coset_t *coset, stabiliser_state_t *state,
                                      mzd_t* temp /* n x n */, mzd_t* temp_transpose /* n x n */, mzp_t* p/* n */, mzp_t *q /* n  */, mzd_t* augmented_matrix /*n x (n+1)*/, mzd_t* augmented_transpose /* n+1 x n */, stabiliser_state_t *new_state){
    int n = STATE_N(state);
    // 1. Do Gaussian elimination on the Z parts of the stabiliser generators.
    mzd_copy(temp,state->z_parts); // copy over the Z parts
    int r = mzd_pluq(temp,p,q,0);
    lwlog_debug("PLUQ done"); // L is n x r, U is r x n.
    // The rref = U Q is a basis for the vector space of the coset. Copy it over
    coset->set_generators->nrows=n;
    mzd_copy(coset->set_generators,temp);
    // only want the upper diagonal part in order to get U.
    make_upper_diagonal(r,coset->set_generators);
    mzd_apply_p_right(coset->set_generators,q);
    coset->set_generators->nrows=r;// Coset generators now equal to UQ^T. Careful: lower rows may be garbage.

    // Now need to calculate L^-1 P^T - the row ops.
    mzd_transpose(temp_transpose,temp); // TODO can you do transpose in place? unsure if worth it. update: looks like you can for small matrices at least
    make_upper_diagonal(r,temp_transpose); // temp_transpose now contains L^T
    fill_rows_with_identity(temp_transpose,r,n); // temp_transpose now contains L^T a unit n x n upper triangular matrix
    mzd_trtri_upper(temp_transpose); // temp_transpose now contains (L^T)^-1.
    // apply P on the left
    mzd_apply_p_left(temp_transpose,p);
    // Now temp_transpose contains P (L^T)^-1 = ( L^-1 P^T )^T. So the transpose of what we wanted.
    // Apply the row_ops to get the new stabilisers
    mzd_transpose(temp,temp_transpose);
    left_multiply(state,temp,new_state);
    state=new_state;
    // Need to find a Pauli which commutes or anticommutes with the X parts of the stabilisers as appropriate.
    // Use gaussian elimination on an augmented matrix to do this.
    // TODO set up augmented matrix here. copy in sign bits as appropriate
    mzd_set_ui(augmented_matrix,0);
    // Figure out how to copy submatrices better
    for(int i = r; i < n; i ++){
      mzd_write_bit(augmented_matrix,i-r,n,mzd_read_bit(state->epsilon_vector,0,i));
      for(int j = 0; j < n; j++){
        mzd_write_bit(augmented_matrix,i-r,j,mzd_read_bit(state->x_parts,i,j));
      }
    }
    mzd_echelonize_naive(augmented_matrix,1);
    // Final column of augmented matrix now contains a solution. transpose it and get it into factor
    mzd_transpose(augmented_transpose,augmented_matrix); // augmented transpose is n+1 x n. final row is factor
    mzd_copy_row(coset->factor_transpose,0,augmented_transpose,n); // or are the remaining columns a basis for the null space? well, we got them from the gaussian elimination anyway. but maybe not necessary to copy them over above then
  }

#ifdef __cplusplus
}
#endif
