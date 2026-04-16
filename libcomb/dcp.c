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
    _mzd_addmul_naive(result->h_z,c1->F,c2->h_x_transpose);
    add_and(result->h_z,c1->F_diagonal,c2->F_diagonal);
    return result;
  } // transposes once. keep track of h_x and h_x transpose? should be able to avoid all transposing when multiplying dcp operators.

  // 2. Turn a DCP into a stabiliser state
  void set_stabiliser_state_to_dcp(stabiliser_state_t* state, dcp_t* dcp){
    lwlog_debug("setting stabiliser state to dcp");
    int n = DCP_N(dcp);
    mzd_copy(state->epsilon_bits,dcp->h_z);
    mzd_set_ui(state->x_parts,1);
    mzd_copy(state->z_parts,dcp->F);
    mzd_copy(state->delta_vector,dcp->F_diagonal); //? think it works
  } // no transposing here either... note that as F is symmetric, and x the identity, state can choose which transposition convention it wants...

  // 3. Carry out measurements on a stabiliser state.
  // First, we'll need to be able to reorder the generators according to some specified row operations.
  void multiply_stabiliser_generators(stabiliser_state_t* state, mzd_t* M /*n x n*/, mzd_t* temp_delta_storage /* n x n*/, mzd_t* temp_a_storage /* n x n*/, stabiliser_state_t* new_state){
    // Performs elementary row operations on the stabiliser generators of a state by multiplying on the left by a matrix M and stores the results in a different state.
    mzd_mul_naive(temp_delta_storage,state->delta_vector,mzd_transpose(state->delta_vector)); // delta delta^T
    mzd_mul_naive(temp_a_storage,state->x_parts,mzd_transpose(state->z_parts)); //SX S^Z^T
    // do the multiplications
    mzd_mul_naive(new_state->x_parts,M,state->x_parts);
    mzd_mul_naive(new_state->z_parts,M,state->z_parts);
    mzd_mul_naive(new_state->delta_vector,M,state->delta_vector);
    mzd_mul_naive(new_state->epsilon_vector,M,state->epsilon_vector);
    // TODO fix the epsilons using the delta and the a storage (important)
  } // lots of transposing.
  // Postselect a stabiliser state on a set of X measurements.
  int postselect_on_measuring_plus_one_x(stabiliser_state_t* state, mzd_t* x_measurement_matrix /* n_x x n*/
                                         ){
    // TODO
  }
  // Postselect a stabiliser state on a set of Z measurements.
  int postselect_on_measuring_plus_one_z(stabiliser_state_t* state, mzd_t* z_measurement_matrix /* n_z x n*/
                                         ){
    // TODO
  }
  // 4. Turn a stabiliser state into a coset.
    void *set_coset_to_stabiliser_state(pauli_z_coset_t *coset, stabiliser_state_t *state,
                                      mzd_t* temp /* n x n */, mzd_t* temp_transpose /* n x n */, mzp_t* p/* n */, mzp_t *q /* n  */, mzd_t* augmented_matrix /*n x (n+1)*/, mzd_t* augmented_transpose /* n+1 x n */, stabiliser_state_t *new_state){
    int n = STABILISER_N(state);
    // 1. Do Gaussian elimination on the Z parts of the stabiliser generators.

    mzd_copy(temp,state->z_parts); // copy over the Z parts

    int r = mzd_pluq(temp,p,q,0);
    lwlog_debug("PLUQ done"); // L is n x r, U is r x n.

    // The rref = U Q is a basis for the vector space of the coset. Copy it over
    coset->set_generators.nrows=n;
    mzd_set_ui(coset->set_generators,0);
    mzd_copy(coset->set_generators,temp);
    // only want the upper diagonal part in order to get U.
    for(int i = 0; i < r; i++){
      mzd_write_bit(coset->set_generators,i,i,1);
      for(int j = 0; j < i ; j++){
        mzd_write_bit(coset->set_generators,i,j,0); // TODO use clear bits on row to speed it up?
      }
    }
    mzd_apply_p_right(coset->set_generators,q);
    coset->set_generators.nrows=r;// Coset generators now equal to UQ^T. Careful: lower rows may be garbage.

    // Now need to calculate L^-1 P^T - the row ops.

    mzd_transpose(temp_transpose,temp); // TODO can you do transpose in place? unsure if worth it. update: looks like you can for small matrices at least
    for(int i = 0; i < n; i++){
      mzd_write_bit(temp_transpose,i,i,1);
      for(int j = 0; j < i ; j++){
        mzd_write_bit(temp_transpose,i,j,0); // TODO use clear bits on row to speed it up?
      }
    } // temp_transpose now contains L^T - a n x n upper triangular matrix. the lower rows have been filled up with the identity.
    mzd_trtri_upper(temp_transpose); //temp_transpose now contains (L^T)^-1.
    // apply P on the left
    mzd_apply_p_left(temp_transpose,p);
    // Now temp_transpose contains P (L^T)^-1 = ( L^-1 P^T )^T. So the transpose of what we wanted.
    // Apply the row_ops to get the new stabilisers
    mzd_transpose(temp,temp_transpose);
    multiply_stabiliser_generators(state,temp,new_state);
    state=new_state;

    // Need to find a Pauli which commutes or anticommutes with the X parts of the stabilisers as appropriate.
    // Use gaussian elimination on an augmented matrix to do this.
    // TODO set up augmented matrix here. copy in sign bits as appropriate
    mzd_set_ui(augmented_matrix,0);
    // Figure out how to copy submatrices better
    for(int i = r; i < n; i ++){
      mzd_write_bit(augmented_matrix,i-r,n,mzd_read_bit(state->sign_vector,0,i));
      for(int j = 0; j < n; j++){
        mzd_write_bit(augmented_matrix,i-r,j,mzd_read_bit(state->x_parts,i,j));
      }
    }
    mzd_echelonize_naive(augmented_matrix,1);
    // Final column of augmented matrix now contains a solution. transpose it and get it into factor
    mzd_transpose(augmented_transpose,augmented_matrix); // augmented transpose is n+1 x n. final row is factor
    mzd_copy_row(coset->factor_transpose,0,augmented_transpose,n); // or are the remaining columns a basis for the null space? well, we got them from the gaussian elimination anyway. but maybe not necessary to copy them over above then
    return coset;
  }

#ifdef __cplusplus
}
#endif
