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


  int _postselect_on_x(stabiliser_state_t* state, mzd_t* x_parts, mzd_t* z_parts, mzd_t* x_measurement_matrix, mzd_t* delta_bits, mzd_t* epsilon_bits){
    int n_x = x_measurement_matrix->nrows;
    int n = x_measurement_matrix->ncols;
    //mzd_t* m_i = mzd_init(1,n); // row vector
    //mzd_t* z_j = mzd_init(1,n); // current stabiliser: row vector
      for(int i = 0; i < n_x; i++){
        lwlog_debug("measuring %d -th operator on state %s",i,print_stabiliser_state(state,1));
        //mzd_copy_row(m_i,0,x_measurement_matrix,i); // get the measurement out
        int anticommuting_stabiliser_found = 0;
        for(int j = 0; j < n && !anticommuting_stabiliser_found; j++){
          //mzd_copy_row(z_j,0,z_parts,j);
          int anticommuting_p = row_product(x_measurement_matrix,i,z_parts,j);
          if(anticommuting_p){ // An anticommuting stabiliser was found
            anticommuting_stabiliser_found=1;
            for(int k = j+1; k < n; k++){
              if(row_product(z_parts,k,x_measurement_matrix,i)){ // if k-th row anticommutes with the measurement operator...
                  add_row(state,j,k); // add row j to k
                }
            }
            // now set the row equal to the measurement operator.
            set_row(x_parts,z_parts,j,x_measurement_matrix,i,epsilon_bits, delta_bits);
          }
        }
        if(!anticommuting_stabiliser_found){
          // if no stabiliser was found, we need to find out whether the measurement operator or its opposite is generated by the stabiliser group. otherwise, we can continue to the next operator.
          expand_state(state); // ensure state has enough room for n+1 stabiliser generators
          set_row(x_parts,z_parts,n,x_measurement_matrix,i, epsilon_bits, delta_bits); // copy measurement operator to last row
          // now do gaussian elimination to ensure that final row is either -I or I. if this is not possible signal an error
          clear_rows(state,x_parts,n);
          clear_rows(state,z_parts,n);

          for(int i = 0; i < n; i++){
            if(mzd_read_bit(x_parts,n,i))printf("warning: x_parts %d not zero",i);
            if(mzd_read_bit(z_parts,n,i))printf("warning: z_parts %d not zero",i);
          }
          BIT identity_sign = mzd_read_bit(epsilon_bits,n,0);
          contract_state(state);
          if(identity_sign) // -I was formed so we postselect. Otherwise, have to continue with other measurements.
            return 1;
        }
       }
      return 0;
  }

  // Postselect a stabiliser state on a set of X measurements.
  int postselect_on_measuring_plus_one_x(stabiliser_state_t* state, mzd_t* x_measurement_matrix /* n_x x n; rows are paulis*/
                                         ){
    return _postselect_on_x(state, state->x_parts, state->z_parts, x_measurement_matrix, state->delta_vector, state->epsilon_vector);
  }



  // Postselect a stabiliser state on a set of Z measurements.
  int postselect_on_measuring_plus_one_z(stabiliser_state_t* state, mzd_t* z_measurement_matrix /* n_z x n*/
                                         ){
    return _postselect_on_x(state, state->z_parts, state->x_parts, z_measurement_matrix, state->delta_vector, state->epsilon_vector);
  }
  // 4. Turn a stabiliser state into a coset. Idea: reduce size of stabiliser state by means of the permutation first? you know it will have zero support on ancilla qubits. what's the point in doing too much pluq. you can do it simply on the reduced data indices.
#define VERIFY_COSET_RREF 1
  int row_empty_p(mzd_t* M, int row){
    for(int col = 0; col < M->ncols; col++)
      if(mzd_read_bit(M,row,col))return 0;
    return 1;
}
  void set_coset_to_stabiliser_state(pauli_z_coset_t *coset, stabiliser_state_t *state){
    int n = STATE_N(state);
    int r = n;
    mzd_t* row_mask = mzd_init(1,n);
    lwlog_debug("State before gaussian elimination:\n %s",print_stabiliser_state(state,1));
    // Gaussian elimination on Z parts - ensuring that last n-r stabilisers have no Z parts
    for(int row = 0; (r>0) && row < n; row++){
      // search for leading entry.
      int leading_column=-1;
      for(int col = 0; (leading_column==-1) && col < n; col++){
        if(mzd_read_bit(state->z_parts,row,col)){
          leading_column = col;
        }
      }
      if(leading_column==-1){
        // row was zero. swap it to the end and try again.
        do {
          r--;
        } while(r > 0 && row_empty_p(state->z_parts,r));
        swap_rows(state,row,r);
        row--; // have to now check the new row again, so decf row
      }
      else{
        // do gaussian elimination on the column.
        for(int row2 = 0; row2 < n; row2++){
          if(row2!=row){
            if(mzd_read_bit(state->z_parts,row2,leading_column)){
              add_row(state,row,row2);
            }
          }
        }
      }
    }

    lwlog_debug("State after gaussian elimination on Z parts: r = %d\n %s",r,print_stabiliser_state(state,1));
    // Gaussian elimination on X parts of last n-r stabilisers
    for(int col = 0; col < n; col++){
      int row_to_add=-1;
      for(int row = r; row < n; row++){
        if(mzd_read_bit(state->x_parts,row,col)){
          mzd_write_bit(row_mask,0,col,1);
          if(row_to_add==-1){ // no row yet found with a 1 in this column; if there are any others, then add this one to them.
            row_to_add=row;
          }
          else{
            add_row(state,row,row_to_add);
          }
        }
      }
      // column should only have one 1 now.
    }
    lwlog_debug("State after gaussian elimination on X parts:\n %s",print_stabiliser_state(state,1));
    // optional: verify state is in rref
    #if VERIFY_COSET_RREF
    // verify that z parts, for all i greater than r, vanish
    for(int row = r; row < n; row++){
      for(int col = 0; col < n; col++){
        if(mzd_read_bit(state->z_parts,row,col)){
          lwlog_err( "state->z_parts[%d]{%d] not zero",row,col);
          exit(1);
        }
      }
    }
    // verify that x parts, for all greater than r, are in rref
    // TODO
    #endif
    // Now read off results and store them in coset
    coset->set_generators->nrows=n; // bring back up to size
    mzd_copy(coset->set_generators,state->z_parts);
    coset->set_generators->nrows=r; // only r of them are valid
    mzd_set_ui(coset->factor_transpose,0);
    for(int i = r; i < n; i++){
      if(mzd_read_bit(state->epsilon_vector,i,0))
        mzd_combine(coset->factor_transpose,0,0,
                    coset->factor_transpose,0,0,
                    state->x_parts,i,0); // if the sign bit is 1, add it to the coset element
    }
    mask_row_vector(coset->factor_transpose,row_mask);
    mzd_free(row_mask);
  }

#ifdef __cplusplus
}
#endif
