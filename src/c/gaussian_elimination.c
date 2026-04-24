// row operations
#include "comb.h"

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

    BIT a2TUa1 = row_product(state->x_parts,sourcerow,state->z_parts,destrow); // this is symmetric under x <-> z for stabiliser states actually
    mzd_row_add(state->x_parts, sourcerow, destrow);
    mzd_row_add(state->z_parts, sourcerow, destrow);
    mzd_row_add(state->delta_vector, sourcerow, destrow);
    mzd_row_add(state->epsilon_vector, sourcerow, destrow);
    // correct epsilon according to dehaene and de moor
    BIT epsilon_bit = mzd_read_bit(state->epsilon_vector,destrow,0);
    mzd_write_bit(state->epsilon_vector,destrow,0,epsilon_bit ^ (delta1 & delta2) ^ a2TUa1);
  }

/* Attempts to form either I or -I from the rows of the extended stabiliser tableau via Gaussian elimination.
 * If LOG_LEVEL is greater than or equal to DEBUG, additional verifications are made to check that this is done
 * and exits otherwise.
 * Returns 0 if I was formed, and 1 if -I was formed.
 * The return value is undefined if neither I or -I were formed.
 */
int form_identity_ge(stabiliser_state_t* state){
  // if in debug mode, check manually
  lwlog_debug("in form_identity_ge:\n %s",print_stabiliser_state(state,1));
  int n = state->x_parts->ncols; // its extended so not a good idea to use STATE_N - its not square
  int row_to_start_eliminating=0;
  for(int col = 0; col < n; col++){
    // Try to eliminate starting from row_to_start_eliminating down. We dont want to add rows that have already been used to eliminate others
    int found_row=-1;
    for(int row = row_to_start_eliminating; row <= n; row++){
      if(mzd_read_bit(state->x_parts, row, col)){
        if(found_row == -1)found_row = row;
        else add_row(state, found_row, row);
      }
    }
    if(found_row == -1){
      // means we have a zero column. ignore it and carry on
    }
    else{
      // swap the found row with the starting row, so that its in rref
      swap_rows(state,found_row,row_to_start_eliminating);
      // now incf row_to_start eliminating so we dont do it again
      row_to_start_eliminating++;
        }
    //lwlog_debug("after trying to eliminate x_%d: \n %s",col,print_stabiliser_state(state,1));
    // same for z
    found_row=-1;
    for(int row = row_to_start_eliminating; row <= n; row++){
      if(mzd_read_bit(state->z_parts, row, col)){
        if(found_row == -1)found_row = row;
        else add_row(state, found_row, row);
      }
    }
    if(found_row == -1){
      // means we have a zero column. ignore it and carry on without incrementing the row to start eliminating
    }
    else{
      // swap the found row with the starting row, so that its in rref
      swap_rows(state,found_row,row_to_start_eliminating);
      // now incf row_to_start eliminating so we dont do it again
      row_to_start_eliminating++;
        }
    //lwlog_debug("after trying to eliminate z_%d: \n %s",col,print_stabiliser_state(state,1));
  }


  #if LOG_LEVEL >= DEBUG
  for(int i = 0; i < n; i++){
    if(mzd_read_bit(state->x_parts, n, i)){
      lwlog_err("warning: x_parts %d not zero",i);
      exit(1);
    }
    if(mzd_read_bit(state->z_parts, n, i)){
      lwlog_err("warning: z_parts %d not zero",i);
      exit(1);
    }
  }
  if(mzd_read_bit(state->delta_vector, n, 0)){
    lwlog_err("error: delta bits %d not zero",n);
    exit(1);
  }
  #endif
  BIT identity_sign = mzd_read_bit(state->epsilon_vector,n,0);
  return identity_sign;
}

#if LOG_LEVEL >= DEBUG
/* Checks to see if the state is in the canonical form required for transformation into a stabiliser state.
 * Returns 1 if it is, and 0 otherwise.
 */
int verify_coset_ge(int r, stabiliser_state_t* state){
  // verify that z parts, for all i greater than r, vanish
  int n = state->x_parts->ncols;
  for(int row = r; row < n; row++){
    for(int col = 0; col < n; col++){
      if(mzd_read_bit(state->z_parts,row,col)){
        lwlog_debug( "state->z_parts[%d]{%d] not zero",row,col);
        return 0;
      }
    }
  }
  // verify that z parts, for i less than r, are full rank (linearly independent)
  mzd_t* coeff = mzd_init(r,n);
  mzp_t* p = mzp_init(r);
  mzp_t* q = mzp_init(n);
  mzd_submatrix(coeff,state->z_parts,0,0,r,n);
  int rank = mzd_pluq(coeff,p,q,0);
  mzp_free(p);
  mzp_free(q);
  mzd_free(coeff);
  if(rank!=r){
    lwlog_debug("state->z_parts not linearly independent: rank = %d, r=%d",rank,r);
    return 0;
  }
  // verify that x parts, for all greater than r, are in rref
  int last_leading_column = -2;
  int started_to_find_zeros_p=0;
  for(int row = r; row < n; row++){
    int leading_column = -1;
    for(int col = 0; (leading_column==-1) && col < n; col++){
      if(mzd_read_bit(state->x_parts,row,col))leading_column=col;
    }
    if(leading_column!=-1){
      if(started_to_find_zeros_p){
        lwlog_debug("non trivial leading column found after a zero row was detected");
        return 0;
      }
      // check to see if its not shared with any other rows
      if(last_leading_column >= leading_column){
        lwlog_debug("leading column not to the right of a previous leading column %d, %d",leading_column, last_leading_column);
        return 0;
      }
      else last_leading_column=leading_column;
      for(int row2 = 0; row2 < n; row2++){
        if(row!=row2){
          if(mzd_read_bit(state->x_parts,row2,leading_column)){
            lwlog_debug( "state->x_parts[%d][%d], state->x_parts[%d]{%d] not zero",row,leading_column,row2,leading_column);
            return 0;
          }
        }
      }
    }
    else{
      started_to_find_zeros_p=1;
    }
  }
  return 1;
}
#endif
