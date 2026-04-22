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

    mzd_row_add(state->x_parts, sourcerow, destrow);
    mzd_row_add(state->z_parts, sourcerow, destrow);
    mzd_row_add(state->delta_vector, sourcerow, destrow);
    mzd_row_add(state->epsilon_vector, sourcerow, destrow);
    // correct epsilon
    BIT a2TUa1 = row_product(state->x_parts,sourcerow,state->z_parts,destrow); // this is symmetric under x <-> z for stabiliser states actually
    BIT epsilon_bit = mzd_read_bit(state->epsilon_vector,destrow,0);
    mzd_write_bit(state->epsilon_vector,destrow,0,epsilon_bit ^ (delta1 & delta2) ^ a2TUa1);
  }
