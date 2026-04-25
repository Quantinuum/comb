#include "comb.h"

void fix_dcp(dcp_t* dcp){
  // shouldn't normally need to be called, except for testing purposes
  int n = DCP_N(dcp);
  for(int i = 0; i < n; i++){
    mzd_write_bit(dcp->F_diagonal,i,0,mzd_read_bit(dcp->F,i,i));
  }
  mzd_transpose(dcp->h_x_transpose,dcp->h_x);
}

BIT row_product(mzd_t* x, int i, mzd_t* z, int j){
  // Forms the inner product of row i of x and row j of z (viewing both rows as vectors over F2).
  word* row1 = mzd_row(x,i);
  word* row2 = mzd_row(z,j);
  wi_t const wide     = x->width;
  word acc = 0;
  for(int i = 0; i < wide; i++){
    word part_of_product = row1[i] & row2[i];
    acc+=__builtin_popcount(part_of_product);
  }
  return acc&1;
}

void expand_state(stabiliser_state_t* state){
  state->x_parts->nrows++;
  state->z_parts->nrows++;
  state->delta_vector->nrows++;
  state->epsilon_vector->nrows++;
}

void contract_state(stabiliser_state_t* state){
  state->x_parts->nrows--;
  state->z_parts->nrows--;
  state->delta_vector->nrows--;
  state->epsilon_vector->nrows--;
}
void set_row(mzd_t* x_parts, mzd_t* z_parts, int j, mzd_t* x_measurement_matrix,
             int i, mzd_t* epsilon_bits, mzd_t* delta_bits){ // moves m_i to row j of x_parts, and clears row j of z_parts
  mzd_copy_row(x_parts,j,x_measurement_matrix,i);
  mzd_row_clear_offset(z_parts,j,0);
  mzd_write_bit(epsilon_bits,j,0,0);
  mzd_write_bit(delta_bits,j,0,0);
}

inline unsigned compress4(unsigned x, unsigned m) { // from hacker's delight
  unsigned mk, mp, mv, t;
  int i;
  x = x & m;           // Clear irrelevant bits.
  mk = ~m << 1;        // We will count 0's to right.
  for (i = 0; i < 5; i++) {
    mp = mk ^ (mk << 1);              // Parallel suffix.
    mp = mp ^ (mp << 2);
    mp = mp ^ (mp << 4);
    mp = mp ^ (mp << 8);
    mp = mp ^ (mp << 16);
    mv = mp & m;                      // Bits to move.
    m = (m ^ mv) | (mv >> (1 << i));    // Compress m.
    t = x & mv;
    x = (x ^ t) | (t >> (1 << i));      // Compress x.
    mk = mk & ~mp;
  }
  return x;
}
void permute_columns(mzd_t* data, unsigned mask){
  for(int i = 0 ; i < data->nrows; i++){
    word *b = mzd_row(data,i);
    // assume small enough so that only have to do the first 64 bits
    b[0]=compress4(b[0],mask);
  }
}
void convert_row(mzd_t* M, int row, comb_vector_t* vector){
  for(int i = 0; i < M->ncols; i++){
    mzd_write_bit(M,row,i,TestBit(vector->data,i)?1:0);
  }
}

void add_and_column_vector(mzd_t* result, mzd_t* F1, mzd_t* F2){
  // assumes that f1, f2, and result are congruous column vectors
  // necessarily linear time.
  int n = result->nrows;
  for(int i = 0; i < n; i++){
    mzd_write_bit(result,i,0,(mzd_read_bit(result,i,0) ^ (mzd_read_bit(F1,i,0) & mzd_read_bit(F2,i,0))));
  }
}
