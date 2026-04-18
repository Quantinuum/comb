#ifdef __cplusplus
extern "C" {
#endif
#include "comb.h"


  void make_upper_diagonal(int r, mzd_t* M){
    // Given M a matrix, clear the lower diagonal parts to ensure that M is an upper diagonal r x n matrix. Useful to obtain U from a PLUQ decomposition.
    for(int i = 0; i < r; i++){
      mzd_write_bit(M,i,i,1);
      for(int j = 0; j < i ; j++){
        mzd_write_bit(M,i,j,0); // TODO use clear bits on row to speed it up? Or maybe a massive and mask
      }
    }
  }

  void fix_dcp(dcp_t* dcp){
    // shouldn't normally need to be called, except for testing purposes
    int n = DCP_N(dcp);
    for(int i = 0; i < n; i++){
      mzd_write_bit(dcp->F_diagonal,i,0,mzd_read_bit(dcp->F,i,i));
    }
    mzd_transpose(dcp->h_x_transpose,dcp->h_x);
  }

  void fill_rows_with_identity(mzd_t* temp_transpose, int r, int n){
    for(int i = r; i < n; i++){ // Fill the r+1 - n th rows with the identity.
      for (int j = 0; j < n; j++){
        mzd_write_bit(temp_transpose,i,j,(i==j)?1:0);
      }
    }
  }

  BIT row_product(mzd_t* x, int i, mzd_t* z, int j){
    // Forms the inner product of row i of x and row j of z (viewing both rows as vectors over F2).
    // x and z are assumed to have the same number of columns.
    // TODO this should admit a very fast implementation using mzd_combine.
    int acc = 0;
    int n = z->ncols;
    for(int k = 0; k < n; k++){
      if(mzd_read_bit(x,i,k) && mzd_read_bit(z,j,k))acc^=1;
    }
    return acc;
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
  void set_row(mzd_t* x_parts, mzd_t* z_parts, int j, mzd_t* x_measurement_matrix,int i,mzd_t* epsilon_bits, mzd_t* delta_bits){ // moves m_i to row j of x_parts, and clears row j of z_parts
    mzd_copy_row(x_parts,j,x_measurement_matrix,i);
    mzd_row_clear_offset(z_parts,j,0);
    mzd_write_bit(epsilon_bits,j,0,0);
    mzd_write_bit(delta_bits,j,0,0); // we are done, so return
  }
  void clear_rows(stabiliser_state_t* state, mzd_t* parts, int n){ // cubic?
    for(int col = 0; col < n; col++){
      // go down the rows until you find one with a non zero entry.
      int row_found_p = -1;
      for(int row = 0; row < n && (row_found_p==-1); row++){
        if(mzd_read_bit(parts,row,col))row_found_p=row;
      }
      if(row_found_p!=-1){
        // if a row was found, eliminate all rows below it with elementary ops.
        // otherwise, can carry on with the loop
        for(int row_to_clear= row_found_p + 1; row_to_clear < n+1 ; row_to_clear++){
          if(mzd_read_bit(parts,row_to_clear,col))add_row(state,row_found_p,row_to_clear);
        }
      }
    }
  }

  uint32_t gray_bit_to_switch(uint32_t n) // allocation free */
     {
       uint32_t n_add_one = n + 1; 
       uint32_t diff = n ^ n_add_one ^ (n >> 1) ^ (n_add_one >> 1);
       return n ? 31 - __builtin_clz(diff) : 0;
     }

     uint32_t mzd_data(mzd_t* mzd){ // allocation free */
       return (uint32_t)mzd->data[0]; 
     } 
  mzd_t *temp_n_n_1,*temp_n_n_2,*temp_n_plus_one_n,*temp_n_n_plus_one;
  void init_temps(int n){
    temp_n_n_1=mzd_init(n,n);
    temp_n_n_2=mzd_init(n,n);
    temp_n_n_plus_one=mzd_init(n,n+1);
    temp_n_plus_one_n=mzd_init(n+1,n);
  }
  void permute_columns(mzd_t* data, unsigned mask){
    for(int i = 0 ; i < data->nrows; i++){
      word *b = mzd_row(data,i);
      // assume small enough so that only have to do the first 64 bits
      b[0]=compress4(b[0],mask);
    }
    
  }
  unsigned compress2(unsigned x, unsigned m) { // from hacker's delight
    unsigned r, s, b;    // Result, shift, mask bit.

    r = 0;
    s = 0;
    do {
      b = m & 1;
      r = r | ((x & b) << s);
      s = s + b;
      x = x >> 1;
      m = m >> 1;
    } while (m != 0);
    return r;
  }
  unsigned compress4(unsigned x, unsigned m) { // from hacker's delight
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
  void convert_row(mzd_t* M, int row, comb_vector_t* vector){
    for(int i = 0; i < M->ncols; i++){
      mzd_write_bit(M,row,i,TestBit(vector->data,i)?1:0);
    }
  }
#ifdef __cplusplus
}
#endif
