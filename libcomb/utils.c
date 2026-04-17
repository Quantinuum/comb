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
  void init_temps(){
    
  }


#ifdef __cplusplus
}
#endif
