#ifdef __cplusplus
extern "C" {
#endif
#include "comb.h"

  void add_and_row_vector(mzd_t* result, mzd_t* F1,mzd_t* F2){
    // Adds F1 & F2 to result. Assumes only one row. must be initialised to right shape
    for(int i = 0; i < result->width; i++){
      result->data[i]^=(F1->data[i] & F2->data[i]);
    }
    result->data[result->width - 1] &= result->high_bitmask; // on the safe side
  }

  void mask_row_vector(mzd_t* result, mzd_t* mask){
    // Adds F1 & F2 to result. Assumes only one row. must be initialised to right shape
    for(int i = 0; i < result->width; i++){
      result->data[i]&=(result->data[i] & mask->data[i]);
    }
    result->data[result->width - 1] &= result->high_bitmask; // on the safe side
  }
  void add_and_column_vector(mzd_t* result, mzd_t* F1, mzd_t* F2){
    // TODO
  }

  void incf_coset(pauli_z_coset_t *coset, int index){
    // Iterates over the coset by adding one of the generators, chosen by a gray code, to factor.
    // !!! won't work with less than 64 bits (reasonable assumption if we're using lookup tables)
    coset->factor_transpose->data[0]^=mzd_row(coset->set_generators,gray_bit_to_switch(index))[0];
  }


#ifdef __cplusplus
}
#endif
