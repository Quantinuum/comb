#ifdef __cplusplus
extern "C" {
#endif
#include "comb.h"


  pauli_z_coset_t *diagonal_clifford_to_coset(clifford_t *clifford){
    int n = CLIFFORD_N(clifford);
    // Work out the stabilisers
    pauli_t** stabilisers=malloc(n*sizeof(pauli_t*));
    for(int i = 0 ; i< n; i++){
      pauli_t* x_i = pauli_x(i,n);
      stabilisers[i] = clifford_conjugate_pauli(0,clifford,x_i);
      free_pauli(x_i);
    }
    // Do Gaussian elimination on the Z parts of the stabilisers by means of a PLE decomposition
    // 1. Fill up a matrix with the Z parts of the stabilisers
    mzd_t* z_parts = mzd_init(n,n);
    for(int i = 0; i< n; i++){
      for(int j = 0; j < n; j++){
        mzd_write_bit(z_parts,i,j,mzd_read_bit(stabilisers[i]->data,j,0));
      }
    }
    // 2. Do the Gaussian elimination by performing a PLUQ decomposition
    mzp_t* p = mzp_init(n);
    mzp_t* q = mzp_init(n);
    int k = mzd_pluq(z_parts,p,q,0);
    mzd_t* U  = mzd_init(n,n);
    mzd_t* L  = mzd_init(n,n);
    // from m4ri examples
    for (rci_t i = 0; i < n; ++i){
      for (rci_t j = 0; j < i && j < n;++j)
        mzd_write_bit(U,i,j, 0);
      for (rci_t j = i + 1; j < n; ++j)
        mzd_write_bit(L,i,j, 0);
      if(i<n)
        mzd_write_bit(U,i,i, 1);
      mzd_write_bit(L,i,i, 1);
    }
    // 3. Work out the new stabilisers by applying the PLUQ decomposition
    pauli_t** new_stabilisers=malloc(k*sizeof(pauli_t*));
    pauli_z_coset_t* coset = init_pauli_z_coset(n,k);

    for(int i = 0; i < k; i++){
      new_stabilisers[i]=init_pauli(n);
      for(int j = 0; j < n; j++){
        mzd_write_bit(new_stabilisers[i]->data,j,0,mzd_read_bit(U,i,j));
      }
      coset->set_generators[i]=new_stabilisers[i];
    }
    return coset;
  }



#ifdef __cplusplus
}
#endif
