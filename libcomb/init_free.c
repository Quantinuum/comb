#ifdef __cplusplus
extern "C" {
#endif
#include "comb.h"


  clifford_t* init_clifford(size_t n){
    clifford_t* tableau = (clifford_t*)malloc(sizeof(clifford_t));
    // follows demoor and dehaene
    tableau->c_bar=mzd_init(1+(2*n),1+(2*n));
    tableau->h_bar=mzd_init(1+(2*n),1);
    for(int i = 0 ;i< 1+ (2*n);i++)mzd_write_bit(tableau->c_bar,i,i,1);
    return tableau;
  }

  pauli_t* init_pauli(size_t n){
    pauli_t* pauli = (pauli_t*)malloc(sizeof(pauli_t));
    pauli->sign_bit=0;
    pauli->data=mzd_init(1+(2*n),1);
    return pauli;
  }
  pauli_z_coset_t* init_pauli_z_coset(size_t n, size_t k){
    pauli_z_coset_t* coset = (pauli_z_coset_t*)malloc(sizeof(pauli_z_coset_t));
    coset->factor=init_pauli(n);
    coset->k=k;
    coset->set_generators=(pauli_t**)malloc(k*sizeof(pauli_t*));
    for(int i = 0;i < k; i++){
      coset->set_generators[i] = init_pauli(n);
    }
    return coset;
  }

  void free_clifford(clifford_t* clifford){
    mzd_free(clifford->c_bar);
    mzd_free(clifford->h_bar);
    free(clifford);
  }
  void free_pauli(pauli_t* pauli){
    mzd_free(pauli->data);
    free(pauli);
  }
  void free_pauli_z_coset(pauli_z_coset_t* coset){
    for(int i = 0 ; i < coset-> k ; i++){
      free_pauli(coset->set_generators[i]);
    }
    free_pauli(coset->factor);
    free(coset);
  }

#ifdef __cplusplus
}
#endif
