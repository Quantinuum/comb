#ifdef __cplusplus
extern "C" {
#endif
#include "comb.h"


  clifford_t* init_clifford(size_t n){
    clifford_t* tableau = (clifford_t*)malloc(sizeof(clifford_t));
    // follows demoor and dehaene
    tableau->c_bar=mzd_init(1+(2*n),1+(2*n));
    tableau->h_bar=mzd_init(1+(2*n),1);
    return tableau;
  }

  pauli_t* init_pauli(size_t n){
    pauli_t* pauli = (pauli_t*)malloc(sizeof(pauli_t));
    pauli->sign_bit=0;
    pauli->data=mzd_init(1+(2*n),1);
    return pauli;
  }
  stabiliser_state_t* init_stabiliser_state(size_t n){
    stabiliser_state_t* state = (stabiliser_state_t*)malloc(sizeof(stabiliser_state_t));
    state->data=mzd_init(n,1+(2*n));
    state->h_bits=mzd_init(n,1);
    return state;
  }
  pauli_z_coset_t* init_pauli_z_coset(size_t n, size_t k){
    pauli_z_coset_t* coset = (pauli_z_coset_t*)malloc(sizeof(pauli_z_coset_t));
    coset->factor=mzd_init(n,1);
    coset->set_generators=mzd_init(n,k);
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
  void free_stabiliser_state(stabiliser_state_t* state){
    mzd_free(state->data);
    mzd_free(state->h_bits);
    free(state);
  }
  void free_pauli_z_coset(pauli_z_coset_t* coset){
    mzd_free(coset->set_generators);
    mzd_free(coset->factor);
    free(coset);
  }


#ifdef __cplusplus
}
#endif
