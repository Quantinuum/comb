#ifdef __cplusplus
extern "C" {
#endif
#include "comb.h"


  dcp_t* init_dcp(size_t n){
    dcp_t* tableau = (dcp_t*)malloc(sizeof(dcp_t));
    // follows demoor and dehaene
    tableau->F=mzd_init(n,n);
    tableau->h_x=mzd_init(n,1);
    tableau->h_z=mzd_init(n,1);
    return tableau;
  }

  void init_dcp_to_size(dcp_t* tableau, size_t n){
    tableau->F=mzd_init(n,n);
    tableau->h_x=mzd_init(n,1);
    tableau->h_z=mzd_init(n,1);
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

  void free_dcp(dcp_t* dcp){
    mzd_free(dcp->F);
    mzd_free(dcp->h_x);
    mzd_free(dcp->h_z);
    free(dcp);
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

  dcp_t* dcp_pauli_x(size_t i, size_t n){
    dcp_t* dcp = init_dcp(n);
    mzd_write_bit(dcp->h_z,i,0,1);
    return dcp;
  }
  dcp_t* dcp_pauli_z(size_t i, size_t n){
    dcp_t* dcp = init_dcp(n);
    mzd_write_bit(dcp->h_x,i,0,1);
    return dcp;
  }
  dcp_t* dcp_cz(size_t i, size_t j,size_t n){
    dcp_t* dcp = init_dcp(n);
    mzd_write_bit(dcp->F,i,j,1);
    mzd_write_bit(dcp->F,j,i,1);
    return dcp;
  }
  dcp_t* dcp_s(size_t i, size_t n){
    dcp_t* dcp = init_dcp(n);
    mzd_write_bit(dcp->F,i,i,1);
    mzd_write_bit(dcp->h_x,i,0,1);
    return dcp;
  }
  dcp_t* dcp_transvection(pauli_t* z_pauli,size_t n){
    // z_pauli assumed to be of z type.
    // Proved in dehaene and demoor in the appendix.
    dcp_t* dcp = init_dcp(n);
    for(int i = 0; i < n; i++){
      for(int j =0 ; j <n ; j++){
        mzd_write_bit(dcp->F,i,j,mzd_read_bit(z_pauli->data,i,0) & mzd_read_bit(z_pauli->data,j,0));
      }
    }
    return dcp;
  }
  void set_dcp_to_transvection(dcp_t* dcp, pauli_t* z_pauli){
    int n = PAULI_N(z_pauli);
    for(int i = 0; i < n; i++){
      for(int j =0 ; j <n ; j++){
        mzd_write_bit(dcp->F,i,j,mzd_read_bit(z_pauli->data,i,0) & mzd_read_bit(z_pauli->data,j,0));
      }
    }
    for(int i =0; i < n; i++){
      mzd_write_bit(dcp->h_x,i,0,0);
      mzd_write_bit(dcp->h_z,i,0,0);
    }
  }
  void set_dcp_to(dcp_t* a, dcp_t *b){
    // assumed to be of same size
    mzd_copy(a->F,b->F);
    mzd_copy(a->h_x,b->h_x);
    mzd_copy(a->h_z,b->h_z);
  }

#ifdef __cplusplus
}
#endif
