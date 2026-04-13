#ifdef __cplusplus
extern "C" {
#endif
#include "comb.h"


  dcp_t* init_dcp(uint32_t n){
    dcp_t* tableau = (dcp_t*)malloc(sizeof(dcp_t));
    // follows demoor and dehaene
    tableau->F=mzd_init(n,n);
    tableau->h_x=mzd_init(n,1);
    tableau->h_z=mzd_init(n,1);
    return tableau;
  }

  void init_dcp_to_size(dcp_t* tableau, uint32_t n){
    tableau->F=mzd_init(n,n);
    tableau->h_x=mzd_init(n,1);
    tableau->h_z=mzd_init(n,1);
  }

  pauli_t* init_pauli(uint32_t n){
    pauli_t* pauli = (pauli_t*)malloc(sizeof(pauli_t));
    pauli->sign_bit=0;
    pauli->data=mzd_init(1+(2*n),1);
    return pauli;
  }
  pauli_z_coset_t* init_pauli_z_coset(uint32_t n, uint32_t k){
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

  dcp_t* dcp_pauli_x(uint32_t i, uint32_t n){
    dcp_t* dcp = init_dcp(n);
    mzd_write_bit(dcp->h_z,i,0,1);
    return dcp;
  }
  dcp_t* dcp_pauli_z(uint32_t i, uint32_t n){
    dcp_t* dcp = init_dcp(n);
    mzd_write_bit(dcp->h_x,i,0,1);
    return dcp;
  }
  dcp_t* dcp_cz(uint32_t i, uint32_t j,uint32_t n){
    dcp_t* dcp = init_dcp(n);
    mzd_write_bit(dcp->F,i,j,1);
    mzd_write_bit(dcp->F,j,i,1);
    return dcp;
  }
  dcp_t* dcp_s(uint32_t i, uint32_t n){
    dcp_t* dcp = init_dcp(n);
    mzd_write_bit(dcp->F,i,i,1);
    mzd_write_bit(dcp->h_x,i,0,1);
    return dcp;
  }
  dcp_t* dcp_transvection(pauli_t* z_pauli,uint32_t n){
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

  comb_dcp_t *comb_calloc_dcp(int count, int n){
    comb_dcp_t* pointer = malloc(count * sizeof(comb_dcp_t));
    for(int i =0; i <count; i++){
      pointer[i].F.m=n;
      pointer[i].F.n=n;
      pointer[i].F.data=comb_alloc_bit_array(n*n);
      pointer[i].h_x.n=n;
      pointer[i].h_x.data=comb_alloc_bit_array(n);
      pointer[i].h_z.n=n;
      pointer[i].h_z.data=comb_alloc_bit_array(n);
    }
    return pointer;
  }
  void comb_cfree_dcp(int count, comb_dcp_t* pointer){
    for(int i =0; i < count; i++){
      comb_free_bit_array(pointer[i].F.data);
      comb_free_bit_array(pointer[i].h_x.data);
      comb_free_bit_array(pointer[i].h_z.data);
    }
    free(pointer);
  }
  void comb_cfree_vector(int count, comb_vector_t* pointer){
    for(int i = 0; i < count; i++){
      comb_free_bit_array(pointer[i].data);
    }
    free(pointer);
  }
  comb_vector_t *comb_calloc_vector(int count, int n){
    comb_vector_t* pointer = (comb_vector_t*)malloc(count * sizeof(comb_vector_t));
    for(int i =0; i < count; i++){
      pointer[i].n=n;
      pointer[i].data=comb_alloc_bit_array(n);
    }
    return pointer;
  }

#ifdef __cplusplus
}
#endif
