#include "comb.h"


  dcp_t* init_dcp(uint32_t n){
    dcp_t* tableau = (dcp_t*)malloc(sizeof(dcp_t));
    // follows demoor and dehaene
    init_dcp_to_size(tableau,n);
    return tableau;
  }
  stabiliser_state_t* init_state(int n){
    stabiliser_state_t* state = (stabiliser_state_t*)malloc(sizeof(stabiliser_state_t));
    state->x_parts=mzd_init(n+1,n); // one more row for gaussian elimination on measurement operator
    state->z_parts=mzd_init(n+1,n);
    state->delta_vector=mzd_init(n+1,1);
    state->epsilon_vector=mzd_init(n+1,1);
    state->x_parts->nrows=n;
    state->z_parts->nrows=n;
    state->delta_vector->nrows=n;
    state->epsilon_vector->nrows=n;
    return state;
  }
  

  void init_dcp_to_size(dcp_t* tableau, uint32_t n){
    tableau->F=mzd_init(n,n);
    tableau->h_x=mzd_init(n,1);
    tableau->h_x_transpose=mzd_init(1,n);
    tableau->h_z=mzd_init(n,1);
    tableau->F_diagonal=mzd_init(n,1);
  }

  pauli_z_coset_t* init_pauli_z_coset(uint32_t n, uint32_t k){
    pauli_z_coset_t* coset = (pauli_z_coset_t*)malloc(sizeof(pauli_z_coset_t));
    coset->factor_transpose=mzd_init(1,n);
    coset->set_generators=mzd_init(n,n);
    return coset;
  }
  void free_dcp_array(uint32_t count, dcp_t* array){
    for(int i = 0; i < count; i++){
      mzd_free(array[i].F);
      mzd_free(array[i].h_x);
      mzd_free(array[i].h_x_transpose);
      mzd_free(array[i].h_z);
      mzd_free(array[i].F_diagonal);
    }
    free(array);
  }
  void free_dcp(dcp_t* dcp){
    mzd_free(dcp->F);
    mzd_free(dcp->h_x);
    mzd_free(dcp->h_x_transpose);
    mzd_free(dcp->h_z);
    mzd_free(dcp->F_diagonal);
    free(dcp);
  }
  void free_pauli_z_coset(pauli_z_coset_t* coset){
    mzd_free(coset->factor_transpose);
    coset->set_generators->nrows=coset->set_generators->ncols; // bring back up to size before freeing
    mzd_free(coset->set_generators);
    free(coset);
  }

  void set_dcp_to(dcp_t* a, dcp_t *b){
    // assumed to be of same size
    mzd_copy(a->F,b->F);
    mzd_copy(a->F_diagonal,b->F_diagonal);
    mzd_copy(a->h_x,b->h_x);
    mzd_copy(a->h_x_transpose,b->h_x_transpose);
    mzd_copy(a->h_z,b->h_z);
  }

  comb_dcp_t *comb_calloc_dcp(int count, int n){
    comb_dcp_t* pointer = (comb_dcp_t*)calloc(count , sizeof(comb_dcp_t));
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
    comb_vector_t* pointer = (comb_vector_t*)calloc(count , sizeof(comb_vector_t));
    for(int i =0; i < count; i++){
      pointer[i].n=n;
      pointer[i].data=comb_alloc_bit_array(n);
    }
    return pointer;
  }
