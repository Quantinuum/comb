#ifdef __cplusplus
extern "C" {
#endif
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "comb.h"



  void convert_permutation(mzd_t* permutation, uint32_t n_z, uint32_t n, uint32_t *rrefs)
  {
    for(int j=0, i = (n - n_z); i < n; i ++,j++){
      int stabiliser_index = j;
      int column_to_eliminate = rrefs[j];
      //mzd_write_bit(permutation,column_to_eliminate,i,1);
      mzd_write_bit(permutation,i,column_to_eliminate,1);
    }
    // now do the other columns
    int current_index=0;
    for(int i = 0; i < n; i++){
      int rref_p = 0;
      for(int j = 0; j < n_z; j++){
        if(rrefs[j]==i)rref_p=1;
      }
      if(!rref_p){ // wasnt an rref; need to add to permutation
        //mzd_write_bit(permutation,i,current_index,1);
        mzd_write_bit(permutation,current_index,i,1);
        current_index++;
      }
    }
  }




  comb_bit_array_t comb_alloc_bit_array(int n){
    //lwlog_debug("allocating bit array of size %d",n);
    return (comb_bit_array_t)malloc(n);
  }
  void comb_free_bit_array(comb_bit_array_t p){
    free(p);
  }


  void comb_set_cvector_aref(comb_vector_t* vector, int n, int i, int value){
    if(value)SetBit(vector[n].data,i);
    else ClearBit(vector[n].data,i);
  }
  void comb_set_cdcp_F_aref(comb_dcp_t* dcp, int n, int i, int j, int value){
    if(value)SetBit(dcp[n].F.data,(i*dcp[n].F.n) + j);
    else ClearBit(dcp[n].F.data,(i*dcp[n].F.n) + j);
  }
  void comb_set_cdcp_h_x_aref(comb_dcp_t* dcp, int n, int i, int value){
    if(value)SetBit(dcp[n].h_x.data,i);
    else ClearBit(dcp[n].h_x.data,i);
  }
  void comb_set_cdcp_h_z_aref(comb_dcp_t* dcp, int n, int i, int value){
    if(value)SetBit(dcp[n].h_z.data,i);
    else ClearBit(dcp[n].h_z.data,i);
  }

  void comb_set_permutation_aref(comb_permutation_t* permutation, int i, int value){
    permutation->indices[i]=value;
  }
  void comb_set_bit_array_aref(comb_bit_array_t array, int i, int value){
    if(value)SetBit(array,i);
    else ClearBit(array,i);
  }

  void convert_dcp(dcp_t* dcp, comb_dcp_t* value)
  {
    int n = value->h_x.n;
    dcp->F=mzd_init(n,n);
    dcp->h_x=mzd_init(n,1);
    dcp->h_z=mzd_init(n,1);
    for(int i = 0; i < n; i++){
      mzd_write_bit(dcp->h_x,i,0,TestBit(value->h_x.data,i)?1:0);
      mzd_write_bit(dcp->h_z,i,0,TestBit(value->h_z.data,i)?1:0);
      for(int j = 0; j < n; j++){
        mzd_write_bit(dcp->F,i,j,TestBit(value->F.data,(i*n) + j)?1:0);
      }
    }
  }
  void convert_pauli_x(pauli_t* pauli, comb_pauli_x_t* value)
  {
    int n = value->n;
    pauli->data=mzd_init(1+(2*n),1);
    pauli->sign_bit=0;
    for(int i = 0; i < n;i++){
      //printf("bit %d=%d\n",i,TestBit(value->data,i));
      if(TestBit(value->data,i)){
        //printf("writing bit %d\n",i);
      mzd_write_bit(pauli->data,i+n,0,1);
      }
    }
  }
  void convert_pauli_z(pauli_t* pauli, comb_pauli_z_t* value){
    int n = value->n;
    pauli->data=mzd_init(1+(2*n),1);
    pauli->sign_bit=0;
    for(int i = 0; i < n;i++){
      //printf("bit %d=%d\n",i,TestBit(value->data,i));
      if(TestBit(value->data,i)){
        //printf("writing bit %d\n",i);
        mzd_write_bit(pauli->data,i,0,1);
      }
    }
  }

#ifdef __cplusplus
}
#endif
