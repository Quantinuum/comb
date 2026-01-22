#ifdef __cplusplus
extern "C" {
#endif
#include "comb.h"


  comb_dcp_tableau_t* tableau_init(size_t n){
    comb_dcp_tableau_t* tableau = (comb_dcp_tableau_t*)
      malloc(sizeof(comb_dcp_tableau_t));
    // follows demoor and dehaene
    tableau->c_bar=mzd_init(1+(2*n),1+(2*n));
    tableau->h_bar=mzd_init(1+(2*n),1);
    return tableau;
  }
#define CLIFFORD_N(t) ((((t)->c_bar->nrows)-1)/2)


  void* tableau_print(comb_dcp_tableau_t* tableau){
    //printf("clifford tableau of size %d:\n",CLIFFORD_N(tableau));
    //mzd_print(tableau->c_bar);
    //printf("\n");
    //mzd_print(tableau->h_bar);
    //printf("\n");
    int n = CLIFFORD_N(tableau);
    for(int i = 0; i < n; i++){
      comb_pauli_x_t* x = pauli_x_init(n);
      comb_pauli_z_t* z = pauli_z_init(n);
      SetBit(x->data,i);
      SetBit(z->data,i);
      print_pauli_x(x,1);
      print_pauli_z(z,0);
      printf(" -> ");
      print_pauli_z(z,1);
    }
  }
  void tableau_free(comb_dcp_tableau_t* tableau){
    mzd_free(tableau->c_bar);
    mzd_free(tableau->h_bar);
    free(tableau);
  }

  comb_pauli_z_t* pauli_z_init(size_t n){
    comb_pauli_z_t* p = (comb_pauli_z_t*)malloc(sizeof(comb_pauli_z_t));
    p->pauli_length = n;
    p->data=init_bit_array(n);
    return p;
  }
  comb_pauli_x_t* pauli_x_init(size_t n){
    comb_pauli_x_t* p = (comb_pauli_x_t*)malloc(sizeof(comb_pauli_x_t));
    p->pauli_length = n;
    p->data=init_bit_array(n);
    return p;
    }

  comb_bit_array_t init_bit_array(size_t n){
    comb_bit_array_t i = (comb_bit_array_t)malloc(1+(n/sizeof(comb_bit_word_t)));
    for(int j = 0; j< (n/sizeof(comb_bit_word_t)); j++)
      i[j]=0;
    return i;
  }

  void print_pauli_x(comb_pauli_x_t* p, char newline){
    for(int i = 0;i<p->pauli_length;i++){
      printf(TestBit(p->data,i)?"X":"I");
    }
    if(newline)printf("\n");
  }
  void print_pauli_z(comb_pauli_z_t* p, char newline){
    for(int i = 0;i<p->pauli_length;i++){
      printf(TestBit(p->data,i)?"Z":"I");
    }
    if(newline)printf("\n");
  }
#ifdef __cplusplus
}
#endif
