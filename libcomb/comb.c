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

  pauli_t* pauli_x(size_t i, size_t n){
    pauli_t* p = init_pauli(n);
    mzd_write_bit(p->data,i + n,0,1);
    return p;
  }
  pauli_t* pauli_z(size_t i,size_t n){
    pauli_t* p = init_pauli(n);
    mzd_write_bit(p->data,i,0,1);
    return p;
  }

  void print_pauli(pauli_t *p, char print_newline_p){
    int phase=0;
    int n = PAULI_N(p);
    for(int i= 0;i<n;i++){
      if(mzd_read_bit(p->data,i,0) && mzd_read_bit(p->data,i+n,0))phase++;
    }

    if(p->sign_bit)phase+=2;
    if(mzd_read_bit(p->data,2*n,0))phase++;
    phase &= 3;
    switch(phase){
    case 0: break;
    case 1: printf("i"); break;
    case 2: printf("-"); break;
    case 3: printf("-i"); break;
    default: printf("error in printing pauli");
    }
    for(int i = 0; i< n; i++){
      char x_p = mzd_read_bit(p->data,i + n,0);
      char z_p = mzd_read_bit(p->data,i ,0);
      printf(z_p?(x_p?"Y":"Z"):(x_p?"X":"I"));
    }
    if(print_newline_p)printf("\n");
  }

  void print_clifford(clifford_t* clifford, char print_newline_p){ // print_newline_p ignored actually
    int n = CLIFFORD_N(clifford);
    for(int i = 0 ;i<n; i++){
      pauli_t* p_x = pauli_x(i,n);
      pauli_t* p_z = pauli_z(i,n);
      pauli_t* result;

      print_pauli(p_x,0);
      printf(" -> ");
      result = clifford_conjugate_pauli(0,clifford,p_x);
      print_pauli(result,1);
      free_pauli(result);

      print_pauli(p_z,0);
      printf(" -> ");
      result = clifford_conjugate_pauli(0,clifford,p_z);
      print_pauli(result,1);
      free_pauli(result);

      free_pauli(p_x);
      free_pauli(p_z);
    }
  }

  pauli_t *clifford_conjugate_pauli(pauli_t* result, clifford_t* clifford, pauli_t *pauli){
    // TODO check, improve efficiency
    int n = PAULI_N(pauli);
    if(!result)result = init_pauli(n);
    // matrix multiplication here
    mzd_mul_naive(result->data,clifford->c_bar,pauli->data);
    mzd_t* h_bar_transpose=mzd_transpose(0,clifford->h_bar);
    mzd_t* h_bar_b_1_bar = mzd_mul_naive(0,h_bar_transpose, pauli->data);
    mzd_t* u = clifford_u_matrix(n);
    mzd_t* uc = mzd_mul_naive(0,u,clifford->c_bar);
    mzd_t* c_transpose = mzd_transpose(0,clifford->c_bar);
    mzd_t* cuc = mzd_mul_naive(0,c_transpose,uc);
    set_to_lows(cuc);
    mzd_t* lows_b1 = mzd_mul_naive(0,cuc,pauli->data);
    mzd_t* p_t = mzd_transpose(0,pauli->data);
    mzd_t* b_lows_b = mzd_mul_naive(0,p_t,lows_b1);
    int sign_bit = pauli->sign_bit + mzd_read_bit(h_bar_b_1_bar,0,0) +  mzd_read_bit(b_lows_b,0,0);
    result->sign_bit = (sign_bit & 1);
    mzd_free(h_bar_transpose);
    mzd_free(h_bar_b_1_bar);
    mzd_free(u);
    mzd_free(uc);
    mzd_free(c_transpose);
    mzd_free(cuc);
    mzd_free(lows_b1);
    mzd_free(p_t);
    mzd_free(b_lows_b);
    return result;
  }
  clifford_t *clifford_mul(clifford_t *result, clifford_t *c2, clifford_t *c1){
    int n = CLIFFORD_N(c1);
    if(!result)result = init_clifford(n);
    mzd_mul_naive(result->c_bar,c2->c_bar,c1->c_bar);
    mzd_t* t1 = c1->h_bar;
    mzd_t* c1_bar_transpose=mzd_transpose(0,c1->c_bar);
    mzd_t* t2= mzd_mul_naive(0,c1_bar_transpose,c2->h_bar);
    mzd_t* c2_bar_transpose=mzd_transpose(0,c2->c_bar);
    mzd_t* u = clifford_u_matrix(n);
    mzd_t* uc2 = mzd_mul_naive(0,u,c2_bar_transpose);
    mzd_t* lows = mzd_mul_naive(0,c2_bar_transpose,uc2);
    set_to_lows(lows);
    mzd_t* lows_c1 = mzd_mul_naive(0,lows,c1->c_bar);
    mzd_t* product = mzd_mul_naive(0,c1_bar_transpose,lows_c1);
    mzd_t* t3=diag(product);
    mzd_t* t1t2 = mzd_add(0,t1,t2);
    mzd_add(result->h_bar,t1t2,t3);
    mzd_free(t1);
    mzd_free(c1_bar_transpose);
    mzd_free(t2);
    mzd_free(c2_bar_transpose);
    mzd_free(u);
    mzd_free(uc2);
    mzd_free(lows);
    mzd_free(lows_c1);
    mzd_free(product);
    mzd_free(t3);
    mzd_free(t1t2);
    return result;
  }
  mzd_t* diag(mzd_t* matrix){
    // returns the diagonal part of a matrix
    int n = matrix->ncols;
    mzd_t* d = mzd_init(n,1);
    for(int i = 0 ; i < n ; i++){
      mzd_write_bit(d,i,0,mzd_read_bit(matrix,i,i));
    }
    return d;
  }

  mzd_t* clifford_u_matrix(int n){
    mzd_t* matrix = mzd_init(1+(2*n),1+(2*n));
    for(int i = 0 ; i < n ; i++){
      mzd_write_bit(matrix,i,i+n,1);
    }
    mzd_write_bit(matrix,2*n,2*n,1);
    return matrix;
  }

  void* set_to_lows(mzd_t* matrix){ // erases the upper triangular part of matrix
    for(int i = 0; i< matrix->nrows; i++){
      for(int j = i; j < matrix->ncols; j++){
        mzd_write_bit(matrix,i,j,0);
      }
    }
  }

  clifford_t *cz(size_t i, size_t j, size_t n){
    clifford_t* c = init_clifford(n);
    mzd_write_bit(c->c_bar,i,j+n,1);
    mzd_write_bit(c->c_bar,j,i+n,1);
    return c;
  }

  clifford_t *cnot(size_t i, size_t j, size_t n){
    clifford_t* c = init_clifford(n);
    mzd_write_bit(c->c_bar,i,j,1);
    mzd_write_bit(c->c_bar,j+n,i+n,1);
    return c;
  }

  pauli_t *pauli_mul(pauli_t *result, pauli_t *p1, pauli_t *p2){
    int n = PAULI_N(p1);
    if(!result)result = init_pauli(n);
    mzd_add(result->data, p1->data, p2->data);
    int d1 = mzd_read_bit(p1->data,2*n,0);
    int d2 = mzd_read_bit(p2->data,2*n,0);
    mzd_t* u = clifford_u_matrix(n);
    mzd_t* p2t = mzd_transpose(0,p2->data);
    mzd_t* au = mzd_mul_naive(0,p2t,u);
    mzd_t* aua = mzd_mul_naive(0,au,p1->data);
    result->sign_bit = (p1->sign_bit) ^ (p2->sign_bit) ^ (d1 & d2) ^ mzd_read_bit(aua,0,0);
    mzd_free(u);
    mzd_free(p2t);
    mzd_free(au);
    mzd_free(aua);
    return result;
  }

#ifdef __cplusplus
}
#endif
