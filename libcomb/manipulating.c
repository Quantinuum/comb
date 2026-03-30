#ifdef __cplusplus
extern "C" {
#endif
#include "comb.h"
  mzd_t* clifford_u_matrix(int n){
    mzd_t* matrix = mzd_init(1+(2*n),1+(2*n));
    for(int i = 0 ; i < n ; i++){
      mzd_write_bit(matrix,i,i+n,1);
    }
    mzd_write_bit(matrix,2*n,2*n,1);
    return matrix;
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

  void add_and(mzd_t* result, mzd_t* F1,mzd_t* F2){
    for(int i = 0; i < F1->nrows; i++){
      mzd_write_bit(result,i,i,mzd_read_bit(result,i,i) ^ (mzd_read_bit(F1,i,i) & mzd_read_bit(F2,i,i)));
    }
  }

  dcp_t *dcp_mul(dcp_t *result, dcp_t *c2, dcp_t *c1){
    int n = DCP_N(c1);
    if(!result)result = init_dcp(n);
    mzd_add(result->F,c2->F,c1->F); // F parts add
    mzd_add(result->h_x,c2->h_x,c1->h_x);// h_x
    // h_z:
    mzd_add(result->h_z,c2->h_z,c1->h_z);
    mzd_addmul_naive(result->h_z,c1->F,c2->h_x);
    add_and(result->h_z,c1->F,c2->F);
    return result;
  }

  void* set_to_lows(mzd_t* matrix){ // erases the upper triangular part of matrix
    for(int i = 0; i< matrix->nrows; i++){
      for(int j = i; j < matrix->ncols; j++){
        mzd_write_bit(matrix,i,j,0);
      }
    }
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




#ifdef __cplusplus
}
#endif
