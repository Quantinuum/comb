#ifdef __cplusplus
extern "C" {
#endif
#include "comb.h"


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


  clifford_t *clifford_mul(clifford_t *result, clifford_t *c2, clifford_t *c1){
    int n = CLIFFORD_N(c1);
    if(!result)result = init_clifford(n);
    mzd_mul_naive(result->c_bar,c2->c_bar,c1->c_bar);
    mzd_t* t1 = c1->h_bar;
    mzd_t* c1_bar_transpose=mzd_transpose(0,c1->c_bar);
    mzd_t* t2= mzd_mul_naive(0,c1_bar_transpose,c2->h_bar);
    mzd_t* c2_bar_transpose=mzd_transpose(0,c2->c_bar);
    mzd_t* u = clifford_u_matrix(n);
    mzd_t* uc2 = mzd_mul_naive(0,u,c2->c_bar);
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

  void* set_to_lows(mzd_t* matrix){ // erases the upper triangular part of matrix
    for(int i = 0; i< matrix->nrows; i++){
      for(int j = i; j < matrix->ncols; j++){
        mzd_write_bit(matrix,i,j,0);
      }
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


  clifford_t *pauli_to_clifford(pauli_t *pauli){
    int n = PAULI_N(pauli);
    clifford_t* clifford = init_clifford(n);
    for(int i = 0; i < n ; i++){
      mzd_write_bit(clifford->h_bar,i+n,0,mzd_read_bit(pauli->data,i,0));
      mzd_write_bit(clifford->h_bar,i,0,mzd_read_bit(pauli->data,i+n,0));
    }
    return clifford;
  }
  clifford_t *pauli_z_transvection(pauli_t *P, int c){
    // does it work? looks like it, they square properly
    // returns e^(i pi/4 c P). c must be a z pauli otherwise answer not defined
    while(c<0)c+=4;
    while(c>=4)c-=4;
    int n = PAULI_N(P);
    switch(c){
    case 0:
      return init_clifford(n);
    case 1:
      clifford_t* c = init_clifford(n);
      for(int i = 0; i < n; i++){
        mzd_write_bit(c->c_bar, 2 * n, i + n, mzd_read_bit(P->data,i,0));
        mzd_write_bit(c->h_bar, i + n, 0, mzd_read_bit(P->data,i,0));
        for(int j = 0; j < n ; j++){
          mzd_write_bit(c->c_bar,i ,j +n ,mzd_read_bit(P->data,i,0) & mzd_read_bit(P->data,j,0));
        }
      }
      return c;
    case 2:
      return pauli_to_clifford(P);
    case 3:
      clifford_t* c1 = pauli_z_transvection(P,1);
      clifford_t* c2 = pauli_z_transvection(P,2);
      clifford_t* cr = clifford_mul(0,c1,c2);
      free_clifford(c1);
      free_clifford(c2);
      return cr;
    }
  }



  pauli_z_coset_t *diagonal_clifford_to_coset(clifford_t *clifford){
    int n = CLIFFORD_N(clifford);
    pauli_t** old_stabilisers=malloc(n*sizeof(pauli_t*));
    for(int i = 0 ; i< n; i++){
      pauli_t* x_i = pauli_x(i,n);
      old_stabilisers[i] = clifford_conjugate_pauli(0,clifford,x_i);
      free_pauli(x_i);
    }
    mzd_t* z_parts = mzd_init(n,n);
    for(int i = 0; i< n; i++){
      for(int j = 0; j < n; j++){
        mzd_write_bit(z_parts,i,j,mzd_read_bit(old_stabilisers[i]->data,j,0));
      }
    }
    mzp_t* p = mzp_init(n);
    mzp_t* q = mzp_init(n);
    int r = mzd_pluq(z_parts,p,q,0);
    mzd_t* u  = mzd_init(n,n);
    mzd_t* l  = mzd_init(n,n);
    // from m4ri examples
    for (rci_t i = 0; i < r; ++i){
      for (rci_t j = 0; j < i; ++j)
        mzd_write_bit(l, i, j, mzd_read_bit(z_parts,i,j));
      for (rci_t j = i + 1; j < n; ++j)
        mzd_write_bit(u, i, j, mzd_read_bit(z_parts,i,j));
    }

    for (rci_t i = r; i < n; ++i)
      for (rci_t j = 0; j < r; ++j)
        mzd_write_bit(l, i, j, mzd_read_bit(z_parts,i,j));
    for (rci_t i = 0; i < r; ++i){
      mzd_write_bit(l,i,i, 1);
      mzd_write_bit(u,i,i, 1);
    }
    mzd_t* rref = u;
    mzd_apply_p_right(rref,q); // warning: cannot later on refer to u
    u=0;

    mzd_t* row_ops_t = mzd_init(n,n);
    // set row_ops_t = L^T to invert
    for(int i = 0 ; i< n; i++){
      
        for(int j = 0 ; j <n; j++){
          mzd_write_bit(row_ops_t,i,j,mzd_read_bit(l,j,i));
        }
        mzd_write_bit(row_ops_t,i,i,1);
      
    }
    mzd_trtri_upper(row_ops_t); // row_ops_t = L^T^-1
    mzd_t* row_ops = mzd_transpose(0,row_ops_t); // row_ops = L^-1 (extended)

    mzd_apply_p_right_trans(row_ops,p); // row_ops = L^-1 P^T
    //mzd_print(row_ops_t);
    //mzd_print(row_ops);

    // Apply the row_ops to get the new stabilisers

    pauli_t** new_stabilisers=malloc(n*sizeof(pauli_t*));
    for(int i = 0 ; i< n; i++){
      new_stabilisers[i]=init_pauli(n);
      for(int j = 0 ; j < n; j++){
        if(mzd_read_bit(row_ops,i,j)){
          new_stabilisers[i]=pauli_mul(0,new_stabilisers[i],old_stabilisers[j]);
        }
      }
    }
    // set up system solving
    mzd_t* coefficient_matrix = mzd_init(n,n);
    mzd_t* coefficient_vector = mzd_init(n,1);
    for(int i = r; i< n; i++){
      mzd_write_bit(coefficient_vector,i-r,0,new_stabilisers[i]->sign_bit);
      for(int j = 0; j< n; j++){
        mzd_write_bit(coefficient_matrix,i-r,j,mzd_read_bit(new_stabilisers[i]->data,j+n,0));
      }
    }
    if(mzd_solve_left(coefficient_matrix,coefficient_vector,0,1)){
      printf("solution not found! \n");
      return 0;
    }


    pauli_z_coset_t* coset = init_pauli_z_coset(n,r);

    for(int i = 0; i < n; i++){
      new_stabilisers[i]=init_pauli(n);
      for(int j = 0; j < n; j++){
        mzd_write_bit(new_stabilisers[i]->data,j,0,mzd_read_bit(rref,i,j));
      }
      coset->set_generators[i]=new_stabilisers[i];
    }
    for(int i = 0; i < n; i++){
      mzd_write_bit(coset->factor->data,i,0,mzd_read_bit(coefficient_vector,i,0));
      }
    return coset;
  }

  clifford_t *czs(size_t n, int k,int* indices){
    //usage e.g. czs(5,4,{{0, 1},{0,2}})
    clifford_t* return_clifford=init_clifford(n);
    for(int i = 0 ; i < k ; i++){
      clifford_t* cze = cz(indices[2*i],indices[(2*i) + 1],n);
      return_clifford=clifford_mul(0,cze,return_clifford);

    }
    return return_clifford;
  }


#ifdef __cplusplus
}
#endif
