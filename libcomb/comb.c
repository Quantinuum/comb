#ifdef __cplusplus
extern "C" {
#endif
#include "comb.h"


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
    printf("old stabilisers:\n");
    for(int i = 0; i < n; i++){
      print_pauli(old_stabilisers[i],1);
    }
    printf("new stabilisers:\n");
    for(int i = 0; i < n; i++){
      print_pauli(new_stabilisers[i],1);
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
    mzd_print(coefficient_matrix);
    mzd_print(coefficient_vector);
    if(0==mzd_solve_left(coefficient_matrix,coefficient_vector,0,1)){
      printf("solution found \n"); // unsure if working.
      mzd_print(coefficient_matrix);
      mzd_print(coefficient_vector);
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
