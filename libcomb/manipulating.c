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
      mzd_write_bit(result,i,0,mzd_read_bit(result,i,0) ^ (mzd_read_bit(F1,i,i) & mzd_read_bit(F2,i,i)));
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


  pauli_t* dcp_image_x(dcp_t *dcp, uint32_t i){
    int n = DCP_N(dcp);
    pauli_t* p_z_parts = init_pauli(n);
    mzd_write_bit(p_z_parts->data,i+n,0,1);
    mzd_write_bit(p_z_parts->data,2*n,0,mzd_read_bit(dcp->F,i,i));
    for(int j = 0;j<n;j++)mzd_write_bit(p_z_parts->data,j,0,mzd_read_bit(dcp->F,i,j));

    p_z_parts->sign_bit=mzd_read_bit(dcp->h_z,i,0); //? think it works
    return p_z_parts;
}
  stabiliser_state_t* dcp_to_stabiliser_state(dcp_t* dcp){
    lwlog_debug("in dcp_to_stabiliser_state");
    int n = DCP_N(dcp);
    stabiliser_state_t* state = (stabiliser_state_t*)malloc(sizeof(stabiliser_state_t));
    state->n=n;
    state->generators=(pauli_t**)malloc(n*sizeof(pauli_t*));
    for(int i = 0 ; i < n; i++){
      state->generators[i]=dcp_image_x(dcp,i);
    }
    lwlog_debug("returning from dcp_to_stabiliser_state");
    return state;
  }

  pauli_t* pauli_x(uint32_t i, uint32_t n){
    pauli_t* p = init_pauli(n);
    mzd_write_bit(p->data,i + n,0,1);
    return p;
  }
  pauli_t* pauli_z(uint32_t i,uint32_t n){
    pauli_t* p = init_pauli(n);
    mzd_write_bit(p->data,i,0,1);
    return p;
  }

  pauli_z_coset_t *diagonal_clifford_to_coset(stabiliser_state_t *state){ // WARNING: dcp will be destroyed after this!
    //print_stabiliser_state(state,1);
    // Gaussian elimination on the Z parts is simply Gaussian elimination on F.
    int n = state->n;
    pauli_t** old_stabilisers= state->generators;
    mzp_t* p = mzp_init(n);
    mzp_t* q = mzp_init(n);
    mzd_t* z_parts = mzd_init(n,n);
    for(int i = 0; i< n; i++){
      for(int j = 0; j < n; j++){
        mzd_write_bit(z_parts,i,j,mzd_read_bit(old_stabilisers[i]->data,j,0));
      }
    }
    int r = mzd_pluq(z_parts,p,q,0);
    lwlog_debug("PLUQ done");
    mzd_t* u  = mzd_init(n,n);
    mzd_t* l  = mzd_init(n,n);

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

    for(int i = 0; i < r; i++){
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
  int pauli_commute_p(pauli_t *p1, pauli_t *p2){
    int commute_p=1;
    int n = PAULI_N(p1);
    for(int i =0 ; i < n; i++){
      if(mzd_read_bit(p1->data,i,0) && mzd_read_bit(p2->data,i+n,0))commute_p=1-commute_p;
      if(mzd_read_bit(p1->data,i+n,0) && mzd_read_bit(p2->data,i,0))commute_p=1-commute_p;
    }
    //print_pauli(p1,0);
    //printf(commute_p?" commutes with ": " anticommutes with ");
    //print_pauli(p2,1);
    return commute_p;
  }

  int postselect_on_measuring_plus_one(stabiliser_state_t* state, pauli_t* p){
    // warning: destructive
       // assume that p is of x type.
    lwlog_debug("measuring state %s ",print_stabiliser_state(state,0));
    //print_pauli(p,1);
    fflush(stdout);
    int n = state->n;
    pauli_t** old_stabilisers=state->generators;
    // case 1: p is in the stabiliser group. in this case, measuring is always +1 (in which case no effect) or -1 in which case always detected, and we return non-zero.
    // case 2: p in the normaliser (commutes with all stabilisers) yet does not lie in the stabiliser group. in this case, p measures some of the data. this case is not possible for stabiliser states, as the normaliser equals the stabiliser group, so we neglect it.
    // case 3: p not in the normaliser. in this case, we take one elt S' which does not commute with P, add P as a stabiliser, and for all S not equal to S', we add S if it commutes with P and PS'S otherwise.
    // seeing as how for us diagonal cliffords enter into the picture solely through the images of the x_i bases might not be an issue. probably want a stabiliser state after all though
    int anticommutes_with_some_pauli_p=0;
    int anti_commuting_index=0;
    for(int i = 0; i < n; i++){
      if(!pauli_commute_p(p,old_stabilisers[i])){
        anti_commuting_index=i;
        anticommutes_with_some_pauli_p=1;
      }
    }
    if(anticommutes_with_some_pauli_p){
      // case 3
      lwlog_debug("stabiliser %s anticommutes",print_pauli(old_stabilisers[anti_commuting_index],0));
      pauli_t* anticommuter = old_stabilisers[anti_commuting_index];
      pauli_t** new_stabilisers=malloc(n*sizeof(pauli_t*));
      for(int i =0 ; i< n; i++){
        if(i== anti_commuting_index){
          new_stabilisers[i]=p;
        }
        else{
          if(pauli_commute_p(p,old_stabilisers[i])){
            new_stabilisers[i]=old_stabilisers[i];
          }
          else new_stabilisers[i]=pauli_mul(0,old_stabilisers[i],anticommuter);
        }
      }
      //  printf("new stabilisers:\n");
      //for(int i =0;i<n;i++){
      //  printf("%d: ",i);
      // print_pauli(new_stabilisers[i],1);
      // }
      // set the old stabilisers equal to the new stabilisers.
      state->generators=new_stabilisers;
      //for(int i = 0; i< n; i++)free_pauli(old_stabilisers[i]);
      free(old_stabilisers);
      return 0;
    }
    else{
      lwlog_debug("every stabiliser commutes, trying to form the identity");
      //case 1
      // try to figure out if you can create -1 from it.
      // can try from the x patterns of the images of the x_i?
      // won't generalise to multiple postselections... TODO this is broken so fix it
      pauli_t* r = p;
      //printf("case 1\n");
      for(int i = 0 ; i < n; i ++){
        if(mzd_read_bit(p->data,i+n,0))
          r = pauli_mul(0,r,old_stabilisers[i]);
      }
      //printf("result:\n");
      lwlog_debug("pauli %s formed",print_pauli(r,0));
      //print_pauli(r,1);
      if(r->sign_bit){
        lwlog_debug("-I so always detected; state vanishes");
        return 1; // -I so always detected; state vanishes
      }
      else{ // state is unchanged
        lwlog_debug("I so state unchanged");
        return 0;
      }
    }
    return 0;
  }

#ifdef __cplusplus
}
#endif
