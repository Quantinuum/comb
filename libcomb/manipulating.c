#ifdef __cplusplus
extern "C" {
#endif
#include "comb.h"
  mzd_t* clifford_u_matrix(int n){
    mzd_t* matrix = mzd_init((2*n),(2*n));
    for(int i = 0 ; i < n ; i++){
      mzd_write_bit(matrix,i,i+n,1);
    }
    // mzd_write_bit(matrix,2*n,2*n,1);
    return matrix;
  }

  pauli_t *pauli_mul(pauli_t *result, pauli_t *p1, pauli_t *p2){ // lots of allocation... try to avoid this
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
    mzd_free(p2t);
    mzd_free(au);
    mzd_free(aua);
    mzd_free(u);
    return result;
  }

  inline void add_and(mzd_t* result, mzd_t* F1,mzd_t* F2){
    // Adds F1 & F2 to result. Assumes only one row. must be initialised to right shape
    for(int i = 0; i < result->width; i++){
      result->data[0]^=(F1->data[0] & F2->data[0]);
    }
    result->data[result->width - 1] &= result->high_bitmask; // on the safe side
  }


  void *incf_coset(pauli_z_coset_t *coset, int index){
    // Iterates over the coset by adding one of the generators, chosen by a gray code, to factor.
    // !!! won't work with less than 64 bits (reasonable assumption if we're using lookup tables)
    coset->factor->data[0]^=mzd_row(coset->set_generators,gray_bit_to_switch(index))[0];
  }


  int pauli_commute_p(pauli_t *p1, pauli_t *p2){ // slow but probably eventually shouldn't need this
    int commute_p=1;
    int n = PAULI_N(p1);
    for(int i =0 ; i < n; i++){
      if(mzd_read_bit(p1->data,i,0) && mzd_read_bit(p2->data,i+n,0))commute_p=1-commute_p;
      if(mzd_read_bit(p1->data,i+n,0) && mzd_read_bit(p2->data,i,0))commute_p=1-commute_p;
    }
    return commute_p;
  }

  int postselect_on_measuring_plus_one_x(stabiliser_state_t* state, mzd_t* x_measurement_matrix /* n_x x n*/
                                         ){
    int n = STABILISER_N(state);
    int n_x = x_measurement_matrix->nrows;

    mzd* temp = mzd_init(1,n);
    // z parts submatrix:
    mzd_t* z_parts = mzd_init(n,n);
    mzd_submatrix(z_parts,state->tableau,0,0,n,n);
    mzd_t* temp_commutations=mzd_init(n_x,n);
  }

  int postselect_on_measuring_plus_one(stabiliser_state_t* state, pauli_t* p){
    // warning: destructive
    int n = STABILISER_N(state);
    pauli_t** old_stabilisers=state->generators;
    // case 1: p is in the stabiliser group. in this case, measuring is always +1 (in which case no effect) or -1 in which case always detected, and we return non-zero.
    // case 2: p in the normaliser (commutes with all stabilisers) yet does not lie in the stabiliser group. in this case, p measures some of the data. this case is not possible for stabiliser states, as the normaliser equals the stabiliser group, so we neglect it.
    // case 3: p not in the normaliser. in this case, we take one elt S' which does not commute with P, add P as a stabiliser, and for all S not equal to S', we add S if it commutes with P and PS'S otherwise.
    // seeing as how for us diagonal cliffords enter into the picture solely through the images of the x_i bases might not be an issue. probably want a stabiliser state after all though
    int anticommutes_with_some_pauli_p=0;
    int anti_commuting_index=-1;
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
        lwlog_debug("processing stabiliser generator %s",print_pauli(old_stabilisers[i],0));
        if(i== anti_commuting_index){
          lwlog_debug("it is the anticommuting operator, so replacing with measurement operator %s",print_pauli(p,0));
          new_stabilisers[i]=init_pauli(n);
          set_to_pauli(new_stabilisers[i],p);
        }
        else{
          if(pauli_commute_p(p,old_stabilisers[i])){
            new_stabilisers[i]=init_pauli(n);
            set_to_pauli(new_stabilisers[i],old_stabilisers[i]);
            lwlog_debug("it commutes, so keeping %s",print_pauli(new_stabilisers[i],0));
          }
          else {
          new_stabilisers[i]=pauli_mul(0,old_stabilisers[i],anticommuter);
          lwlog_debug("it anticommutes, so multiplying it by the anticommuter: %s",print_pauli(new_stabilisers[i],0));
          }
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
      mzd_t* A = mzd_init(2*n,n);
      mzd_t* B = mzd_init(2*n,1);
      for(int i = 0; i < (2*n); i++){
        for(int j = 0; j < n; j++){
          mzd_write_bit(A,i,j,mzd_read_bit(old_stabilisers[j]->data,i,0));
        }
        mzd_write_bit(B,i,0,mzd_read_bit(p->data,i,0));
      }
      int solution = mzd_solve_left(A,B,0,0);
      if(solution==-1)lwlog_err("warning: solution was -1! no solution found");
      pauli_t* r = p;
      //printf("case 1\n");
      for(int i = 0 ; i < n; i ++){
        if(mzd_read_bit(B,i,0))
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
