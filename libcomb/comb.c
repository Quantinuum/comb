#ifdef __cplusplus
extern "C" {
#endif
#include "comb.h"



  int postselect_diagonal_clifford(clifford_t* result, clifford_t* clifford, pauli_t *p){
    // assume that p is of x type.
    int n = CLIFFORD_N(clifford);
    pauli_t** old_stabilisers=malloc(n*sizeof(pauli_t*));
    for(int i = 0 ; i< n; i++){
      pauli_t* x_i = pauli_x(i,n);
      old_stabilisers[i] = clifford_conjugate_pauli(0,clifford,x_i);
      free_pauli(x_i);
      print_pauli(old_stabilisers[i],1);
    }
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
      printf("case 3\n");
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
      printf("new stabilisers:\n");
      for(int i =0;i<n;i++){
        printf("%d: ",i);
        print_pauli(new_stabilisers[i],1);
      }
      // set the images of the x_i in result to the new_stabilisers
      for(int i = 0; i< n; i++){
        mzd_write_bit(result->h_bar,i,0,new_stabilisers[i]->sign_bit);
        mzd_write_bit(result->h_bar,i+n,0,0);
        for(int j = 0; j < 2*n; j++){
          mzd_write_bit(result->c_bar,i+n,j,mzd_read_bit(new_stabilisers[i]->data,j,0));
          mzd_write_bit(result->c_bar,i+n,j,0);
        }
      }
      return 0;
    }
    else{
      //case 1
      // try to figure out if you can create -1 from it.
      // can try from the x patterns of the images of the x_i?
      // won't generalise to multiple postselections...
      pauli_t* r = p;
      printf("case 1\n");
      for(int i = 0 ; i < n; i ++){
        if(mzd_read_bit(p->data,i+n,0))
          r = pauli_mul(0,r,old_stabilisers[i]);
      }
      printf("result:\n");
      print_pauli(r,1);
      if(r->sign_bit){
        return 1; // -I so always detected; state vanishes
      }
      else{ // set result = input clifford
        result->c_bar=clifford->c_bar;
        result->h_bar=clifford->h_bar;
        return 0;
      }
    }


    return 0;
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

  int postselect_pauli_z_coset(pauli_z_coset_t* result, pauli_z_coset_t* coset, int k, pauli_t **measurement_operators){
    return 0;
  }
#ifdef __cplusplus
}
#endif
