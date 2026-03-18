#ifdef __cplusplus
extern "C" {
#endif
#include "comb.h"



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

  void print_pauli_z_coset(pauli_z_coset_t *coset, char print_newline_p){
    print_pauli(coset->factor,0);
    printf("<");
    for(int i= 0;i<coset->k;i++){
      print_pauli(coset->set_generators[i],0);
      printf(i+1==coset->k?"":", ");
    }
    printf(">");
    if(print_newline_p)printf("\n");
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


#ifdef __cplusplus
}
#endif
