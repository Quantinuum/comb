#ifdef __cplusplus
extern "C" {
#endif
#include "comb.h"



  void print_dcp(dcp_t* dcp, char print_newline_p){ // print_newline_p ignored actually
    int n = DCP_N(dcp);
    for(int i = 0 ;i<n; i++){
      pauli_t* p_x = pauli_x(i,n);
      pauli_t* p_z = pauli_z(i,n);

      print_pauli(p_x,0);
      printf(" -> ");

      pauli_t* p_z_parts = init_pauli(n);

      for(int j = 0;j<n;j++)mzd_write_bit(p_z_parts->data,j,0,mzd_read_bit(dcp->F,i,j));
      p_z_parts->sign_bit=mzd_read_bit(dcp->h_x,i,0);
      pauli_t* result = init_pauli(n);
      pauli_mul(result,p_x,p_z_parts);
      print_pauli(result,1);

      print_pauli(p_z,0);
      printf(" -> ");
      if(mzd_read_bit(dcp->h_z,i,0))printf("-");
      print_pauli(p_z,1);

      free_pauli(p_z_parts);
      free_pauli(p_x);
      free_pauli(p_z);
    }
    if(print_newline_p)printf("\n");
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
