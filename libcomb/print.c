#ifdef __cplusplus
extern "C" {
#endif
#include "comb.h"

  #define PRINT_SIZE 1024
  char print_buffer[PRINT_SIZE];

  char* print_dcp(dcp_t* dcp, char print_newline_p){ // print_newline_p ignored actually
    int n = DCP_N(dcp);
    int counter = 0;
    for(int i = 0 ;i<n; i++){
      pauli_t* p_x = pauli_x(i,n);
      pauli_t* p_z = pauli_z(i,n);

      print_pauli(p_x,0);
      counter += snprintf(print_buffer+counter, PRINT_SIZE-counter," -> ");
      pauli_t* p_z_parts=dcp_image_x(dcp,i);
      print_pauli(p_z_parts,1);

      print_pauli(p_z,0);
      counter += snprintf(print_buffer+counter, PRINT_SIZE-counter," -> ");
      if(mzd_read_bit(dcp->h_z,i,0))counter += snprintf(print_buffer+counter, PRINT_SIZE-counter,"-");
      print_pauli(p_z,1);

      free_pauli(p_z_parts);
      free_pauli(p_x);
      free_pauli(p_z);
    }
    if(print_newline_p)counter += snprintf(print_buffer+counter, PRINT_SIZE-counter,"\n");
    return print_buffer;
  }

  char* print_stabiliser_state(stabiliser_state_t *coset, char print_newline_p){
    int counter = 0;
    counter += snprintf(print_buffer+counter, PRINT_SIZE-counter,"{");
    for(int i = 0 ; i< coset->n; i++){
      print_pauli(coset->generators[i],0);
        if(i +1 != coset->n)counter += snprintf(print_buffer+counter, PRINT_SIZE-counter,", ");
    }
      counter += snprintf(print_buffer+counter, PRINT_SIZE-counter,print_newline_p?"}\n":"}");
      return print_buffer;
  }

  char* print_pauli_z_coset(pauli_z_coset_t *coset, char print_newline_p){
    int counter = 0;
    print_pauli(coset->factor,0);
    counter += snprintf(print_buffer+counter, PRINT_SIZE-counter,"<");
    for(int i= 0;i<coset->k;i++){
      print_pauli(coset->set_generators[i],0);
      counter += snprintf(print_buffer+counter, PRINT_SIZE-counter,i+1==coset->k?"":", ");
    }
    counter += snprintf(print_buffer+counter, PRINT_SIZE-counter,">");
    if(print_newline_p)counter += snprintf(print_buffer+counter, PRINT_SIZE-counter,"\n");
    return print_buffer;
  }

  char* print_pauli(pauli_t *p, char print_newline_p){
    int phase=0;
    int counter = 0;
    int n = PAULI_N(p);
    for(int i= 0;i<n;i++){
      if(mzd_read_bit(p->data,i,0) && mzd_read_bit(p->data,i+n,0))phase++;
    }

    if(p->sign_bit)phase+=2;
    if(mzd_read_bit(p->data,2*n,0))phase++;
    phase &= 3;
    switch(phase){
    case 0: break;
    case 1: counter += snprintf(print_buffer+counter, PRINT_SIZE-counter,"i"); break;
    case 2: counter += snprintf(print_buffer+counter, PRINT_SIZE-counter,"-"); break;
    case 3: counter += snprintf(print_buffer+counter, PRINT_SIZE-counter,"-i"); break;
    default: counter += snprintf(print_buffer+counter, PRINT_SIZE-counter,"error in printing pauli");
    }
    for(int i = 0; i< n; i++){
      char x_p = mzd_read_bit(p->data,i + n,0);
      char z_p = mzd_read_bit(p->data,i ,0);
      counter+=snprintf(print_buffer+counter, PRINT_SIZE-counter,z_p?(x_p?"Y":"Z"):(x_p?"X":"I"));
    }
    if(print_newline_p)counter+=snprintf(print_buffer+counter, PRINT_SIZE-counter,"\n");
    return print_buffer;
  }
  int do_print_pauli(pauli_t* p, char print_newline_p, int counter){
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
    case 1: counter += snprintf(print_buffer+counter, PRINT_SIZE-counter,"i"); break;
    case 2: counter += snprintf(print_buffer+counter, PRINT_SIZE-counter,"-"); break;
    case 3: counter += snprintf(print_buffer+counter, PRINT_SIZE-counter,"-i"); break;
    default: counter += snprintf(print_buffer+counter, PRINT_SIZE-counter,"error in printing pauli");
    }
    for(int i = 0; i< n; i++){
      char x_p = mzd_read_bit(p->data,i + n,0);
      char z_p = mzd_read_bit(p->data,i ,0);
      counter+=snprintf(print_buffer+counter, PRINT_SIZE-counter,z_p?(x_p?"Y":"Z"):(x_p?"X":"I"));
    }
    if(print_newline_p)counter+=snprintf(print_buffer+counter, PRINT_SIZE-counter,"\n");
    return counter;
  }

#ifdef __cplusplus
}
#endif
