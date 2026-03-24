#include "lwlog.h"
#include "comb.h"

int main(int argv, char** args){
  int i1[] = {0,1};
  int i2[] = {0, 1, 1, 2, 2, 0};
  int i3[]={1, 0 ,2, 0, 3, 0, 4, 0 ,4, 1, 4, 2 ,4, 3} ;
  clifford_t* cs[] = {czs(2,1,i1),
                      czs(3,3,i2),
                      czs(5,7,i3)};
  for(int i =0 ; i< 3; i++){
     clifford_t* c = cs[i];
    printf("%d: ",i);
    print_clifford(c,1);
    mzd_print(c->c_bar);
    mzd_print(c->h_bar);
    pauli_z_coset_t* coset = diagonal_clifford_to_coset(c);
    if(coset)print_pauli_z_coset(coset,1);
  }

  return 0;
}

/*
 *
 mzd_t *A = mzd_from_str(4, 4, args[1]);
 mzd_print(A);
 //printf("%d",__M4RI_STRASSEN_MUL_CUTOFF);
 lwlog_emerg("This a emerge log.");
 lwlog_alert("This a alert log.");
 lwlog_crit("This a crit log.");
 lwlog_err("This a err log.");
 lwlog_warning("This a warning log.");
 lwlog_notice("This a notice log.%d %d",2, 3);
 lwlog_info("This a info log.");
 lwlog_debug("This a debug log.");
 mzp_t *P = mzp_init(7);
 mzp_print(P);
 printf("\n");
 return 0;
 */
