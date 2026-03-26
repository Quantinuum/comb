#include "lwlog.h"
#include "comb.h"

int main(int argv, char** args){
  int i1[] = {0,1};
  int i2[] = {0, 1, 1, 2, 2, 0};
  int i3[]={1, 0 ,2, 0, 3, 0, 4, 0 ,4, 1, 4, 2 ,4, 3} ;
  clifford_t* cs[] = {czs(3,1,i1),
                      czs(3,3,i2),
                      czs(5,7,i3)};
  clifford_t* rs[] = {czs(3,1,i1),
                      czs(3,3,i2),
                      czs(5,7,i3)};
  pauli_t* ps[] = {pauli_x(2,3),pauli_x(1,3),pauli_x(2,5)};
  for(int i =0 ; i< 3; i++){
    clifford_t* c = cs[i];
    clifford_t* r = rs[i];
    pauli_t* p = ps[i];
    printf("%d: ",i);
    print_clifford(c,1);
    print_pauli(p,1);
    pauli_z_coset_t* coset = diagonal_clifford_to_coset(c);
    print_pauli_z_coset(coset,1);
    if(!postselect_diagonal_clifford(r,c,p))print_clifford(r,1);
    printf("after postselect %d\n",i);
  }
  printf("about to return 0 \n");
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
