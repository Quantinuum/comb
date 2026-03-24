#include "lwlog.h"
#include "comb.h"

int main(int argv, char** args){
  clifford_t* cs[] = {pauli_z_transvection(pauli_z(0,1),1),
                     cz(0,1,2),
                     clifford_mul(0,pauli_z_transvection(pauli_z(1,2),1),
                                  pauli_z_transvection(pauli_z(0,2),1)),
                     clifford_mul(0,cz(0,1,3),cz(1,2,3)),
                     clifford_mul(0,clifford_mul(0,cz(0,1,3),cz(1,2,3)),cz(2,0,3))};
  for(int i =0 ; i< 5; i++){
     clifford_t* c = cs[i];
    printf("%d: ",i);
    print_clifford(c,1);
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
