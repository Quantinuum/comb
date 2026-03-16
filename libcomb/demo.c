#include "lwlog.h"
#include "comb.h"

int main(int argv, char** args){
  //comb_dcp_tableau_t* tableau = tableau_init(3);
  //tableau_print(tableau);
  //tableau_free(tableau);
  pauli_t* p = pauli_x(3, 5);
  print_pauli(p, 1);
  mzd_write_bit(p->data, 7, 0, 1);
  print_pauli(p, 1);
  print_pauli(pauli_z(2, 6), 1);
  clifford_t* c =clifford_mul(0, cz(0,2,4),cnot(0,2,4));
  print_clifford(c,0);
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
