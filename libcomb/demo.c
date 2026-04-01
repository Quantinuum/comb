#include "lwlog.h"
#include "comb.h"

int main(int argv, char** args){
  dcp_t* dcp = dcp_pauli_z(0,1);//dcp_mul(0, dcp_cz(0,1,3), dcp_mul(0, dcp_cz(0,2,3), dcp_cz(1,2,3) ));
  print_dcp(dcp,1);
  stabiliser_state_t* state = dcp_to_stabiliser_state(dcp);
  print_stabiliser_state(state,1);
  printf("postselect: %d\n",postselect_on_measuring_plus_one(state,pauli_x(0,1)));
  printf("after postselection\n");
  print_stabiliser_state(state,1);
  pauli_z_coset_t* coset = diagonal_clifford_to_coset(state);
  print_pauli_z_coset(coset,1);
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
