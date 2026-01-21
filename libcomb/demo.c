#include <stdlib.h>
#include <m4ri/m4ri.h>
#include "lwlog.h"
#include "comb.h"

int main(int argv, char** args){
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
}
