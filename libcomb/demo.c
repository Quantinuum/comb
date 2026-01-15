#include <stdlib.h>
#include <m4ri/m4ri.h>

int main(char** args){
	mzd_t *A = mzd_from_str(4, 4, "1000010000100001");
	mzd_print(A);
	//printf("%d",__M4RI_STRASSEN_MUL_CUTOFF);
	return 0;
}
