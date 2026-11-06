#ifndef COMB_RUN_FILE_H
#define COMB_RUN_FILE_H

#include "comb_cuda.h"

/* Owner of a legacy Comb .run problem.  Zero-initialize before loading; free
 * releases every allocation, including allocations made by a failed load. */
typedef struct {
  comb_problem_t problem;
} comb_run_file_t;

comb_cuda_status_t comb_run_file_load(const char *path, comb_run_file_t *run);
void comb_run_file_free(comb_run_file_t *run);

#endif
