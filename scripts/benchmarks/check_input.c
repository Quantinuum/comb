/* Validate an unchanged serialized example before benchmarking. Compile with
 * COMB_BENCHMARK_M4RI_CPU to additionally compare the full main CPU decision.
 * Reference limits are explicit: pilots are never presented as full scans. */
#include "comb_cuda_internal.h"
#include "comb_rank.h"
#include "comb_reference.h"
#include "comb_run_file.h"

#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char **argv) {
  comb_run_file_t run = {0};
  comb_rank_plan_t *plan = NULL;
  comb_reference_context_t *reference = NULL;
  comb_cuda_result_t result = COMB_CUDA_RESULT_INITIALIZER;
  void *blob = NULL;
  uint64_t blob_size = 0, total = 0, limit = 0;
  int cpu = 0, code = EXIT_FAILURE;
  comb_cuda_status_t status;
  if (argc < 2) {
    fprintf(stderr, "Usage: check_input INPUT.run [REFERENCE_COUNT] [--cpu]\n");
    return EXIT_FAILURE;
  }
  if (argc >= 3) limit = strtoull(argv[2], NULL, 10);
  if (argc >= 4 && strcmp(argv[3], "--cpu") == 0) cpu = 1;
  status = comb_run_file_load(argv[1], &run);
  if (status != COMB_CUDA_STATUS_SUCCESS) goto done;
  status = comb_cuda_flat_problem_pack(&run.problem, &blob, &blob_size);
  if (status != COMB_CUDA_STATUS_SUCCESS) goto done;
  status = comb_cuda_flat_problem_validate_science(blob, blob_size);
  if (status != COMB_CUDA_STATUS_SUCCESS) goto done;
  if (comb_rank_plan_create(run.problem.fault_location_count,
                           run.problem.fault_location_partitions,
                           run.problem.s, &plan) != COMB_RANK_OK ||
      comb_rank_plan_total(plan, &total) != COMB_RANK_OK) goto done;
  printf("input=%s n_all=%u n=%u s=%u total=%" PRIu64 " flat_bytes=%" PRIu64 "\n",
         argv[1], run.problem.n_all, run.problem.n, run.problem.s, total, blob_size);
  if (limit != 0) {
    if (limit > total) limit = total;
    status = comb_reference_context_create(&run.problem, &reference);
    if (status != COMB_CUDA_STATUS_SUCCESS) goto done;
    status = comb_reference_verify_range(reference, 0, limit,
                COMB_CUDA_WITNESS_DETERMINISTIC, &result, NULL);
    if (status != COMB_CUDA_STATUS_SUCCESS) goto done;
    printf("reference_checked=%" PRIu64 " reference_full=%s code=%d bad_rank=%" PRIu64 "\n",
           result.checked_combination_count, limit == total ? "yes" : "no",
           result.verification_code, result.minimum_bad_rank);
    if (result.verification_code != ALL_FAULTS_GOOD) goto done;
  }
  if (cpu) {
#ifdef COMB_BENCHMARK_M4RI_CPU
    const comb_problem_t *p = &run.problem;
    uint32_t *faults = calloc(p->s, sizeof(*faults));
    comb_code_t actual;
    if (faults == NULL) goto done;
    actual = comb_verify_css_magic_state(
        p->n_all, p->n, p->n_x, p->n_z, p->k,
        (comb_pauli_x_t *)p->x_stabilisers, (comb_bit_array_t)p->x_lookup_table,
        (comb_pauli_z_t *)p->z_stabilisers, (comb_bit_array_t)p->z_lookup_table,
        (comb_dcp_t *)p->dcp_stabilisers, p->fault_location_count,
        (uint32_t *)p->fault_location_partitions, p->s, faults,
        p->output_error_count, (comb_dcp_t *)p->output_errors,
        p->x_measurement_count, (comb_pauli_x_t *)p->x_measurement_paulis,
        p->z_measurement_count, (comb_pauli_z_t *)p->z_measurement_paulis, 1U);
    free(faults);
    printf("main_cpu_full=yes code=%d\n", actual);
    if (actual != ALL_FAULTS_GOOD) goto done;
#else
    fputs("--cpu requires a build with COMB_BENCHMARK_M4RI_CPU\n", stderr);
    goto done;
#endif
  }
  code = EXIT_SUCCESS;
done:
  if (code != EXIT_SUCCESS) fprintf(stderr, "validation failed, status=%d\n", status);
  comb_reference_context_destroy(reference);
  comb_rank_plan_destroy(plan);
  comb_cuda_flat_problem_free(blob);
  comb_run_file_free(&run);
  return code;
}
