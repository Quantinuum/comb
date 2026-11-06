#ifndef COMB_REFERENCE_H
#define COMB_REFERENCE_H

/*
 * Dependency-free sequential scientific reference implementation for Comb.
 *
 * This path consumes the same normalized FlatProblem as the CUDA backend but
 * intentionally retains the pre-batched measurement sequence.  It is an
 * independent differential implementation.  The merged main M4RI CPU path is
 * the compatibility source of truth; this path does not establish equivalence
 * with its batched-Z transform by itself.
 */

#include <stdint.h>

#include "comb_cuda.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct comb_reference_context comb_reference_context_t;

/* Validate, normalize, and own a complete copy of problem. */
comb_cuda_status_t comb_reference_context_create(
    const comb_problem_t *problem, comb_reference_context_t **context);

void comb_reference_context_destroy(comb_reference_context_t *context);

/*
 * Verify every rank in [rank_begin, rank_end).  The reference implementation
 * always completes the range and returns its lowest bad rank; FAST_ANY is
 * accepted for API compatibility but does not weaken deterministic ordering.
 */
comb_cuda_status_t comb_reference_verify_range(
    comb_reference_context_t *context, uint64_t rank_begin,
    uint64_t rank_end, comb_cuda_witness_mode_t mode,
    comb_cuda_result_t *result, comb_cuda_witness_t *witness);

comb_cuda_status_t comb_reference_verify(
    comb_reference_context_t *context, comb_cuda_witness_mode_t mode,
    comb_cuda_result_t *result, comb_cuda_witness_t *witness);

#ifdef __cplusplus
}
#endif

#endif /* COMB_REFERENCE_H */
