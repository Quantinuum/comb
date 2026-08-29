#ifndef COMB_RANK_H
#define COMB_RANK_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Status codes returned by the checked fault-combination rank API.
 */
typedef enum {
  COMB_RANK_OK = 0,
  COMB_RANK_INVALID_ARGUMENT = 1,
  COMB_RANK_OVERFLOW = 2,
  COMB_RANK_OUT_OF_RANGE = 3,
  COMB_RANK_ALLOCATION_FAILED = 4
} comb_rank_status_t;

/**
 * Immutable lookup plan for weighted fault-combination rank and unrank.
 *
 * The structure is opaque so callers cannot invalidate the precomputed
 * dynamic-programming table. A plan owns a copy of the partition sizes.
 */
typedef struct comb_rank_plan comb_rank_plan_t;

/**
 * Construct a checked weighted rank plan.
 *
 * Locations are numbered 0 through location_count - 1. Each partition size
 * must be positive, and faults_per_combination must not exceed
 * location_count. The mathematically useful empty combination is supported:
 * when both values are zero, partition_sizes may be NULL and the total is 1.
 *
 * On every failure, *out_plan is set to NULL. COMB_RANK_OVERFLOW is returned
 * when either a weighted count exceeds UINT64_MAX or the required lookup
 * storage cannot be addressed with size_t.
 */
comb_rank_status_t comb_rank_plan_create(
    uint32_t location_count, const uint32_t *partition_sizes,
    uint32_t faults_per_combination, comb_rank_plan_t **out_plan);

/** Release a plan created by comb_rank_plan_create. NULL is accepted. */
void comb_rank_plan_destroy(comb_rank_plan_t *plan);

/** Return the total number of weighted combinations represented by plan. */
comb_rank_status_t comb_rank_plan_total(const comb_rank_plan_t *plan,
                                        uint64_t *out_total);

/**
 * Rank one weighted fault tuple in the canonical global order.
 *
 * locations must contain faults_per_combination strictly increasing location
 * indices. offsets[i] is a zero-based fault offset in locations[i]. Location
 * subsets are ordered colexicographically. Within one subset, offsets use
 * mixed-radix order with offsets[0] varying fastest.
 */
comb_rank_status_t comb_rank_plan_rank(const comb_rank_plan_t *plan,
                                       const uint32_t *locations,
                                       const uint32_t *offsets,
                                       uint64_t *out_rank);

/**
 * Invert comb_rank_plan_rank for a rank in [0, total).
 *
 * The output arrays must each have room for faults_per_combination entries.
 * For the empty-combination plan they may be NULL.
 */
comb_rank_status_t comb_rank_plan_unrank(const comb_rank_plan_t *plan,
                                         uint64_t rank,
                                         uint32_t *out_locations,
                                         uint32_t *out_offsets);

/** Return a stable human-readable name for a rank status code. */
const char *comb_rank_status_string(comb_rank_status_t status);

#ifdef __cplusplus
}
#endif

#endif
