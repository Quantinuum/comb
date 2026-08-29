#include "comb_rank.h"

#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

struct comb_rank_plan {
  uint32_t location_count;
  uint32_t faults_per_combination;
  size_t dp_columns;
  uint32_t *partition_sizes;
  uint64_t *weighted_counts;
  uint64_t total;
};

static int checked_size_multiply(size_t left, size_t right, size_t *result) {
  if (left != 0 && right > SIZE_MAX / left) return 0;
  *result = left * right;
  return 1;
}

static int checked_u64_add(uint64_t left, uint64_t right, uint64_t *result) {
  if (right > UINT64_MAX - left) return 0;
  *result = left + right;
  return 1;
}

static int checked_u64_multiply(uint64_t left, uint64_t right,
                                uint64_t *result) {
  if (left != 0 && right > UINT64_MAX / left) return 0;
  *result = left * right;
  return 1;
}

static size_t dp_index(const comb_rank_plan_t *plan, size_t locations,
                       size_t faults) {
  return locations * plan->dp_columns + faults;
}

static uint64_t weighted_count(const comb_rank_plan_t *plan, size_t locations,
                               size_t faults) {
  return plan->weighted_counts[dp_index(plan, locations, faults)];
}

static comb_rank_status_t subset_start_and_size(
    const comb_rank_plan_t *plan, const uint32_t *locations,
    uint64_t *out_start, uint64_t *out_size) {
  size_t selected_count = (size_t)plan->faults_per_combination;
  uint64_t start = 0;
  uint64_t higher_partition_product = 1;

  /*
   * Expanding the weighted colex recurrence gives
   *
   *   start(c_0,...,c_{s-1}) = W[c_{s-1},s]
   *       + p[c_{s-1}] W[c_{s-2},s-1] + ... .
   *
   * Iterating from the highest selected location avoids recursion whose depth
   * would otherwise depend on the total location count.
   */
  for (size_t j = selected_count; j > 0; --j) {
    size_t location = (size_t)locations[j - 1];
    uint64_t term;
    uint64_t new_start;
    uint64_t new_product;

    if (!checked_u64_multiply(higher_partition_product,
                              weighted_count(plan, location, j), &term) ||
        !checked_u64_add(start, term, &new_start) ||
        !checked_u64_multiply(
            higher_partition_product,
            (uint64_t)plan->partition_sizes[location], &new_product)) {
      return COMB_RANK_OVERFLOW;
    }
    start = new_start;
    higher_partition_product = new_product;
  }

  *out_start = start;
  *out_size = higher_partition_product;
  return COMB_RANK_OK;
}

comb_rank_status_t comb_rank_plan_create(
    uint32_t location_count, const uint32_t *partition_sizes,
    uint32_t faults_per_combination, comb_rank_plan_t **out_plan) {
  comb_rank_plan_t *plan = NULL;
  size_t location_rows;
  size_t fault_columns;
  size_t cell_count;
  size_t table_bytes;
  size_t partition_bytes = 0;

  if (out_plan == NULL) return COMB_RANK_INVALID_ARGUMENT;
  *out_plan = NULL;

  if (faults_per_combination > location_count ||
      (location_count > 0 && partition_sizes == NULL)) {
    return COMB_RANK_INVALID_ARGUMENT;
  }

  if ((size_t)location_count == SIZE_MAX ||
      (size_t)faults_per_combination == SIZE_MAX) {
    return COMB_RANK_OVERFLOW;
  }

  location_rows = (size_t)location_count + 1;
  fault_columns = (size_t)faults_per_combination + 1;
  if (!checked_size_multiply(location_rows, fault_columns, &cell_count) ||
      !checked_size_multiply(cell_count, sizeof(uint64_t), &table_bytes) ||
      !checked_size_multiply((size_t)location_count, sizeof(uint32_t),
                             &partition_bytes)) {
    return COMB_RANK_OVERFLOW;
  }

  for (size_t location = 0; location < (size_t)location_count; ++location) {
    if (partition_sizes[location] == 0) return COMB_RANK_INVALID_ARGUMENT;
  }

  plan = (comb_rank_plan_t *)calloc(1, sizeof(*plan));
  if (plan == NULL) return COMB_RANK_ALLOCATION_FAILED;

  if (partition_bytes != 0) {
    plan->partition_sizes = (uint32_t *)malloc(partition_bytes);
    if (plan->partition_sizes == NULL) {
      comb_rank_plan_destroy(plan);
      return COMB_RANK_ALLOCATION_FAILED;
    }
    memcpy(plan->partition_sizes, partition_sizes, partition_bytes);
  }

  plan->weighted_counts = (uint64_t *)calloc(1, table_bytes);
  if (plan->weighted_counts == NULL) {
    comb_rank_plan_destroy(plan);
    return COMB_RANK_ALLOCATION_FAILED;
  }

  plan->location_count = location_count;
  plan->faults_per_combination = faults_per_combination;
  plan->dp_columns = fault_columns;

  for (size_t locations = 0; locations < location_rows; ++locations) {
    plan->weighted_counts[dp_index(plan, locations, 0)] = 1;
  }

  for (size_t locations = 1; locations < location_rows; ++locations) {
    size_t max_faults = locations;
    if (max_faults > (size_t)faults_per_combination) {
      max_faults = (size_t)faults_per_combination;
    }
    for (size_t faults = 1; faults <= max_faults; ++faults) {
      uint64_t included;
      uint64_t count;
      uint64_t excluded = weighted_count(plan, locations - 1, faults);
      uint64_t prefix = weighted_count(plan, locations - 1, faults - 1);
      uint64_t partition =
          (uint64_t)plan->partition_sizes[locations - 1];

      if (!checked_u64_multiply(partition, prefix, &included) ||
          !checked_u64_add(excluded, included, &count)) {
        comb_rank_plan_destroy(plan);
        return COMB_RANK_OVERFLOW;
      }
      plan->weighted_counts[dp_index(plan, locations, faults)] = count;
    }
  }

  plan->total = weighted_count(plan, (size_t)location_count,
                               (size_t)faults_per_combination);
  *out_plan = plan;
  return COMB_RANK_OK;
}

void comb_rank_plan_destroy(comb_rank_plan_t *plan) {
  if (plan == NULL) return;
  free(plan->weighted_counts);
  free(plan->partition_sizes);
  free(plan);
}

comb_rank_status_t comb_rank_plan_total(const comb_rank_plan_t *plan,
                                        uint64_t *out_total) {
  if (plan == NULL || out_total == NULL) return COMB_RANK_INVALID_ARGUMENT;
  *out_total = plan->total;
  return COMB_RANK_OK;
}

comb_rank_status_t comb_rank_plan_rank(const comb_rank_plan_t *plan,
                                       const uint32_t *locations,
                                       const uint32_t *offsets,
                                       uint64_t *out_rank) {
  uint64_t subset_start;
  uint64_t subset_size;
  uint64_t local_rank = 0;
  uint64_t stride = 1;
  uint64_t rank;
  size_t selected_count;

  if (plan == NULL || out_rank == NULL) return COMB_RANK_INVALID_ARGUMENT;
  selected_count = (size_t)plan->faults_per_combination;
  if (selected_count != 0 && (locations == NULL || offsets == NULL)) {
    return COMB_RANK_INVALID_ARGUMENT;
  }

  for (size_t i = 0; i < selected_count; ++i) {
    uint32_t location = locations[i];
    uint64_t term;
    uint64_t new_local_rank;
    uint64_t new_stride;

    if (location >= plan->location_count ||
        (i != 0 && locations[i - 1] >= location)) {
      return COMB_RANK_INVALID_ARGUMENT;
    }
    if (offsets[i] >= plan->partition_sizes[location]) {
      return COMB_RANK_OUT_OF_RANGE;
    }

    if (!checked_u64_multiply(stride, (uint64_t)offsets[i], &term) ||
        !checked_u64_add(local_rank, term, &new_local_rank) ||
        !checked_u64_multiply(
            stride, (uint64_t)plan->partition_sizes[location], &new_stride)) {
      return COMB_RANK_OVERFLOW;
    }
    local_rank = new_local_rank;
    stride = new_stride;
  }

  comb_rank_status_t status =
      subset_start_and_size(plan, locations, &subset_start, &subset_size);
  if (status != COMB_RANK_OK) return status;
  if (stride != subset_size || local_rank >= subset_size) {
    return COMB_RANK_INVALID_ARGUMENT;
  }
  if (!checked_u64_add(subset_start, local_rank, &rank)) {
    return COMB_RANK_OVERFLOW;
  }
  if (rank >= plan->total) return COMB_RANK_OUT_OF_RANGE;

  *out_rank = rank;
  return COMB_RANK_OK;
}

comb_rank_status_t comb_rank_plan_unrank(const comb_rank_plan_t *plan,
                                         uint64_t rank,
                                         uint32_t *out_locations,
                                         uint32_t *out_offsets) {
  size_t remaining_locations;
  size_t remaining_faults;
  uint64_t membership_rank;
  uint64_t subset_start;
  uint64_t subset_size;
  uint64_t local_rank;
  comb_rank_status_t status;

  if (plan == NULL) return COMB_RANK_INVALID_ARGUMENT;
  if (rank >= plan->total) return COMB_RANK_OUT_OF_RANGE;
  if (plan->faults_per_combination != 0 &&
      (out_locations == NULL || out_offsets == NULL)) {
    return COMB_RANK_INVALID_ARGUMENT;
  }

  remaining_locations = (size_t)plan->location_count;
  remaining_faults = (size_t)plan->faults_per_combination;
  membership_rank = rank;

  while (remaining_faults > 0) {
    uint64_t excluded_count;
    uint64_t partition;

    if (remaining_locations == 0) return COMB_RANK_INVALID_ARGUMENT;
    excluded_count = weighted_count(plan, remaining_locations - 1,
                                    remaining_faults);
    if (membership_rank < excluded_count) {
      --remaining_locations;
      continue;
    }

    partition =
        (uint64_t)plan->partition_sizes[remaining_locations - 1];
    membership_rank = (membership_rank - excluded_count) / partition;
    out_locations[remaining_faults - 1] =
        (uint32_t)(remaining_locations - 1);
    --remaining_faults;
    --remaining_locations;
  }

  status = subset_start_and_size(plan, out_locations, &subset_start,
                                 &subset_size);
  if (status != COMB_RANK_OK) return status;
  if (rank < subset_start || rank - subset_start >= subset_size) {
    return COMB_RANK_INVALID_ARGUMENT;
  }

  local_rank = rank - subset_start;
  for (size_t i = 0; i < (size_t)plan->faults_per_combination; ++i) {
    uint64_t partition =
        (uint64_t)plan->partition_sizes[out_locations[i]];
    out_offsets[i] = (uint32_t)(local_rank % partition);
    local_rank /= partition;
  }
  if (local_rank != 0) return COMB_RANK_INVALID_ARGUMENT;
  return COMB_RANK_OK;
}

const char *comb_rank_status_string(comb_rank_status_t status) {
  switch (status) {
    case COMB_RANK_OK:
      return "ok";
    case COMB_RANK_INVALID_ARGUMENT:
      return "invalid argument";
    case COMB_RANK_OVERFLOW:
      return "overflow";
    case COMB_RANK_OUT_OF_RANGE:
      return "out of range";
    case COMB_RANK_ALLOCATION_FAILED:
      return "allocation failed";
    default:
      return "unknown rank status";
  }
}
