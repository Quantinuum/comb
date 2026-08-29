#include "comb_rank.h"

#ifdef NDEBUG
#undef NDEBUG
#endif
#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

static unsigned popcount_u32(uint32_t value) {
  unsigned count = 0;
  while (value != 0) {
    count += value & 1U;
    value >>= 1;
  }
  return count;
}

/*
 * Increasing integer masks enumerate equal-popcount subsets in colex order:
 * the largest element in the symmetric difference determines the order. This
 * intentionally simple oracle does not use the production DP recurrence.
 */
static void check_exhaustive(uint32_t location_count,
                             const uint32_t *partitions, uint32_t faults) {
  comb_rank_plan_t *plan = NULL;
  uint64_t total = 0;
  uint64_t expected_rank = 0;
  uint32_t locations[16] = {0};
  uint32_t offsets[16] = {0};
  uint32_t decoded_locations[16] = {0};
  uint32_t decoded_offsets[16] = {0};

  assert(location_count <= 16);
  assert(comb_rank_plan_create(location_count, partitions, faults, &plan) ==
         COMB_RANK_OK);
  assert(plan != NULL);
  assert(comb_rank_plan_total(plan, &total) == COMB_RANK_OK);

  for (uint32_t mask = 0; mask < (UINT32_C(1) << location_count); ++mask) {
    uint64_t subset_size = 1;
    size_t selected = 0;

    if (popcount_u32(mask) != faults) continue;
    for (uint32_t location = 0; location < location_count; ++location) {
      if ((mask & (UINT32_C(1) << location)) != 0) {
        locations[selected++] = location;
        subset_size *= partitions[location];
      }
    }
    assert(selected == faults);

    for (uint64_t local = 0; local < subset_size; ++local) {
      uint64_t quotient = local;
      uint64_t actual_rank = UINT64_MAX;

      for (size_t i = 0; i < selected; ++i) {
        offsets[i] = (uint32_t)(quotient % partitions[locations[i]]);
        quotient /= partitions[locations[i]];
      }
      assert(quotient == 0);

      assert(comb_rank_plan_rank(plan, locations, offsets, &actual_rank) ==
             COMB_RANK_OK);
      assert(actual_rank == expected_rank);

      memset(decoded_locations, 0xff, sizeof(decoded_locations));
      memset(decoded_offsets, 0xff, sizeof(decoded_offsets));
      assert(comb_rank_plan_unrank(plan, expected_rank, decoded_locations,
                                   decoded_offsets) == COMB_RANK_OK);
      assert(memcmp(decoded_locations, locations,
                    selected * sizeof(locations[0])) == 0);
      assert(memcmp(decoded_offsets, offsets,
                    selected * sizeof(offsets[0])) == 0);
      ++expected_rank;
    }
  }

  assert(expected_rank == total);
  assert(comb_rank_plan_unrank(plan, total, decoded_locations,
                               decoded_offsets) == COMB_RANK_OUT_OF_RANGE);
  comb_rank_plan_destroy(plan);
}

static void check_empty_combination(void) {
  comb_rank_plan_t *plan = NULL;
  uint64_t total = 0;
  uint64_t rank = UINT64_MAX;

  assert(comb_rank_plan_create(0, NULL, 0, &plan) == COMB_RANK_OK);
  assert(comb_rank_plan_total(plan, &total) == COMB_RANK_OK);
  assert(total == 1);
  assert(comb_rank_plan_rank(plan, NULL, NULL, &rank) == COMB_RANK_OK);
  assert(rank == 0);
  assert(comb_rank_plan_unrank(plan, 0, NULL, NULL) == COMB_RANK_OK);
  assert(comb_rank_plan_unrank(plan, 1, NULL, NULL) ==
         COMB_RANK_OUT_OF_RANGE);
  comb_rank_plan_destroy(plan);
}

static void check_small_partition_family(void) {
  uint32_t partitions[5] = {0};

  /* Exhaust every assignment p_i in {1,2,3}, every s, and every tuple. */
  for (uint32_t location_count = 1; location_count <= 5; ++location_count) {
    uint32_t assignment_count = 1;
    for (uint32_t i = 0; i < location_count; ++i) assignment_count *= 3;

    for (uint32_t assignment = 0; assignment < assignment_count;
         ++assignment) {
      uint32_t encoded = assignment;
      for (uint32_t i = 0; i < location_count; ++i) {
        partitions[i] = encoded % 3 + 1;
        encoded /= 3;
      }
      for (uint32_t faults = 0; faults <= location_count; ++faults) {
        check_exhaustive(location_count, partitions, faults);
      }
    }
  }
}

static void check_invalid_inputs(void) {
  const uint32_t partitions[] = {2, 3, 4};
  const uint32_t zero_partition[] = {2, 0, 4};
  comb_rank_plan_t *plan = (comb_rank_plan_t *)(uintptr_t)1;
  uint64_t value = 0;
  uint32_t locations[] = {0, 2};
  uint32_t offsets[] = {0, 0};

  assert(comb_rank_plan_create(3, partitions, 2, NULL) ==
         COMB_RANK_INVALID_ARGUMENT);
  assert(comb_rank_plan_create(3, NULL, 2, &plan) ==
         COMB_RANK_INVALID_ARGUMENT);
  assert(plan == NULL);
  plan = (comb_rank_plan_t *)(uintptr_t)1;
  assert(comb_rank_plan_create(3, partitions, 4, &plan) ==
         COMB_RANK_INVALID_ARGUMENT);
  assert(plan == NULL);
  plan = (comb_rank_plan_t *)(uintptr_t)1;
  assert(comb_rank_plan_create(3, zero_partition, 2, &plan) ==
         COMB_RANK_INVALID_ARGUMENT);
  assert(plan == NULL);

  assert(comb_rank_plan_create(3, partitions, 2, &plan) == COMB_RANK_OK);
  assert(comb_rank_plan_total(NULL, &value) == COMB_RANK_INVALID_ARGUMENT);
  assert(comb_rank_plan_total(plan, NULL) == COMB_RANK_INVALID_ARGUMENT);
  assert(comb_rank_plan_rank(NULL, locations, offsets, &value) ==
         COMB_RANK_INVALID_ARGUMENT);
  assert(comb_rank_plan_rank(plan, NULL, offsets, &value) ==
         COMB_RANK_INVALID_ARGUMENT);
  assert(comb_rank_plan_rank(plan, locations, NULL, &value) ==
         COMB_RANK_INVALID_ARGUMENT);
  assert(comb_rank_plan_rank(plan, locations, offsets, NULL) ==
         COMB_RANK_INVALID_ARGUMENT);
  assert(comb_rank_plan_unrank(NULL, 0, locations, offsets) ==
         COMB_RANK_INVALID_ARGUMENT);
  assert(comb_rank_plan_unrank(plan, 0, NULL, offsets) ==
         COMB_RANK_INVALID_ARGUMENT);
  assert(comb_rank_plan_unrank(plan, 0, locations, NULL) ==
         COMB_RANK_INVALID_ARGUMENT);

  locations[0] = 1;
  locations[1] = 1;
  assert(comb_rank_plan_rank(plan, locations, offsets, &value) ==
         COMB_RANK_INVALID_ARGUMENT);
  locations[0] = 2;
  locations[1] = 1;
  assert(comb_rank_plan_rank(plan, locations, offsets, &value) ==
         COMB_RANK_INVALID_ARGUMENT);
  locations[0] = 0;
  locations[1] = 3;
  assert(comb_rank_plan_rank(plan, locations, offsets, &value) ==
         COMB_RANK_INVALID_ARGUMENT);
  locations[1] = 2;
  offsets[1] = 4;
  assert(comb_rank_plan_rank(plan, locations, offsets, &value) ==
         COMB_RANK_OUT_OF_RANGE);

  comb_rank_plan_destroy(plan);
  comb_rank_plan_destroy(NULL);
}

static void check_overflow_and_boundary(void) {
  const uint32_t maximum_partitions[] = {UINT32_MAX, UINT32_MAX,
                                         UINT32_MAX};
  comb_rank_plan_t *plan = NULL;
  uint64_t total = 0;
  uint64_t rank = 0;
  uint32_t locations[] = {0, 1};
  uint32_t offsets[] = {UINT32_MAX - 1, UINT32_MAX - 1};
  const uint64_t expected =
      (uint64_t)UINT32_MAX * (uint64_t)UINT32_MAX;

  assert(comb_rank_plan_create(3, maximum_partitions, 3, &plan) ==
         COMB_RANK_OVERFLOW);
  assert(plan == NULL);

  assert(comb_rank_plan_create(2, maximum_partitions, 2, &plan) ==
         COMB_RANK_OK);
  assert(comb_rank_plan_total(plan, &total) == COMB_RANK_OK);
  assert(total == expected);
  assert(comb_rank_plan_rank(plan, locations, offsets, &rank) ==
         COMB_RANK_OK);
  assert(rank == total - 1);
  offsets[0] = offsets[1] = 0;
  assert(comb_rank_plan_unrank(plan, total - 1, locations, offsets) ==
         COMB_RANK_OK);
  assert(offsets[0] == UINT32_MAX - 1);
  assert(offsets[1] == UINT32_MAX - 1);
  comb_rank_plan_destroy(plan);

  assert(strcmp(comb_rank_status_string(COMB_RANK_OK), "ok") == 0);
  assert(strcmp(comb_rank_status_string((comb_rank_status_t)999),
                "unknown rank status") == 0);
}

int main(void) {
  const uint32_t singleton[] = {1};
  const uint32_t s_one[] = {2, 1, 3, 2};
  const uint32_t heterogeneous[] = {2, 3, 1, 4, 2};
  const uint32_t unit[] = {1, 1, 1, 1, 1, 1};

  check_empty_combination();
  check_exhaustive(1, singleton, 1);
  check_exhaustive(4, s_one, 1);
  check_exhaustive(5, heterogeneous, 2);
  check_exhaustive(5, heterogeneous, 3);
  check_exhaustive(6, unit, 3);
  check_small_partition_family();
  check_invalid_inputs();
  check_overflow_and_boundary();

  puts("comb_rank: all tests passed");
  return 0;
}
