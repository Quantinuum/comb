#include "comb.h"

#include <m4ri/m4ri_config.h>
#include <assert.h>
#include <omp.h>
#include <stdio.h>
#include <string.h>

#if __M4RI_ENABLE_MZD_CACHE || __M4RI_ENABLE_MMC
#error "OpenMP checks require M4RI configured with --enable-thread-safe"
#endif

// C(128,2) = 8128 pairs, split into eleven chunks by the current scheduler.
enum { N = 3, LOCATIONS = 128, S = 2, CHUNK_PAIRS = 769, CHUNKS = 11 };

static int has_x_error(unsigned location, unsigned offset, int mixed) {
  return mixed && location >= 96U && ((location + offset) % 2U == 0U);
}

int main(int argc, char **argv) {
  assert(argc == 3);
  const int formatting = strcmp(argv[1], "format") == 0;
  const int all_bad = strcmp(argv[1], "bad") == 0;
  const int mixed = strcmp(argv[1], "mixed") == 0;
  assert(formatting || all_bad || mixed || strcmp(argv[1], "good") == 0);
  const int requested_threads = atoi(argv[2]);
  int actual_threads = 0;
  #pragma omp parallel
  {
    #pragma omp single
    actual_threads = omp_get_num_threads();
  }
  assert(actual_threads == requested_threads);
  assert(actual_threads <= CHUNKS);

  if (formatting) {
    #pragma omp parallel
    {
      const int worker = omp_get_thread_num();
      mzd_t *row = mzd_init(1, CHUNKS);
      mzd_write_bit(row, 0, worker, 1);
      char expected[CHUNKS + 1];
      memset(expected, 'I', CHUNKS);
      expected[worker] = 'X';
      expected[CHUNKS] = '\0';
      char *rendered = print_row_vector_as_pauli(row, 0, "X");
      // Other workers must not overwrite the returned string.
      #pragma omp barrier
      assert(strcmp(rendered, expected) == 0);
      mzd_free(row);
    }
    printf("openmp-format threads=%d passed\n", actual_threads);
    return EXIT_SUCCESS;
  }

  uint32_t partitions[LOCATIONS], prefixes[LOCATIONS + 1] = {0};
  uint32_t bad[S] = {UINT32_MAX, UINT32_MAX};
  uint32_t expected_bad[S] = {UINT32_MAX, UINT32_MAX};
  for (unsigned location = 0; location < LOCATIONS; ++location) {
    partitions[location] = mixed ? 1U + location % 3U : 1U;
    prefixes[location + 1] = prefixes[location] + partitions[location];
  }
  const uint32_t error_count = prefixes[LOCATIONS];
  comb_pauli_x_t *xs = comb_calloc_vector(1, N);
  comb_pauli_z_t *zs = comb_calloc_vector(1, N);
  comb_dcp_t *logical = comb_calloc_dcp(1, N);
  comb_dcp_t *errors = comb_calloc_dcp((int)error_count, N);
  comb_bit_array_t xl = comb_alloc_bit_array(4);
  comb_bit_array_t zl = comb_alloc_bit_array(4);
  assert(xs && zs && logical && errors && xl && zl);
  comb_set_cvector_aref(xs, 0, 1, 1); /* X_1 */
  comb_set_cvector_aref(zs, 0, 0, 1); /* Z_0; logical qubit is 2. */
  for (int bit = 0; bit < 4; ++bit) {
    comb_set_bit_array_aref(xl, bit, !all_bad && (!mixed || bit == 0));
    comb_set_bit_array_aref(zl, bit, !all_bad);
  }
  for (unsigned location = 0; location < LOCATIONS; ++location)
    for (unsigned offset = 0; offset < partitions[location]; ++offset)
      comb_set_cdcp_h_x_aref(errors, (int)(prefixes[location] + offset), 2,
                            has_x_error(location, offset, mixed));

  // Count expected checks in colex order without using Comb's iterator.
  // In the mixed case, exactly one X_2 error makes a pair bad.
  unsigned long long expected_count = 0;
  unsigned pair_rank = 0, chunk_count = 0;
  int chunk_bad = 0, any_bad = 0;
  for (unsigned second = 1; second < LOCATIONS; ++second) {
    for (unsigned first = 0; first < second; ++first, ++pair_rank) {
      if (pair_rank % CHUNK_PAIRS == 0) {
        ++chunk_count;
        chunk_bad = 0;
      }
      for (unsigned b = 0; b < partitions[second] && !chunk_bad; ++b) {
        for (unsigned a = 0; a < partitions[first] && !chunk_bad; ++a) {
          ++expected_count;
          if (all_bad || (has_x_error(first, a, mixed) ^
                          has_x_error(second, b, mixed))) {
            chunk_bad = any_bad = 1;
            expected_bad[0] = prefixes[first] + a;
            expected_bad[1] = prefixes[second] + b;
          }
        }
      }
    }
  }
  assert(pair_rank == 8128U && chunk_count == CHUNKS);
  assert(mixed || expected_count == (all_bad ? 11ULL : 8128ULL));
  if (all_bad) assert(expected_bad[0] == 64U && expected_bad[1] == 124U);

  const comb_code_t code = comb_verify_css_magic_state(
      N, N, 1, 1, 1, xs, xl, zs, zl, logical, LOCATIONS, partitions, S,
      bad, error_count, errors, 0, NULL, 0, NULL, 999);
  assert(code == (any_bad ? BAD_FAULT_FOUND : ALL_FAULTS_GOOD));
  assert(bad[0] == expected_bad[0] && bad[1] == expected_bad[1]);
  printf("openmp-test mode=%s threads=%d expected_count=%llu witness=%u,%u\n",
         argv[1], actual_threads, expected_count, bad[0], bad[1]);
  comb_cfree_vector(1, xs);
  comb_cfree_vector(1, zs);
  comb_cfree_dcp(1, logical);
  comb_cfree_dcp((int)error_count, errors);
  comb_free_bit_array(xl);
  comb_free_bit_array(zl);
  return EXIT_SUCCESS;
}
