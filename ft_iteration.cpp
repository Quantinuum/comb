#include "ft_iteration.h"
#include "comb.hpp"
#include <iostream>
#include <vector>

int ft_iteration_iterate(int s, int *m, int k, fault_t *propagated_faults,
                         int *bad_fault_array, ft_iteration_evaluator evaluator,
                         ft_iteration_updater updater, int update_stride,
                         unsigned n_worker_threads) {
  if (updater) {
    std::cout << "Updater function not supported: ignoring." << std::endl;
  }
  (void)update_stride;
  std::vector<std::vector<fault_t>> faults(s);
  unsigned p = 0;
  for (int i = 0; i < s; i++) {
    std::vector<fault_t> bucket(m[i]);
    for (int j = 0; j < m[i]; j++) {
      bucket[j] = propagated_faults[p];
      p++;
    }
    faults[i] = bucket;
  }
  std::optional<std::vector<std::pair<unsigned, unsigned>>> bad_faults =
      evaluate(faults, k, evaluator, n_worker_threads);
  if (bad_faults) {
    std::vector<unsigned> r(s);
    r[0] = 0;
    for (int i = 0; i + 1 < s; i++) {
      r[i + 1] = r[i] + m[i];
    }
    for (int i = 0; i < k; i++) {
      std::pair<unsigned, unsigned> bucket_and_index = bad_faults->at(i);
      bad_fault_array[i] = r[bucket_and_index.first] + bucket_and_index.second;
    }
    return BAD_FAULT_FOUND;
  } else {
    return ALL_FAULTS_GOOD;
  }
}
