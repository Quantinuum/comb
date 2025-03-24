#include "comb.hpp"
#include <cassert>
#include <cstdint>
#include <iostream>
#include <stdint.h>
#include <vector>

static int g(uint64_t x) { return x % 3 == 0; }
static int g1(uint64_t x) { return 0; }

static void run_test(const std::vector<std::vector<uint64_t>> &faults,
                     unsigned k, eval_fn g, unsigned n_worker_threads = 8) {
  std::optional<std::vector<std::pair<unsigned, unsigned>>> bad_faults =
      evaluate(faults, k, g, n_worker_threads);
  std::cout << "============" << std::endl;
  if (bad_faults) {
    std::cout << "bad fault combination found: [";
    uint64_t x = 0;
    for (const std::pair<unsigned, unsigned> &bucket_and_index : *bad_faults) {
      unsigned bucket = bucket_and_index.first;
      unsigned index = bucket_and_index.second;
      std::cout << "(" << bucket << "," << index << "), ";
      x ^= faults[bucket][index];
    }
    std::cout << "]" << std::endl;
    assert(g(x) != 0);
  } else {
    std::cout << "no bad fault combinations" << std::endl;
  }
}

int main() {
  run_test({{1, 1, 2}, {3}, {5, 8}}, 2, g);
  run_test({{1, 1, 2}, {5}, {8, 13}}, 2, g, 1);
  run_test({{1, 1, 2}, {3}, {5, 8}}, 2, g1);
  run_test({{1, 2, 3, 4}, {5, 6}, {}, {7, 8, 9}}, 3, g, 24);
  run_test({{1, 2, 4, 5}}, 1, g);
  return 0;
}
