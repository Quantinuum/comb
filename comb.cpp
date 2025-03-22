#include "comb.hpp"
#include <iostream>
#include <optional>
#include <vector>

/// Class to generate all combinations of k distinct integers in the range
/// [0,n).
class CombsGenerator {
public:
  CombsGenerator(unsigned n, unsigned k)
      : n(n), k(k), initial(true), indices(k) {
    for (unsigned i = 0; i < k; i++) {
      indices[i] = i;
    }
  }
  // Return the next combination, or nullopt if we're finished.
  std::optional<std::vector<unsigned>> next() {
    if (initial) {
      initial = false;
      return indices;
    } else if (update()) {
      return indices;
    } else {
      return std::nullopt;
    }
  }

private:
  // Try to update indices. Return true if successful, false if not possible.
  bool update() {
    // 1. Find the max i s.t. indices[i] != n - k + i.
    bool found = false;
    unsigned d = n - k;
    unsigned i = k;
    while (i > 0) {
      i--;
      if (indices[i] != d + i) {
        found = true;
        break;
      }
    }
    // 2. If none found, we're finished.
    if (!found) {
      return false;
    }
    // 3. Increment indices[i] and set all later indices to the min possible.
    unsigned v = indices[i] + 1;
    unsigned t = k - i;
    for (unsigned j = 0; j < t; j++) {
      indices[i + j] = v + j;
    }
    return true;
  }
  unsigned n;
  unsigned k;
  bool initial;
  std::vector<unsigned> indices;
};

std::optional<std::vector<std::pair<unsigned, unsigned>>>
evaluate(const std::vector<std::vector<fault_t>> faults, unsigned k, eval_fn g,
         unsigned n_worker_threads) {
  unsigned n = faults.size();
  CombsGenerator cg(n, k);
  // TODO
  while (true) {
    std::optional<std::vector<unsigned>> combo = cg.next();
    if (!combo) {
      break;
    }
    std::cout << "(";
    for (unsigned i = 0; i < k; i++) {
      std::cout << combo->at(i) << ",";
    }
    std::cout << ")" << std::endl;
  }
  return std::nullopt;
}
