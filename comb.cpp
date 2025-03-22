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

/// Class to generate all (i0, i1, ..., i_{k-1}) where 0<=i_j<s_j.
class ProdGenerator {
public:
  ProdGenerator(const std::vector<unsigned> &sizes)
      : k(sizes.size()), sizes(sizes), initial(true), indices(k) {
    for (unsigned i = 0; i < k; i++) {
      indices[i] = 0;
    }
  }
  // Return the next vector of indices, or nullopt if we're finished.
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
    // 1. Find the max i s.t. indices[i] != sizes[i] - 1
    bool found = false;
    unsigned i = k;
    while (i > 0) {
      i--;
      if (indices[i] + 1 != sizes[i]) {
        found = true;
        break;
      }
    }
    // 2. If none found, we're finished.
    if (!found) {
      return false;
    }
    // 3. Increment indices[i] and set all later indices to 0.
    indices[i]++;
    for (unsigned j = i + 1; j < k; j++) {
      indices[j] = 0;
    }
    return true;
  }
  unsigned k;
  std::vector<unsigned> sizes;
  bool initial;
  std::vector<unsigned> indices;
};

/// Search for bad combinations of faults from the baskets indicated by their
/// indices in `faults`.
/// If a bad combination is found it is written to result (otherwise nullopt
/// is).
void search_for_bad_combs(const std::vector<std::vector<fault_t>> &faults,
                          const std::vector<unsigned> &baskets,
                          const eval_fn &g,
                          std::optional<std::vector<unsigned>> *result) {
  unsigned k = baskets.size();
  std::vector<unsigned> sizes(k);
  for (unsigned i = 0; i < k; i++) {
    sizes[i] = faults[baskets[i]].size();
  }
  ProdGenerator pg(sizes);
  *result = std::nullopt;
  while (true) {
    std::optional<std::vector<unsigned>> indices = pg.next();
    if (indices) {
      fault_t x = 0;
      for (unsigned i = 0; i < k; i++) {
        x ^= faults[baskets[i]][indices->at(i)];
      }
      if (g(x)) {
        *result = *indices;
        return;
      }
    } else {
      return;
    }
  }
}

std::optional<std::vector<std::pair<unsigned, unsigned>>>
evaluate(const std::vector<std::vector<fault_t>> &faults, unsigned k,
         const eval_fn &g, unsigned n_worker_threads) {
  // TODO Check all vectors in `faults` are non-empty
  unsigned n = faults.size();
  CombsGenerator cg(n, k);
  std::vector<std::optional<std::vector<unsigned>>> results(n_worker_threads);

  // TODO
  return std::nullopt;
}
