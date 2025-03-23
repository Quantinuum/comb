#include "comb.hpp"
#include <cstdint>
#include <future>
#include <optional>
#include <vector>

using namespace std::chrono_literals;

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
/// All elements of `sizes` must be positive.
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
/// If a bad combination is found it is written to `result` and the function
/// returns `true`. Otherwise the function returns `false`.
static bool
search_for_bad_combs(const std::vector<std::vector<uint64_t>> &faults,
                     const std::vector<unsigned> &baskets, const eval_fn &g,
                     std::vector<std::pair<unsigned, unsigned>> *result) {
  unsigned k = baskets.size();
  std::vector<unsigned> sizes(k);
  for (unsigned i = 0; i < k; i++) {
    unsigned basket_sz = faults[baskets[i]].size();
    if (basket_sz == 0) {
      return false;
    }
    sizes[i] = basket_sz;
  }
  ProdGenerator pg(sizes);
  while (true) {
    std::optional<std::vector<unsigned>> indices = pg.next();
    if (indices) {
      uint64_t x = 0;
      for (unsigned i = 0; i < k; i++) {
        x ^= faults[baskets[i]][indices->at(i)];
      }
      if (g(x)) {
        for (unsigned i = 0; i < k; i++) {
          (*result)[i] = {baskets[i], indices->at(i)};
        }
        return true;
      }
    } else {
      return false;
    }
  }
}

std::optional<std::vector<std::pair<unsigned, unsigned>>>
evaluate(const std::vector<std::vector<uint64_t>> &faults, unsigned k,
         const eval_fn &g, unsigned n_worker_threads) {
  unsigned n = faults.size();
  CombsGenerator cg(n, k);
  std::vector<std::vector<std::pair<unsigned, unsigned>>> results(
      n_worker_threads);
  for (unsigned j = 0; j < n_worker_threads; j++) {
    results[j] = std::vector<std::pair<unsigned, unsigned>>(k);
  }
  std::vector<std::future<bool>> futures(n_worker_threads);
  std::optional<std::vector<unsigned>> indices = cg.next();
  std::vector<bool> busy(n_worker_threads);
  while (true) {
    // Poll all existing threads for completion. For each completed thread: if
    // results were found, return them; otherwise empty the slot.
    for (unsigned j = 0; j < n_worker_threads; j++) {
      if (busy[j] && futures[j].wait_for(0ms) == std::future_status::ready) {
        if (futures[j].get()) {
          return results[j];
        } else {
          busy[j] = false;
        }
      }
    }
    if (indices.has_value()) {
      // If any thread slot is empty, start a new thread in it with indices, and
      // update indices.
      for (unsigned j = 0; j < n_worker_threads; j++) {
        if (!busy[j]) {
          std::vector<std::pair<unsigned, unsigned>> *result =
              &results.data()[j];
          futures[j] =
              std::async(std::launch::async, [faults, indices, g, result] {
                return search_for_bad_combs(faults, *indices, g, result);
              });
          busy[j] = true;
          indices = cg.next();
          break;
        }
      }
    } else {
      // If all thread slots are free, return.
      bool all_free = true;
      for (unsigned j = 0; j < n_worker_threads; j++) {
        if (busy[j]) {
          all_free = false;
          break;
        }
      }
      if (all_free) {
        return std::nullopt;
      }
    }
  }
}
