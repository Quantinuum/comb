#pragma once

#include <cstdint>
#include <optional>
#include <vector>

typedef int (*eval_fn)(uint64_t);

/// Look for bad combinations of `k` faults from different buckets.
/// `faults` is the vector of buckets. The combination x_0, ..., x_{k-1} of
/// faults is considered bad iff g(x_0 ^ ... ^ x_{k-1}) != 0.
/// If a bad combination is found, retuen it as a vector of pairs where each
/// pair consists of the bucket index and the index within the bucket.
std::optional<std::vector<std::pair<unsigned, unsigned>>>
evaluate(const std::vector<std::vector<uint64_t>> &faults, unsigned k,
         const eval_fn &g, unsigned n_worker_threads = 8);
