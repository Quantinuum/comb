#include <cstdint>
#include <optional>
#include <vector>

typedef uint64_t fault_t;

typedef int (*eval_fn)(fault_t);

std::optional<std::vector<std::pair<unsigned, unsigned>>>
evaluate(const std::vector<std::vector<fault_t>> faults, unsigned k, eval_fn g,
         unsigned n_worker_threads = 8);
