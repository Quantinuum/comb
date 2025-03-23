#include "comb.hpp"
#include <cstdint>
#include <iostream>
#include <vector>

int g(uint64_t x) { return x % 3 == 0; }

int main() {
  std::vector<std::vector<uint64_t>> faults{{1, 1, 2}, {3}, {5, 8}};
  std::optional<std::vector<std::pair<unsigned, unsigned>>> bad_faults =
      evaluate(faults, 2, g);
  if (bad_faults) {
    std::cout << "bad fault combination found: [";
    for (const std::pair<unsigned, unsigned> &bucket_and_index : *bad_faults) {
      std::cout << "(" << bucket_and_index.first << ","
                << bucket_and_index.second << "), ";
    }
    std::cout << "]" << std::endl;
  } else {
    std::cout << "no bad fault combinations" << std::endl;
  }
  return 0;
}
