#include "ft_iteration.h"
#include <stdio.h>

int g(fault_t x) { return x % 3 == 0; }

int main() {
  int s = 3;
  int m[3] = {3, 1, 2};
  int k = 2;
  fault_t propagated_faults[6] = {1, 1, 2, 3, 5, 8};
  int bad_fault_array[2];
  int x = ft_iteration_iterate(s, m, k, propagated_faults, bad_fault_array, g,
                               0, 5);
  switch (x) {
  case ALL_FAULTS_GOOD:
    printf("no bad fault combinations\n");
    break;
  case BAD_FAULT_FOUND:
    printf("bad fault combination found: [%d, %d]\n", bad_fault_array[0],
           bad_fault_array[1]);
    break;
  default:
    printf("ft_iteration_iterate returned %d\n", x);
  }
  return 0;
}
