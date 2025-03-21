#include <stdint.h>

#define ALL_FAULTS_GOOD 0
#define BAD_FAULT_FOUND 1
#define PROBLEM_TOO_LARGE 2
#define INTERNAL_ERROR 3
#define OTHER_ERROR 4

typedef uint64_t fault_t;

typedef int (*ft_iteration_evaluator)(fault_t output_fault);

typedef void (*ft_iteration_updater)(int fault_combination_count);

int ft_iteration_iterate(int s, int *m, int k, fault_t *propagated_faults,
                         int *bad_fault_array, ft_iteration_evaluator evaluator,
                         ft_iteration_updater updater, int update_stride);
