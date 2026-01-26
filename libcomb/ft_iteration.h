// Fault tolerant verification API

#include <stdint.h>

#define ALL_FAULTS_GOOD 0
#define BAD_FAULT_FOUND 1
#define PROBLEM_TOO_LARGE 2
#define INTERNAL_ERROR 3
#define OTHER_ERROR 4

// 64 bits ought to be enough for anybody...


typedef uint64_t fault_t;

/**
   A function that takes an output fault - the result of propagating a
   combination of faults through to the end of the circuit and returns 0 if and
   only if the output fault is good (low weight, not detected etc.)
*/
typedef int (*ft_iteration_evaluator)(fault_t output_fault);

/**
   A function that is called with the number of fault combinations that have
   occurred since last update. ft_iteration_iterate may return before the
   updater is called with all information. Called from within the main thread.
*/
typedef void (*ft_iteration_updater)(int fault_combination_count);

/**
   Iterates over all combinations of exactly `k' faults, such that no two of
   them belong to the same fault location, computes their sum and verifies if
   they are valid according to some user-supplied function. `s' names the number
   of fault locations. `m' is an array such that m[i] denotes the number of
   faults at the i-th location. `m' must have enough room for s elements. `k'
   denotes the number of faults that are combined to make a given fault
   combination. `propagated_faults' is an array such that propagated_fault[i] is
   the Pauli that results at the end of the circuit from propagating the i-th
   fault, where faults are labelled lexicographically according to their
   location in `m`. `propagated_faults' must have enough room for m[0] + ... +
   m[s-1] elements. `bad_fault_array' stores the indices of the faults that make
   up a bad combination if one is found. `bad_fault_array' must have enough room
   for k elements. `evaluator' is the function that is called to evaluate each
   output fault combination. `updater' is the function that is called
   periodically to notify the program that a certain number of fault
   combinations have been evaluated. It may be NULL if none is needed.
   `update_stride' indicates the minimum number of faults that must be processed
   between successive invocations of `updater'. Ignored if `updater' is NULL.
*/
int ft_iteration_iterate(int s, int *m, int k, fault_t *propagated_faults,
                         fault_t *bad_fault_array,
                         ft_iteration_evaluator evaluator,
                         ft_iteration_updater updater, int update_stride);
