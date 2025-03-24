# comb

A parallelized implementation of `ft_iteration_iterate()`.

See https://github.com/CQCL/ft-iterator for more details. The signature of
`ft_iteration_iterate()` defined here differs from the one there only in that:

* the `bad_fault_array` argument is an `int*` instead of a `fault_t*` (this is
  more correct as it stores indices not faults);
* there is an extra argument `n_worker_threads` specifying the maximum number of
  worker threads used by the library.

Note also that the `updater` and `update_stride` arguments are ignored.
