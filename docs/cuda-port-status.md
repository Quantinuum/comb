# CUDA port implementation status

This file distinguishes implemented infrastructure from scientific behavior
that has passed the plan's acceptance gates.  It is not a substitute for the
design in `cuda-port-plan.md` or the acceptance procedure in
`cuda-validation.md`.

## Implemented

- The legacy `comb_verify_css_magic_state` signature and `thread_count`
  behavior are unchanged.
- The additive, versioned CUDA context/range/result/witness API is defined.
- Structural and scientific input validation plus a pointer-free, tail-masked
  `uint64_t` flat problem blob are implemented.  The contract checks CSS
  dimensions and commutation, stabiliser RREF, symmetric DCP matrices, the
  supported batched-Z dual relation, checked weighted-colex counts, and
  rank/unrank tables before upload.
- Dependency-free dense GF(2), phase-aware tableau, DCP, postselection, and
  coset primitives have executable deterministic and randomized tests.
- A complete sequential reference verifier executes the pre-batched
  measurement contract and returns the lowest bad global rank.  A separate
  M4RI-gated test compares it with the legacy batched verifier.
- The legacy fault loop now composes later faults on the left, locking the
  paper's `E_s ... E_1` order.  A dependency-free two-location fixture uses
  noncommuting DCPs whose lookup decision changes if that order is reversed.
  Residual coset support outside the code lookup domain is rejected instead of
  reproducing the legacy missing-measurement undefined behavior.
- CUDA context creation uploads the immutable blob.  Device code contains
  weighted unranking, packed XOR/parity/phase primitives, atomic-min witness
  support, a race-free persistent scheduling scaffold, and distinct one-word
  and multiword dispatch points.
- The MPI layer broadcasts one validated flat blob, maps one rank to one
  physical GPU, schedules fine cyclic chunks inside invariant epochs, reduces
  deterministic or fast witnesses, and supports hash-checked checkpoint and
  resume.  The core CUDA library remains MPI-independent.
- A single-rank mock CTest drives the production MPI scheduler through
  deterministic and fast witnesses, cyclic chunking, checkpoint/resume,
  caller-witness reuse, and the physical-device identity hook.  It does not
  stand in for multi-rank validation.
- CMake/CTest is the primary optional CPU/CUDA/MPI build.  The Makefile remains
  a compatibility wrapper, generated Doxygen output is under `build/docs`, and
  `make clean` preserves this tracked `docs` directory.

## Deliberately not activated

For every nonempty interval, `comb_cuda_verify_range` currently returns
`COMB_CUDA_STATUS_NOT_IMPLEMENTED`.  It does not launch the scheduler scaffold
and cannot return a scientific decision.  Empty intervals return the valid
vacuous `ALL_FAULTS_GOOD` result.

Activation requires all of the following, in order:

1. Execute the M4RI legacy/sequential differential CTest and investigate every
   mismatch, especially DCP multiplication order and batched-Z equivalence.
2. Import at least one user-supplied production known-good and one known-bad
   problem with expected outcomes.
3. Extend the flat preprocessing contract with the validated pivot,
   permutation, batched measurement-transform, and row-operation data required
   by the device pipeline.
4. Implement the complete per-candidate CUDA path and obtain exact CPU/GPU
   decision and lowest-witness agreement before enabling nonempty ranges.
5. Pass the A100 Compute Sanitizer, multi-GPU/world-size, and multi-node gates
   in `cuda-validation.md`.

## Validation completed on this host

- Release CMake/CTest with the M4RI-dependent legacy target disabled.
- Strict C11 warnings-as-errors builds for rank, dense primitives, flat
  packing, the CUDA-unavailable stub, the sequential verifier, and the MPI
  orchestration translation unit.
- AddressSanitizer and UndefinedBehaviorSanitizer runs for all six executable
  host-side contract suites available when MPI C is installed.
- Host-side CUDA syntax parsing with a mock runtime header.

## Open environment and input gates

- No production corpus or expected production outcomes were available.
- This macOS/Arm64 host has no CUDA compiler, NVIDIA driver, or A100, so no
  `nvcc`, device, profiler, or Compute Sanitizer run is claimed.
- The pinned M4RI bootstrap could not run because `autoreconf` is absent, so
  the M4RI-gated legacy/sequential differential test has not executed here.
- The installed Open MPI and PMIx components on this host cannot launch a
  reliable multi-rank job.  The MPI source was compiled strictly, but runtime
  collectives still require validation on the target Slurm system.
