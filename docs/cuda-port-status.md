# CUDA and MPI implementation status

Updated 2026-11-06. The MPI source studied for this port is
[`mpi` at 8dd5ce73f80994b84b28b759b17459924c77d630](https://github.com/Quantinuum/comb/tree/8dd5ce73f80994b84b28b759b17459924c77d630).
This records implementation and observed validation separately.

## Implemented

- The existing CPU ABI remains available. CUDA and MPI are optional CMake
  targets, and CUDA/MPI can be built without M4RI (`COMB_BUILD_CPU=OFF`).
- `comb-cuda-mpi` reads the MPI branch's existing native-endian `.run` files.
  Rank zero loads once; the checked pointer-free problem is broadcast to all
  ranks and uploaded once per GPU. The decoder rejects malformed dimensions,
  truncated input and trailing bytes.
- The complete scientific candidate pipeline executes inside the CUDA kernel:
  weighted-colex fault unranking, ordered DCP composition, logical variants,
  X-syndrome rejection, CSS X lookup, phase-aware Z/X postselection, tableau
  elimination and universal membership of the resulting Z coset.
- Faults compose as `E1 E2 ... Es`; selected logical DCPs left-multiply in
  ascending index order. Both match the MPI branch's scientific loop.
- Stabiliser pivots, free columns and measurement sequences are prepared on
  the host and uploaded as immutable data. The GPU uses the sequential
  postselection formulation, with the validated dual-measurement relation;
  it does not copy the CPU's optimized batched-Z/PLUQ transform.
- GPU threads evaluate independent candidates, with deterministic pivot
  selection inside each candidate. Workspaces are allocated once per context,
  disjoint between threads and reused. The default scratch budget is the
  smaller of 256 MiB and one quarter of free GPU memory; callers can set it.
  No device allocation or CPU fallback occurs in the candidate loop.
- A persistent kernel claims ranks using a saturating atomic counter.
  Deterministic range calls complete their entire interval and return its
  smallest bad rank. Fast-any calls may stop after any valid bad candidate.
  Runtime errors remain distinct from scientific decisions.
- MPI assigns fine cyclic chunks inside common bounded epochs. Status,
  counts, minimum rank and the selected tuple are reduced. One rank owns one
  physical GPU; duplicate PCI identities within a node are rejected.
  Checkpoint/resume retains input-hash and epoch-boundary validation.
- GPU functions are shared with host differential tests. This checks their
  mathematics independently of CUDA scheduling. Real-device tests use the
  public CUDA API and are built only when CUDA is enabled.

## Validation on this development host

This macOS/Arm64 host provides Clang, M4RI and Open MPI, but has no NVIDIA
GPU, driver or CUDA toolkit.

- CPU-enabled CMake/CTest passes 18/18 tests, including the added CLI
  scientific test. The earlier dependency-free suite passed 16/16; the final
  native build also passes all 17 common host tests. These include rank, dense, flat, stub and reference tests and the
  expanded CPU/reference differential.
- The host candidate differential checks every rank in 116 fixtures, including
  heterogeneous partitions, `s=1,2,3`, noncommuting DCPs, CSS stabilisers,
  logical variants, measured ancillas and 65-qubit multiword inputs.
- Shared GPU math primitives pass deterministic/randomized comparisons with
  the independent dense implementation at widths 1, 63, 64, 65 and 129.
- The expanded M4RI CPU/reference differential includes measured-ancilla
  fixtures, checking the sequential formulation against CPU batched-Z logic.
- Production MPI orchestration passes one-, two- and three-process tests with
  test backends, including the actual candidate functions executing on the
  host against the independent scientific reference. These exercise real MPI
  collectives, but do not execute CUDA.
- Strict host compiler checks and AddressSanitizer/UndefinedBehaviorSanitizer
  runs pass for GPU mathematics, the complete host candidate pipeline and
  run-file parsing. Clang's CUDA frontend checks host/device source and emits
  sm_80 PTX with temporary runtime declarations. That fixture does not provide
  real CUDA runtime headers or linking, and is not an NVIDIA toolkit build.

## Native validation on duvel3

CUDA 11.8.89 and GCC 8.4 compile and link the real CUDA/MPI targets on `duvel3`.
The host has three NVIDIA A30 GPUs (24 GiB each, `sm_80`) and driver
550.163.01. CMake 3.31.10 and MPICH 5.0.2 are installed only in an isolated
validation environment, using uv; no system packages were changed.

- The final CMake/CTest suite passes 24/24 tests, with no skips. It includes
  the 116-fixture device differential and scheduling checks, native scientific
  MPI tests at 1/2/3 ranks, and the actual `.run` CLI at those GPU counts.
- CLI tests check good/bad decisions, minimum witnesses, exact checked counts,
  checkpoint/resume and fast-any behavior against fixed synthetic fixtures.
- `cuobjdump` confirms native `sm_80` and embedded PTX. A separate device run
  with `CUDA_FORCE_PTX_JIT=1` passes the full device differential.
- Compute Sanitizer memcheck, initcheck, racecheck and synccheck pass on
  the final build. Three-rank MPI memcheck reports zero errors and zero leaks
  on every rank. Raw logs are in
  [the native validation record](validation/duvel3-2026-11-06/README.md).
- The GCC build exposed an always-false allocation comparison; the MPI layer
  now checks wire-width counts before narrowing. Removing unnecessary CUDA
  device linking preserves the configured PTX image in the shared library.

## Remaining validation and optimization

No production input corpus or performance baseline was supplied. Native
correctness evidence currently covers synthetic fixtures on one node with
three A30s; it does not establish production outcomes, multi-node behavior,
other GPU architectures or a speedup. Use `cuda-validation.md` for the
remaining acceptance matrix.

This is an initial candidate-parallel implementation. Warp-cooperative
elimination, a batched-Z GPU transform and performance tuning remain future
work, to be justified by actual GPU profiles. The benchmark record stays
empty until measurements are available.
