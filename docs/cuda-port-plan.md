# CUDA and distributed-GPU port of Comb

## Summary

The CUDA port will parallelize over globally ranked fault combinations. Each
combination's complete DCP, postselection, Gaussian-elimination, and coset
pipeline remains on one GPU. Candidate-local elimination is cooperative within
a warp or thread block; an elimination is never split across GPUs. Immutable
problem data is replicated on each GPU, and workers communicate only status,
counts, and bad-witness ranks.

This decomposition preserves the scientific contract in the accompanying
paper: the multiplication order is \(E=E_s\cdots E_1\), tableau row operations
track phase, coherent terms are cancelled before weight tests, the
logical-stabilizer search is existential, and coset checking is universal. The
primary reference is the user-supplied
`CNOT___T_for_magic_state_preparation-10.pdf`; the working copy is intentionally
kept outside this repository.

The unmerged `multithreaded` and `chunked` branches are research references
only. Their weighted partitioning, termination, MPI, and witness logic are not
correct enough to merge.

## Interfaces and build

- Preserve `comb_verify_css_magic_state` and the existing CPU ABI unchanged.
  In particular, do not repurpose `thread_count`.
- Add an opaque CUDA context API:
  - `comb_cuda_context_create(problem, options, &context)` validates,
    preprocesses, packs, and uploads one problem.
  - `comb_cuda_verify_range(context, rank_begin, rank_end, mode, result,
    witness)` checks an explicit rank interval.
  - `comb_cuda_verify(...)` checks the complete range.
  - `comb_cuda_context_destroy(...)` releases all host and device resources.
- Group existing inputs into a versioned, non-owning `comb_problem_t`. Options
  select the CUDA device, epoch size, scratch limit, and deterministic or
  fast-any-witness mode.
- Return the verification code, total and checked combination counts, minimum
  bad rank, decoded fault tuple, and setup/verification timings.
- Keep the legacy external bit representation intact. Normalize it into an
  internal, pointer-free `FlatProblem` containing dense, tail-masked `uint64_t`
  rows. Never reinterpret M4RI objects or the legacy eight-bits-per-word macros
  on the GPU.
- Define one canonical global order: location subsets follow the current colex
  ordering; fault offsets within a subset use mixed-radix order with the first
  tuple element varying fastest. Implement checked weighted rank/unrank using
  dynamic programming over partition sizes. Return `PROBLEM_TOO_LARGE` if
  counts exceed `UINT64_MAX` or lookup tables cannot be addressed safely.
- Use CMake/CTest as the primary build, with optional CUDA and MPI targets, and
  retain Makefile compatibility wrappers. Build A100-native `sm_80` code plus
  `compute_80` PTX and include `-lineinfo`. See NVIDIA's
  [GPU compute-capability table](https://developer.nvidia.com/cuda/gpus) and
  [CUDA compilation documentation](https://docs.nvidia.com/cuda/cuda-programming-guide/02-basics/nvcc.html).

Generated Doxygen output belongs under `build/docs`; the tracked `docs/`
directory is reserved for source documentation such as this plan.

## Implementation phases

### 1. Correctness contract and CPU oracle

- Add executable tests with assertions; the current diagnostic test programs
  are not a usable oracle.
- Retain a test-only sequential postselection path and prove it equivalent to
  `main`'s batched-Z transform before porting that transform.
- Lock DCP multiplication order, phase-aware row-add semantics, pivot
  selection, coset factor construction, lookup addressing, unsupported-input
  behavior, and witness ordering.
- Add comprehensive dimension, shift, rank, partition, table-length, and
  overflow validation.
- Import the user-supplied production corpus with at least one known-good and
  one known-bad problem.

### 2. Correct single-A100 implementation

- Keep validation and candidate-independent M4RI preprocessing on the CPU,
  performed once per context.
- Store output errors, stabilizers, measurements, lookup tables, PLUQ
  permutation, and row-operation schedules as packed structure-of-arrays device
  data.
- Implement the general multiword kernel first: one cooperative block per
  candidate, global reusable scratch, and a shared-memory pivot-row tile.
- Execute the complete candidate pipeline in order: weighted unranking, DCP
  composition, logical representatives, cheap detection/lookups, batched-Z
  transform, X postselection, coset reduction, and coset lookup scan.
- Use custom XOR, AND, parity, and `__popcll` primitives. cuBLAS is not suitable
  for GF(2) arithmetic with nonlinear Pauli phase tracking.

### 3. Optimized single-A100 path

- Add an `n_all <= 64` fast path using one warp per candidate and multiple
  warps per block when shared-memory occupancy permits.
- Use lowest-index deterministic pivots, parallel destination-row XORs, and
  cooperative parity reductions.
- Parallelize coset elements across lanes instead of preserving the CPU's
  serial Gray-code dependency.
- Process bounded contiguous epochs through a persistent device work counter so
  scratch scales with resident teams, not total combinations.
- Add NVTX stage instrumentation. Consider BMMA, lookup caching, fusion, or
  further preprocessing only when Nsight profiles identify the corresponding
  bottleneck. Follow the
  [CUDA Best Practices Guide](https://docs.nvidia.com/cuda/cuda-c-best-practices-guide/index.html).

### 4. Multiple A100s on one node

- Use one MPI rank per GPU under Slurm. Map local rank to one visible device and
  verify that assignments are unique.
- Replicate `FlatProblem` on every GPU. Divide each contiguous epoch into fine
  cyclic chunks to reduce variation between cheap detected candidates and
  expensive coset candidates.
- In deterministic mode, complete the epoch and use host-side
  `MPI_Allreduce(MIN)` to select its lowest bad rank. Earlier epochs are already
  complete, so this is the globally earliest witness.
- In fast mode, propagate a valid bad-witness flag through periodic nonblocking
  reductions and stop as soon as practical.
- Do not put peer-to-peer transfers, NVLink, or NCCL on the hot path; the
  algorithm has no shared scientific state between candidates. This follows
  NVIDIA's [multi-GPU programming model](https://docs.nvidia.com/cuda/cuda-programming-guide/03-advanced/multi-gpu-systems.html).

### 5. Tested multi-node prototype

- Keep the core CUDA library MPI-independent. A separate Slurm/MPI driver owns
  `MPI_Init` and `MPI_Finalize`.
- Rank zero validates and broadcasts one versioned flat-problem blob; every
  rank uploads its local copy.
- Retain identical global ranks and epoch boundaries at every world size.
- Reduce only minimum bad rank, checked count, and rank-wide error status.
- Add checkpoint/resume using the problem hash and next unprocessed epoch.
- Start with static cyclic chunks. Add deterministic work redistribution only
  if measured rank imbalance exceeds 10%.
- Keep CUDA communication and MPI control phases separate. NVIDIA documents
  progress and deadlock hazards when mixing NCCL and CUDA-aware MPI in its
  [NCCL/MPI guidance](https://docs.nvidia.com/deeplearning/nccl/user-guide/docs/mpi.html),
  so NCCL and GPUDirect RDMA are excluded from the initial prototype.

## Verification and performance gates

- Exhaustively compare CPU and GPU for valid instances with `n_all <= 5`.
- Run at least 10,000 seeded randomized cases for DCP multiplication,
  phase-aware row addition, postselection, Gaussian elimination, coset
  extraction, and lookup behavior.
- Cover heterogeneous fault partitions, empty and singleton ranges, `s=1`, the
  maximum supported `s`, count overflow, and exact rank/unrank agreement with
  the serial iterator.
- Exercise widths 31, 32, 33, 63, 64, 65, and larger multiword tableaus; test
  tail-bit masking and lookup-index limits.
- Require exact decision and lowest-witness agreement between CPU, one GPU,
  1/2/4/8 GPUs, and different MPI world sizes. Fast mode may return a different
  witness, but it must independently reproduce as bad on the CPU.
- Run Compute Sanitizer `memcheck`, `racecheck`, `initcheck`, and `synccheck`
  with zero findings; see NVIDIA's
  [Compute Sanitizer documentation](https://docs.nvidia.com/compute-sanitizer/ComputeSanitizer/index.html).
- Benchmark:
  - The paper's 824-fault and 5,619,552/12,662,112 pair workloads.
  - The supplied production corpus.
  - Seeded distance-seven `s=3` stress cases with varied detection and coset
    ranks.
  - CPU, one A100, 2/4/8 A100s, then 1/2/4 nodes.
- Record combinations/s, end-to-end speedup, setup/H2D/reduction fractions,
  stage hit rates, load imbalance, scratch usage, occupancy, register spills,
  cache behavior, and scaling efficiency.
- Establish CPU and one-A100 baselines before optimization, then freeze numeric
  promotion gates in a benchmark document. Accept an optimization only if it
  improves median large-workload throughput without any correctness or
  small-workload regression.

## Assumptions and boundaries

- The user will provide representative production files and expected outcomes
  before the correctness phase closes.
- The target is Linux, Slurm, MPI, and A100 GPUs. Record exact CUDA, MPI, GPU
  topology, and network details from the target node before benchmarking.
- General multiword tableaus are required, with a specialized one-word fast
  path. Physically unaddressable lookup tables remain rejected explicitly.
- CSS-code verification remains in scope. Fault propagation, lookup-table
  generation, non-CSS generalization, and changes to the paper's scientific
  criterion are out of scope.
- The CPU implementation remains the supported fallback and differential
  oracle.
- The live NVIDIA skill catalog contained no closer skill specifically for this
  classical CUDA-C verifier; the architecture therefore follows installed CUDA
  guidance and current primary NVIDIA documentation.
