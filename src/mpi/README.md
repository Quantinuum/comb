# Comb CUDA MPI orchestration

This directory contains the MPI layer for the candidate-parallel CUDA
verifier.  It deliberately does not add MPI calls to the CUDA library.  Rank
zero packs the non-owning `comb_problem_t` once, broadcasts the versioned flat
blob in chunks no larger than `INT_MAX`, and every rank validates that blob and
creates its own `comb_cuda_context_t`.

The `comb-cuda-mpi` executable reads the same native-endian `.run` files as
[the MPI branch](https://github.com/Quantinuum/comb/tree/8dd5ce7).
Rank zero reads the file and broadcasts the immutable problem; other ranks do
not need access to the input path. The decoder checks dimensions, truncation,
and trailing bytes before any GPU work. It preserves the historic external
bit representation (eight logical bits per `uint64_t`) before the flat packer
converts it to dense device words. Native-endian `.run` files must be generated
on a host with the same byte order.

```bash
cmake -S . -B build-cuda -DCOMB_BUILD_CPU=OFF \
  -DCOMB_ENABLE_CUDA=ON -DCOMB_ENABLE_MPI=ON \
  -DCMAKE_CUDA_ARCHITECTURES=80
cmake --build build-cuda -j
mpiexec -n 3 build-cuda/comb-cuda-mpi input.run \
  --epoch-size 1048576 --chunk-size 16384
```

Use the architecture of the target GPUs in place of `80` for other NVIDIA
hardware. `--checkpoint PATH --resume` resumes a prior checkpoint;
`--checkpoint PATH` starts a new run and saves completed epochs. Optional
`--fast-any` permits earlier termination, and `--scratch-bytes N` bounds the
per-GPU working storage. Exit codes preserve `comb_code_t`: zero means all
faults are good, one means a bad witness was found, and other values report an
input or runtime failure.

An embedding application may also construct a `comb_problem_t` directly:

```c
MPI_Init(&argc, &argv);
MPI_Comm_rank(MPI_COMM_WORLD, &rank);

/* Application-supplied loading: only rank zero needs a populated problem. */
const comb_problem_t *root_problem = rank == 0 ? &problem : NULL;
comb_cuda_mpi_options_t options = COMB_CUDA_MPI_OPTIONS_INITIALIZER;
comb_cuda_mpi_result_t result = COMB_CUDA_MPI_RESULT_INITIALIZER;
comb_cuda_mpi_status_t status = comb_cuda_mpi_verify_comm(
    MPI_COMM_WORLD, root_problem, &options, &result, NULL);

MPI_Finalize();
```

`comb_cuda_mpi_run` is a convenience wrapper for applications that already
have the root problem available before MPI initialization.  It initializes
MPI with `MPI_THREAD_FUNNELED` only when necessary and finalizes MPI only when
it performed that initialization itself.

## Scheduling and witness semantics

All ranks use rank zero's normalized options and therefore have identical
contiguous epoch boundaries.  Fine chunks within an epoch are assigned
cyclically: chunk `j` belongs to world rank `j % world_size`.

Deterministic mode calls `comb_cuda_verify_range` on every assigned chunk and
requires that each call reports the full chunk length as checked, including a
chunk containing a bad candidate.  Each rank retains its minimum bad rank;
`MPI_Allreduce(MIN)` after the complete epoch selects the globally first bad
rank.  Because every earlier epoch was completely good, that rank is also the
globally earliest witness.  A CUDA backend that short-circuits deterministic
ranges violates this orchestration contract and is rejected as an internal
error.

Fast-any mode processes a configurable number of chunks per rank between
`MPI_Iallreduce` polls.  A pending reduction overlaps the following work group,
so at most one polling group of extra work is expected after a witness is
visible globally.  A final blocking reduction selects one valid rank from the
work actually completed.  The result need not be the globally earliest bad
rank, and must be reproduced independently with the merged `main` CPU
implementation, which is the compatibility oracle.

Both modes sum checked counts and reduce CUDA status across all ranks.  Only
the owner of the selected rank contributes the decoded fault tuple; that tuple
is then broadcast from the owner and returned through rank zero's witness.

## Rank-to-GPU mapping

The driver creates an MPI shared-memory communicator for each node.

- If every process sees the node's full `CUDA_VISIBLE_DEVICES` list, CUDA
  ordinal `local_rank` is selected.
- If Slurm gives every process one isolated visible GPU, each process selects
  CUDA ordinal zero.
- If `CUDA_VISIBLE_DEVICES` is absent, ordinal `local_rank` is selected.

After selecting the runtime ordinal, the layer asks CUDA for its physical PCI
bus identifier.  Those identifiers are gathered within the node and duplicate
assignments are rejected.  This remains correct when Slurm isolates one GPU per
task and every rank consequently sees the same local ordinal and token (`0`).
The caller must leave `options.cuda_options.cuda_device` at `-1`; the MPI layer
owns this mapping.  Context creation is the final check that the selected CUDA
ordinal exists.

## Checkpoint format

When `checkpoint_path` is non-null, rank zero writes a checkpoint after every
completely good epoch.  It writes and `fsync`s a temporary file in the same
directory, atomically renames it over the checkpoint path, and then `fsync`s
the containing directory so the new name is durable.  Other ranks never access
the file.

The portable text format is:

```text
COMB_CUDA_MPI_CHECKPOINT 1
problem_hash 0123456789abcdef
epoch_size 1048576
next_epoch 2097152
```

`problem_hash` is a stable 64-bit FNV-1a-derived hash of the exact flat blob.
It is intended to reject accidental resume against different scientific input,
not to provide cryptographic authentication.  Resume requires an exact version,
blob hash, and epoch-size match; `next_epoch` must be an epoch boundary (or the
total count).  A changed world size is allowed because global ranks, chunks,
and epoch boundaries remain unchanged.

## Initial boundaries

- Immutable problem data is replicated on every GPU.  The hot path exchanges
  only scalar status, counts, minimum ranks, and the selected fault tuple.
- NCCL, peer-to-peer copies, CUDA-aware MPI, and GPUDirect RDMA are not used.
- Static cyclic chunks are the initial policy.  Deterministic redistribution
  belongs in a later measured optimization only if rank imbalance exceeds the
  plan's threshold.
- The flat blob is a native binary representation.  All nodes in one run must
  use a compatible host ABI and byte order, as expected for a homogeneous A100
  Slurm partition.

See `scripts/slurm/comb-cuda-mpi.sbatch` for a one-rank-per-A100 launch
template. For example:

```bash
COMB_MPI_EXECUTABLE="$PWD/build-cuda/comb-cuda-mpi" \
  sbatch scripts/slurm/comb-cuda-mpi.sbatch input.run --checkpoint progress.chk
```

## Relationship to the MPI branch

The source studied is MPI commit `8dd5ce7` (`MPI`). It uses a rank-zero work
coordinator and OpenMP worker threads, requiring `MPI_THREAD_MULTIPLE`.
Workers request CPU chunks and abort the whole MPI job on the first bad
candidate. The CUDA implementation runs disjoint candidate ranges on each
rank's GPU, uses `MPI_THREAD_FUNNELED`, and returns a collective result and a
decoded witness through the normal API. Complete deterministic epochs make
the minimum witness rank independent of process scheduling.

The scientific computation belongs to the CUDA backend. The MPI mock tests
exercise scheduling, count reduction, witness ownership, device placement,
and error propagation; they do not validate GPU arithmetic. Native GPU
validation must compare real CUDA decisions and witnesses with the CPU
reference and run Compute Sanitizer on the target system.

## Native GPU tests

On a host with three allocated GPUs, enable the additional native tests:

```bash
cmake -S . -B build-cuda -DCOMB_BUILD_CPU=OFF \
  -DCOMB_ENABLE_CUDA=ON -DCOMB_ENABLE_MPI=ON \
  -DCOMB_BUILD_CUDA_MPI_TESTS=ON
cmake --build build-cuda --parallel
CUDA_VISIBLE_DEVICES=0,1,2 ctest --test-dir build-cuda --output-on-failure
```

`COMB_CUDA_MPI_TEST_RANKS` defaults to `1;2;3`; set it to `1` for one allocated
GPU. The native tests compare real CUDA decisions, minimum witnesses and
counts with the reference and run the actual `.run` launcher through good,
bad, checkpoint/resume and fast-any cases. These opt-in tests require their
requested GPUs and fail if those devices are unavailable. Host mock tests
remain separately labeled.

The initial native validation on `duvel3` used three A30 GPUs, CUDA 11.8 and
MPICH 5.0.2. All 24 tests passed. See
[the validation record](../../docs/validation/duvel3-2026-11-06/README.md) for
commands, sanitizer evidence and remaining production/performance limits.
