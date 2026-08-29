# CUDA validation runbook

This runbook is the acceptance procedure for the CUDA and MPI backends.  A
successful host-only build is not evidence that the scientific CUDA pipeline
is correct.  The GPU rows below must be run on the target NVIDIA A100 system.

The current implementation is deliberately pre-activation: nonempty CUDA
ranges return `COMB_CUDA_STATUS_NOT_IMPLEMENTED`.  Do not run or publish CUDA
performance numbers until the activation gates in `cuda-port-status.md` are
closed and that guard is removed.

## Record the target

Capture this information beside every baseline:

```sh
date --iso-8601=seconds
git rev-parse HEAD
nvcc --version
nvidia-smi --query-gpu=index,uuid,name,compute_cap,memory.total,driver_version --format=csv
nvidia-smi topo -m
mpirun --version
ompi_info --parsable --all
```

Also record the Linux distribution, CUDA toolkit and driver versions, MPI
implementation, Slurm version, CPU model, node count, A100 count per node,
MIG state, GPU and NIC affinity, and interconnect topology.  Do not compare
runs made with different clocks, MIG layouts, or process binding as if they
were the same baseline.

## Configure and build

Install M4RI under `m4ri/` with the compatibility target, then configure the
primary build:

```sh
make m4ri
cmake -S . -B build -G Ninja \
  -DCMAKE_BUILD_TYPE=RelWithDebInfo \
  -DCOMB_ENABLE_CUDA=ON \
  -DCOMB_ENABLE_MPI=ON
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

The CUDA target must contain native `sm_80` code, forward-compatible
`compute_80` PTX, and line information.  Confirm the built image before a
benchmark:

```sh
cuobjdump --list-elf build/libcomb_cuda.so
cuobjdump --dump-ptx build/libcomb_cuda.so
```

## Correctness matrix

All deterministic configurations must agree on the decision, checked count,
minimum bad global rank, and decoded fault tuple.

| Gate | Required cases |
| --- | --- |
| Rank order | Heterogeneous partitions, `s=1`, maximum supported `s`, empty and singleton rank intervals, overflow |
| Packed widths | 31, 32, 33, 63, 64, 65, 127, and a larger production width |
| Primitive differential | At least 10,000 fixed-seed cases per DCP multiply, phase-aware row add, postselection, elimination, coset extraction, and lookup operation |
| End-to-end exhaustive | Every valid generated instance with `n_all <= 5` |
| Distributed deterministic | CPU, one GPU, and 1/2/4/8 GPUs at more than one MPI world size |
| Fast witness | Every returned witness rechecked by the CPU oracle |
| Production corpus | At least one supplied known-good and one supplied known-bad problem |

The production-corpus gate remains open until the corpus and expected outcomes
are supplied.  Tests must not substitute synthetic outcomes for that evidence.

## Compute Sanitizer

Run a small but stage-complete workload under every tool.  Zero findings are
required:

```sh
compute-sanitizer --tool memcheck  --error-exitcode 1 <test-command>
compute-sanitizer --tool racecheck --error-exitcode 1 <test-command>
compute-sanitizer --tool initcheck --error-exitcode 1 <test-command>
compute-sanitizer --tool synccheck --error-exitcode 1 <test-command>
```

Repeat multi-GPU tests with one MPI rank per GPU.  Sanitizer invocations should
use small ranges so diagnostic overhead does not obscure hangs.

## Profiling and benchmarks

Establish the CPU and one-A100 baselines before changing kernel layout:

```sh
nsys profile --trace=cuda,nvtx,osrt --stats=true -o comb-nsys <test-command>
ncu --set full --target-processes all -o comb-ncu <test-command>
```

Record combinations per second, end-to-end speedup, setup/H2D/reduction time,
stage hit rates, per-rank checked counts and time, scratch bytes, occupancy,
register spills, and cache metrics.  Add a numeric optimization gate to
`cuda-benchmarks.md` only after repeatable CPU and one-A100 baselines exist.

## Distributed invariants

- Map the node-local MPI rank to one CUDA device and reject duplicate device
  assignments.
- Use identical global ranks and epoch boundaries at every world size.
- Complete an epoch before deterministic `MPI_MIN` witness reduction.
- Checkpoint only the next wholly unprocessed epoch, together with the flat
  problem hash.
- Do not put NCCL, peer memory, or CUDA-aware MPI in the candidate hot path.
