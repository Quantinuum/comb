# CUDA validation runbook

This runbook is the acceptance procedure for the CUDA and MPI backends.  A
successful host-only build is not evidence that the scientific CUDA pipeline
is correct.  The GPU rows below require a target NVIDIA system. The initial A30
results are recorded in `validation/duvel3-2026-11-06/README.md`.

Nonempty ranges launch the complete candidate kernel. Native CUDA compilation
and single-/multi-GPU synthetic correctness checks pass on three A30s. The
checks below describe how to reproduce that evidence and extend it to a
production corpus or another system before reporting performance. The
implementation uses one thread per candidate; cooperative row elimination and the CPU batched-Z transform are not implemented on GPU.

## Record the target

Capture this information beside every baseline:

```sh
date --iso-8601=seconds
git rev-parse HEAD
nvcc --version
nvidia-smi --query-gpu=index,uuid,name,compute_cap,memory.total,driver_version --format=csv
nvidia-smi topo -m
mpirun --version
ompi_info --parsable --all  # Open MPI only
```

Also record the Linux distribution, CUDA toolkit and driver versions, MPI
implementation, Slurm version, CPU model, node count, GPU count per node,
MIG state, GPU and NIC affinity, and interconnect topology.  Do not compare
runs made with different clocks, MIG layouts, or process binding as if they
were the same baseline.

## Configure and build

For the M4RI CPU comparison, install M4RI under `m4ri/` with the compatibility
target, then configure the primary build. For a CUDA/MPI-only build, omit
`make m4ri` and add `-DCOMB_BUILD_CPU=OFF`:

```sh
make m4ri
cmake -S . -B build -G Ninja \
  -DCMAKE_BUILD_TYPE=RelWithDebInfo \
  -DCOMB_ENABLE_CUDA=ON \
  -DCOMB_ENABLE_MPI=ON \
  -DCOMB_BUILD_CUDA_MPI_TESTS=ON
cmake --build build --parallel
ctest --test-dir build --output-on-failure
ctest --test-dir build -R comb.cuda_device --output-on-failure
mpiexec -n 3 build/comb-cuda-mpi path/to/comb-example.run
```

The default Ampere CUDA target contains native `sm_80` code, forward-compatible
`compute_80` PTX, and line information. When changing an existing build from
separable to whole-program CUDA compilation, rebuild with `--clean-first`
to discard old relocatable device objects. Confirm the built image before a
benchmark:

```sh
cuobjdump --list-elf build/libcomb_cuda.so
cuobjdump --dump-ptx build/libcomb_cuda.so
```

## Correctness matrix

All configurations must agree with the merged `main` CPU implementation on
the scientific decision and on each candidate's classification.  That
implementation is the compatibility oracle for the algorithm.  The main host
loop does not define a globally minimum witness when a problem spans multiple
chunks: later bad chunks can overwrite `bad_fault_array`, even with one OpenMP
thread.  Deterministic CUDA range calls must agree with the sequential reference on
checked count, minimum bad global rank and decoded fault tuple. MPI completes
the first bad epoch and stops, so its minimum witness must match the
reference while its checked count covers only the completed epoch prefix.
For a bad rank r, that prefix ends at min(total, (floor(r / epoch_size) + 1)
* epoch_size); an all-good run reaches total. Subtract the resumed prefix
when checking newly completed work. Direct witness equality with main is
required only for explicit one-chunk fixtures.

| Gate | Required cases |
| --- | --- |
| Output-error fold | Noncommuting `s=2` and `s>2` cases that distinguish `E1 E2 ... Es` from the reversed fold |
| Rank order | Heterogeneous partitions, `s=1`, maximum supported `s`, empty and singleton rank intervals, overflow |
| Packed widths | 31, 32, 33, 63, 64, 65, 127, and a larger production width |
| Primitive differential | At least 10,000 fixed-seed cases per DCP multiply, phase-aware row add, postselection, elimination, coset extraction, and lookup operation |
| End-to-end exhaustive | Every valid generated instance with `n_all <= 5` |
| Distributed deterministic | Main decision agreement plus sequential-reference counts and lowest witnesses across available GPU counts (1/2/3 tested on duvel3); extend to larger allocations and multiple nodes before making scaling claims |
| Fast witness | Every returned witness rechecked with the main candidate logic and the sequential reference |
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
