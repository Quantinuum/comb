# Native CUDA/MPI validation on duvel3

Date: 2026-11-06 (JST). Validation used the working tree on `codex/cuda`, based
on commit `4e8ff181a1b024d4a5821dee187b20d6634ed509`. It includes uncommitted
changes. `source-sha256.txt` identifies the source actually compiled;
`final-build-metadata.txt` identifies the resulting binaries. All 55 captured
C/CUDA/header/CMake file hashes match the local working tree. Documentation
was updated after the source snapshot, so its hashes may differ.

## System and isolated workspace

- Host: `duvel3`, Linux x86_64, three NVIDIA A30 GPUs, 24 GiB each, SM80.
- NVIDIA driver: 550.163.01; CUDA toolkit: 11.8.89; GCC: 8.4.0.
- CMake: 3.31.10; MPICH: 5.0.2; tooling installed using uv in a private
  environment. No system packages or existing workloads were modified.
- Remote workspace:
  `/home/erinaldi/codex-validation/comb-cuda-20261106-lsTaEl`.
- Tool prefix: `tooling/env/bin` inside that workspace. The built executable
  remains at `build-gpu/comb-cuda-mpi`.
- Tests used `CUDA_VISIBLE_DEVICES=1,2,0` to place the first ranks on the two
  initially idle GPUs. Three-rank runs also used GPU 0, with small workspaces.

See `environment.txt`, `final-build-metadata.txt` and
`tooling-requirements.txt` for captured versions and configuration.

## Results

| Check | Result | Evidence |
| --- | --- | --- |
| Native CMake build with nvcc | Passed | `cmake-configure.log`, `cmake-build.log` |
| Full CTest suite | 24/24 passed, no skips | `ctest-native.log`, `ctest-details.log` |
| Device scientific differential | Every rank in 116 synthetic fixtures; scheduling and scratch-limit cases passed | `ctest-details.log` |
| Native scientific MPI | 1, 2 and 3 GPUs; decisions/witnesses match the reference, counts match the completed-epoch expectation | `ctest-details.log` |
| Native `.run` CLI | 1, 2 and 3 GPUs; good/bad, checkpoint/resume and fast-any passed | `ctest-details.log` |
| Native code and embedded PTX | Both present in the final shared library | `binary-inspection.txt` |
| Forced PTX JIT | Complete device differential passed | `ptx-jit.log` |
| Compute Sanitizer memcheck | Zero errors | `cuda-memcheck.log` |
| Compute Sanitizer initcheck | Zero errors | `cuda-initcheck.log` |
| Compute Sanitizer racecheck | Zero hazards/errors/warnings | `cuda-racecheck.log` |
| Compute Sanitizer synccheck | Zero errors | `cuda-synccheck.log` |
| Three-GPU MPI memcheck | Zero errors and zero leaks on all three ranks | `cuda-mpi-memcheck.log` |

The shared GPU mathematical functions also pass host differential and sanitizer
tests. Separately, the M4RI CPU/reference differential passes 24 synthetic
fixtures on the Mac, including measured ancillas. M4RI is not required for the
native CUDA/MPI-only build below.

The native build exposed a GCC warning in an allocation guard. The check now
keeps counts at their wire width before narrowing. Binary inspection exposed
missing embedded PTX after unnecessary CUDA device linking. The single CUDA
translation unit now uses whole-program compilation; the clean rebuild,
complete native suite and forced PTX run all pass.

## Reproduce on duvel3

From the recorded remote workspace:

```bash
export PATH="$PWD/tooling/env/bin:/usr/local/cuda-11.8/bin:$PATH"
cmake -S . -B build-gpu \
  -DCMAKE_BUILD_TYPE=RelWithDebInfo \
  -DCOMB_BUILD_CPU=OFF \
  -DCOMB_ENABLE_CUDA=ON -DCOMB_ENABLE_MPI=ON \
  -DCOMB_BUILD_CUDA_MPI_TESTS=ON \
  -DCOMB_CUDA_MPI_TEST_RANKS='1;2;3' \
  -DCMAKE_CUDA_COMPILER=/usr/local/cuda-11.8/bin/nvcc \
  -DCUDAToolkit_ROOT=/usr/local/cuda-11.8 \
  -DMPI_C_COMPILER="$PWD/tooling/env/bin/mpicc"
cmake --build build-gpu --clean-first --parallel 4
CUDA_VISIBLE_DEVICES=1,2,0 ctest --test-dir build-gpu --output-on-failure -j1
CUDA_VISIBLE_DEVICES=1 CUDA_FORCE_PTX_JIT=1 ./build-gpu/comb_cuda_device_test
```

The default architecture list is `80-real;80-virtual`. A clean rebuild is
needed when switching an existing build between relocatable-device and
whole-program compilation, so that stale device objects are discarded.

```bash
for sanitizer in memcheck initcheck racecheck synccheck; do
  CUDA_VISIBLE_DEVICES=1 compute-sanitizer --tool "$sanitizer" \
    --error-exitcode 99 ./build-gpu/comb_cuda_device_test < /dev/null
done
CUDA_VISIBLE_DEVICES=1,2,0 mpiexec -n 3 compute-sanitizer \
  --tool memcheck --leak-check full --error-exitcode 99 \
  ./build-gpu/comb_cuda_mpi_device_test < /dev/null
```

To run an existing production input:

```bash
CUDA_VISIBLE_DEVICES=1,2,0 mpiexec -n 3 ./build-gpu/comb-cuda-mpi input.run \
  --checkpoint progress.checkpoint
```

Exit status 1 means a bad witness was found; status 0 means all checked faults
were good. Add `--resume` to continue from a matching checkpoint.

## Limits

These are live native results on one node with three A30s. The fixtures are
synthetic; no user-supplied production corpus was available. They do not
establish production correctness, multi-node behavior, behavior on another
GPU architecture or a performance advantage. No benchmark or speedup is
reported. Racecheck mostly targets shared-memory hazards; this kernel uses
private candidate workspaces and global atomics, so its clean result alone
does not establish global scheduling correctness. Exact distributed outcomes
are checked separately by the native MPI tests.
