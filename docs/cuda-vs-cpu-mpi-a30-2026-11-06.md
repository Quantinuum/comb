# CUDA versus native CPU MPI on two A30 GPUs

Live measurements on `duvel3`, 2026-11-06. For the 5,781,600-pair workload,
two A30s were **3.05× faster than two CPU workers**, while the 31-worker
CPU configuration was **2.11× faster than the GPUs**. Performance depends
strongly on workload and scheduling.

Compared the current `codex/cuda` working tree with the **unmodified native
`mpi` branch at `8dd5ce7`**. Hardware: two NVIDIA A30 24 GiB GPUs (physical
GPUs **1 and 2 only**) and two Xeon Gold 6326 CPUs, 32 physical cores total.
Both builds use `-O3`; host code uses GCC 8.4 with `-march=native`, CUDA uses
nvcc 11.8/SM80, and both use MPICH 5.0.2. The CPU baseline uses thread-safe
M4RI 20240729.

Times below are **median complete-process wall times from three repetitions**
after one warmup per configuration. Profile order rotated between repetitions.
They include MPI startup, loading, preprocessing, CUDA initialization/transfers,
verification and teardown. All 48 primary runs returned `ALL_FAULTS_GOOD` and
checked the exact full combination count on the same unchanged verification inputs.

| Workload | Qubits incl. ancillas | Combinations | CPU: 3 MPI ranks | CPU: 32 MPI ranks | CUDA: 2 A30s |
| --- | ---: | ---: | ---: | ---: | ---: |
| Steane, single faults | 15 | 944 | 0.162 s | 0.640 s | 0.774 s |
| Colour code, single faults | 58 | 2,712 | 0.837 s | 1.710 s | 0.851 s |
| Golay, single faults | 143 | 7,420 | 11.981 s | 14.897 s | 2.149 s |
| Distance-5 optimised, fault pairs | 62 | 5,781,600 | 34.534 s | 5.357 s | 11.321 s |

CPU configurations use one coordinator plus 2 or 31 compute workers, with
`OMP_NUM_THREADS=1` and disjoint physical-core binding. CUDA uses two MPI ranks
with default epoch/chunk sizes of 1,048,576/16,384. **Each single-fault case
fits in one CPU chunk and one CUDA chunk**, so only one worker/GPU executes
faults despite the larger allocations. The pair case used all 31 CPU workers
and both GPUs: GPU throughput was 510,718 combinations/s versus
1,079,259 for the whole-host CPU job. Its measured ranges were
11.301–11.339 s (GPU), 34.182–34.612 s (2 CPU workers), and 5.284–5.420 s
(31 CPU workers).

A separate, later [GPU follow-up](benchmarks/duvel3-a30-2026-11-06/balanced-results.json) used `chunk_size = ceil(N/2)` to engage both
A30s on each small case. Medians were **0.765 s** (Steane),
**0.783 s** (colour code), and
**2.015 s** (Golay); all 12 runs passed. These
noninterleaved follow-up results are not substituted into the primary table.

The default Golay job was 5.58× faster than the smaller CPU job, while Steane
was slower and colour-code timings were comparable within observed variation.
These are implementation/end-to-end results, not isolated GPU-kernel timings
or evidence of general two-GPU scaling. Only four all-good verification inputs
and three repetitions were tested. The original CPU branch's MPI shutdown
warnings remain in the logs. Neither implementation was changed for these
benchmark runs.

[Individual timings and ranges](benchmarks/duvel3-a30-2026-11-06/results.json),
[input hashes](benchmarks/duvel3-a30-2026-11-06/inputs/inputs.json),
and [commands/configuration](benchmarks/duvel3-a30-2026-11-06/comparison.json) are retained.
Raw timing logs are kept privately outside version control.
Reproduce in the recorded remote `benchmark` directory with
`python3 run_comparison.py comparison.json --output-directory results-new`.
Circuit inputs are kept outside version control. Obtain the named inputs
separately, check them against the recorded hashes, and update the input paths
in your comparison configuration. Historical logs and metadata retain the
paths used for the measurements.
