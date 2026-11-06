# CUDA benchmark record

This file intentionally contains no invented performance thresholds.  Fill it
from repeatable measurements on the target A100 system, then review and freeze
the promotion gates before optimizing the implementation.

## Baseline metadata

| Field | Value |
| --- | --- |
| Git commit | pending |
| Date | pending |
| Corpus version/hash | pending |
| CUDA toolkit / driver | pending |
| MPI / Slurm | pending |
| Nodes and GPUs per node | pending |
| GPU model, memory, clocks, MIG | pending |
| CPU and affinity | pending |
| GPU/NIC topology | pending |
| Build type and compiler flags | pending |

## Required workloads

| Workload | Main CPU | 1 A100 | 2 A100 | 4 A100 | 8 A100 | 2 nodes | 4 nodes |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Paper: 824 faults | pending | pending | pending | pending | pending | pending | pending |
| Paper: 5,619,552 pairs | pending | pending | pending | pending | pending | pending | pending |
| Paper: 12,662,112 pairs | pending | pending | pending | pending | pending | pending | pending |
| Production known-good | pending | pending | pending | pending | pending | pending | pending |
| Production known-bad | pending | pending | pending | pending | pending | pending | pending |
| Seeded distance-seven `s=3` | pending | pending | pending | pending | pending | pending | pending |

For each cell retain raw repetitions and report the median, not only the best
run.  Record setup and transfer separately from verification throughput.

## Promotion gate

The numeric gate is pending the first CPU and one-A100 baselines.  Once frozen,
an optimization is promoted only when it:

1. preserves all deterministic correctness and sanitizer gates;
2. improves median throughput on the agreed large-workload set;
3. causes no agreed small-workload regression; and
4. has a profile-supported explanation rather than a timing-only guess.

Any later change to the gate must be recorded with its rationale and a new
baseline commit.
