#!/usr/bin/env bash
# Example only: adapt account, partition, module names, and task count to the
# target cluster.  This requests one MPI rank and one A100 for each task.
#SBATCH --job-name=comb-cuda-mpi
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=4
#SBATCH --gpus-per-node=a100:4
#SBATCH --gpus-per-task=1
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
#SBATCH --output=comb-cuda-mpi-%j.out
#SBATCH --error=comb-cuda-mpi-%j.err

set -euo pipefail

# Pass the serialized .run input file and any CLI options to sbatch.
# Override COMB_MPI_EXECUTABLE when the build directory is elsewhere.
COMB_MPI_EXECUTABLE="${COMB_MPI_EXECUTABLE:-./build-cuda/comb-cuda-mpi}"

# Many clusters use pmix; replace this with the MPI plugin configured locally.
# GPU binding should give each task a distinct CUDA_VISIBLE_DEVICES assignment.
srun --mpi=pmix --gpu-bind=single:1 "${COMB_MPI_EXECUTABLE}" "$@"
