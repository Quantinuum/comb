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

# The repository intentionally has no generic production-problem parser.
# Point this at an application that constructs comb_problem_t on MPI rank zero
# and calls comb_cuda_mpi_verify_comm.
: "${COMB_MPI_EXECUTABLE:?Set COMB_MPI_EXECUTABLE to the embedding executable}"

# Many clusters use pmix; replace this with the MPI plugin configured locally.
# GPU binding should give each task a distinct CUDA_VISIBLE_DEVICES assignment.
srun --mpi=pmix --gpu-bind=single:1 "${COMB_MPI_EXECUTABLE}" "$@"
