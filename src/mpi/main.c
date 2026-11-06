#include "comb_cuda_mpi.h"
#include "comb_run_file.h"

#include <errno.h>
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void usage(FILE *stream) {
  fprintf(stream,
          "Usage: comb-cuda-mpi input.run [--fast-any] [--epoch-size N]\n"
          "       [--chunk-size N] [--scratch-bytes N] [--checkpoint PATH] [--resume]\n"
          "Run with mpiexec or srun, with one MPI rank per NVIDIA GPU.\n");
}

static int read_number(const char *text, uint64_t *value) {
  char *end = NULL;
  unsigned long long parsed;
  if (text == NULL || text[0] == '\0' || text[0] == '-') return 0;
  errno = 0;
  parsed = strtoull(text, &end, 10);
  if (errno != 0 || end == text || *end != '\0' || parsed == 0) return 0;
  *value = (uint64_t)parsed;
  return 1;
}

static int parse_options(int argc, char **argv,
                         comb_cuda_mpi_options_t *options) {
  for (int i = 2; i < argc; ++i) {
    if (strcmp(argv[i], "--fast-any") == 0) {
      options->witness_mode = COMB_CUDA_WITNESS_FAST_ANY;
    } else if (strcmp(argv[i], "--resume") == 0) {
      options->resume_from_checkpoint = 1;
    } else if (strcmp(argv[i], "--checkpoint") == 0) {
      if (++i == argc || argv[i][0] == '\0') return 0;
      options->checkpoint_path = argv[i];
    } else {
      uint64_t *destination = NULL;
      if (strcmp(argv[i], "--epoch-size") == 0)
        destination = &options->epoch_size;
      else if (strcmp(argv[i], "--chunk-size") == 0)
        destination = &options->chunk_size;
      else if (strcmp(argv[i], "--scratch-bytes") == 0)
        destination = &options->cuda_options.scratch_limit_bytes;
      if (destination == NULL || ++i == argc ||
          !read_number(argv[i], destination)) return 0;
    }
  }
  return !options->resume_from_checkpoint || options->checkpoint_path != NULL;
}

int main(int argc, char **argv) {
  int rank = -1;
  int provided = 0;
  int input_status = 0;
  int exit_code = OTHER_ERROR;
  comb_run_file_t run = {0};
  comb_cuda_mpi_options_t options = COMB_CUDA_MPI_OPTIONS_INITIALIZER;
  comb_cuda_mpi_result_t result = COMB_CUDA_MPI_RESULT_INITIALIZER;
  comb_cuda_witness_t witness = COMB_CUDA_WITNESS_INITIALIZER;
  uint32_t *faults = NULL;

  if (MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided) != MPI_SUCCESS)
    return OTHER_ERROR;
  if (MPI_Comm_rank(MPI_COMM_WORLD, &rank) != MPI_SUCCESS) {
    MPI_Finalize();
    return OTHER_ERROR;
  }
  if (provided < MPI_THREAD_FUNNELED) {
    if (rank == 0) fprintf(stderr, "MPI_THREAD_FUNNELED is required.\n");
    MPI_Finalize();
    return OTHER_ERROR;
  }
  if (rank == 0) {
    if (argc == 2 && strcmp(argv[1], "--help") == 0) {
      usage(stdout);
      input_status = -1;
    } else if (argc < 2 || !parse_options(argc, argv, &options)) {
      usage(stderr);
      input_status = INVALID_INPUT;
    } else {
      const comb_cuda_status_t status = comb_run_file_load(argv[1], &run);
      if (status != COMB_CUDA_STATUS_SUCCESS) {
        fprintf(stderr, "Cannot read %s: %s\n", argv[1],
                comb_cuda_status_string(status));
        input_status = INVALID_INPUT;
      } else if (run.problem.s != 0) {
        faults = (uint32_t *)calloc(run.problem.s, sizeof(*faults));
        if (faults == NULL) input_status = OTHER_ERROR;
        witness.fault_capacity = run.problem.s;
        witness.fault_indices = faults;
      }
    }
  }
  if (MPI_Bcast(&input_status, 1, MPI_INT, 0, MPI_COMM_WORLD) != MPI_SUCCESS)
    goto done;
  if (input_status != 0) {
    exit_code = input_status == -1 ? 0 : input_status;
    goto done;
  }
  {
    const comb_cuda_mpi_status_t status = comb_cuda_mpi_verify_comm(
        MPI_COMM_WORLD, rank == 0 ? &run.problem : NULL,
        rank == 0 ? &options : NULL, &result, rank == 0 ? &witness : NULL);
    if (status == COMB_CUDA_MPI_STATUS_SUCCESS) {
      exit_code = (int)result.verification_code;
      if (rank == 0) {
        printf("code = %d\nchecked = %" PRIu64 " / %" PRIu64 "\n",
               result.verification_code, result.checked_combination_count,
               result.total_combination_count);
        if (result.verification_code == BAD_FAULT_FOUND) {
          printf("bad_rank = %" PRIu64 "\n", witness.rank);
          for (uint32_t i = 0; i < witness.fault_count; ++i)
            printf("bad_fault_array[%" PRIu32 "]=%" PRIu32 "\n", i, faults[i]);
        }
      }
    } else if (rank == 0) {
      fprintf(stderr, "MPI verification failed: %s; CUDA: %s (MPI code %d)\n",
              comb_cuda_mpi_status_string(status),
              comb_cuda_status_string(result.cuda_status), result.mpi_error_code);
    }
  }
done:
  free(faults);
  comb_run_file_free(&run);
  if (MPI_Finalize() != MPI_SUCCESS) return OTHER_ERROR;
  return exit_code;
}
