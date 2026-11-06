#define _POSIX_C_SOURCE 200809L
#include "comb_run_file.h"
#include "comb_cuda_internal.h"

#include <assert.h>

#undef assert
#define assert(condition) do { if (!(condition)) { \
  fprintf(stderr, "check failed at %s:%d: %s\n", __FILE__, __LINE__, #condition); \
  exit(EXIT_FAILURE); } } while (0)
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

static void u32(FILE *stream, uint32_t value) {
  assert(fwrite(&value, sizeof(value), 1, stream) == 1);
}

static void bits(FILE *stream, uint32_t count, uint32_t set) {
  uint64_t *external = calloc(1 + count / 8, sizeof(*external));
  assert(external != NULL);
  if (set < count) external[set / 8] |= UINT64_C(1) << (set % 8);
  assert(fwrite(external, 4, 1 + count / 4, stream) == 1 + count / 4);
  free(external);
}

static void vector(FILE *stream, uint32_t width, uint32_t set) {
  u32(stream, width);
  bits(stream, width, set);
}

static void dcp(FILE *stream, uint32_t width, uint32_t set) {
  u32(stream, width);
  u32(stream, width);
  bits(stream, width * width, UINT32_MAX);
  vector(stream, width, set);
  vector(stream, width, UINT32_MAX);
}

static void write_fixture(FILE *stream, uint32_t vector_width) {
  const uint32_t header[] = {9, 2, 1, 1, 0, 1, 1, 1, 0, 0};
  for (size_t i = 0; i < sizeof(header) / sizeof(header[0]); ++i)
    u32(stream, header[i]);
  vector(stream, vector_width, 0);
  bits(stream, 2, 0);
  vector(stream, 2, 1);
  bits(stream, 2, 0);
  u32(stream, 1);
  dcp(stream, 9, 8);
}

int main(void) {
  char path[] = "/tmp/comb-run-file-XXXXXX";
  int descriptor = mkstemp(path);
  FILE *stream;
  comb_run_file_t run = {0};
  void *blob = NULL;
  uint64_t blob_size = 0;
  assert(descriptor >= 0);
  stream = fdopen(descriptor, "wb");
  assert(stream != NULL);
  write_fixture(stream, 2);
  assert(fclose(stream) == 0);
  assert(comb_run_file_load(path, &run) == COMB_CUDA_STATUS_SUCCESS);
  assert(run.problem.n_all == 9 && run.problem.n == 2);
  assert(run.problem.output_errors[0].h_x.data[1] == 1);
  assert(run.problem.fault_location_partitions[0] == 1);
  assert(comb_cuda_flat_problem_pack(&run.problem, &blob, &blob_size) ==
         COMB_CUDA_STATUS_SUCCESS);
  {
    const comb_cuda_flat_problem_t *flat = blob;
    const uint64_t *hx = (const uint64_t *)((const char *)blob + flat->output_dcp_hx_offset);
    assert(hx[0] == (UINT64_C(1) << 8));
  }
  comb_cuda_flat_problem_free(blob);
  comb_run_file_free(&run);

  /* Every truncated prefix must fail cleanly, including partial DCP rows. */
  {
    unsigned char encoded[512];
    size_t length;
    stream = fopen(path, "rb");
    assert(stream != NULL);
    length = fread(encoded, 1, sizeof(encoded), stream);
    assert(length > 0 && length < sizeof(encoded));
    assert(fclose(stream) == 0);
    for (size_t prefix = 0; prefix < length; ++prefix) {
      stream = fopen(path, "wb");
      assert(stream != NULL);
      assert(fwrite(encoded, 1, prefix, stream) == prefix);
      assert(fclose(stream) == 0);
      assert(comb_run_file_load(path, &run) != COMB_CUDA_STATUS_SUCCESS);
      comb_run_file_free(&run);
    }
    stream = fopen(path, "wb");
    assert(stream != NULL);
    assert(fwrite(encoded, 1, length, stream) == length);
    assert(fputc(0, stream) == 0);
    assert(fclose(stream) == 0);
    assert(comb_run_file_load(path, &run) == COMB_CUDA_STATUS_INVALID_PROBLEM);
  }

  stream = fopen(path, "wb");
  assert(stream != NULL);
  write_fixture(stream, 3);
  assert(fclose(stream) == 0);
  assert(comb_run_file_load(path, &run) == COMB_CUDA_STATUS_INVALID_PROBLEM);
  comb_run_file_free(&run);

  stream = fopen(path, "wb");
  assert(stream != NULL);
  u32(stream, UINT32_MAX);
  assert(fclose(stream) == 0);
  assert(comb_run_file_load(path, &run) == COMB_CUDA_STATUS_INVALID_PROBLEM);
  assert(unlink(path) == 0);
  assert(comb_run_file_load(path, &run) == COMB_CUDA_STATUS_INVALID_PROBLEM);
  puts("Legacy .run parser: compatible bit layout and invalid-input checks passed");
  return 0;
}
