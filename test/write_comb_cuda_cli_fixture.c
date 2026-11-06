/* Produce native-endian legacy .run files for the actual MPI executable.
 * The bad fixture's only bad candidate is rank 3 (output-error index 3).
 * With epoch_size=2, a checkpoint records the preceding all-good prefix. */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void u32(FILE *stream, uint32_t value) {
  if (fwrite(&value, sizeof(value), 1, stream) != 1) {
    fputs("cannot write CLI scientific fixture\n", stderr);
    exit(EXIT_FAILURE);
  }
}

static void vector(FILE *stream, uint32_t value) {
  u32(stream, 1);
  /* The historic serializer writes 1 + bit_count / 4 uint32_t words. */
  u32(stream, value);
}

static void dcp(FILE *stream, uint32_t hx) {
  u32(stream, 1);
  u32(stream, 1);
  u32(stream, 0);
  vector(stream, hx);
  vector(stream, 0);
}

static void fixture(const char *directory, const char *name, uint32_t lookup) {
  static const uint32_t header[] = {1, 1, 0, 0, 1, 4, 1, 4, 0, 0};
  const size_t length = strlen(directory) + strlen(name) + 2;
  char *path = malloc(length);
  FILE *stream;
  if (path == NULL) exit(EXIT_FAILURE);
  snprintf(path, length, "%s/%s", directory, name);
  stream = fopen(path, "wb");
  if (stream == NULL) {
    fprintf(stderr, "cannot open fixture %s\n", path);
    free(path);
    exit(EXIT_FAILURE);
  }
  free(path);
  for (size_t i = 0; i < sizeof(header) / sizeof(header[0]); ++i)
    u32(stream, header[i]);
  u32(stream, lookup);
  u32(stream, 3); /* Both canonical Z errors are accepted. */
  /* Z is an independent logical DCP stabiliser on this one-qubit code. */
  u32(stream, 1);
  u32(stream, 1);
  u32(stream, 0);
  vector(stream, 0);
  vector(stream, 1);
  for (unsigned i = 0; i < 4; ++i) u32(stream, 1);
  for (unsigned i = 0; i < 4; ++i) dcp(stream, i == 3 ? 1 : 0);
  if (fclose(stream) != 0) exit(EXIT_FAILURE);
}

int main(int argc, char **argv) {
  if (argc != 2) {
    fputs("Usage: write_comb_cuda_cli_fixture DIRECTORY\n", stderr);
    return EXIT_FAILURE;
  }
  fixture(argv[1], "good.run", 3);
  fixture(argv[1], "bad.run", 1);
  return EXIT_SUCCESS;
}
