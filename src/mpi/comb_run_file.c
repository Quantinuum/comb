#include "comb_run_file.h"

#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* The historic writer stores native-endian uint32_t words.  Its external bit
 * representation has eight logical bits in each uint64_t, and it writes
 * 1 + logical_bit_count / 4 uint32_t words (including padding).  Preserve that
 * representation here; the flat CUDA packer subsequently makes dense words. */
typedef struct {
  FILE *stream;
  uint64_t remaining;
  comb_cuda_status_t status;
} reader_t;

static int read_u32(reader_t *reader, uint32_t *value) {
  if (reader->status != COMB_CUDA_STATUS_SUCCESS) return 0;
  if (reader->remaining < sizeof(*value) ||
      fread(value, sizeof(*value), 1, reader->stream) != 1) {
    reader->status = COMB_CUDA_STATUS_INVALID_PROBLEM;
    return 0;
  }
  reader->remaining -= sizeof(*value);
  return 1;
}

static void *allocate(reader_t *reader, uint64_t count, size_t size) {
  void *result;
  if (reader->status != COMB_CUDA_STATUS_SUCCESS || count == 0) return NULL;
  if (count > SIZE_MAX / size) {
    reader->status = COMB_CUDA_STATUS_PROBLEM_TOO_LARGE;
    return NULL;
  }
  result = calloc((size_t)count, size);
  if (result == NULL) reader->status = COMB_CUDA_STATUS_OUT_OF_MEMORY;
  return result;
}

static comb_bit_word_t *read_bits(reader_t *reader, uint64_t bits) {
  const uint64_t serialized_words = 1 + bits / 4;
  const uint64_t external_words = 1 + bits / 8;
  const uint64_t byte_count = serialized_words * sizeof(uint32_t);
  comb_bit_word_t *result;
  if (reader->status != COMB_CUDA_STATUS_SUCCESS) return NULL;
  /* Original serializer uses uint32_t dimensions and int iteration counters.
   * Reject sizes it cannot have represented without overflow. */
  if (bits > INT_MAX || byte_count > reader->remaining) {
    reader->status = bits > INT_MAX ? COMB_CUDA_STATUS_PROBLEM_TOO_LARGE
                                  : COMB_CUDA_STATUS_INVALID_PROBLEM;
    return NULL;
  }
  result = (comb_bit_word_t *)allocate(reader, external_words, sizeof(*result));
  if (result != NULL &&
      fread(result, 1, (size_t)byte_count, reader->stream) != byte_count) {
    reader->status = COMB_CUDA_STATUS_INVALID_PROBLEM;
  }
  if (reader->status == COMB_CUDA_STATUS_SUCCESS) reader->remaining -= byte_count;
  return result;
}

static comb_vector_t *read_vectors(reader_t *reader, uint32_t count,
                                    uint32_t width) {
  comb_vector_t *vectors;
  /* Each vector has at least a dimension and one data word. */
  if ((uint64_t)count > reader->remaining / 8) {
    reader->status = COMB_CUDA_STATUS_INVALID_PROBLEM;
    return NULL;
  }
  vectors = (comb_vector_t *)allocate(reader, count, sizeof(*vectors));
  for (uint32_t i = 0; vectors != NULL && i < count &&
                       reader->status == COMB_CUDA_STATUS_SUCCESS; ++i) {
    if (!read_u32(reader, &vectors[i].n)) break;
    if (vectors[i].n != width) {
      reader->status = COMB_CUDA_STATUS_INVALID_PROBLEM;
      break;
    }
    vectors[i].data = read_bits(reader, width);
  }
  return vectors;
}

static comb_dcp_t *read_dcps(reader_t *reader, uint32_t count, uint32_t width) {
  comb_dcp_t *dcps;
  /* Matrix dimensions, vector dimensions and at least three data words. */
  if ((uint64_t)count > reader->remaining / 28) {
    reader->status = COMB_CUDA_STATUS_INVALID_PROBLEM;
    return NULL;
  }
  dcps = (comb_dcp_t *)allocate(reader, count, sizeof(*dcps));
  for (uint32_t i = 0; dcps != NULL && i < count &&
                       reader->status == COMB_CUDA_STATUS_SUCCESS; ++i) {
    comb_dcp_t *dcp = &dcps[i];
    if (!read_u32(reader, &dcp->F.m) || !read_u32(reader, &dcp->F.n)) break;
    if (dcp->F.m != width || dcp->F.n != width) {
      reader->status = COMB_CUDA_STATUS_INVALID_PROBLEM;
      break;
    }
    dcp->F.data = read_bits(reader, (uint64_t)width * width);
    if (!read_u32(reader, &dcp->h_x.n)) break;
    if (dcp->h_x.n != width) {
      reader->status = COMB_CUDA_STATUS_INVALID_PROBLEM;
      break;
    }
    dcp->h_x.data = read_bits(reader, width);
    if (!read_u32(reader, &dcp->h_z.n)) break;
    if (dcp->h_z.n != width) {
      reader->status = COMB_CUDA_STATUS_INVALID_PROBLEM;
      break;
    }
    dcp->h_z.data = read_bits(reader, width);
  }
  return dcps;
}

static void free_vectors(const comb_vector_t *vectors, uint32_t count) {
  if (vectors == NULL) return;
  for (uint32_t i = 0; i < count; ++i) free(vectors[i].data);
  free((void *)vectors);
}

static void free_dcps(const comb_dcp_t *dcps, uint32_t count) {
  if (dcps == NULL) return;
  for (uint32_t i = 0; i < count; ++i) {
    free(dcps[i].F.data);
    free(dcps[i].h_x.data);
    free(dcps[i].h_z.data);
  }
  free((void *)dcps);
}

void comb_run_file_free(comb_run_file_t *run) {
  comb_problem_t *p;
  if (run == NULL) return;
  p = &run->problem;
  free_vectors(p->x_stabilisers, p->n_x);
  free_vectors(p->z_stabilisers, p->n_z);
  free((void *)p->x_lookup_table);
  free((void *)p->z_lookup_table);
  free_dcps(p->dcp_stabilisers, p->k);
  free_dcps(p->output_errors, p->output_error_count);
  free_vectors(p->x_measurement_paulis, p->x_measurement_count);
  free_vectors(p->z_measurement_paulis, p->z_measurement_count);
  free((void *)p->fault_location_partitions);
  memset(run, 0, sizeof(*run));
}

comb_cuda_status_t comb_run_file_load(const char *path, comb_run_file_t *run) {
  comb_problem_t *p;
  long length;
  reader_t reader = {NULL, 0, COMB_CUDA_STATUS_SUCCESS};
  if (path == NULL || run == NULL) return COMB_CUDA_STATUS_INVALID_ARGUMENT;
  /* The caller must release a prior load before loading a new file. */
  memset(run, 0, sizeof(*run));
  p = &run->problem;
  p->struct_size = sizeof(*p);
  p->abi_version = COMB_CUDA_ABI_VERSION;
  reader.stream = fopen(path, "rb");
  if (reader.stream == NULL) return COMB_CUDA_STATUS_INVALID_PROBLEM;
  if (fseek(reader.stream, 0, SEEK_END) != 0 ||
      (length = ftell(reader.stream)) < 0 ||
      fseek(reader.stream, 0, SEEK_SET) != 0) {
    fclose(reader.stream);
    return COMB_CUDA_STATUS_INVALID_PROBLEM;
  }
  reader.remaining = (uint64_t)length;
  if (!read_u32(&reader, &p->n_all) || !read_u32(&reader, &p->n) ||
      !read_u32(&reader, &p->n_x) || !read_u32(&reader, &p->n_z) ||
      !read_u32(&reader, &p->k) ||
      !read_u32(&reader, &p->fault_location_count) ||
      !read_u32(&reader, &p->s) ||
      !read_u32(&reader, &p->output_error_count) ||
      !read_u32(&reader, &p->x_measurement_count) ||
      !read_u32(&reader, &p->z_measurement_count)) goto done;
  if (p->n_all == 0 || p->n > p->n_all || p->n_x > p->n ||
      p->n_z > p->n || p->n - p->n_x >= 31 || p->n - p->n_z >= 31 ||
      p->s > p->fault_location_count || p->k > p->n ||
      (uint64_t)p->n_all * p->n_all > INT_MAX) {
    reader.status = COMB_CUDA_STATUS_INVALID_PROBLEM;
    goto done;
  }
  p->x_stabilisers = read_vectors(&reader, p->n_x, p->n);
  p->x_lookup_bit_count = UINT64_C(1) << (p->n - p->n_x);
  p->x_lookup_table = read_bits(&reader, p->x_lookup_bit_count);
  p->z_stabilisers = read_vectors(&reader, p->n_z, p->n);
  p->z_lookup_bit_count = UINT64_C(1) << (p->n - p->n_z);
  p->z_lookup_table = read_bits(&reader, p->z_lookup_bit_count);
  p->dcp_stabilisers = read_dcps(&reader, p->k, p->n_all);
  if ((uint64_t)p->fault_location_count > reader.remaining / sizeof(uint32_t)) {
    reader.status = COMB_CUDA_STATUS_INVALID_PROBLEM;
  } else {
    uint32_t *partitions = (uint32_t *)allocate(
        &reader, p->fault_location_count, sizeof(*partitions));
    p->fault_location_partitions = partitions;
    for (uint32_t i = 0; partitions != NULL && i < p->fault_location_count;
         ++i) {
      if (!read_u32(&reader, &partitions[i])) break;
    }
  }
  p->output_errors = read_dcps(&reader, p->output_error_count, p->n_all);
  p->x_measurement_paulis = read_vectors(
      &reader, p->x_measurement_count, p->n_all);
  p->z_measurement_paulis = read_vectors(
      &reader, p->z_measurement_count, p->n_all);
  if (reader.status == COMB_CUDA_STATUS_SUCCESS && reader.remaining != 0)
    reader.status = COMB_CUDA_STATUS_INVALID_PROBLEM;
done:
  if (fclose(reader.stream) != 0 && reader.status == COMB_CUDA_STATUS_SUCCESS)
    reader.status = COMB_CUDA_STATUS_INVALID_PROBLEM;
  if (reader.status != COMB_CUDA_STATUS_SUCCESS) comb_run_file_free(run);
  return reader.status;
}
