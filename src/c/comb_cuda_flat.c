#include "comb_cuda_internal.h"

#include <limits.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#define FLAT_ALIGNMENT UINT64_C(8)

typedef struct {
  unsigned char *data;
  uint64_t size;
} flat_builder_t;

static int add_u64(uint64_t left, uint64_t right, uint64_t *result) {
  if (result == NULL || right > UINT64_MAX - left) return 0;
  *result = left + right;
  return 1;
}

static int multiply_u64(uint64_t left, uint64_t right, uint64_t *result) {
  if (result == NULL || (left != 0 && right > UINT64_MAX / left)) return 0;
  *result = left * right;
  return 1;
}

static int align_u64(uint64_t value, uint64_t *aligned) {
  const uint64_t mask = FLAT_ALIGNMENT - 1;
  if (value > UINT64_MAX - mask) return 0;
  *aligned = (value + mask) & ~mask;
  return 1;
}

static uint32_t dense_word_count(uint32_t bit_count) {
  return (uint32_t)(((uint64_t)bit_count + 63) / 64);
}

static uint64_t tail_mask_u64(uint64_t bit_count) {
  const uint32_t tail = (uint32_t)(bit_count & UINT64_C(63));
  return tail == 0 ? UINT64_MAX : (UINT64_C(1) << tail) - 1;
}

static uint32_t legacy_bit(const comb_bit_word_t *data, uint64_t bit_index) {
  return (uint32_t)((data[bit_index / sizeof(comb_bit_word_t)] >>
                     (bit_index % sizeof(comb_bit_word_t))) &
                    UINT64_C(1));
}

static void pack_legacy_bits(const comb_bit_word_t *source,
                             uint64_t source_bit_begin, uint64_t bit_count,
                             uint64_t *destination) {
  const uint64_t word_count = (bit_count + 63) / 64;
  memset(destination, 0, (size_t)(word_count * sizeof(uint64_t)));
  for (uint64_t bit = 0; bit < bit_count; ++bit) {
    if (legacy_bit(source, source_bit_begin + bit) != 0) {
      destination[bit / 64] |= UINT64_C(1) << (bit & UINT64_C(63));
    }
  }
}

static comb_cuda_status_t builder_init(flat_builder_t *builder) {
  builder->size = sizeof(comb_cuda_flat_problem_t);
  builder->data = (unsigned char *)calloc((size_t)builder->size, 1);
  return builder->data == NULL ? COMB_CUDA_STATUS_OUT_OF_MEMORY
                               : COMB_CUDA_STATUS_SUCCESS;
}

static comb_cuda_status_t builder_append(flat_builder_t *builder,
                                         uint64_t byte_count,
                                         uint64_t *offset) {
  uint64_t aligned = 0;
  uint64_t end = 0;
  if (builder == NULL || offset == NULL ||
      !align_u64(builder->size, &aligned) ||
      !add_u64(aligned, byte_count, &end) || end > (uint64_t)SIZE_MAX) {
    return COMB_CUDA_STATUS_PROBLEM_TOO_LARGE;
  }
  if (byte_count == 0) {
    *offset = 0;
    return COMB_CUDA_STATUS_SUCCESS;
  }
  unsigned char *resized = (unsigned char *)realloc(builder->data, (size_t)end);
  if (resized == NULL) return COMB_CUDA_STATUS_OUT_OF_MEMORY;
  builder->data = resized;
  memset(builder->data + builder->size, 0, (size_t)(end - builder->size));
  builder->size = end;
  *offset = aligned;
  return COMB_CUDA_STATUS_SUCCESS;
}

static comb_cuda_status_t validate_dcp_array(const comb_dcp_t *dcps,
                                             uint32_t count,
                                             uint32_t n_all) {
  if (count != 0 && dcps == NULL) return COMB_CUDA_STATUS_INVALID_PROBLEM;
  for (uint32_t index = 0; index < count; ++index) {
    const comb_dcp_t *dcp = dcps + index;
    if (dcp->F.m != n_all || dcp->F.n != n_all || dcp->F.data == NULL ||
        dcp->h_x.n != n_all || dcp->h_x.data == NULL ||
        dcp->h_z.n != n_all || dcp->h_z.data == NULL) {
      return COMB_CUDA_STATUS_INVALID_PROBLEM;
    }
  }
  return COMB_CUDA_STATUS_SUCCESS;
}

static comb_cuda_status_t validate_vector_array(const comb_vector_t *vectors,
                                                uint32_t count,
                                                uint32_t width) {
  if (count != 0 && vectors == NULL) return COMB_CUDA_STATUS_INVALID_PROBLEM;
  for (uint32_t index = 0; index < count; ++index) {
    if (vectors[index].n != width || vectors[index].data == NULL) {
      return COMB_CUDA_STATUS_INVALID_PROBLEM;
    }
  }
  return COMB_CUDA_STATUS_SUCCESS;
}

static comb_cuda_status_t validate_problem(const comb_problem_t *problem,
                                           uint64_t *partition_sum) {
  if (problem == NULL || partition_sum == NULL ||
      problem->struct_size < sizeof(comb_problem_t) ||
      problem->abi_version != COMB_CUDA_ABI_VERSION) {
    return COMB_CUDA_STATUS_INVALID_ARGUMENT;
  }
  if (problem->n_all == 0 || problem->n == 0 || problem->n > problem->n_all ||
      problem->n_x > problem->n || problem->n_z > problem->n ||
      problem->k > problem->n_all || problem->fault_location_count == 0 ||
      problem->fault_location_count == UINT32_MAX || problem->s == 0 ||
      problem->s == UINT32_MAX || problem->s > problem->fault_location_count ||
      problem->output_error_count == 0) {
    return COMB_CUDA_STATUS_INVALID_PROBLEM;
  }
  const uint32_t x_lookup_exponent = problem->n - problem->n_x;
  const uint32_t z_lookup_exponent = problem->n - problem->n_z;
  if (x_lookup_exponent >= 64 || z_lookup_exponent >= 64) {
    return COMB_CUDA_STATUS_PROBLEM_TOO_LARGE;
  }
  if (problem->x_lookup_bit_count != (UINT64_C(1) << x_lookup_exponent) ||
      problem->z_lookup_bit_count != (UINT64_C(1) << z_lookup_exponent) ||
      problem->x_lookup_table == NULL || problem->z_lookup_table == NULL ||
      problem->fault_location_partitions == NULL) {
    return COMB_CUDA_STATUS_INVALID_PROBLEM;
  }

  comb_cuda_status_t status = validate_vector_array(
      problem->x_stabilisers, problem->n_x, problem->n);
  if (status != COMB_CUDA_STATUS_SUCCESS) return status;
  status = validate_vector_array(problem->z_stabilisers, problem->n_z,
                                 problem->n);
  if (status != COMB_CUDA_STATUS_SUCCESS) return status;
  status = validate_vector_array(problem->x_measurement_paulis,
                                 problem->x_measurement_count, problem->n_all);
  if (status != COMB_CUDA_STATUS_SUCCESS) return status;
  status = validate_vector_array(problem->z_measurement_paulis,
                                 problem->z_measurement_count, problem->n_all);
  if (status != COMB_CUDA_STATUS_SUCCESS) return status;
  status =
      validate_dcp_array(problem->dcp_stabilisers, problem->k, problem->n_all);
  if (status != COMB_CUDA_STATUS_SUCCESS) return status;
  status = validate_dcp_array(problem->output_errors,
                              problem->output_error_count, problem->n_all);
  if (status != COMB_CUDA_STATUS_SUCCESS) return status;

  uint64_t sum = 0;
  for (uint32_t location = 0; location < problem->fault_location_count;
       ++location) {
    const uint32_t partition = problem->fault_location_partitions[location];
    if (partition == 0) return COMB_CUDA_STATUS_INVALID_PROBLEM;
    if (!add_u64(sum, partition, &sum)) {
      return COMB_CUDA_STATUS_PROBLEM_TOO_LARGE;
    }
  }
  if (sum > problem->output_error_count) {
    return COMB_CUDA_STATUS_INVALID_PROBLEM;
  }
  *partition_sum = sum;
  return COMB_CUDA_STATUS_SUCCESS;
}

static comb_cuda_status_t append_raw_array(flat_builder_t *builder,
                                           const void *source,
                                           uint64_t element_count,
                                           uint64_t element_size,
                                           uint64_t *offset) {
  uint64_t byte_count = 0;
  if (!multiply_u64(element_count, element_size, &byte_count)) {
    return COMB_CUDA_STATUS_PROBLEM_TOO_LARGE;
  }
  comb_cuda_status_t status = builder_append(builder, byte_count, offset);
  if (status != COMB_CUDA_STATUS_SUCCESS || byte_count == 0) return status;
  memcpy(builder->data + *offset, source, (size_t)byte_count);
  return COMB_CUDA_STATUS_SUCCESS;
}

static comb_cuda_status_t append_vector_rows(flat_builder_t *builder,
                                             const comb_vector_t *vectors,
                                             uint32_t row_count,
                                             uint32_t bit_count,
                                             uint64_t *offset) {
  const uint64_t words = dense_word_count(bit_count);
  uint64_t element_count = 0;
  uint64_t byte_count = 0;
  if (!multiply_u64(row_count, words, &element_count) ||
      !multiply_u64(element_count, sizeof(uint64_t), &byte_count)) {
    return COMB_CUDA_STATUS_PROBLEM_TOO_LARGE;
  }
  comb_cuda_status_t status = builder_append(builder, byte_count, offset);
  if (status != COMB_CUDA_STATUS_SUCCESS || row_count == 0) return status;
  uint64_t *destination = (uint64_t *)(builder->data + *offset);
  for (uint32_t row = 0; row < row_count; ++row) {
    pack_legacy_bits(vectors[row].data, 0, bit_count,
                     destination + (uint64_t)row * words);
  }
  return COMB_CUDA_STATUS_SUCCESS;
}

static comb_cuda_status_t append_lookup(flat_builder_t *builder,
                                        const comb_bit_word_t *lookup,
                                        uint64_t bit_count, uint64_t *offset) {
  const uint64_t words = (bit_count + 63) / 64;
  uint64_t byte_count = 0;
  if (!multiply_u64(words, sizeof(uint64_t), &byte_count)) {
    return COMB_CUDA_STATUS_PROBLEM_TOO_LARGE;
  }
  comb_cuda_status_t status = builder_append(builder, byte_count, offset);
  if (status != COMB_CUDA_STATUS_SUCCESS) return status;
  pack_legacy_bits(lookup, 0, bit_count,
                   (uint64_t *)(builder->data + *offset));
  return COMB_CUDA_STATUS_SUCCESS;
}

static comb_cuda_status_t append_dcp_soa(
    flat_builder_t *builder, const comb_dcp_t *dcps, uint32_t count,
    uint32_t n_all, uint64_t *f_offset, uint64_t *hx_offset,
    uint64_t *hz_offset) {
  const uint64_t words = dense_word_count(n_all);
  uint64_t matrix_rows = 0;
  uint64_t matrix_words = 0;
  uint64_t vector_words = 0;
  uint64_t matrix_bytes = 0;
  uint64_t vector_bytes = 0;
  if (!multiply_u64(count, n_all, &matrix_rows) ||
      !multiply_u64(matrix_rows, words, &matrix_words) ||
      !multiply_u64(count, words, &vector_words) ||
      !multiply_u64(matrix_words, sizeof(uint64_t), &matrix_bytes) ||
      !multiply_u64(vector_words, sizeof(uint64_t), &vector_bytes)) {
    return COMB_CUDA_STATUS_PROBLEM_TOO_LARGE;
  }
  comb_cuda_status_t status = builder_append(builder, matrix_bytes, f_offset);
  if (status != COMB_CUDA_STATUS_SUCCESS) return status;
  status = builder_append(builder, vector_bytes, hx_offset);
  if (status != COMB_CUDA_STATUS_SUCCESS) return status;
  status = builder_append(builder, vector_bytes, hz_offset);
  if (status != COMB_CUDA_STATUS_SUCCESS || count == 0) return status;

  uint64_t *flat_f = (uint64_t *)(builder->data + *f_offset);
  uint64_t *flat_hx = (uint64_t *)(builder->data + *hx_offset);
  uint64_t *flat_hz = (uint64_t *)(builder->data + *hz_offset);
  for (uint32_t index = 0; index < count; ++index) {
    for (uint32_t row = 0; row < n_all; ++row) {
      pack_legacy_bits(dcps[index].F.data, (uint64_t)row * n_all, n_all,
                       flat_f + ((uint64_t)index * n_all + row) * words);
    }
    pack_legacy_bits(dcps[index].h_x.data, 0, n_all,
                     flat_hx + (uint64_t)index * words);
    pack_legacy_bits(dcps[index].h_z.data, 0, n_all,
                     flat_hz + (uint64_t)index * words);
  }
  return COMB_CUDA_STATUS_SUCCESS;
}

static comb_cuda_status_t append_rank_dp(
    flat_builder_t *builder, const uint32_t *partitions,
    uint32_t location_count, uint32_t s, uint64_t *offset,
    uint64_t *total_combination_count) {
  const uint32_t rows = location_count + 1;
  const uint32_t columns = s + 1;
  uint64_t element_count = 0;
  uint64_t byte_count = 0;
  if (!multiply_u64(rows, columns, &element_count) ||
      !multiply_u64(element_count, sizeof(uint64_t), &byte_count)) {
    return COMB_CUDA_STATUS_PROBLEM_TOO_LARGE;
  }
  comb_cuda_status_t status = builder_append(builder, byte_count, offset);
  if (status != COMB_CUDA_STATUS_SUCCESS) return status;
  uint64_t *dp = (uint64_t *)(builder->data + *offset);
  dp[0] = 1;
  for (uint32_t m = 1; m <= location_count; ++m) {
    uint64_t *current = dp + (uint64_t)m * columns;
    const uint64_t *previous = current - columns;
    current[0] = 1;
    const uint32_t maximum_r = m < s ? m : s;
    for (uint32_t r = 1; r <= maximum_r; ++r) {
      uint64_t term = 0;
      if (!multiply_u64(previous[r - 1], partitions[m - 1], &term) ||
          !add_u64(previous[r], term, current + r)) {
        return COMB_CUDA_STATUS_PROBLEM_TOO_LARGE;
      }
    }
  }
  *total_combination_count = dp[(uint64_t)location_count * columns + s];
  return COMB_CUDA_STATUS_SUCCESS;
}

comb_cuda_status_t comb_cuda_flat_problem_pack(const comb_problem_t *problem,
                                               void **blob,
                                               uint64_t *blob_size) {
  if (blob == NULL || blob_size == NULL) {
    return COMB_CUDA_STATUS_INVALID_ARGUMENT;
  }
  *blob = NULL;
  *blob_size = 0;
  uint64_t partition_sum = 0;
  comb_cuda_status_t status = validate_problem(problem, &partition_sum);
  if (status != COMB_CUDA_STATUS_SUCCESS) return status;

  flat_builder_t builder = {0};
  status = builder_init(&builder);
  if (status != COMB_CUDA_STATUS_SUCCESS) return status;
  comb_cuda_flat_problem_t header;
  memset(&header, 0, sizeof(header));
  header.magic = COMB_CUDA_FLAT_MAGIC;
  header.version = COMB_CUDA_FLAT_VERSION;
  header.header_size = sizeof(header);
  header.n_all = problem->n_all;
  header.n = problem->n;
  header.n_x = problem->n_x;
  header.n_z = problem->n_z;
  header.k = problem->k;
  header.words_all = dense_word_count(problem->n_all);
  header.words_code = dense_word_count(problem->n);
  header.fault_location_count = problem->fault_location_count;
  header.s = problem->s;
  header.rank_dp_rows = problem->fault_location_count + 1;
  header.rank_dp_columns = problem->s + 1;
  header.output_error_count = problem->output_error_count;
  header.x_measurement_count = problem->x_measurement_count;
  header.z_measurement_count = problem->z_measurement_count;
  header.x_lookup_bit_count = problem->x_lookup_bit_count;
  header.z_lookup_bit_count = problem->z_lookup_bit_count;

#define RETURN_PACK_ERROR(call)       \
  do {                                \
    status = (call);                  \
    if (status != COMB_CUDA_STATUS_SUCCESS) { \
      free(builder.data);             \
      return status;                  \
    }                                 \
  } while (0)

  RETURN_PACK_ERROR(append_raw_array(
      &builder, problem->fault_location_partitions,
      problem->fault_location_count, sizeof(uint32_t),
      &header.fault_partitions_offset));

  uint64_t prefix_count = (uint64_t)problem->fault_location_count + 1;
  uint64_t prefix_bytes = 0;
  if (!multiply_u64(prefix_count, sizeof(uint64_t), &prefix_bytes)) {
    free(builder.data);
    return COMB_CUDA_STATUS_PROBLEM_TOO_LARGE;
  }
  RETURN_PACK_ERROR(
      builder_append(&builder, prefix_bytes, &header.fault_prefixes_offset));
  uint64_t *prefixes =
      (uint64_t *)(builder.data + header.fault_prefixes_offset);
  prefixes[0] = 0;
  for (uint32_t location = 0; location < problem->fault_location_count;
       ++location) {
    prefixes[location + 1] =
        prefixes[location] + problem->fault_location_partitions[location];
  }
  if (prefixes[problem->fault_location_count] != partition_sum) {
    free(builder.data);
    return COMB_CUDA_STATUS_INTERNAL_ERROR;
  }

  RETURN_PACK_ERROR(append_rank_dp(
      &builder, problem->fault_location_partitions,
      problem->fault_location_count, problem->s, &header.rank_dp_offset,
      &header.total_combination_count));
  RETURN_PACK_ERROR(append_vector_rows(
      &builder, problem->x_stabilisers, problem->n_x, problem->n,
      &header.x_stabilisers_offset));
  RETURN_PACK_ERROR(append_vector_rows(
      &builder, problem->z_stabilisers, problem->n_z, problem->n,
      &header.z_stabilisers_offset));
  RETURN_PACK_ERROR(append_lookup(&builder, problem->x_lookup_table,
                                  problem->x_lookup_bit_count,
                                  &header.x_lookup_offset));
  RETURN_PACK_ERROR(append_lookup(&builder, problem->z_lookup_table,
                                  problem->z_lookup_bit_count,
                                  &header.z_lookup_offset));
  RETURN_PACK_ERROR(append_dcp_soa(
      &builder, problem->dcp_stabilisers, problem->k, problem->n_all,
      &header.logical_dcp_f_offset, &header.logical_dcp_hx_offset,
      &header.logical_dcp_hz_offset));
  RETURN_PACK_ERROR(append_dcp_soa(
      &builder, problem->output_errors, problem->output_error_count,
      problem->n_all, &header.output_dcp_f_offset,
      &header.output_dcp_hx_offset, &header.output_dcp_hz_offset));
  RETURN_PACK_ERROR(append_vector_rows(
      &builder, problem->x_measurement_paulis, problem->x_measurement_count,
      problem->n_all, &header.x_measurements_offset));
  RETURN_PACK_ERROR(append_vector_rows(
      &builder, problem->z_measurement_paulis, problem->z_measurement_count,
      problem->n_all, &header.z_measurements_offset));

#undef RETURN_PACK_ERROR

  header.blob_size = builder.size;
  memcpy(builder.data, &header, sizeof(header));
  *blob = builder.data;
  *blob_size = builder.size;
  return COMB_CUDA_STATUS_SUCCESS;
}

static int consume_section(uint64_t offset, uint64_t element_count,
                           uint64_t element_size, uint64_t blob_size,
                           uint64_t *cursor) {
  uint64_t byte_count = 0;
  uint64_t aligned = 0;
  uint64_t end = 0;
  if (!multiply_u64(element_count, element_size, &byte_count)) return 0;
  if (byte_count == 0) return offset == 0;
  if (!align_u64(*cursor, &aligned) || offset != aligned ||
      !add_u64(aligned, byte_count, &end) || end > blob_size) {
    return 0;
  }
  *cursor = end;
  return 1;
}

static int rows_tail_masked(const unsigned char *blob, uint64_t offset,
                            uint64_t row_count, uint32_t bit_count) {
  if (row_count == 0) return offset == 0;
  const uint32_t words = dense_word_count(bit_count);
  const uint64_t mask = tail_mask_u64(bit_count);
  if (mask == UINT64_MAX) return 1;
  const uint64_t *rows = (const uint64_t *)(blob + offset);
  for (uint64_t row = 0; row < row_count; ++row) {
    if ((rows[row * words + words - 1] & ~mask) != 0) return 0;
  }
  return 1;
}

static unsigned flat_word_bit(const uint64_t *row, uint32_t bit) {
  return (unsigned)((row[bit / 64U] >> (bit % 64U)) & UINT64_C(1));
}

static unsigned flat_row_dot(const uint64_t *left, const uint64_t *right,
                             uint32_t words) {
  unsigned parity = 0U;
  uint32_t word;
  for (word = 0U; word < words; ++word) {
#if defined(__GNUC__) || defined(__clang__)
    parity ^= (unsigned)__builtin_parityll(left[word] & right[word]);
#else
    uint64_t value = left[word] & right[word];
    value ^= value >> 32U;
    value ^= value >> 16U;
    value ^= value >> 8U;
    value ^= value >> 4U;
    value &= UINT64_C(0xf);
    parity ^= (unsigned)((UINT64_C(0x6996) >> value) & UINT64_C(1));
#endif
  }
  return parity & 1U;
}

static int flat_rows_are_rref(const uint64_t *rows, uint32_t row_count,
                              uint32_t width, uint32_t words) {
  uint32_t previous_pivot = 0U;
  unsigned have_previous = 0U;
  uint32_t row;
  for (row = 0U; row < row_count; ++row) {
    const uint64_t *current = rows + (uint64_t)row * words;
    uint32_t pivot = width;
    uint32_t column;
    uint32_t other;
    for (column = 0U; column < width; ++column) {
      if (flat_word_bit(current, column) != 0U) {
        pivot = column;
        break;
      }
    }
    if (pivot == width ||
        (have_previous != 0U && pivot <= previous_pivot)) {
      return 0;
    }
    for (other = 0U; other < row_count; ++other) {
      if (other != row &&
          flat_word_bit(rows + (uint64_t)other * words, pivot) != 0U) {
        return 0;
      }
    }
    previous_pivot = pivot;
    have_previous = 1U;
  }
  return 1;
}

static int flat_dcp_matrices_are_symmetric(const unsigned char *bytes,
                                           uint64_t offset, uint32_t count,
                                           uint32_t width, uint32_t words) {
  const uint64_t *matrices = (const uint64_t *)(bytes + offset);
  uint32_t index;
  for (index = 0U; index < count; ++index) {
    const uint64_t *matrix =
        matrices + (uint64_t)index * width * words;
    uint32_t row;
    for (row = 0U; row < width; ++row) {
      uint32_t column;
      for (column = row + 1U; column < width; ++column) {
        if (flat_word_bit(matrix + (uint64_t)row * words, column) !=
            flat_word_bit(matrix + (uint64_t)column * words, row)) {
          return 0;
        }
      }
    }
  }
  return 1;
}

comb_cuda_status_t comb_cuda_flat_problem_validate(const void *blob,
                                                   uint64_t blob_size) {
  if (blob == NULL || ((uintptr_t)blob & (FLAT_ALIGNMENT - 1)) != 0 ||
      blob_size < sizeof(comb_cuda_flat_problem_t) ||
      blob_size > (uint64_t)SIZE_MAX) {
    return COMB_CUDA_STATUS_INVALID_ARGUMENT;
  }
  const unsigned char *bytes = (const unsigned char *)blob;
  const comb_cuda_flat_problem_t *header =
      (const comb_cuda_flat_problem_t *)blob;
  if (header->magic != COMB_CUDA_FLAT_MAGIC ||
      header->version != COMB_CUDA_FLAT_VERSION ||
      header->header_size != sizeof(comb_cuda_flat_problem_t) ||
      header->blob_size != blob_size || header->n_all == 0 || header->n == 0 ||
      header->n > header->n_all || header->n_x > header->n ||
      header->n_z > header->n || header->k > header->n_all ||
      header->fault_location_count == 0 ||
      header->fault_location_count == UINT32_MAX || header->s == 0 ||
      header->s == UINT32_MAX || header->s > header->fault_location_count ||
      header->rank_dp_rows != header->fault_location_count + 1 ||
      header->rank_dp_columns != header->s + 1 ||
      header->output_error_count == 0 ||
      header->words_all != dense_word_count(header->n_all) ||
      header->words_code != dense_word_count(header->n)) {
    return COMB_CUDA_STATUS_INVALID_PROBLEM;
  }
  const uint32_t x_lookup_exponent = header->n - header->n_x;
  const uint32_t z_lookup_exponent = header->n - header->n_z;
  if (x_lookup_exponent >= 64 || z_lookup_exponent >= 64) {
    return COMB_CUDA_STATUS_PROBLEM_TOO_LARGE;
  }
  if (header->x_lookup_bit_count != (UINT64_C(1) << x_lookup_exponent) ||
      header->z_lookup_bit_count != (UINT64_C(1) << z_lookup_exponent)) {
    return COMB_CUDA_STATUS_INVALID_PROBLEM;
  }

  uint64_t cursor = sizeof(comb_cuda_flat_problem_t);
  uint64_t rank_dp_elements = 0;
  uint64_t x_stabiliser_words = 0;
  uint64_t z_stabiliser_words = 0;
  uint64_t logical_matrix_rows = 0;
  uint64_t logical_matrix_words = 0;
  uint64_t logical_vector_words = 0;
  uint64_t output_matrix_rows = 0;
  uint64_t output_matrix_words = 0;
  uint64_t output_vector_words = 0;
  uint64_t x_measurement_words = 0;
  uint64_t z_measurement_words = 0;
  if (!multiply_u64(header->rank_dp_rows, header->rank_dp_columns,
                    &rank_dp_elements) ||
      !multiply_u64(header->n_x, header->words_code, &x_stabiliser_words) ||
      !multiply_u64(header->n_z, header->words_code, &z_stabiliser_words) ||
      !multiply_u64(header->k, header->n_all, &logical_matrix_rows) ||
      !multiply_u64(logical_matrix_rows, header->words_all,
                    &logical_matrix_words) ||
      !multiply_u64(header->k, header->words_all, &logical_vector_words) ||
      !multiply_u64(header->output_error_count, header->n_all,
                    &output_matrix_rows) ||
      !multiply_u64(output_matrix_rows, header->words_all,
                    &output_matrix_words) ||
      !multiply_u64(header->output_error_count, header->words_all,
                    &output_vector_words) ||
      !multiply_u64(header->x_measurement_count, header->words_all,
                    &x_measurement_words) ||
      !multiply_u64(header->z_measurement_count, header->words_all,
                    &z_measurement_words)) {
    return COMB_CUDA_STATUS_PROBLEM_TOO_LARGE;
  }
  const uint64_t x_lookup_words = (header->x_lookup_bit_count + 63) / 64;
  const uint64_t z_lookup_words = (header->z_lookup_bit_count + 63) / 64;
  if (!consume_section(header->fault_partitions_offset,
                       header->fault_location_count, sizeof(uint32_t),
                       blob_size, &cursor) ||
      !consume_section(header->fault_prefixes_offset,
                       (uint64_t)header->fault_location_count + 1,
                       sizeof(uint64_t), blob_size, &cursor) ||
      !consume_section(header->rank_dp_offset, rank_dp_elements,
                       sizeof(uint64_t), blob_size, &cursor) ||
      !consume_section(header->x_stabilisers_offset, x_stabiliser_words,
                       sizeof(uint64_t), blob_size, &cursor) ||
      !consume_section(header->z_stabilisers_offset, z_stabiliser_words,
                       sizeof(uint64_t), blob_size, &cursor) ||
      !consume_section(header->x_lookup_offset, x_lookup_words,
                       sizeof(uint64_t), blob_size, &cursor) ||
      !consume_section(header->z_lookup_offset, z_lookup_words,
                       sizeof(uint64_t), blob_size, &cursor) ||
      !consume_section(header->logical_dcp_f_offset, logical_matrix_words,
                       sizeof(uint64_t), blob_size, &cursor) ||
      !consume_section(header->logical_dcp_hx_offset, logical_vector_words,
                       sizeof(uint64_t), blob_size, &cursor) ||
      !consume_section(header->logical_dcp_hz_offset, logical_vector_words,
                       sizeof(uint64_t), blob_size, &cursor) ||
      !consume_section(header->output_dcp_f_offset, output_matrix_words,
                       sizeof(uint64_t), blob_size, &cursor) ||
      !consume_section(header->output_dcp_hx_offset, output_vector_words,
                       sizeof(uint64_t), blob_size, &cursor) ||
      !consume_section(header->output_dcp_hz_offset, output_vector_words,
                       sizeof(uint64_t), blob_size, &cursor) ||
      !consume_section(header->x_measurements_offset, x_measurement_words,
                       sizeof(uint64_t), blob_size, &cursor) ||
      !consume_section(header->z_measurements_offset, z_measurement_words,
                       sizeof(uint64_t), blob_size, &cursor) ||
      cursor != blob_size) {
    return COMB_CUDA_STATUS_INVALID_PROBLEM;
  }

  if (!rows_tail_masked(bytes, header->x_stabilisers_offset, header->n_x,
                        header->n) ||
      !rows_tail_masked(bytes, header->z_stabilisers_offset, header->n_z,
                        header->n) ||
      !rows_tail_masked(bytes, header->logical_dcp_f_offset,
                        logical_matrix_rows, header->n_all) ||
      !rows_tail_masked(bytes, header->logical_dcp_hx_offset, header->k,
                        header->n_all) ||
      !rows_tail_masked(bytes, header->logical_dcp_hz_offset, header->k,
                        header->n_all) ||
      !rows_tail_masked(bytes, header->output_dcp_f_offset,
                        output_matrix_rows, header->n_all) ||
      !rows_tail_masked(bytes, header->output_dcp_hx_offset,
                        header->output_error_count, header->n_all) ||
      !rows_tail_masked(bytes, header->output_dcp_hz_offset,
                        header->output_error_count, header->n_all) ||
      !rows_tail_masked(bytes, header->x_measurements_offset,
                        header->x_measurement_count, header->n_all) ||
      !rows_tail_masked(bytes, header->z_measurements_offset,
                        header->z_measurement_count, header->n_all)) {
    return COMB_CUDA_STATUS_INVALID_PROBLEM;
  }
  const uint64_t *x_lookup =
      (const uint64_t *)(bytes + header->x_lookup_offset);
  const uint64_t *z_lookup =
      (const uint64_t *)(bytes + header->z_lookup_offset);
  if ((x_lookup[x_lookup_words - 1] &
       ~tail_mask_u64(header->x_lookup_bit_count)) != 0 ||
      (z_lookup[z_lookup_words - 1] &
       ~tail_mask_u64(header->z_lookup_bit_count)) != 0) {
    return COMB_CUDA_STATUS_INVALID_PROBLEM;
  }

  const uint32_t *partitions =
      (const uint32_t *)(bytes + header->fault_partitions_offset);
  const uint64_t *prefixes =
      (const uint64_t *)(bytes + header->fault_prefixes_offset);
  const uint64_t *dp = (const uint64_t *)(bytes + header->rank_dp_offset);
  if (prefixes[0] != 0 || dp[0] != 1) {
    return COMB_CUDA_STATUS_INVALID_PROBLEM;
  }
  for (uint32_t r = 1; r <= header->s; ++r) {
    if (dp[r] != 0) return COMB_CUDA_STATUS_INVALID_PROBLEM;
  }
  for (uint32_t m = 1; m <= header->fault_location_count; ++m) {
    if (partitions[m - 1] == 0 ||
        prefixes[m - 1] > UINT64_MAX - partitions[m - 1] ||
        prefixes[m] != prefixes[m - 1] + partitions[m - 1]) {
      return COMB_CUDA_STATUS_INVALID_PROBLEM;
    }
    const uint64_t *previous = dp + (uint64_t)(m - 1) * header->rank_dp_columns;
    const uint64_t *current = dp + (uint64_t)m * header->rank_dp_columns;
    if (current[0] != 1) return COMB_CUDA_STATUS_INVALID_PROBLEM;
    for (uint32_t r = 1; r <= header->s; ++r) {
      uint64_t expected = 0;
      if (r <= m) {
        uint64_t term = 0;
        if (!multiply_u64(previous[r - 1], partitions[m - 1], &term) ||
            !add_u64(previous[r], term, &expected)) {
          return COMB_CUDA_STATUS_PROBLEM_TOO_LARGE;
        }
      }
      if (current[r] != expected) return COMB_CUDA_STATUS_INVALID_PROBLEM;
    }
  }
  if (prefixes[header->fault_location_count] > header->output_error_count ||
      dp[(uint64_t)header->fault_location_count * header->rank_dp_columns +
         header->s] != header->total_combination_count) {
    return COMB_CUDA_STATUS_INVALID_PROBLEM;
  }
  return COMB_CUDA_STATUS_SUCCESS;
}

comb_cuda_status_t comb_cuda_flat_problem_validate_science(
    const void *blob, uint64_t blob_size) {
  comb_cuda_status_t status =
      comb_cuda_flat_problem_validate(blob, blob_size);
  const unsigned char *bytes;
  const comb_cuda_flat_problem_t *header;
  const uint64_t *x_rows;
  const uint64_t *z_rows;
  const uint64_t *z_measurements;
  uint32_t x;
  uint32_t z;
  if (status != COMB_CUDA_STATUS_SUCCESS) return status;

  bytes = (const unsigned char *)blob;
  header = (const comb_cuda_flat_problem_t *)blob;
  if (header->k >= 64U) return COMB_CUDA_STATUS_PROBLEM_TOO_LARGE;
  if (header->n_x > header->n - header->n_z ||
      header->k != header->n - header->n_x - header->n_z ||
      header->z_measurement_count > header->n_all - header->n_z) {
    return COMB_CUDA_STATUS_INVALID_PROBLEM;
  }

  x_rows = (const uint64_t *)(bytes + header->x_stabilisers_offset);
  z_rows = (const uint64_t *)(bytes + header->z_stabilisers_offset);
  if (!flat_rows_are_rref(x_rows, header->n_x, header->n,
                          header->words_code) ||
      !flat_rows_are_rref(z_rows, header->n_z, header->n,
                          header->words_code)) {
    return COMB_CUDA_STATUS_INVALID_PROBLEM;
  }
  for (x = 0U; x < header->n_x; ++x) {
    for (z = 0U; z < header->n_z; ++z) {
      if (flat_row_dot(x_rows + (uint64_t)x * header->words_code,
                       z_rows + (uint64_t)z * header->words_code,
                       header->words_code) != 0U) {
        return COMB_CUDA_STATUS_INVALID_PROBLEM;
      }
    }
  }

  if (!flat_dcp_matrices_are_symmetric(
          bytes, header->logical_dcp_f_offset, header->k, header->n_all,
          header->words_all) ||
      !flat_dcp_matrices_are_symmetric(
          bytes, header->output_dcp_f_offset, header->output_error_count,
          header->n_all, header->words_all)) {
    return COMB_CUDA_STATUS_INVALID_PROBLEM;
  }

  /* The retained batched-Z construction pairs C_z = [M; S_z] with
   * C_x = [M; E_pivot].  Require C_z C_x^T = I explicitly: the measurement
   * Gram matrix is identity, measurements commute with S_z, and no
   * measurement occupies an S_z pivot.  This rejects the unsupported general
   * parity inputs called out by the legacy implementation. */
  z_measurements =
      (const uint64_t *)(bytes + header->z_measurements_offset);
  for (z = 0U; z < header->z_measurement_count; ++z) {
    const uint64_t *measurement =
        z_measurements + (uint64_t)z * header->words_all;
    uint32_t other;
    for (other = 0U; other < header->z_measurement_count; ++other) {
      const uint64_t *other_measurement =
          z_measurements + (uint64_t)other * header->words_all;
      const unsigned expected = other == z ? 1U : 0U;
      if (flat_row_dot(measurement, other_measurement, header->words_all) !=
          expected) {
        return COMB_CUDA_STATUS_INVALID_PROBLEM;
      }
    }
    for (other = 0U; other < header->n_z; ++other) {
      const uint64_t *stabiliser =
          z_rows + (uint64_t)other * header->words_code;
      uint32_t pivot = 0U;
      while (pivot < header->n && flat_word_bit(stabiliser, pivot) == 0U) {
        ++pivot;
      }
      if (pivot == header->n ||
          flat_row_dot(measurement, stabiliser, header->words_code) != 0U ||
          flat_word_bit(measurement, pivot) != 0U) {
        return COMB_CUDA_STATUS_INVALID_PROBLEM;
      }
    }
  }
  return COMB_CUDA_STATUS_SUCCESS;
}

void comb_cuda_flat_problem_free(void *blob) { free(blob); }
