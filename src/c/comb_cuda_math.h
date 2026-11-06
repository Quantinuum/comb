#ifndef COMB_CUDA_MATH_H
#define COMB_CUDA_MATH_H

/* Allocation-free packed stabiliser operations shared by CUDA and host
 * differential tests. One thread owns each view and its scratch storage.
 * The caller validates dimensions and zero tail padding before execution.
 * Tableau storage has at least qubits + 1 rows; all strides equal words.
 */
#include <stdint.h>

#if defined(__CUDACC__) || defined(__CUDA__)
#define COMB_GPU_HD __host__ __device__
#else
#define COMB_GPU_HD
#endif

typedef struct {
  uint32_t qubits, words, rows;
  uint64_t *x, *z;
  uint8_t *delta, *epsilon;
} comb_gpu_tableau_t;

typedef struct {
  uint32_t qubits, words;
  uint64_t *f, *hx, *hz;
} comb_gpu_dcp_t;

COMB_GPU_HD static inline unsigned comb_gpu_bit(const uint64_t *a,
                                                uint64_t bit) {
  return (unsigned)((a[bit / 64U] >> (bit % 64U)) & UINT64_C(1));
}

COMB_GPU_HD static inline void comb_gpu_set(uint64_t *a, uint32_t bit,
                                            unsigned value) {
  const uint64_t mask = UINT64_C(1) << (bit % 64U);
  if (value != 0U)
    a[bit / 64U] |= mask;
  else
    a[bit / 64U] &= ~mask;
}

COMB_GPU_HD static inline unsigned comb_gpu_parity(uint64_t a) {
  a ^= a >> 32U;
  a ^= a >> 16U;
  a ^= a >> 8U;
  a ^= a >> 4U;
  return (unsigned)((UINT64_C(0x6996) >> (a & UINT64_C(15))) & 1U);
}

COMB_GPU_HD static inline unsigned comb_gpu_dot(const uint64_t *a,
                                                const uint64_t *b,
                                                uint32_t words) {
  unsigned parity = 0U;
  for (uint32_t w = 0U; w < words; ++w) parity ^= comb_gpu_parity(a[w] & b[w]);
  return parity;
}

COMB_GPU_HD static inline void comb_gpu_clear_row(comb_gpu_tableau_t *t,
                                                  uint32_t row) {
  for (uint32_t w = 0U; w < t->words; ++w) {
    t->x[(uint64_t)row * t->words + w] = 0U;
    t->z[(uint64_t)row * t->words + w] = 0U;
  }
  t->delta[row] = 0U;
  t->epsilon[row] = 0U;
}

COMB_GPU_HD static inline void comb_gpu_row_add(comb_gpu_tableau_t *t,
                                                uint32_t source,
                                                uint32_t destination) {
  const uint64_t s = (uint64_t)source * t->words;
  const uint64_t d = (uint64_t)destination * t->words;
  const uint8_t sd = t->delta[source], dd = t->delta[destination];
  const uint8_t se = t->epsilon[source], de = t->epsilon[destination];
  const unsigned phase = comb_gpu_dot(t->x + s, t->z + d, t->words);
  for (uint32_t w = 0U; w < t->words; ++w) {
    t->x[d + w] ^= t->x[s + w];
    t->z[d + w] ^= t->z[s + w];
  }
  t->delta[destination] = (uint8_t)(sd ^ dd);
  t->epsilon[destination] = (uint8_t)(se ^ de ^ (sd & dd) ^ phase);
}

COMB_GPU_HD static inline void comb_gpu_row_swap(comb_gpu_tableau_t *t,
                                                 uint32_t a, uint32_t b) {
  for (uint32_t w = 0U; w < t->words; ++w) {
    uint64_t tmp = t->x[(uint64_t)a * t->words + w];
    t->x[(uint64_t)a * t->words + w] = t->x[(uint64_t)b * t->words + w];
    t->x[(uint64_t)b * t->words + w] = tmp;
    tmp = t->z[(uint64_t)a * t->words + w];
    t->z[(uint64_t)a * t->words + w] = t->z[(uint64_t)b * t->words + w];
    t->z[(uint64_t)b * t->words + w] = tmp;
  }
  uint8_t tmp = t->delta[a];
  t->delta[a] = t->delta[b];
  t->delta[b] = tmp;
  tmp = t->epsilon[a];
  t->epsilon[a] = t->epsilon[b];
  t->epsilon[b] = tmp;
}

/* 0 success, 1 not generated, 2 inconsistent. Destructive, matching dense GE.
 */
COMB_GPU_HD static inline int comb_gpu_form_identity(comb_gpu_tableau_t *t,
                                                     unsigned *negative) {
  uint32_t start = 0U;
  for (uint32_t c = 0U; c < t->qubits; ++c) {
    for (unsigned part = 0U; part < 2U; ++part) {
      const uint64_t *matrix = part == 0U ? t->x : t->z;
      uint32_t found = t->rows;
      for (uint32_t r = start; r < t->rows; ++r) {
        if (comb_gpu_bit(matrix + (uint64_t)r * t->words, c) != 0U) {
          if (found == t->rows)
            found = r;
          else
            comb_gpu_row_add(t, found, r);
        }
      }
      if (found != t->rows) comb_gpu_row_swap(t, found, start++);
    }
  }
  for (uint32_t w = 0U; w < t->words; ++w)
    if (t->x[(uint64_t)t->qubits * t->words + w] != 0U ||
        t->z[(uint64_t)t->qubits * t->words + w] != 0U)
      return 1;
  if (t->delta[t->qubits] != 0U) return 2;
  *negative = t->epsilon[t->qubits];
  return 0;
}

/* basis 0 = X, 1 = Z; measurements are packed rows of t->words words. */
COMB_GPU_HD static inline int comb_gpu_postselect(comb_gpu_tableau_t *t,
                                                  const uint64_t *measurements,
                                                  uint32_t count,
                                                  unsigned basis,
                                                  unsigned *rejected) {
  *rejected = 0U;
  for (uint32_t m = 0U; m < count; ++m) {
    const uint64_t *measurement = measurements + (uint64_t)m * t->words;
    const uint64_t *opposite = basis == 0U ? t->z : t->x;
    uint32_t pivot = t->qubits;
    for (uint32_t r = 0U; r < t->qubits; ++r)
      if (comb_gpu_dot(measurement, opposite + (uint64_t)r * t->words,
                       t->words) != 0U) {
        pivot = r;
        break;
      }
    if (pivot != t->qubits) {
      for (uint32_t r = pivot + 1U; r < t->qubits; ++r)
        if (comb_gpu_dot(measurement, opposite + (uint64_t)r * t->words,
                         t->words) != 0U)
          comb_gpu_row_add(t, pivot, r);
    }
    comb_gpu_clear_row(t, pivot);
    uint64_t *part = basis == 0U ? t->x : t->z;
    for (uint32_t w = 0U; w < t->words; ++w)
      part[(uint64_t)pivot * t->words + w] = measurement[w];
    if (pivot == t->qubits) {
      unsigned negative = 0U;
      t->rows = t->qubits + 1U;
      const int status = comb_gpu_form_identity(t, &negative);
      t->rows = t->qubits;
      if (status != 0) return status;
      if (negative != 0U) {
        *rejected = 1U;
        return 0;
      }
    }
  }
  return 0;
}

COMB_GPU_HD static inline uint32_t comb_gpu_reduce_coset(
    comb_gpu_tableau_t *t) {
  uint32_t start = 0U, z_rank = 0U;
  for (unsigned part = 0U; part < 2U; ++part) {
    const uint64_t *matrix = part == 0U ? t->z : t->x;
    for (uint32_t c = 0U; c < t->qubits; ++c) {
      uint32_t found = t->rows;
      for (uint32_t r = start; r < t->rows; ++r) {
        if (comb_gpu_bit(matrix + (uint64_t)r * t->words, c) != 0U) {
          if (found == t->rows)
            found = r;
          else
            comb_gpu_row_add(t, found, r);
        }
      }
      if (found != t->rows) comb_gpu_row_swap(t, found, start++);
    }
    if (part == 0U) z_rank = start;
  }
  for (uint32_t r = 0U; r < t->rows; ++r) {
    uint32_t leading = t->qubits;
    for (uint32_t c = 0U; c < t->qubits; ++c)
      if (comb_gpu_bit(t->x + (uint64_t)r * t->words, c) != 0U) {
        leading = c;
        break;
      }
    if (leading != t->qubits)
      for (uint32_t above = r; above != 0U;) {
        --above;
        if (comb_gpu_bit(t->x + (uint64_t)above * t->words, leading) != 0U)
          comb_gpu_row_add(t, r, above);
      }
  }
  return z_rank;
}

/* Output must not alias either input. Compute c2 * c1, preserving phase order.
 */
COMB_GPU_HD static inline void comb_gpu_dcp_multiply(comb_gpu_dcp_t *out,
                                                     const comb_gpu_dcp_t *c2,
                                                     const comb_gpu_dcp_t *c1) {
  for (uint32_t w = 0U; w < out->words; ++w) {
    out->hx[w] = c2->hx[w] ^ c1->hx[w];
    out->hz[w] = c2->hz[w] ^ c1->hz[w];
  }
  for (uint32_t r = 0U; r < out->qubits; ++r) {
    const uint64_t offset = (uint64_t)r * out->words;
    const unsigned phase =
        comb_gpu_dot(c1->f + offset, c2->hx, out->words) ^
        (comb_gpu_bit(c1->f + offset, r) & comb_gpu_bit(c2->f + offset, r));
    if (phase != 0U) comb_gpu_set(out->hz, r, comb_gpu_bit(out->hz, r) ^ 1U);
    for (uint32_t w = 0U; w < out->words; ++w)
      out->f[offset + w] = c2->f[offset + w] ^ c1->f[offset + w];
  }
}

COMB_GPU_HD static inline void comb_gpu_tableau_from_dcp(
    comb_gpu_tableau_t *t, const comb_gpu_dcp_t *dcp) {
  t->rows = t->qubits;
  for (uint32_t r = 0U; r <= t->qubits; ++r) comb_gpu_clear_row(t, r);
  for (uint32_t r = 0U; r < t->qubits; ++r) {
    const uint64_t offset = (uint64_t)r * t->words;
    for (uint32_t w = 0U; w < t->words; ++w)
      t->z[offset + w] = dcp->f[offset + w];
    comb_gpu_set(t->x + offset, r, 1U);
    t->delta[r] = (uint8_t)comb_gpu_bit(dcp->f + offset, r);
    t->epsilon[r] = (uint8_t)comb_gpu_bit(dcp->hz, r);
  }
}

COMB_GPU_HD static inline uint64_t comb_gpu_compress(const uint64_t *row,
                                                     const uint32_t *columns,
                                                     uint32_t count) {
  uint64_t result = 0U;
  for (uint32_t c = 0U; c < count; ++c)
    result |= (uint64_t)comb_gpu_bit(row, columns[c]) << c;
  return result;
}

/* Reduce then check every member of the Z coset. generators has qubits slots.
 * Returns 0 success, 2 invalid scientific state, 3 unrepresentable coset.
 */
COMB_GPU_HD static inline int comb_gpu_coset_good(
    comb_gpu_tableau_t *t, const uint32_t *free_columns, uint32_t free_count,
    const uint64_t *lookup, uint64_t lookup_bits, uint64_t *generators,
    unsigned *good) {
  *good = 0U;
  const uint32_t rank = comb_gpu_reduce_coset(t);
  if (rank >= 64U || free_count >= 64U) return 3;
  for (uint32_t r = 0U; r < rank; ++r) {
    const uint64_t *row = t->z + (uint64_t)r * t->words;
    for (uint32_t c = 0U; c < t->qubits; ++c) {
      if (comb_gpu_bit(row, c) != 0U) {
        unsigned found = 0U;
        for (uint32_t j = 0U; j < free_count; ++j)
          if (free_columns[j] == c) {
            found = 1U;
            break;
          }
        if (found == 0U) return 2;
      }
    }
    generators[r] = comb_gpu_compress(row, free_columns, free_count);
  }
  uint64_t factor = 0U;
  for (uint32_t r = rank; r < t->qubits; ++r) {
    if (t->epsilon[r] == 0U) continue;
    uint32_t leading = t->qubits;
    for (uint32_t c = 0U; c < t->qubits; ++c)
      if (comb_gpu_bit(t->x + (uint64_t)r * t->words, c) != 0U) {
        leading = c;
        break;
      }
    unsigned found = 0U;
    for (uint32_t j = 0U; j < free_count; ++j)
      if (free_columns[j] == leading) {
        factor |= UINT64_C(1) << j;
        found = 1U;
        break;
      }
    if (found == 0U) return 2;
  }
  /* Gray-code enumeration changes one generator at a time. */
  uint64_t reduced = factor;
  const uint64_t count = UINT64_C(1) << rank;
  for (uint64_t i = 0U; i < count; ++i) {
    if (reduced >= lookup_bits) return 2;
    if (((lookup[reduced / 64U] >> (reduced % 64U)) & 1U) == 0U) return 0;
    if (i + 1U < count) {
      uint64_t next = i + 1U;
      uint32_t changed = 0U;
      while ((next & 1U) == 0U) {
        ++changed;
        next >>= 1U;
      }
      reduced ^= generators[changed];
    }
  }
  *good = 1U;
  return 0;
}

#endif /* COMB_CUDA_MATH_H */
