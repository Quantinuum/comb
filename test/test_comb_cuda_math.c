#include "comb_cuda_math.h"
#include "comb_dense.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define CHECK(value) do { if (!(value)) { \
  fprintf(stderr, "%s:%d: %s\n", __FILE__, __LINE__, #value); exit(1); \
} } while (0)

static uint64_t random_state = UINT64_C(0xda942042e4dd58b5);
static uint64_t random_word(void) {
  random_state ^= random_state >> 12U;
  random_state ^= random_state << 25U;
  random_state ^= random_state >> 27U;
  return random_state * UINT64_C(2685821657736338717);
}

static comb_gpu_tableau_t make_gpu_tableau(uint32_t n) {
  comb_gpu_tableau_t t;
  t.qubits = n; t.words = (n + 63U) / 64U; t.rows = n;
  t.x = (uint64_t *)calloc((size_t)(n + 1U) * t.words, sizeof(uint64_t));
  t.z = (uint64_t *)calloc((size_t)(n + 1U) * t.words, sizeof(uint64_t));
  t.delta = (uint8_t *)calloc(n + 1U, 1U);
  t.epsilon = (uint8_t *)calloc(n + 1U, 1U);
  CHECK(t.x && t.z && t.delta && t.epsilon);
  return t;
}

static void free_gpu_tableau(comb_gpu_tableau_t *t) {
  free(t->x); free(t->z); free(t->delta); free(t->epsilon);
}

static comb_gpu_dcp_t make_gpu_dcp(uint32_t n) {
  comb_gpu_dcp_t d;
  d.qubits = n; d.words = (n + 63U) / 64U;
  d.f = (uint64_t *)calloc((size_t)n * d.words, sizeof(uint64_t));
  d.hx = (uint64_t *)calloc(d.words, sizeof(uint64_t));
  d.hz = (uint64_t *)calloc(d.words, sizeof(uint64_t));
  CHECK(d.f && d.hx && d.hz);
  return d;
}

static void free_gpu_dcp(comb_gpu_dcp_t *d) {
  free(d->f); free(d->hx); free(d->hz);
}

static comb_gpu_dcp_t dense_dcp_view(comb_dense_dcp_t *d) {
  comb_gpu_dcp_t g;
  g.qubits = (uint32_t)d->qubits; g.words = (uint32_t)d->word_count;
  g.f = d->f.words; g.hx = d->h_x; g.hz = d->h_z;
  return g;
}

static void random_dcp(comb_dense_dcp_t *d) {
  for (uint32_t r = 0U; r < d->qubits; ++r) {
    for (uint32_t c = r; c < d->qubits; ++c) {
      const unsigned value = (unsigned)(random_word() & 1U);
      CHECK(comb_dense_matrix_set(&d->f, r, c, value) == COMB_DENSE_OK);
      CHECK(comb_dense_matrix_set(&d->f, c, r, value) == COMB_DENSE_OK);
    }
    comb_gpu_set(d->h_x, r, (unsigned)(random_word() & 1U));
    comb_gpu_set(d->h_z, r, (unsigned)(random_word() & 1U));
  }
}

static void dense_from_gpu(const comb_gpu_tableau_t *g,
                           comb_dense_tableau_t *d) {
  const size_t words = (size_t)(g->qubits + 1U) * g->words;
  d->row_count = g->rows;
  memcpy(d->x.words, g->x, words * sizeof(uint64_t));
  memcpy(d->z.words, g->z, words * sizeof(uint64_t));
  memcpy(d->delta, g->delta, g->qubits + 1U);
  memcpy(d->epsilon, g->epsilon, g->qubits + 1U);
}

static void compare_tableaux(const comb_gpu_tableau_t *g,
                             const comb_dense_tableau_t *d) {
  const size_t bytes = (size_t)(g->qubits + 1U) * g->words * sizeof(uint64_t);
  CHECK(g->rows == d->row_count);
  CHECK(memcmp(g->x, d->x.words, bytes) == 0);
  CHECK(memcmp(g->z, d->z.words, bytes) == 0);
  CHECK(memcmp(g->delta, d->delta, g->qubits + 1U) == 0);
  CHECK(memcmp(g->epsilon, d->epsilon, g->qubits + 1U) == 0);
}

static void test_width(uint32_t n) {
  comb_dense_dcp_t a, b, product;
  comb_dense_tableau_t dense;
  CHECK(comb_dense_dcp_init(&a, n) == COMB_DENSE_OK);
  CHECK(comb_dense_dcp_init(&b, n) == COMB_DENSE_OK);
  CHECK(comb_dense_dcp_init(&product, n) == COMB_DENSE_OK);
  CHECK(comb_dense_tableau_init(&dense, n, n + 1U, n) == COMB_DENSE_OK);
  comb_gpu_dcp_t gpu_product = make_gpu_dcp(n);
  comb_gpu_tableau_t gpu = make_gpu_tableau(n);
  comb_dense_matrix_t measurements;
  CHECK(comb_dense_matrix_init(&measurements, 4U, n) == COMB_DENSE_OK);

  for (unsigned iteration = 0U; iteration < 16U; ++iteration) {
    random_dcp(&a); random_dcp(&b);
    comb_gpu_dcp_t av = dense_dcp_view(&a), bv = dense_dcp_view(&b);
    comb_gpu_dcp_multiply(&gpu_product, &av, &bv);
    CHECK(comb_dense_dcp_multiply(&product, &a, &b) == COMB_DENSE_OK);
    CHECK(memcmp(gpu_product.f, product.f.words,
                 (size_t)n * av.words * sizeof(uint64_t)) == 0);
    CHECK(memcmp(gpu_product.hx, product.h_x,
                 (size_t)av.words * sizeof(uint64_t)) == 0);
    CHECK(memcmp(gpu_product.hz, product.h_z,
                 (size_t)av.words * sizeof(uint64_t)) == 0);

    comb_gpu_tableau_from_dcp(&gpu, &gpu_product);
    dense_from_gpu(&gpu, &dense);
    for (unsigned operation = 0U; operation < 32U; ++operation) {
      const uint32_t source = (uint32_t)(random_word() % n);
      const uint32_t destination = (uint32_t)(random_word() % n);
      if (source == destination) continue;
      if ((operation & 1U) != 0U) {
        comb_gpu_row_add(&gpu, source, destination);
        CHECK(comb_dense_tableau_row_add(&dense, source, &dense,
                                         destination) == COMB_DENSE_OK);
      } else {
        comb_gpu_row_swap(&gpu, source, destination);
        CHECK(comb_dense_tableau_row_swap(&dense, source,
                                          destination) == COMB_DENSE_OK);
      }
      compare_tableaux(&gpu, &dense);
    }
    /* An appended stabiliser product must reduce to signed identity. */
    comb_gpu_clear_row(&gpu, n); gpu.rows = n + 1U;
    for (uint32_t r = 0U; r < n; ++r)
      if ((random_word() & 1U) != 0U) comb_gpu_row_add(&gpu, r, n);
    const unsigned expected_negative = iteration & 1U;
    gpu.epsilon[n] ^= (uint8_t)expected_negative;
    dense_from_gpu(&gpu, &dense);
    unsigned gpu_negative = 2U, dense_negative = 2U;
    CHECK(comb_gpu_form_identity(&gpu, &gpu_negative) == 0);
    CHECK(comb_dense_tableau_form_identity(&dense,
                                          &dense_negative) == COMB_DENSE_OK);
    CHECK(gpu_negative == expected_negative && gpu_negative == dense_negative);
    compare_tableaux(&gpu, &dense);

    /* Random sequential measurements exercise pivot replacement and forced
     * outcomes, including destructive GE for deterministic measurements. */
    comb_gpu_tableau_from_dcp(&gpu, &gpu_product);
    dense_from_gpu(&gpu, &dense);
    for (unsigned basis = 0U; basis < 2U; ++basis) {
      for (uint32_t r = 0U; r < 4U; ++r)
        for (uint32_t c = 0U; c < n; ++c)
          CHECK(comb_dense_matrix_set(&measurements, r, c,
                   (unsigned)(random_word() & 1U)) == COMB_DENSE_OK);
      unsigned gpu_rejected = 2U, dense_rejected = 2U;
      CHECK(comb_gpu_postselect(&gpu, measurements.words, 4U, basis,
                                &gpu_rejected) == 0);
      CHECK(comb_dense_tableau_postselect(&dense, &measurements,
          basis == 0U ? COMB_DENSE_MEASURE_X : COMB_DENSE_MEASURE_Z,
          &dense_rejected) == COMB_DENSE_OK);
      CHECK(gpu_rejected == dense_rejected);
      compare_tableaux(&gpu, &dense);
      if (gpu_rejected != 0U) break;
    }

    const uint32_t gpu_rank = comb_gpu_reduce_coset(&gpu);
    size_t dense_rank = SIZE_MAX;
    CHECK(comb_dense_tableau_reduce_coset(&dense, &dense_rank) == COMB_DENSE_OK);
    CHECK(gpu_rank == dense_rank);
    compare_tableaux(&gpu, &dense);
  }
  comb_dense_matrix_destroy(&measurements);
  comb_dense_tableau_destroy(&dense);
  comb_dense_dcp_destroy(&a); comb_dense_dcp_destroy(&b);
  comb_dense_dcp_destroy(&product);
  free_gpu_dcp(&gpu_product); free_gpu_tableau(&gpu);
}

static void test_universal_coset(void) {
  const uint32_t columns[3] = {0U, 1U, 2U};
  comb_gpu_tableau_t t = make_gpu_tableau(65U);
  uint64_t generators[65];
  for (unsigned phase = 0U; phase < 2U; ++phase) {
    for (unsigned removed = 0U; removed <= 4U; ++removed) {
      for (uint32_t r = 0U; r <= t.qubits; ++r) comb_gpu_clear_row(&t, r);
      for (uint32_t r = 0U; r < t.qubits; ++r)
        comb_gpu_set((r < 2U ? t.z : t.x) + (uint64_t)r * t.words, r, 1U);
      t.epsilon[2] = (uint8_t)phase;
      /* Scramble the basis without changing its stabiliser group. */
      for (unsigned i = 0U; i < 64U; ++i) {
        const uint32_t a = (uint32_t)(random_word() % t.qubits);
        const uint32_t b = (uint32_t)(random_word() % t.qubits);
        if (a != b) comb_gpu_row_add(&t, a, b);
        comb_gpu_row_swap(&t, a, b);
      }
      uint64_t lookup = UINT64_C(0xff);
      if (removed < 4U) lookup &= ~(UINT64_C(1) << (phase * 4U + removed));
      unsigned good = 2U;
      CHECK(comb_gpu_coset_good(&t, columns, 3U, &lookup, 8U,
                                generators, &good) == 0);
      CHECK(good == (removed == 4U ? 1U : 0U));
    }
  }
  /* A coset with 64 independent generators cannot be enumerated in uint64. */
  for (uint32_t r = 0U; r <= t.qubits; ++r) comb_gpu_clear_row(&t, r);
  for (uint32_t r = 0U; r < t.qubits; ++r)
    comb_gpu_set(t.z + (uint64_t)r * t.words, r, 1U);
  uint64_t lookup = UINT64_C(0xff);
  unsigned good = 2U;
  CHECK(comb_gpu_coset_good(&t, columns, 3U, &lookup, 8U,
                            generators, &good) == 3);
  free_gpu_tableau(&t);
}

static void test_invalid_states(void) {
  comb_gpu_tableau_t t = make_gpu_tableau(1U);
  unsigned negative = 2U;
  /* The requested operator is outside the generated stabiliser group. */
  t.rows = 2U;
  comb_gpu_set(t.x, 0U, 1U);
  comb_gpu_set(t.z + t.words, 0U, 1U);
  CHECK(comb_gpu_form_identity(&t, &negative) == 1);
  comb_gpu_clear_row(&t, 0U); comb_gpu_clear_row(&t, 1U);
  comb_gpu_set(t.x, 0U, 1U);
  t.delta[1] = 1U;
  CHECK(comb_gpu_form_identity(&t, &negative) == 2);
  free_gpu_tableau(&t);

  t = make_gpu_tableau(4U);
  const uint32_t columns[] = {0U, 1U};
  uint64_t lookup = UINT64_C(15), generators[4];
  unsigned good = 2U;
  for (uint32_t r = 0U; r < t.qubits; ++r)
    comb_gpu_set(t.x + (uint64_t)r * t.words, r, 1U);
  t.epsilon[3] = 1U;
  CHECK(comb_gpu_coset_good(&t, columns, 2U, &lookup, 4U,
                            generators, &good) == 2);
  for (uint32_t r = 0U; r <= t.qubits; ++r) comb_gpu_clear_row(&t, r);
  for (uint32_t r = 0U; r < t.qubits; ++r)
    comb_gpu_set((r == 3U ? t.z : t.x) + (uint64_t)r * t.words, r, 1U);
  CHECK(comb_gpu_coset_good(&t, columns, 2U, &lookup, 4U,
                            generators, &good) == 2);
  free_gpu_tableau(&t);
}

int main(void) {
  const uint32_t widths[] = {1U, 63U, 64U, 65U, 129U};
  for (size_t i = 0U; i < sizeof(widths) / sizeof(widths[0]); ++i)
    test_width(widths[i]);
  test_universal_coset();
  test_invalid_states();
  puts("CUDA arithmetic host differentials passed (1, 63, 64, 65, 129 qubits)");
  return 0;
}
