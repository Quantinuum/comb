#include "comb.h"

#define PRINT_SIZE 1000000
// Each worker needs its own buffer when logging is enabled.
static _Thread_local char print_buffer[PRINT_SIZE];

int do_print_pauli(mzd_t *x_parts, mzd_t *z_parts, mzd_t *delta_vector,
                   mzd_t *epsilon_vector, int row_number, char print_newline_p,
                   int counter) {
  int phase = 0;
  int n = x_parts->ncols;
  for (int i = 0; i < n; i++) {
    if (mzd_read_bit(x_parts, row_number, i) &&
        mzd_read_bit(z_parts, row_number, i))
      phase++;
  }

  if (mzd_read_bit(epsilon_vector, row_number, 0)) phase += 2;
  if (mzd_read_bit(delta_vector, row_number, 0)) phase++;
  switch (phase & 3) {
    case 0:
      counter += snprintf(print_buffer + counter, PRINT_SIZE - counter, "  ");
      break;
    case 1:
      counter += snprintf(print_buffer + counter, PRINT_SIZE - counter, " i");
      break;
    case 2:
      counter += snprintf(print_buffer + counter, PRINT_SIZE - counter, " -");
      break;
    case 3:
      counter += snprintf(print_buffer + counter, PRINT_SIZE - counter, "-i");
      break;
    default:
      counter += snprintf(print_buffer + counter, PRINT_SIZE - counter,
                          "error in printing pauli");
  }
  for (int i = 0; i < n; i++) {
    char x_p = mzd_read_bit(x_parts, row_number, i);
    char z_p = mzd_read_bit(z_parts, row_number, i);
    counter += snprintf(print_buffer + counter, PRINT_SIZE - counter,
                        z_p ? (x_p ? "Y" : "Z") : (x_p ? "X" : "I"));
  }
  if (print_newline_p)
    counter += snprintf(print_buffer + counter, PRINT_SIZE - counter, "\n");
  return counter;
}
char *print_row_vector_as_pauli(mzd_t *row_vector, int row,
                                const char *character) {
  int counter = 0;
  const char *identity = "I";
  for (int i = 0; i < row_vector->ncols; i++) {
    counter +=
        snprintf(print_buffer + counter, PRINT_SIZE - counter, "%s",
                 mzd_read_bit(row_vector, row, i) ? character : identity);
  }
  return print_buffer;
}
char *print_dcp(dcp_t *dcp, char print_newline_p) {
  int n = DCP_N(dcp);
  int counter = 0;
  mzd_t *zero = mzd_init(n, n);
  mzd_t *identity = mzd_init(n, n);
  mzd_t *zero_vector = mzd_init(n, 1);
  mzd_set_ui(zero, 0);
  mzd_set_ui(zero_vector, 0);
  mzd_set_ui(identity, 1);
  for (int i = 0; i < n; i++) {
    counter =
        do_print_pauli(identity, zero, zero_vector, zero_vector, i, 0, counter);
    counter += snprintf(print_buffer + counter, PRINT_SIZE - counter, "  ->");
    counter = do_print_pauli(identity, dcp->F, dcp->F_diagonal, dcp->h_z, i, 1,
                             counter);
    counter =
        do_print_pauli(zero, identity, zero_vector, zero_vector, i, 0, counter);
    counter += snprintf(print_buffer + counter, PRINT_SIZE - counter, "  ->");
    counter =
        do_print_pauli(zero, identity, zero_vector, dcp->h_x, i, 1, counter);
  }
  mzd_free(zero);
  mzd_free(zero_vector);
  mzd_free(identity);
  if (print_newline_p)
    counter += snprintf(print_buffer + counter, PRINT_SIZE - counter, "\n");
  return print_buffer;
}

char *print_stabiliser_state(stabiliser_state_t *state, char print_newline_p) {
  int counter = 0;
  int n = STATE_N(state);
  counter += snprintf(print_buffer + counter, PRINT_SIZE - counter, "{");
  for (int i = 0; i < state->x_parts->nrows; i++) {
    counter =
        do_print_pauli(state->x_parts, state->z_parts, state->delta_vector,
                       state->epsilon_vector, i, 0, counter);
    if (i + 1 != n)
      counter += snprintf(print_buffer + counter, PRINT_SIZE - counter, ", ");
  }
  counter += snprintf(print_buffer + counter, PRINT_SIZE - counter,
                      print_newline_p ? "}\n" : "}");
  return print_buffer;
}

char *print_pauli_z_coset(pauli_z_coset_t *coset, char print_newline_p) {
  int counter = 0;
  int n = COSET_N(coset);
  int k = COSET_K(coset);
  for (int i = 0; i < n; i++) {
    counter +=
        snprintf(print_buffer + counter, PRINT_SIZE - counter,
                 mzd_read_bit(coset->factor_transpose, 0, i) ? "Z" : "I");
  }
  counter += snprintf(print_buffer + counter, PRINT_SIZE - counter, "<");
  for (int i = 0; i < k; i++) {
    for (int j = 0; j < n; j++) {
      counter +=
          snprintf(print_buffer + counter, PRINT_SIZE - counter,
                   mzd_read_bit(coset->set_generators, i, j) ? "Z" : "I");
    }
    if (i != (k - 1))
      counter += snprintf(print_buffer + counter, PRINT_SIZE - counter, ", ");
  }
  counter += snprintf(print_buffer + counter, PRINT_SIZE - counter, ">");
  if (print_newline_p)
    counter += snprintf(print_buffer + counter, PRINT_SIZE - counter, "\n");
  return print_buffer;
}
