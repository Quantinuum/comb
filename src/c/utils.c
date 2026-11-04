#include "comb.h"

void fix_dcp(dcp_t *dcp) {
  // shouldn't normally need to be called, except for testing purposes
  int n = DCP_N(dcp);
  for (int i = 0; i < n; i++) {
    mzd_write_bit(dcp->F_diagonal, i, 0, mzd_read_bit(dcp->F, i, i));
  }
  mzd_transpose(dcp->h_x_transpose, dcp->h_x);
}

word _popcnt(word x){
  word r = 0;
  lwlog_debug("testing %d",x);
  for(int i = 0 ; i < sizeof(x)*8; i++){
    if(x & 1){
      lwlog_debug("incf r, now %d",r);
      r ++;
    }
    x>>=1;
    lwlog_debug("x is now %d",x);
  }
  return r;
}
inline word popcnt(word x){
  uint8_t* arr = (uint8_t*)&x;
  word r = 0;
  for(int i = 0; i <8; i++){
    r+=__builtin_popcount(arr[i]); // simd popcount? do a whole simd one at a time
  }
  return r;
}

BIT row_product(mzd_t *x, int i, mzd_t *z, int j) {
  // Forms the inner product of row i of x and row j of z (viewing both rows as
  // vectors over F2).
  word acc = 0;
  lwlog_debug("in dot product");
  word *row1 = mzd_row(x, i);
  word *row2 = mzd_row(z, j);
  word mask_end = x->high_bitmask | z->high_bitmask;
  wi_t width = x->width - 1;
  for (int k = 0; k < width; k++) {
    word x_part = row1[k];
    word z_part = row2[k];
    lwlog_debug("x part = %d, z part = %d",x_part,z_part);
    word part_of_product = x_part & z_part;
     acc = acc ^ (popcnt(part_of_product) & 1);
    }
    word x_part = row1[width] & mask_end;
    word z_part = row2[width] & mask_end;
    word part_of_product = x_part & z_part;
    lwlog_debug("x part = %d, z part = %d, part of product = %d, popcount = %d ",x_part,z_part,part_of_product, popcnt(part_of_product));
    acc = acc ^ (popcnt(part_of_product) & 1);
    #if LOG_LEVEL >= DEBUG
    word acc2 = 0;
    int w = x->ncols;
    // slow impl:
    for (int k = 0; k < w; k++) {
      if(mzd_read_bit(x,i,k) && mzd_read_bit(z,j,k))acc2=1-acc2;
    }
    if(acc2 != acc){
      printf("x:\n");
      mzd_print(x);
      printf("z:\n");
      mzd_print(z);
      lwlog_err("dot products not equal: acc = %d should be %d, i = %d, j = %d",acc,acc2,i,j);
      exit(1);
    }
#endif
    return acc;
}

void expand_state(stabiliser_state_t *state) {
  state->x_parts->nrows++;
  state->z_parts->nrows++;
  state->delta_vector->nrows++;
  state->epsilon_vector->nrows++;
}

void contract_state(stabiliser_state_t *state) {
  state->x_parts->nrows--;
  state->z_parts->nrows--;
  state->delta_vector->nrows--;
  state->epsilon_vector->nrows--;
}
void clear_state(stabiliser_state_t *state) {
  // clears the state. assumes it hasn't been expanded.
  expand_state(state);
  mzd_set_ui(state->x_parts, 0);
  mzd_set_ui(state->z_parts, 0);
  mzd_set_ui(state->epsilon_vector, 0);
  mzd_set_ui(state->delta_vector, 0);
  contract_state(state);
}
void set_row(mzd_t *x_parts, mzd_t *z_parts, int j, mzd_t *x_measurement_matrix,
             int i, mzd_t *epsilon_bits,
             mzd_t *delta_bits) {  // moves m_i to row j of x_parts, and clears
                                   // row j of z_parts
  mzd_row_clear_offset(x_parts, j, 0);
  mzd_copy_row(x_parts, j, x_measurement_matrix, i);
  mzd_row_clear_offset(z_parts, j, 0);
  mzd_write_bit(epsilon_bits, j, 0, 0);
  mzd_write_bit(delta_bits, j, 0, 0);
}

inline unsigned compress4(unsigned x, unsigned m) {  // from hacker's delight
  unsigned mk, mp, mv, t;
  int i;
  x = x & m;     // Clear irrelevant bits.
  mk = ~m << 1;  // We will count 0's to right.
  for (i = 0; i < 5; i++) {
    mp = mk ^ (mk << 1);  // Parallel suffix.
    mp = mp ^ (mp << 2);
    mp = mp ^ (mp << 4);
    mp = mp ^ (mp << 8);
    mp = mp ^ (mp << 16);
    mv = mp & m;                      // Bits to move.
    m = (m ^ mv) | (mv >> (1 << i));  // Compress m.
    t = x & mv;
    x = (x ^ t) | (t >> (1 << i));  // Compress x.
    mk = mk & ~mp;
  }
  return x;
}
void permute_columns(mzd_t *data, unsigned mask) {
  for (int i = 0; i < data->nrows; i++) {
    word *b = mzd_row(
        data,
        i);  // assume small enough so that only have to do the first 64 bits
    b[0] = compress4(b[0], mask);
  }
}
void convert_row(mzd_t *M, int row, comb_vector_t *vector) {
  for (int i = 0; i < M->ncols; i++) {
    mzd_write_bit(M, row, i, TestBit(vector->data, i) ? 1 : 0);
  }
}

void add_and_column_vector(mzd_t *result, mzd_t *F1, mzd_t *F2) {
  // assumes that f1, f2, and result are congruous column vectors. necessarily
  // linear time.
  for (int i = 0; i < result->nrows; i++) {
    mzd_write_bit(result, i, 0,
                  (mzd_read_bit(result, i, 0) ^
                   (mzd_read_bit(F1, i, 0) & mzd_read_bit(F2, i, 0))));
  }
}


int stabilisers_rref_p(comb_vector_t *operators, uint32_t n, uint32_t m,
                       uint32_t *coefficients, unsigned *permutation_mask) {
  // ensures the operators are in rref. n operators on m qubits.
  int current_leading_index = -1;
  for (int row = 0; row < n; row++) {
    comb_vector_t *current_operator = &(operators[row]);
    int new_leading_index = m;  // new_leading_index = n means that the row
                                // contains no entries at all.
    for (int column = 0; column < m && (new_leading_index == m); column++) {
      if (TestBit(current_operator->data, column)) new_leading_index = column;
    }
    if (new_leading_index < current_leading_index &&
        !((new_leading_index = m) && (current_leading_index == m)))
      return 0;  // new leading index must be to the right of the old one
    // verify if new leading index has no entries in other rows
    for (int row_2 = 0; row_2 < n; row_2++) {
      if (row_2 != row) {
        if (TestBit(operators[row_2].data, new_leading_index))
          return 0;  // leading index shares entries in other rows
      }
    }
    current_leading_index = new_leading_index;
    coefficients[row] = current_leading_index;
    *permutation_mask |= (1 << current_leading_index);
    lwlog_info("rref: coefficients[%d]=%d", row, current_leading_index);
  }
  *permutation_mask = ~*permutation_mask;
  return 1;
}
