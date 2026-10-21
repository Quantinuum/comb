#include "comb.h"

// Attempts to increment the bit vector v in place. Returns 1 if this is not
// possible, 0 otherwise.
int incf_bit_vector(comb_bit_array_t v, int n) {
  for (int m = (n - 1); m >= 0; m--) {
    if (TestBit(v, m)) {  // m the index of the rightmost 1
      for (int i = (m - 1); i >= 0; i--) {
        if (TestBit(v, i) == 0) {
          SetBit(v, i);
          ClearBit(v, i + 1);
          int onecount = 0;
          for (int k = (1 + i); k < n; k++) {
            if (TestBit(v, k)) onecount++;
          }
          for (int k = i + 2; k < (n - onecount); k++) ClearBit(v, k);
          for (int k = (n - onecount); k < n; k++) SetBit(v, k);
          return 0;
        }
      }
    }
  }
  return 1;
}
// Increments the fault combination n_0, ... , n_(k-1) and returns 0. if not
// possible it returns 1.
int ft_iteration_combination_incf(int k, uint32_t *sizes, int n,
                                  int *first_balls, int *last_balls,
                                  int *current_combination) {
  // k number of balls to take - size of current_combination
  // n number of sizes/baskets
  for (int i = 0; i < k; i++) {
    if (last_balls[n - k + i] !=
        current_combination[i]) {  // this means we are not at the final fault
                                   // combination yet
      for (int j = 0; j < n; j++) {
        // move from left to right, if you encounter a basket with a ball under
        // it, incf that ball and set baskets to the left
        for (int i = 0; i < k; i++) {
          int last_ball_p = 0;
          for (int p = 0; p < n; p++) {
            if (last_balls[p] == current_combination[i]) last_ball_p = 1;
          }
          if (!last_ball_p) {
            for (int j = 0; j < i; j++) {
              int basket_ref = -1;
              for (int r = 0; r < n; r++) {
                if ((first_balls[r] <= current_combination[j]) &&
                    (current_combination[j] <= last_balls[r])) {
                  basket_ref = r;
                  break;
                }
              }
              current_combination[j] = first_balls[basket_ref];
            }
            current_combination[i]++;
            return 1;  // incf was successful
          }
        }
        // we know that every ball is a last one.
        // this means that seeing if the y-th basket is occupied is equivalent
        // to checking if last_ball[y] is in current_combination.

        comb_bit_word_t v[MAX_SIZE_COUNT];  // v the bit string which has the y-th entry
                   // as 1 if the y-th basket is occupied.
        for (int i = 0; i < n; i++) ClearBit(v, i);
        for (int i = 0; i < n; i++) {
          for (int j = 0; j < k; j++) {
            if (last_balls[i] == current_combination[j]) SetBit(v, (n - 1) - i);
          }
        }

        if (incf_bit_vector(v, n)) {
          lwlog_err("incf bit vector returned 1; should not be here\n");
          return 0;
        }

        int j2 = 0;
        for (int i = 0; i < n; i++) {
          if (TestBit(v, (n - 1) - i)) {
            current_combination[j2++] = first_balls[i];
          }
        }
        return 1;
      }
    }
  }
  return 0;  // we had the last element so return 0 to indicate there is no
             // successor
}

// https://stackoverflow.com/questions/24294192/computing-the-binomial-coefficient-in-c
unsigned long long combi(int n,int k)
{
    if(n==0 || n < k)return 0;
  unsigned long long ans=1;
  k=k>n-k?n-k:k;
  int j=1;
  for(;j<=k;j++,n--)
      ans = ans/j*n + ans%j*n/j;
  return ans;
}
void unrank(int s, int* combination, unsigned long long m){
  // sets the combination in combination to the m-th combination ranked lexicographically.
  lwlog_debug("Unranking %llu m:\n",m);
  for(int k = s-1; k >=0; k--){
    unsigned long long c_k = 0;
    while(combi(c_k+1,k+1) <= m){
      c_k++;
    }
    combination[k]=c_k;
    lwlog_debug("combination[%d] = %d\n",k,c_k);
    m-=combi(c_k,k+1);
  //printf("m is %llu:\n",m);
  }
}
unsigned long long rank(int s, int* combination){
  unsigned long long n = 0;
  for(int i = 0; i < s; i++){
    lwlog_debug("combination[%d] = %d m:\n",i,combination[i]);
    n+=combi(combination[i],i+1);
  }
  lwlog_debug("Ranked to %llu:\n",n);
  return n;
}
