#include "libsspar/util.h"

#include <inttypes.h>
#include <stdlib.h>

/**
 * computes the average of an array of integers
 */
unsigned int average_ui(int count, const unsigned int* array) {
  uint64_t sum = 0;
  int i = 0;
  for (i = 0;  i < count;  i ++) {
    sum += array[i];
  }
  return sum / count;
}

/**
 * computes the average of an array of integers
 */
int average_int(int count, const int* array) {
  int64_t sum = 0;
  int i = 0;
  for (i = 0;  i < count;  i ++) {
    sum += array[i];
  }
  return sum / count;
}

/**
 * computes the average of an array of uint64_t
 */
uint64_t average_u64(int count, const uint64_t* array) {
  uint64_t sum = 0;
  int i = 0;
  for (i = 0;  i < count;  i ++) {
    sum += array[i];
  }
  return sum / count;
}

/**
 * Write a gnuplot data file
 */
void write_gnuplot_datfile_u64(const char* filename,
                               const uint64_t* x,
                               const uint64_t* y,
                               const int sample_count) {
  int i;
  FILE* f;
  f = fopen(filename, "w");
  for (i = 0;  i < sample_count;  i ++) {
    fprintf(f, "%"PRIu64", %"PRIu64"\n", x[i], y[i]);
  }
  fclose(f);
}

