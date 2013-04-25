/**
 * @file util.h
 * Utility functions.
 */
#pragma once
#ifndef UTIL__INCLUDED
#define UTIL__INCLUDED

#include <stdarg.h>  // for vprintf
#include <gmp.h>
#include <inttypes.h>
#include <stdint.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>

#include "liboptarith/math64.h"

#define min(a,b) ((a)<(b)?(a):(b))
#define max(a,b) ((a)>(b)?(a):(b))

// FIB_PRIME/2^64 is approximately (sqrt(5)-1)/2 and prime.
#define FIB_PRIME ((uint64_t)11400714819323198549ULL)

// 0 is fatal
// 1 is progress
// 2 is group members used
// 3 is operations
#define verbose_level 0

// To prevent implicit declaration error in cprintf.
// TODO: There should be a better way than including this line.
__GMP_DECLSPEC int gmp_vprintf (const char *, va_list);

/// Print if level is <= verbose_level.
static inline int cprintf(const int level, const char* fmt, ...) {
  int res = 0;
  va_list args;
  if (level <= verbose_level) {
    va_start(args, fmt);
    res = gmp_vprintf(fmt, args);
    va_end(args);
  }
  return res;
}

/// Gives the time from system on in nanoseconds
static inline uint64_t current_nanos(void) {
#ifdef __linux__
  struct timespec res;
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &res);
  return (res.tv_sec * 1000000000ULL) + res.tv_nsec;
#else
  struct timeval tv;
  gettimeofday(&tv, 0);
  return ((uint64_t)tv.tv_sec * 1000000ULL + (uint64_t)tv.tv_usec) * 1000;
#endif
}

/// maps v \in [l1,h1] -> [l2,h2]
static inline int affine_map_i(int v, int l1, int h1, int l2, int h2) {
  return l2 + (v-l1)*(h2-l2)/(h1-l1);
}

/// maps v \in [l1,h1] -> [l2,h2]
static inline double affine_map_d(double v, double l1, double h1, double l2, double h2) {
  return l2 + (v-l1)*(h2-l2)/(h1-l1);
}

/// maps v \in [l1,h1] -> [l2,h2]
static inline int64_t affine_map_s64(int64_t v, int64_t l1, int64_t h1, int64_t l2, int64_t h2) {
  return l2 + (v-l1)*(h2-l2)/(h1-l1);
}

/// maps v \in [l1,h1] -> [l2,h2]
static inline uint64_t affine_map_u64(uint64_t v, uint64_t l1, uint64_t h1, uint64_t l2, uint64_t h2) {
    return l2 + (v-l1)*(h2-l2)/(h1-l1);
}

/// Compute the average value in the array.
int average_int(int count, const int* array);

/// Compute the average value in the array.
unsigned int average_ui(int count, const unsigned int* array);

/// Compute the average value in the array.
uint64_t average_u64(int count, const uint64_t* array);

/// Bog down the CPU by performing sqrts of random 64bit integers
static inline void full_cpu_load(uint64_t secs) {
  static uint64_t s;
  uint64_t nanos = secs * 1000000000ULL;
  uint64_t t;
  uint64_t start = current_nanos();
  while (current_nanos() < start + nanos) {
    // take the sqrt of a random 128bit integer
    t = rand_u64();
    s += sqrt_u64(t);
  }
}

/// Write a gnuplot data file
void write_gnuplot_datfile_u64(const char* filename,
                               const uint64_t* x,
                               const uint64_t* y,
                               const int sample_count);

#endif  // UTIL__INCLUDED

