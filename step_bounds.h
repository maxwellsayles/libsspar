#pragma once
#ifndef STEP_BOUNDS__INCLUDED
#define STEP_BOUNDS__INCLUDED

#define STEP_BOUNDS_MIN_INDEX 0
#define STEP_BOUNDS_MAX_INDEX 175

typedef struct {
  int bound;
  int odd_primorial_index;
  int step_count;
} step_bound_t;

extern step_bound_t step_bounds[176];

#endif  // STEP_BOUNDS__INCLUDED

