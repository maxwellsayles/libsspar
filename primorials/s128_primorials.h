/**
 * @file s128_primorials.h
 * 
 * Factored 2,3 representations of primorials optimized for s128_qforms.
 */

#pragma once
#ifndef S128_PRIMORIALS__INCLUDED
#define S128_PRIMORIALS__INCLUDED

#include <stdint.h>

#include "liboptarith/closest_23.h"

extern const factored_two_three_term16_t* s128_primorial_terms[];
extern const uint16_t s128_primorial_term_counts[];

#endif

