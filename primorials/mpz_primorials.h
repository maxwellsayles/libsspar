/**
 * @file mpz_primorials.h
 * 
 * Factored 2,3 representations of primorials optimized for mpz_qforms.
 */

#pragma once
#ifndef MPZ_PRIMORIALS__INCLUDED
#define MPZ_PRIMORIALS__INCLUDED

#include <stdint.h>

#include "liboptarith/closest_23.h"

extern const factored_two_three_term16_t* mpz_primorial_terms[];
extern const uint16_t mpz_primorial_term_counts[];

#endif

