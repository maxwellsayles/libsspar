/**
 * @file sspar.h
 * Interface for the SuperSPAR factoring algorithm.
 */
#pragma once
#ifndef SSPAR__INCLUDED
#define SSPAR__INCLUDED

#include <gmp.h>
#include <stdint.h>

#include "liboptarith/closest_23.h"
#include "liboptarith/group.h"
#include "liboptarith/group_pow.h"
#include "libqform/mpz_qform.h"
#include "libqform/qform_group.h"
#include "libqform/s64_qform.h"
#include "libqform/s128_qform.h"
#include "libsspar/open_addr_hash.h"

/**
 * The structure that maintains the data for the SuperSPAR algorithm.
 */
typedef struct {
  mpz_t D;  // Discriminant
  mpz_t t;  // Temporary

  // Active group elements.
  // These are polymorphic based on discriminant.
  qform_group_t* qform_group;
  qform_t* initform;
  qform_t* search;
  qform_t* giant;
  qform_t* current;
  qform_t* temp_form;

  // Concrete group elements for each discriminant group.
  s64_qform_group_t qform_group_s64;
  s64_qform_t initform_s64;
  s64_qform_t search_s64;
  s64_qform_t giant_s64;
  s64_qform_t current_s64;
  s64_qform_t temp_form_s64;
  s128_qform_group_t qform_group_s128;
  s128_qform_t initform_s128;
  s128_qform_t search_s128;
  s128_qform_t giant_s128;
  s128_qform_t current_s128;
  s128_qform_t temp_form_s128;
  mpz_qform_group_t qform_group_mpz;
  mpz_qform_t initform_mpz;
  mpz_qform_t search_mpz;
  mpz_qform_t giant_mpz;
  mpz_qform_t current_mpz;
  mpz_qform_t temp_form_mpz;

  // Powering
  group_pow_t* pow;  // Polymorphic
  group_pow_t pow_s64;
  group_pow_t pow_s128;
  group_pow_t pow_mpz;

  // The order of the ideal once found.
  // NOTE: This is limited to 32bits as the hashtable used has 32bit keys.
  uint32_t order;

  // Bound of 2^logh of the even part of the order.
  int logh;  // h_D = sqrt(D)*log(D)/pi

  // Primorial index used for big exponentiation.
  int primorial_index;
  const factored_two_three_term16_t* primorial_terms;
  int primorial_term_count;

  // Take coprime baby-steps \le step_bound and take giant-steps
  // that are a multiple of step_bound.
  int step_bound;

  // The number of coprime baby-steps \le step_bound not including the
  // identity.
  int step_count;

  // Coprimes.  There must be at least step_count number of coprime_deltas.
  int first_coprime;
  const uint8_t* coprime_deltas;

  // Cached difference between coprimes.
  // Since all deltas are even, the largest delta is (gap_capacity-1)*2.
  qform_t** gap;  // Polymorphic.
  s64_qform_t* gap_s64;
  s128_qform_t* gap_s128;
  mpz_qform_t* gap_mpz;
  int gap_capacity;
  int max_gap;

  // Table used by search phase
  hash_table_t table;
} sspar_t;

#define SSPAR_MIN_PRIMORIAL_INDEX  0
#define SSPAR_MAX_PRIMORIAL_INDEX  512

#define SSPAR_ADDITIONAL_QFORMS_AUTO (-1)

/**
 * Initialize the SuperSPAR data structure.
 */
void sspar_init(sspar_t* this);

/**
 * Release any memory used by the SuperSPAR data structure.
 */
void sspar_clear(sspar_t* this);

/**
 * Attempts to factor the integer N.
 * Make sure to call sspar_setup() with appropriate parameters before
 * calling sspar_factor_raw().
 * @param d (out) A divisor of N.
 * @param N The integer to be factored.
 * @param primorial_index The nth primorial to use for exponentiation.
 * @param step_bound Take coprime baby-steps \le step_bound and take giant
 *                   steps that are a multiple of step_bound.
 * @param step_count The number of coprime baby-steps \le step_bound.
 * @param step_first_coprime The first coprime to step_bound.
 * @param step_coprime_deltas Deltas for consecutive coprimes to step_bound.
 * @param additional_qforms The number of qforms to try once a multiple of
 *                          the order of the first qform is known. Use
 *                          SSPAR_ADDITIONAL_QFORMS_AUTO to let the
 *                          algorithm select the value based on the size
 *                          of N.
 * @return 1 if a non-trivial factor of N is found, 0 otherwise.
 */
int sspar_factor_raw(sspar_t* this,
                     mpz_t d,
                     const mpz_t N,
                     const int primorial_index,
                     const int step_bound,
                     const int step_count,
                     const int step_first_coprime,
                     const uint8_t* step_coprime_deltas,
		     const int additional_qforms);

/**
 * Attempts to factor the integer N.
 * @param N The integer to be factored.
 * @param d A divisor of N.
 * @return 1 If a non-trivial factor of N is found, 0 otherwise.
 */
int sspar_factor(sspar_t* this, mpz_t d, const mpz_t N);

#endif

