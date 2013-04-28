/**
 * SuperSPAR factoring implementation.
 *
 * For discriminants sized [32, 40..80], for multipliers [1, 2, 3, 5, 6, 7, 10],
 * and for the first 5 prime forms that existed started with 11, we computed
 * a total of 65,389,135 ideal orders and the factorization of the order.
 *
 * In expectation, we should not need to pick more than three prime forms
 * before a successful factorization. Therefore, if we haven't successfully
 * factorized for a given multiplier within three attempts, we switch multipliers.
 *
 * Furthermore, we noticed that the factorization of each prime form for
 * a given discriminant and multiplier tend to only differ by a few small factors
 * at most, and these are likely to be removed by the initial primorial
 * exponentiation.  Therefore, once we determine the order of a prime form,
 * we use that order on two successive prime forms, if necessary.
 */
#include "libsspar/sspar.h"

#include <assert.h>
#include <gmp.h>
#include <inttypes.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "liboptarith/closest_23.h"
#include "liboptarith/group.h"
#include "liboptarith/math32.h"
#include "liboptarith/math_mpz.h"
#include "liboptarith/primes.h"
#include "liboptarith/primorial.h"
#include "liboptarith/square_free.h"
#include "libqform/mpz_qform.h"
#include "libqform/qform_group.h"
#include "libqform/s64_qform.h"
#include "libqform/s128_qform.h"
#include "libsspar/coprime_deltas.h"
#include "libsspar/open_addr_hash.h"
#include "libsspar/primorials/mpz_primorials.h"
#include "libsspar/primorials/s64_primorials.h"
#include "libsspar/primorials/s128_primorials.h"
#include "libsspar/step_bounds.h"
#include "libsspar/util.h"

// TODO: tidy up these values when tighter bounds are known.
#define FIRST_PRIME_INDEX        4  // First prime form is >= 11.
#define TABLE_SCALAR             3
#define TABLE_MAX_ENTRIES        65536
#define TABLE_MAX_MATCHES        16
#define MAX_GAP                  22  // Empirically determined.

#define SPECIAL_MULTIPLIERS_COUNT       7
const int special_multipliers_1mod4[] = { 6, 10, 3, 1, 7, 2, 5 };
const int special_multipliers_3mod4[] = { 1, 10, 3, 2, 5, 7, 6 };

/// Default values for an instance of SuperSPAR.
void sspar_init(sspar_t* this) {
  mpz_init(this->D);
  mpz_init(this->t);

  // Initialize groups.
  s64_qform_group_init(&this->qform_group_s64);
  s64_qform_init(&this->qform_group_s64, &this->initform_s64);
  s64_qform_init(&this->qform_group_s64, &this->search_s64);
  s64_qform_init(&this->qform_group_s64, &this->giant_s64);
  s64_qform_init(&this->qform_group_s64, &this->current_s64);
  s64_qform_init(&this->qform_group_s64, &this->temp_form_s64);
  group_pow_init(&this->pow_s64, &this->qform_group_s64.desc.group);

  s128_qform_group_init(&this->qform_group_s128);
  s128_qform_init(&this->qform_group_s128, &this->initform_s128);
  s128_qform_init(&this->qform_group_s128, &this->search_s128);
  s128_qform_init(&this->qform_group_s128, &this->giant_s128);
  s128_qform_init(&this->qform_group_s128, &this->current_s128);
  s128_qform_init(&this->qform_group_s128, &this->temp_form_s128);
  group_pow_init(&this->pow_s128, &this->qform_group_s128.desc.group);

  mpz_qform_group_init(&this->qform_group_mpz);
  mpz_qform_init(&this->qform_group_mpz, &this->initform_mpz);
  mpz_qform_init(&this->qform_group_mpz, &this->search_mpz);
  mpz_qform_init(&this->qform_group_mpz, &this->giant_mpz);
  mpz_qform_init(&this->qform_group_mpz, &this->current_mpz);
  mpz_qform_init(&this->qform_group_mpz, &this->temp_form_mpz);
  group_pow_init(&this->pow_mpz, &this->qform_group_mpz.desc.group);
    
  // Orders.
  this->order = 0;
    
  // 2^logh bound.
  this->logh = 0;

  // Primorial.
  this->primorial_index = 0;

  // Create coprime delta table.
  this->first_coprime = 0;
  this->coprime_deltas = 0;

  // Cached difference between coprimes.
  this->gap_capacity = (MAX_GAP>>1) + 1;
  this->gap = (qform_t**)malloc(sizeof(qform_t*) * this->gap_capacity);
  this->gap_s64 = (s64_qform_t*)group_elem_array_alloc(
      &this->qform_group_s64.desc.group, this->gap_capacity);
  this->gap_s128 = (s128_qform_t*)group_elem_array_alloc(
      &this->qform_group_s128.desc.group, this->gap_capacity);
  this->gap_mpz = (mpz_qform_t*)group_elem_array_alloc(
      &this->qform_group_mpz.desc.group, this->gap_capacity);
  this->max_gap = 0;
    
  // Number of primorial steps.
  this->step_bound = 0;
    
  // Initialize the table for the search phase.
  hash_table_init(&this->table, TABLE_MAX_ENTRIES);
}

/// Release SSPAR resources.
void sspar_clear(sspar_t* this) {
  hash_table_clear(&this->table);

  free(this->gap);

  group_pow_clear(&this->pow_s64);
  group_elem_array_free(&this->qform_group_s64.desc.group,
			this->gap_s64, this->gap_capacity);
  s64_qform_clear(&this->qform_group_s64, &this->initform_s64);
  s64_qform_clear(&this->qform_group_s64, &this->search_s64);
  s64_qform_clear(&this->qform_group_s64, &this->giant_s64);
  s64_qform_clear(&this->qform_group_s64, &this->current_s64);
  s64_qform_clear(&this->qform_group_s64, &this->temp_form_s64);
  s64_qform_group_clear(&this->qform_group_s64);

  group_pow_clear(&this->pow_s128);
  group_elem_array_free(&this->qform_group_s128.desc.group,
			this->gap_s128, this->gap_capacity);
  s128_qform_clear(&this->qform_group_s128, &this->initform_s128);
  s128_qform_clear(&this->qform_group_s128, &this->search_s128);
  s128_qform_clear(&this->qform_group_s128, &this->giant_s128);
  s128_qform_clear(&this->qform_group_s128, &this->current_s128);
  s128_qform_clear(&this->qform_group_s128, &this->temp_form_s128);
  s128_qform_group_clear(&this->qform_group_s128);

  group_pow_clear(&this->pow_mpz);
  group_elem_array_free(&this->qform_group_mpz.desc.group,
			this->gap_mpz, this->gap_capacity);
  mpz_qform_clear(&this->qform_group_mpz, &this->initform_mpz);
  mpz_qform_clear(&this->qform_group_mpz, &this->search_mpz);
  mpz_qform_clear(&this->qform_group_mpz, &this->giant_mpz);
  mpz_qform_clear(&this->qform_group_mpz, &this->current_mpz);
  mpz_qform_clear(&this->qform_group_mpz, &this->temp_form_mpz);
  mpz_qform_group_clear(&this->qform_group_mpz);

  mpz_clear(this->D);
  mpz_clear(this->t);
}

/**
 * Configure the internal parameters for using SuperSPAR.
 * @param primorial_index The nth primorial to use for exponentiation.
 * @param step_bound Take coprime baby-steps \le step_bound and take
 *                   giant-steps that are a multiple of step_bound.
 */
static void sspar_setup(sspar_t* this,
			const int primorial_index,
			const int step_bound,
			const int step_count,
			const int step_first_coprime,
			const uint8_t* step_coprime_deltas) {
  assert(primorial_index >= SSPAR_MIN_PRIMORIAL_INDEX &&
	 primorial_index <= SSPAR_MAX_PRIMORIAL_INDEX);
  assert(step_count > 0);
  assert(step_first_coprime > 0);
  assert(step_coprime_deltas != 0);
    
  this->primorial_index = primorial_index;
    
  this->step_bound = step_bound;
  this->step_count = step_count;

  this->first_coprime = step_first_coprime;
  this->coprime_deltas = step_coprime_deltas;    
}

/**
 * Hashes the form and checks for a collision.
 * Returns 1 if we found the order, 0 otherwise.
 * 
 * This takes advantage of inversions, i.e. e-x and e+x
 * where x is a previously cached form.
 */
static int test_and_hash_qform(sspar_t* this,
			       const uint32_t e,
			       const qform_t* form) {
  group_t* group = &this->qform_group->group;
  uint32_t matches[TABLE_MAX_MATCHES];
  uint32_t hash = 0;
  int32_t d = 0;
  int k = 0;
  int j = 0;

  if (verbose_level >= 3) {
    cprintf(3, "exp=%"PRIu32" form=", e);
    group->print(group, form);
    cprintf(3, "\n");
  }

  // TODO: We can now hash 0 values, so this code can change to not
  //       explicitly check for identity.
  if (group->is_id(group, form)) {
    this->order = e;
    return 1;
  }

  // Look up the value.
  hash = group->hash32(group, form);
  k = hash_table_lookup(&this->table, &matches[0], TABLE_MAX_MATCHES, hash);
    
  // TODO: The current hash will return either 0 or 1 value.
  // TODO: Consider not hashing both a form and its inverse.  Then
  //       when we get a lookup, we don't need the extra eponentiation
  //       to check if its actually the order.  We'll just assume that it is.
  // For each match, check to see if one of them is a multiple of the order.
  for (j = 0;  j < k;  j ++) {
    cprintf(3, "Testing match %"PRIu32"\n", matches[j]);

    // Check if e-matches[j] is a multiple of the order.
    d = e - matches[j];
    assert(d > 0);
    qform_pow_u32(this->pow, this->temp_form, this->search, d);
    if (group->is_id(group, this->temp_form)) {
      this->order = d;
      cprintf(3, "Found the order: %"PRIu32"\n", this->order);
      return 1;
    }
    
    // Check if e+matches[j] is a multiple of the order.
    d = e + matches[j];
    assert(d > 0);
    qform_pow_u32(this->pow, this->temp_form, this->search, d);
    if (group->is_id(group, this->temp_form)) {
      this->order = d;
      cprintf(3, "Found the order: %"PRIu32"\n", this->order);
      return 1;
    }
  }
  
  // Hash the value.
  hash_table_insert(&this->table, e, hash);
  return 0;
}

/**
 * Expects this->search to have the form whose order we want to find.
 * Uses this->current, this->gap, this->table, and this->giant.
 *
 * If the order is found, 1 is returned and this->order is the order
 * of the element. Otherwise 0 is returned and this->order = 0.
 */
static int sspar_search(sspar_t* this) {
  group_t* group = &this->qform_group->group;
  int coprime_index = 0;  // also the number of steps taken
  int coprime = 0; 
  int coprime_delta = 0;
  int big_step = 0;  // giant step size
  int last_gap = 2;
    
  assert(this->step_bound % 2 == 0);
  if (verbose_level >= 2) {
    cprintf(2, "Searching with: ");
    group->print(group, this->search);
    cprintf(2, "\n");
  }
    
  // Cache gap 0 and 2.
  group->set_id(group, this->gap[0]);
  group->square(group, this->gap[2>>1], this->search);

  // Reset hash table for baby steps.
  cprintf(2, "Clearing hash table.\n");
  hash_table_reset(&this->table, TABLE_SCALAR * this->step_count);

  // Store the initial form.
  hash_table_insert(&this->table, 1, group->hash32(group, this->search));

  // Exponentiate to first coprime.
  coprime = this->first_coprime;
  qform_pow_u32(this->pow, this->current, this->search, coprime);
  if (test_and_hash_qform(this, coprime, this->current)) {
    // We found the order.
    return 1;
  }

  // Walk coprime baby-steps.
  // NOTE: We took our first step above, so we iterate coprime_index
  //       to one less than this->step_count.
  for (; coprime_index < this->step_count - 1; coprime_index ++) {
    // Load the current delta.
    coprime_delta = this->coprime_deltas[coprime_index];
    coprime += coprime_delta;

    // Make sure we have enough cached gaps.
    assert(coprime_delta <= MAX_GAP);
    while (last_gap < coprime_delta) {
      last_gap += 2;
      group->compose(group,
		     this->gap[last_gap>>1],
		     this->gap[(last_gap-2)>>1],
		     this->gap[2>>1]);
    }
    assert(last_gap >= coprime_delta);
        
    // Exponentiate to the next coprime.
    // Use cached gap (gap is always even).
    group->compose(group,
		   this->current,
		   this->current,
		   this->gap[coprime_delta>>1]);

    if (test_and_hash_qform(this, coprime, this->current)) {
      return 1;
    }
  }
  assert(coprime <= this->step_bound);
  assert(coprime + this->coprime_deltas[coprime_index] > this->step_bound);
    
  // Walk coprime to step_bounds.
  coprime_delta = this->step_bound - coprime;
  assert(coprime_delta % 2 == 1);
  group->compose(group, this->current, this->current, this->search);
  coprime_delta--;
  coprime++;
  if (test_and_hash_qform(this, coprime, this->current)) {
    return 1;
  }
  while (coprime_delta > 0) {  // move by at most last_gap each step
    int delta = min(coprime_delta, last_gap);
    group->compose(group, this->current,
		   this->current,
		   this->gap[delta>>1]);
    coprime_delta -= delta;
    coprime += delta;
    if (test_and_hash_qform(this, coprime, this->current)) {
      return 1;
    }
  }
  assert(coprime == this->step_bound);

  // Giant step size is twice the step bound.
  group->square(group, this->current, this->current);
  coprime += coprime;
  if (test_and_hash_qform(this, coprime, this->current)) {
    return 1;
  }
  group->set(group, this->giant, this->current);
  big_step = coprime;
    
  // Walk giant-steps.
  coprime_index ++;  // This accounts for the first baby-step.
  for (; coprime_index > 0; coprime_index --) {        
    coprime += big_step;
    group->compose(group, this->current, this->current, this->giant);
    if (test_and_hash_qform(this, coprime, this->current)) {
      return 1;
    }
  }

  // we failed to find the order
  this->order = 0;
  return 0;
}

/**
 * Find the next prime ideal that doesn't divide the discriminant.
 * Takes a prime index as input and returns the prime index of the prime
 * ideal as output.
 * TODO: Any prime divisors are non-trivial divisors.
 */
static int next_primeform(qform_group_t* qgroup,
			  qform_t* form,
			  int prime_index) {
  group_t* group = &qgroup->group;
  prime_index = qform_next_primeform(qgroup, form, prime_index);
  assert(prime_index < prime_list_count);
  if (verbose_level >= 2) {
    cprintf(2, "Using primeform: ");
    group->print(group, form);
    cprintf(2, "\n");
  }
  return prime_index;
}

/**
 * Exponentiates this->init_form by primorial[this->primorial_index]
 */
static inline void exponentiate_primorial(sspar_t* this) {
  group_pow_factored23(this->pow, this->initform, this->initform,
		       this->primorial_terms, this->primorial_term_count);
}

/**
 * Takes this->initform and exponentiates it by this->order.
 * Then squares this until an ambigous form is found and attempts to factor N.
 * @param d The non-trivial factor if one is found.
 * @return 1 if a non-trivial factor was found.
 */
static int test_ambiguous_form(sspar_t* this, mpz_t d, const mpz_t N) {
  int i;
  qform_group_t* qgroup = this->qform_group;
  group_t* group = &qgroup->group;
    
  // Exponentiate initform by the odd part of order.
  group_pow_naf_r2l_u32(this->pow,
			this->temp_form, this->initform, this->order);

  // If this is the identity, then the initial form has odd order,
  // and the identity will be the only ambiguous form.
  if (group->is_id(group, this->temp_form)) {
    cprintf(2, "The primeform has an odd order and so only a trivial factorization.\n");
    return 0;
  }
    
  // Since we change primeforms, we aren't guaranteed that we
  // know the odd part of the order anymore.  Therefore, we must
  // limit the number of squares performed, since the remaining
  // order might be odd.
  for (i = 0;
       i < this->logh && !qgroup->is_ambiguous(qgroup, this->temp_form);
       i ++) {
    group->square(group, this->temp_form, this->temp_form);
  }

  // Try to split the ambiguous form.
  if (qgroup->split_ambiguous(qgroup, d, N, this->temp_form)) {
    if (verbose_level >= 2) {
      cprintf(2, "Successful with ", N);
      group->print(group, this->temp_form);
      cprintf(2, "\n");
    }
    return 1;
  } else {
    if (verbose_level >= 2) {
      cprintf(2, "Failed with ", N);
      group->print(group, this->temp_form);
      cprintf(2, "\n");
    }
    return 0;
  }
}

/**
 * Set all polymorphic variables and the group discriminant.
 * NOTE: Requires this->D and this->primorial_index to be set.
 */
static void set_discriminant(sspar_t* this) {
  int i;
  cprintf(2, "Discriminant: %Zd\n", this->D);
  size_t logD = mpz_sizeinbase(this->D, 2);
  if (logD <= s64_qform_group_max_bits) {
    // Use s64 implementations.
    this->qform_group          = &this->qform_group_s64.desc;
    this->initform             = &this->initform_s64;
    this->search               = &this->search_s64;
    this->giant                = &this->giant_s64;
    this->current              = &this->current_s64;
    this->temp_form            = &this->temp_form_s64;
    this->pow                  = &this->pow_s64;
    this->primorial_terms      = s64_primorial_terms[this->primorial_index];
    this->primorial_term_count =
        s64_primorial_term_counts[this->primorial_index];
    for (i = 0; i < this->gap_capacity; i++) {
      this->gap[i] = &this->gap_s64[i];
    }
  } else if (logD <= s128_qform_group_max_bits) {
    // Use s128 implementations.
    this->qform_group          = &this->qform_group_s128.desc;
    this->initform             = &this->initform_s128;
    this->search               = &this->search_s128;
    this->giant                = &this->giant_s128;
    this->current              = &this->current_s128;
    this->temp_form            = &this->temp_form_s128;
    this->pow                  = &this->pow_s128;
    this->primorial_terms      = s128_primorial_terms[this->primorial_index];
    this->primorial_term_count =
        s128_primorial_term_counts[this->primorial_index];
    for (i = 0; i < this->gap_capacity; i++) {
      this->gap[i] = &this->gap_s128[i];
    }
  } else {
    // Use MPZ implementations.
    this->qform_group          = &this->qform_group_mpz.desc;
    this->initform             = &this->initform_mpz;
    this->search               = &this->search_mpz;
    this->giant                = &this->giant_mpz;
    this->current              = &this->current_mpz;
    this->temp_form            = &this->temp_form_mpz;
    this->pow                  = &this->pow_mpz;
    this->primorial_terms      = mpz_primorial_terms[this->primorial_index];
    this->primorial_term_count =
        mpz_primorial_term_counts[this->primorial_index];
    for (i = 0; i < this->gap_capacity; i++) {
      this->gap[i] = &this->gap_mpz[i];
    }
  }
  this->qform_group->set_discriminant(this->qform_group, this->D);
}

// Use additional_qforms unless it is equal to
// SSPAR_ADDITIONAL_QFORMS_AUTO, in which case, compute the value
// from the size of N.
static int determine_additional_qforms(const mpz_t N,
				       const int additional_qforms) {
  assert(additional_qforms_ >= 0 ||
	 additional_qforms_ == SSPAR_ADDITIONAL_QFORMS_AUTO);
  // Indexes correspond to 32, 40, 48, 56, 64, 72, and 80 bit inputs.
  static const int lookup[7] = {4, 6, 8, 9, 10, 11, 12};
  if (additional_qforms == SSPAR_ADDITIONAL_QFORMS_AUTO) {
    int i = (mpz_sizeinbase(N, 2) - 32 + 7) / 8;
    if (i < 0) return 4;
    if (i > 6) return i * 2 - 1;
    return lookup[i];
  }
  return additional_qforms;
}

/**
 * Attempt to find a factor of N.
 * 
 * This works by trying to find the order of a prime ideal that doesn't divide
 * -kN.  If we are not successful finding the order, we switch to the next
 * multiplier.  Otherwise, attempt to factor.  If we aren't successful,
 * try powering 2 other prime ideals.
 */
int sspar_factor_raw(sspar_t* this,
                     mpz_t d,
                     const mpz_t N,
                     const int primorial_index,
                     const int step_bound,
                     const int step_count,
                     const int step_first_coprime,
                     const uint8_t* step_coprime_deltas,
		     const int additional_qforms_) {
  assert(additional_qforms_ >= 0 ||
	 additional_qforms_ == SSPAR_ADDITIONAL_QFORMS_AUTO);
  int multiplier_index = 0;
  int prime_index = 0;
  unsigned int k = 1;
  int i = 0;
  int found_order = 0;
  int additional_qforms = determine_additional_qforms(N, additional_qforms_);
    
  cprintf(1, "Attempting to factor %Zd\n", N);
  sspar_setup(this,
	      primorial_index,
	      step_bound,
	      step_count,
	      step_first_coprime,
	      step_coprime_deltas);
    
  while (multiplier_index < square_free_count) {
    // Only use multipliers that are square-free and relatively prime to N.
    // Experiments showed that some multipliers work better than others
    // depending on N mod 4. So early on, we pick preferred multipliers.
    // After a while, we just use regular square free multipliers.
    if (multiplier_index < SPECIAL_MULTIPLIERS_COUNT) {
      k = (N->_mp_d[0]&3) == 1
	? special_multipliers_1mod4[multiplier_index]
	: special_multipliers_3mod4[multiplier_index];
    } else {
      k = square_free[multiplier_index];
    }
    cprintf(2, "Using multiplier %u\n", k);

    // Verify that k is not a divisor of N
    if ((k > 1) && (mpz_cmp_ui(N, k) != 0) && mpz_divisible_ui_p(N, k)) {
      // k is a non-trivial divisor of n
      mpz_set_ui(d, k);
      return 1;
    }

    // Set up discriminant for this multiplier
    mpz_mul_ui(this->D, N, k);
    if ((this->D->_mp_d[0] & 3) != 3) {
      // D != 3 (mod 4)
      mpz_mul_2exp(this->D, this->D, 2);
    }
    mpz_neg(this->D, this->D);
    set_discriminant(this);

    // NOTE: These are only valid after set_discriminant()
    qform_group_t* qgroup = this->qform_group;
    group_t* group = &qgroup->group;

    // Find logh such that 2^logh >= h_D
    // given the bound h_D < sqrt(|D|)*log(|D|)/pi
    this->logh = mpz_sizeinbase(this->D, 2);
    this->logh = (this->logh>>1) + numbits_s32(this->logh) - 1;
    if (this->logh < 1) this->logh = 1;
    cprintf(3, "logh=%d\n", this->logh);

    // Start with the first prime ideal >= 11.
    prime_index = next_primeform(this->qform_group,
				 this->initform,
				 FIRST_PRIME_INDEX);

    // Big exponentiation by product of odd prime powers (power primorial)
    exponentiate_primorial(this);

    // Now exponentiate by 2^logh
    group->square(group, this->search, this->initform);
    for (i = 0;  i < this->logh;  i ++) {
      group->square(group, this->search, this->search);
    }

    // Check if 'search' is the identity.
    this->order = 1;
    found_order = group->is_id(group, this->search);
    
    if (!found_order) {
      // Do a bounded primorial step search for the order.
      found_order = sspar_search(this);
    }

    if (!found_order) {
      // Search failed.  Switch to the next multiplier.
      if (verbose_level >= 2) {
	cprintf(2, "Failed to find the order of ");
	group->print(group, this->search);
	cprintf(2, " within a bound of %d.\n", this->step_bound);
      }
      this->order = 1;
      multiplier_index ++;
      continue;
    }

    // Use the odd order of initform to find an ambiguous form.
    this->order >>= lsb_u32(this->order);

    // Check if this leads to a factorization.
    if (test_ambiguous_form(this, d, N)) {
      return 1;
    }

    // Try additional qforms.
    for (i = additional_qforms; i; i--) {
      // Try the next prime ideal with the same order.
      prime_index = next_primeform(qgroup,
				   this->initform,
				   prime_index + 1);

      // Big exponentiation by product of odd prime powers (power primorial)
      exponentiate_primorial(this);

      // Check if this leads to a factorization.
      if (test_ambiguous_form(this, d, N)) {
	return 1;
      }
    }

    // We've tried the maximum number of different prime forms
    // and none of them lead to a factorization. We give up on
    // this multiplier.
    multiplier_index ++;
    // Check if we ran out of square free integers.
    if (multiplier_index >= square_free_count) {
      cprintf(1, "We failed to find a factor because we ran out of square free integers.\n");
      cprintf(1, "Consider increasing the list.\n");
    }
  }

  // We failed to factor the integer.
  return 0;
}

/**
 * Attempts to factor the integer N.
 * @param N The integer to be factored.
 * @param d A divisor of N.
 * @return 1 if a non-trivial factor of N is found, 0 otherwise.
 */
int sspar_factor(sspar_t* this, mpz_t d, const mpz_t N) {
  // Indexes are for 16, 18, ..., 100 bit numbers.
  // These values were determined empirically
  // using sspar_search/search_range.
  static const int primorial_indexes[43] = {
    2, 2, 2, 2, 2, 3, 4, 4,
    4, 4, 6, 6, 6, 8, 8, 10,
    14, 14, 18, 21, 23, 26, 34, 35,
    39, 52, 53, 63, 68, 86, 90, 104,
    109, 143, 148, 160, 173, 243, 247,
    258, 283, 378, 411};
  static const int steps_indexes[43] = {
    11, 4, 4, 5, 5, 7, 9, 11, 14, 16, 19,
    23, 23, 23, 29, 29, 29, 33, 33, 33, 33,
    37, 40, 40, 40, 47, 47, 47, 50, 52, 52, 59,
    60, 64, 64, 64, 64, 73, 73, 74, 77, 77, 85};
  
  const size_t logn = mpz_sizeinbase(N, 2);
  int i = ((logn + 1) >> 1) - 8;
  i = max(i, 0);
  i = min(i, 42);
  assert(i >=0 && i < 43);

  const int primorial = primorial_indexes[i];
  const int steps = steps_indexes[i];
  const step_bound_t* bound = &step_bounds[steps];
  //  printf("primorial=%d steps=%d\n", primorial, steps);
  return sspar_factor_raw(this, d, N, primorial,
			  bound->bound,
			  bound->step_count,
			  coprime_first_step[bound->odd_primorial_index],
			  coprime_steps[bound->odd_primorial_index],
			  SSPAR_ADDITIONAL_QFORMS_AUTO);
}
