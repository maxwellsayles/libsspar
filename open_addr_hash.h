/**
 * @file hash_table.h
 *
 * Open addressing hash tables. Loosely based on the implementation in
 * Python. We only store the hash value and not the reference object,
 * so if two objects hash to the same value, only the latter one is
 * remembered.  A special key/value pair signifies an empty slot and
 * in the unlikely event that this pair is generated, the effect is to
 * delete any existing entry or to do nothing.
 *
 * This class is intended to be used with a baby-step/giant-step search
 * that is tolerant of missing table entries.
 *
 * Written by Maxwell Sayles, 2013.
 */
#pragma once
#ifndef OPEN_ADDR_HASH__INCLUDED
#define OPEN_ADDR_HASH__INCLUDED

#include <stdint.h>

#include "liboptarith/math32.h"

typedef struct {
  uint32_t k;
  uint32_t v;
} htab_entry_t;

typedef struct {
  htab_entry_t* tab;
  uint32_t mask;
  uint32_t max_entries;
} hash_table_t;

// Choose a key/value pair that is unlikely to occur.
#define HASH_KEY_EMPTY 0xAAAAAAAA
#define HASH_VAL_EMPTY 0x55555555

/// Reset a hash table to its initial state (i.e. no members).
/// NOTE: entries is rounded up to the nearest power of 2
static inline void hash_table_reset(hash_table_t* hash, const int entries) {
  const uint32_t size = ceil_pow2_u32(entries);
  assert(size <= hash->max_entries);
  hash->mask = size - 1;
  uint32_t i;
  for (i = 0; i < size; i++) {
    hash->tab[i].k = HASH_KEY_EMPTY;
    hash->tab[i].v = HASH_VAL_EMPTY;
  }
}

/// Allocate resources for use by the hash table.
/// NOTE: max_entries is rounded up to the nearest power of 2.
static inline void hash_table_init(hash_table_t* hash,
				   const uint32_t max_entries) {
  hash->max_entries = ceil_pow2_u32(max_entries);
  hash->tab = malloc(sizeof(htab_entry_t) * hash->max_entries);
}

/// Releases the resources held by the hash table.
/// It is named 'clear' for consistency with GMP.
/// NOTE: The hash table will not be usable after this call.
static inline void hash_table_clear(hash_table_t* hash) {
  free(hash->tab);
  hash->tab = 0;
  hash->mask = 0;
  hash->max_entries = 0;
}

// Private function.
static inline uint32_t _hash_table_find_slot(hash_table_t* hash,
					     const uint32_t key) {
  uint32_t i = key & hash->mask;
  uint32_t perturb = key;
  while ((hash->tab[i].k != key) &&
	 (hash->tab[i].k != HASH_KEY_EMPTY || hash->tab[i].v != HASH_VAL_EMPTY)) {
    i = (5 * i + 1 + perturb) & hash->mask;
    perturb >>= 5;
  }
  return i;
}

/// Insert an element into a table.
/// NOTE: datum and key cannot be VAL_EMPTY and KEY_EMPTY respectively,
///       (although the probability of this is quite low).
static inline void hash_table_insert(hash_table_t* hash,
				     const uint32_t datum,
				     const uint32_t key) {
  const uint32_t i = _hash_table_find_slot(hash, key);
  hash->tab[i].k = key;
  hash->tab[i].v = datum;
}

/// Returns the number of entries placed in data.
/// data_size is assumed to be >= 2. There is no error checking done.
static inline int hash_table_lookup(hash_table_t* hash,
				    uint32_t* data, const int data_size,
				    const uint32_t key) {
  const uint32_t i = _hash_table_find_slot(hash, key);
  if (hash->tab[i].k == HASH_KEY_EMPTY &&
      hash->tab[i].v == HASH_VAL_EMPTY) {
    return 0;
  }
  *data = hash->tab[i].v;
  return 1;
}

#endif  // OPEN_ADDR_HASH__INCLUDED

