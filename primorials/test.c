/**
 * @file test.c
 * 
 * This program verifies the 2,3 representations of primorials.
 */

#include <gmp.h>
#include <inttypes.h>
#include <stdio.h>

#include "liboptarith/math_mpz.h"
#include "liboptarith/primorial.h"
#include "mpz_primorials.h"
#include "s64_primorials.h"
#include "s128_primorials.h"

#define MIN_PRIMORIAL 4
#define MAX_PRIMORIAL 2047
#define N (MAX_PRIMORIAL-MIN_PRIMORIAL+1)

void two_three_value(mpz_t result,
		     const factored_two_three_term16_t* terms,
		     const int term_count) {
  int i;
  uint16_t a;
  uint16_t b;
  mpz_t t;
  mpz_init(t);
  mpz_set_ui(result, 0);
  for (i = 0; i < term_count; i++) {
    a = terms[i].a;
    b = terms[i].b;
    
    // Test the high bit of b to determine if we add/subtract 3*b.
    if (b&(1<<15)) {
      b &= ~(1<<15);
      mpz_ui_pow_ui(t, 3, b);
      mpz_sub(result, result, t);
    } else {
      mpz_ui_pow_ui(t, 3, b);
      mpz_add(result, result, t);
    }

    // Multiply by 2^a.
    mpz_set_ui(t, 1);
    mpz_mul_2exp(t, t, a);
    mpz_mul(result, result, t);
  }
  mpz_clear(t);
}

// Pass in primorials so that the 0 index is primorial 0 (not primorial 4).
void verify(mpz_t* primorials,
	    const factored_two_three_term16_t* terms[],
	    const uint16_t counts[]) {
  mpz_t t;
  int i;
  mpz_init(t);
  for (i = 0; i < N; i++) {
    two_three_value(t, terms[i+MIN_PRIMORIAL], counts[i+MIN_PRIMORIAL]);
    if (mpz_cmp(primorials[i], t) != 0) {
      gmp_printf("Expected %Zd but got %Zd\n", primorials[i], t);
      exit(0);
    }
  }
  mpz_clear(t);
}

int main(int argc, char** argv) {
  // Compute the full set of primorials.
  mpz_t* primorials = mpz_primorials(N, 11);
  int i;
  for (i = 0; i < N; i++) {
    // Multiply in 3*3*3*5*5*7*7
    mpz_mul_ui(primorials[i], primorials[i], 3*3*3*5*5*7*7);
  }

  // Verify primorials.
  verify(primorials,
	 s64_primorial_terms,
	 s64_primorial_term_counts);
  verify(primorials,
	 s128_primorial_terms,
	 s128_primorial_term_counts);
  verify(primorials,
	 mpz_primorial_terms,
	 mpz_primorial_term_counts);

  mpz_clear_array(primorials, N);

  printf("All tests passed.\n");
  return 0;
}


