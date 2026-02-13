//
// Created by tom on 2/5/26.
//
#include "../include/poly_utils.h"
#include <flint/flint.h>
#include <flint/fmpq.h>
#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>

// #define DEBUG 0

int count_sign_variations(fmpz_poly_t poly) {
  slong degree = fmpz_poly_degree(poly);

  int var = 0;
  for (int i = 0; i < degree; i++) {
    int s0 = fmpz_sgn(fmpz_poly_get_coeff_ptr(poly, i));
    int s1 = fmpz_sgn(fmpz_poly_get_coeff_ptr(poly, i + 1));
    /*
    slong c0 = fmpz_poly_get_coeff_si(poly, i);
    slong c1 = fmpz_poly_get_coeff_si(poly, i+1);
    */
    if (s0 * s1 < 0) {
      var++;
    }
  }
  return var;
}

void reverse_coeffs(fmpz_poly_t out, fmpz_poly_t poly) {
  slong deg = fmpz_poly_degree(poly);

  if (deg < 0) { // zero polynomial
    fmpz_poly_zero(out);
    return;
  }

  // Ensure out has enough length
  fmpz_poly_fit_length(out, deg + 1);
  fmpz_poly_zero(out);

  fmpz_t c;
  fmpz_init(c);

  for (slong i = 0; i <= deg; i++) {
    fmpz_poly_get_coeff_fmpz(c, poly, i);
    fmpz_poly_set_coeff_fmpz(out, deg - i, c);
  }

  fmpz_clear(c);
  _fmpz_poly_normalise(out);
}

void shift_in_proportions_by_k(fmpz_poly_t outPoly, fmpz_poly_t poly, int k) {
  slong degree = fmpz_poly_degree(poly);
  if (k == 0) {
    return;
  }

  int direction = (k > 0) ? 1 : 0;
  k = abs(k);

  int a = (direction) ? 0 : degree;
  int b = (direction) ? 1 : (-1);

  for (int i = 0; i < degree + 1; i++) {
    fmpz_t c;
    fmpz_poly_get_coeff_fmpz(c, poly, i);
    fmpz_mul_2exp(c, c, a + b * (k * i));
    fmpz_poly_set_coeff_fmpz(outPoly, i, c);
  }

  _fmpz_poly_normalise(outPoly);
}

void cauchy_bound(fmpq_t bound, fmpz_poly_t poly) {
  if (fmpz_poly_is_zero(poly)) {
    fmpq_zero(bound);
    return;
  }
  slong degree = fmpz_poly_degree(poly);
  fmpz_t currMax;
  fmpz_init(currMax);
  fmpz_poly_get_coeff_fmpz(currMax, poly, 0);
  fmpz_abs(currMax, currMax);
  for (int i = 0; i < degree; i++) {
    fmpz_t tmp;
    fmpz_poly_get_coeff_fmpz(tmp, poly, i);
    fmpz_abs(tmp, tmp);
    // printf("current evaluated coeff: %d: ", i);
    // fmpz_print(tmp);
    // printf("\n");
    if (fmpz_cmp(currMax, tmp) < 0) {
      fmpz_set(currMax, tmp);
    }
  }
  // fmpz_t maxDegreeCoeff = fmpz_poly_lead(poly);
  fmpz_t maxDegreeCoeff;
  fmpz_poly_get_coeff_fmpz(maxDegreeCoeff, poly, degree);

  fmpz_abs(maxDegreeCoeff, maxDegreeCoeff);

  fmpq_t currMaxQ;
  fmpz_t one;
  fmpz_one(one);
  fmpq_set_fmpz_frac(currMaxQ, currMax, one);
  fmpq_div_fmpz(bound, currMaxQ, maxDegreeCoeff);
  // fmpq_div_fmpz(bound, currMax, bound);
  fmpq_add_si(bound, bound, 1);
  fmpz_clear(currMax);
}

/**
  Random dense FLINT integer polynomials (fmpz_poly_t)

  - Generates a dense polynomial of EXACT degree `deg` (i.e., leading coeff ≠ 0)
  - Coefficients are random signed integers with up to `bits` bits (roughly
  uniform)
  - Requires FLINT headers and linking with -lflint (and dependencies)

  Example usage:
    flint_rand_t state; flint_randinit(state);
    fmpz_poly_t f; fmpz_poly_init(f);
    random_dense_fmpz_poly(f, state, 20, 64);   // deg=20, 64-bit-ish coeffs
    fmpz_poly_print(f); printf("\n");
    fmpz_poly_clear(f); flint_randclear(state);
*/
void random_dense_fmpz_poly(fmpz_poly_t poly, flint_rand_t state, slong deg,
                            flint_bitcnt_t bits) {
  if (deg < 0) {
    fmpz_poly_zero(poly);
    return;
  }

  fmpz_t temp_coeff;
  fmpz_init(temp_coeff);

  /* Fill all coefficients */
  for (slong i = 0; i <= deg; i++) {
    fmpz_randbits(temp_coeff, state, bits);

    if (n_randint(state, 2))
      fmpz_neg(temp_coeff, temp_coeff);

    /* Ensure leading coefficient is nonzero for exact degree */
    if (i == deg && fmpz_is_zero(temp_coeff)) {
      fmpz_set_si(temp_coeff, (n_randint(state, 2) ? 1 : -1));
    }

    /* This is the safe way to set coefficients */
    fmpz_poly_set_coeff_fmpz(poly, i, temp_coeff);
  }

  fmpz_clear(temp_coeff);
}