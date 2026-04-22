//
// Created by tom on 2/5/26.
//
#include "../include/poly_utils.h"
#include <flint/flint.h>
#include <flint/fmpq.h>
#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>
#include <omp.h>

int count_sign_variations(fmpz_poly_t poly) {
  slong degree = fmpz_poly_degree(poly);
  if (degree < 0)
    return 0;
  int coeff_sign = 0, variations = 0, current_sign = 0;
  for (int i = 0; i < degree + 1; i++) {
    coeff_sign = fmpz_sgn(fmpz_poly_get_coeff_ptr(poly, i));
    if (current_sign == 0 && coeff_sign != 0) {
      current_sign = coeff_sign;
      continue;
    }
    if (current_sign * coeff_sign < 0) {
      variations++;
      current_sign = coeff_sign;
    }
  }
  return variations;
}

void shift_in_proportions_by_k(fmpz_poly_t outPoly, fmpz_poly_t poly, int k) {
  fmpz_t c;
  fmpz_init(c);
  slong degree = fmpz_poly_degree(poly);
  if (k == 0) {
    return;
  }

  int direction = (k > 0) ? 1 : 0;
  k = abs(k);

  int a = (direction) ? 0 : (degree * k);
  int b = (direction) ? 1 : (-1);

  for (int i = 0; i < degree + 1; i++) {
    fmpz_poly_get_coeff_fmpz(c, poly, i);
    fmpz_mul_2exp(c, c, a + b * (k * i));
    fmpz_poly_set_coeff_fmpz(outPoly, i, c);
  }
  fmpz_clear(c);
}

void cauchy_bound(fmpq_t bound, fmpz_poly_t in_poly) {
  fmpz_t height, lead;
  fmpz_init(height);
  fmpz_init(lead);

  fmpz_poly_height(height, in_poly);
  fmpz_poly_get_coeff_fmpz(lead, in_poly, fmpz_poly_degree(in_poly));
  fmpz_abs(lead, lead);
  fmpq_set_fmpz_frac(bound, height, lead);
  fmpq_add_ui(bound, bound, 1);

  fmpz_clear(height);
  fmpz_clear(lead);
}

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

slong fmpq_clog(fmpq_t n, int b) {
  slong b1 = fmpz_clog_ui(fmpq_numref(n), b);
  slong b2 = fmpz_flog_ui(fmpq_denref(n), b);
  return b1 - b2;
}

void neg_varchange(fmpz_poly_t out_poly, fmpz_poly_t in_poly) {
  fmpz_t c;
  fmpz_init(c);
  slong degree = fmpz_poly_degree(in_poly);
  for (int i = 0; i < (degree + 1); i++) {
    fmpz_poly_get_coeff_fmpz(c, in_poly, i);
    if ((i % 2) == 1) {
      fmpz_mul_si(c, c, -1);
    }
    fmpz_poly_set_coeff_fmpz(out_poly, i, c);
  }
  fmpz_clear(c);
}

void fmpz_poly_split_coeffs(fmpz_poly_t P1, fmpz_poly_t P2,
                            const fmpz_poly_t f) {
  slong n = fmpz_poly_length(f);

  fmpz_poly_fit_length(P1, n);
  fmpz_poly_fit_length(P2, n);

  _fmpz_poly_set_length(P1, n);
  _fmpz_poly_set_length(P2, n);

  for (slong i = 0; i < n; i++) {
    const fmpz *c = fmpz_poly_get_coeff_ptr(f, i);

    // P1 = floor(c / 2)
    fmpz_fdiv_q_2exp(P1->coeffs + i, c, 1);

    // P2 = c - P1
    fmpz_sub(P2->coeffs + i, c, P1->coeffs + i);
  }

  _fmpz_poly_normalise(P1);
  _fmpz_poly_normalise(P2);
}

void fmpz_poly_taylor_shift_split_parallel(fmpz_poly_t res, const fmpz_poly_t f,
                                           const fmpz_t a) {
  fmpz_poly_t P1, P2;
  fmpz_poly_init(P1);
  fmpz_poly_init(P2);

  fmpz_poly_split_coeffs(P1, P2, f);

  fmpz_poly_t R1, R2;
  fmpz_poly_init(R1);
  fmpz_poly_init(R2);

// parallel execution
#pragma omp parallel sections
  {
#pragma omp section
    {
      fmpz_poly_taylor_shift(R1, P1, a);
    }

#pragma omp section
    {
      fmpz_poly_taylor_shift(R2, P2, a);
    }
  }

  // combine
  fmpz_poly_add(res, R1, R2);

  fmpz_poly_clear(P1);
  fmpz_poly_clear(P2);
  fmpz_poly_clear(R1);
  fmpz_poly_clear(R2);
}