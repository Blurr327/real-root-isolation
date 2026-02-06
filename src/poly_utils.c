//
// Created by tom on 2/5/26.
//
#include "../include/poly_utils.h"
#include <flint/fmpz.h>
#include <flint/flint.h>


int count_sign_variations(fmpz_poly_t poly) {
    slong degree = fmpz_poly_degree(poly);

    int var = 0;
    for (int i = 0; i < degree-1; i++) {
        slong c0 = fmpz_poly_get_coeff_si(poly, i);
        slong c1 = fmpz_poly_get_coeff_si(poly, i+1);
        if (c0 * c1 < 0) {
            var++;
        }
    }
    return var;
}

void reverse_coeffs(fmpz_poly_t outPoly, fmpz_poly_t poly) {
    slong degree = fmpz_poly_degree(poly);
    fmpz_t *oldCoeffs = malloc(degree * sizeof(fmpz_t));
    for (int i = 0; i < degree; i++) {
        fmpz_poly_get_coeff_fmpz(oldCoeffs[i], poly, i);
    }
    for (int i = 0; i < degree; i++) {
        fmpz_poly_set_coeff_fmpz(outPoly, i, oldCoeffs[degree-i]);
    }
    free(oldCoeffs);
}

void shift_in_proportions_by_k(fmpz_poly_t poly, int k) {
    slong degree = fmpz_poly_degree(poly);
    int direction = (k > 0) ? 1 : 0;
    if (direction) {
        for (int i = 0; i < degree; i++) {
            fmpz_t c;
            fmpz_poly_get_coeff_fmpz(c, poly, i);
            fmpz_mul_2exp(c, c, k*i);
            fmpz_poly_set_coeff_fmpz(poly, i, c);
        }
    } else {
        for (int i = 0; i < degree; i++) {
            fmpz_t c;
            fmpz_poly_get_coeff_fmpz(c, poly, i);
            fmpz_fdiv_q_2exp(c, c, k*i);
            fmpz_poly_set_coeff_fmpz(poly, i, c);
        }
    }
}

void cauchy_bound(fmpz_t bound, fmpz_poly_t poly) {
    slong degree = fmpz_poly_degree(poly);
    fmpz_t *coeffs = malloc(degree * sizeof(fmpz_t));
    for (int i = 0; i < degree; i++) {
        fmpz_poly_get_coeff_fmpz(coeffs[i], poly, i);
    }
    fmpz_t currMin;
    fmpz_poly_get_coeff_fmpz(currMin, poly, 0);
    for (int i = 0; i < degree-1; i++) {
        fmpz_t tmp;
        fmpz_poly_get_coeff_fmpz(tmp, poly, i);
        if (fmpz_cmp(currMin, tmp) < 0) {
            fmpz_set(currMin, tmp);
        }
    }
    fmpz_poly_get_coeff_fmpz(bound, poly, degree);
    fmpz_fdiv_q_ui(bound, bound, *currMin);
    fmpz_add_si(bound, bound, 1);
}


/**
  Random dense FLINT integer polynomials (fmpz_poly_t)

  - Generates a dense polynomial of EXACT degree `deg` (i.e., leading coeff ≠ 0)
  - Coefficients are random signed integers with up to `bits` bits (roughly uniform)
  - Requires FLINT headers and linking with -lflint (and dependencies)

  Example usage:
    flint_rand_t state; flint_randinit(state);
    fmpz_poly_t f; fmpz_poly_init(f);
    random_dense_fmpz_poly(f, state, 20, 64);   // deg=20, 64-bit-ish coeffs
    fmpz_poly_print(f); printf("\n");
    fmpz_poly_clear(f); flint_randclear(state);
*/
void random_dense_fmpz_poly(fmpz_poly_t poly, flint_rand_t state, slong deg, flint_bitcnt_t bits)
{
    if (deg < 0)
    {
        fmpz_poly_zero(poly);
        return;
    }

    /* Allocate space for deg+1 coefficients (dense) */
    fmpz_poly_fit_length(poly, deg + 1);

    /* Fill all coefficients randomly */
    for (slong i = 0; i <= deg; i++)
    {
        /* Random signed integer with up to `bits` bits */
        fmpz_randbits(fmpz_poly_get_coeff_ptr(poly, i), state, bits);

        /* Randomize sign */
        if (n_randint(state, 2))
            fmpz_neg(fmpz_poly_get_coeff_ptr(poly, i), fmpz_poly_get_coeff_ptr(poly, i));
    }

    /* Ensure leading coefficient is nonzero so degree is exactly `deg` */
    if (fmpz_is_zero(fmpz_poly_get_coeff_ptr(poly, deg)))
    {
        /* Set to 1 (or -1) deterministically from RNG */
        fmpz_set_ui(fmpz_poly_get_coeff_ptr(poly, deg), 1);
        if (n_randint(state, 2))
            fmpz_neg(fmpz_poly_get_coeff_ptr(poly, deg), fmpz_poly_get_coeff_ptr(poly, deg));
    }

    /* Normalize internal length/degree */
    _fmpz_poly_set_length(poly, deg + 1);
    _fmpz_poly_normalise(poly);
}



