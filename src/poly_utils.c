//
// Created by tom on 2/5/26.
//
#include "../include/poly_utils.h"
#include <flint/fmpz.h>


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
    free(&oldCoeffs);
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
    fmpz_fdiv_q_ui(bound, bound, currMin);
}
